import argparse
import ast
import json
from pathlib import Path
import random
import string

import gget
import ArangoDB as adb
import pandas as pd
import scanpy as sc

from CellOntology import parse_term


ALPHABET = string.ascii_lowercase + string.digits


def get_uuid():
    """Get an eight character random string.

    Parameters
    ----------
    None

    Returns
    -------
    An eight character random string.
    """
    return "".join(random.choices(ALPHABET, k=8))


# TODO: Use HGNC web service and compare
# See: https://www.genenames.org/help/rest/#!/#tocAnchor-1-1-2
def get_gene_name_to_ids_map():
    """Query BioMart to map gene names to ids.

    Parameters
    ----------
    None

    Returns
    -------
    nm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene id
    """
    nm2ids = sc.queries.biomart_annotations(
        "hsapiens", ["ensembl_gene_id", "external_gene_name"], use_cache=True
    ).set_index("external_gene_name")

    return nm2ids


def map_gene_name_to_ids(name, nm2ids):
    """Map a gene name to a gene id list.

    Parameters
    ----------
    name : str
        Gene name
    nm2ids : pd.DataFrame
        DataFrame indexed by name containing id

    Returns
    -------
    list
        Gene ids
    """
    if name in nm2ids.index:
        id = nm2ids.loc[name, "ensembl_gene_id"]
        if isinstance(id, pd.core.series.Series):
            return id.to_list()
        else:
            return [id]
    else:
        print(f"Could not find ids for name: {name}")
        return []


def load_mdata(mdata_dirname, mdata_filename):
    """Load manually curated data.

    Parameters
    ----------
    mdata_dirname : Path
        Name of directory containing manually curated data
    mdata_filename : str
        Name of Excel file containing manually curated data

    Returns
    -------
    mdata : pd.DataFrame
        DataFrame containing manually curated data
    """
    mdata = pd.read_excel(mdata_dirname / mdata_filename, header=1, skiprows=0)
    mdata["uuid"] = [get_uuid() for idx in mdata.index]
    return mdata


def init_gdata(mdata, nm2ids):
    """Initialize gget data.

    Parameters
    ----------
    mdata : pd.DataFrame
        DataFrame containing manually curated data

    Returns
    -------
    gdata : dict
        Dictionary containing cell types, marker names, and marker ids
        for each marker name all for each source publication, and all
        unique marker ids
    """
    # Collect cell types, and marker names and ids for each source
    gdata = {}
    gdata["source"] = {}
    gdata["marker_ids"] = set()
    for source in ["HLCA", "CellRef"]:
        gdata["source"][source] = {}

        # Iterate through all rows for each source
        for _, row in mdata.iterrows():

            # Collect cell types
            if "cell_types" not in gdata["source"][source]:
                gdata["source"][source]["cell_types"] = {}
            if pd.isna(row[source + "_cellset"]):
                continue
            cell_type = row[source + "_cellset"]

            # Collect marker names
            if pd.isna(row[source + "_NSForestMarkers"]):
                continue
            marker_names = ast.literal_eval(row[source + "_NSForestMarkers"])
            if cell_type not in gdata["source"][source]["cell_types"]:
                gdata["source"][source]["cell_types"][cell_type] = {}
            if "marker_names" not in gdata["source"][source]["cell_types"][cell_type]:
                gdata["source"][source]["cell_types"][cell_type]["marker_names"] = {}
            for marker_name in marker_names:

                # Collect marker id
                marker_ids = map_gene_name_to_ids(marker_name, nm2ids)
                if (
                    marker_name
                    not in gdata["source"][source]["cell_types"][cell_type][
                        "marker_names"
                    ]
                ):
                    gdata["source"][source]["cell_types"][cell_type]["marker_names"][
                        marker_name
                    ] = {}
                gdata["source"][source]["cell_types"][cell_type]["marker_names"][
                    marker_name
                ]["marker_ids"] = marker_ids
                gdata["marker_ids"] |= set(marker_ids)

    # Marker ids are more convenient as a list
    gdata["marker_ids"] = list(gdata["marker_ids"])

    return gdata


def load_gdata(mdata, gdata_dirname, gdata_filename, nm2ids):
    """Call gget for each gene id, collecting all specified resources.

    Parameters
    ----------
    mdata : pd.DataFrame
        DataFrame containing manually curated data
    gdata_dirname : Path
        Name of directory containing gget data
    gdata_filename : str
        Name of JSON file containing gget data

    Returns
    -------
    gdata : dict
        Dictionary containing HLCA cell types, marker names, marker
        ids for each marker name, all specified resources for each
        marker id, and all unique marker ids
    """
    # Create and dump, or load gget data
    if not (gdata_dirname / gdata_filename).exists():
        # Create gget data
        gdata = init_gdata(mdata, nm2ids)
        resources = [
            "diseases",
            "drugs",
            "interactions",
            "pharmacogenetics",
            "tractability",
            "expression",
            "depmap",
        ]

        # Collect gene ids
        gdata["gene_ids"] = {}
        for id in gdata["marker_ids"]:
            gdata["gene_ids"][id] = {}

            # Collect resources
            gdata["gene_ids"][id]["resources"] = {}
            for resource in resources:
                try:
                    gdata["gene_ids"][id]["resources"][resource] = gget.opentargets(
                        id, resource=resource, json=True, verbose=True
                    )
                except Exception as exc:
                    print(f"Could not gget resource: {resource} for id: {id}")
                    gdata["gene_ids"][id]["resources"][resource] = {}

        # Dump gget data
        with open(gdata_dirname / gdata_filename, "w") as fp:
            json.dump(gdata, fp, indent=4)

    else:

        # Load gget data
        with open(gdata_dirname / gdata_filename, "r") as fp:
            gdata = json.load(fp)

    return gdata


def get_graph(db_name, graph_name):
    """Get an graph from an ArangoDB database.

    Parameters
    ----------
    db_name : str
        Name of ArangoDB database
    graph_name : str
        Name of ArangoDB database graph

    Returns
    -------
    adb_graph : arango.graph.Graph
        ArangoDB database graph
    """
    db = adb.create_or_get_database(db_name)
    adb_graph = adb.create_or_get_graph(db, graph_name)
    return adb_graph


def init_collections(adb_graph):
    """Define and create vertex and edge collections.

    Parameters
    ----------
    adb_graph : arango.graph.Graph
        ArangoDB database graph

    Returns
    -------
    vertex_collections : dict
        A dictionary with vertex name keys containing
        arango.collection.VertexCollection instance values
    edge_collections : dict
        A dictionary with edge name keys containing
        arango.collection.EdgeCollection instance values
    """
    # Define and create vertex collections
    vertex_collections = {}
    vertex_names = [
        "CL",  # Equivalent to cell_set_cls
        "anatomic_structure_cls",
        "gene_cls",
        "transcript_cls",
        "publication_ind",
        "publication_cls",
        "cell_set_ind",
        "transcript_ind",
        "biomarker_combination_ind",
        "biomarker_combination_cls",
        "disease_cls",
        "drug_product_cls",
    ]
    for vertex_name in vertex_names:
        collection = adb.create_or_get_vertex_collection(adb_graph, vertex_name)
        vertex_collections[vertex_name] = collection

    # Define and create edge collections
    edge_collections = {}
    vertex_name_pairs = [
        ("CL", "anatomic_structure_cls"),
        ("CL", "gene_cls"),
        ("gene_cls", "transcript_cls"),
        ("publication_ind", "publication_cls"),
        ("cell_set_ind", "CL"),
        ("transcript_ind", "transcript_cls"),
        ("biomarker_combination_ind", "biomarker_combination_cls"),
        ("cell_set_ind", "publication_ind"),
        ("cell_set_ind", "transcript_ind"),
        ("transcript_ind", "biomarker_combination_ind"),
        ("biomarker_combination_ind", "cell_set_ind"),
        ("CL", "CL"),
        ("gene_cls", "disease_cls"),
        ("drug_product_cls", "gene_cls"),
        ("drug_product_cls", "disease_cls"),
    ]
    for vertex_name_pair in vertex_name_pairs:
        from_vertex = vertex_name_pair[0]
        to_vertex = vertex_name_pair[1]
        collection, edge_name = adb.create_or_get_edge_collection(
            adb_graph, from_vertex, to_vertex
        )
        edge_collections[edge_name] = collection

    return vertex_collections, edge_collections


def insert_individual_vertices(vertex_collections):
    """Insert the one anatomic structure class, two publication
    individual, one publication class, and one biomarker combination
    class vertices.

    Parameters
    ----------
    vertex_collections : dict
        A dictionary with vertex name keys containing
        arango.collection.VertexCollection instance values

    Returns
    -------
    anatomic_structure_cls_vertex : dict
        The ArangoDB anatomic structure vertex document
    publication_hlca_vertex : dict
        The ArangoDB HLCA publication vertex document
    publication_cellref_vertex,
        The ArangoDB CellRef publication vertex document
    publication_cls_vertex,
        The ArangoDB publication class vertex document
    biomarker_combination_cls_vertex
        The ArangoDB biomarker combination class vertex document
    """
    # Anatomic structure vertex
    _key = "Organ_12345"
    if not vertex_collections["anatomic_structure_cls"].has(_key):
        vertex_collections["anatomic_structure_cls"].insert(
            {"_key": _key, "label": "lung"}
        )
    anatomic_structure_cls_vertex = vertex_collections["anatomic_structure_cls"].get(
        _key
    )

    # Publication individual and class vertices
    _key = "Sikkema-et-al-2023-Nat-Med"
    if not vertex_collections["publication_ind"].has(_key):
        vertex_collections["publication_ind"].insert(
            {
                "_key": _key,
                "label": "HLCA",
                "citation": "Sikkema, L., Ramírez-Suástegui, C., Strobl, D.C. et al. An integrated cell atlas of the lung in health and disease. Nat Med 29, 1563–1577 (2023).",
                "DOI": "https://doi.org/10.1038/s41591-023-02327-2",
            }
        )
    publication_hlca_vertex = vertex_collections["publication_ind"].get(_key)
    _key = "Guo-et-al-2023-Nat-Commun"
    if not vertex_collections["publication_ind"].has(_key):
        vertex_collections["publication_ind"].insert(
            {
                "_key": _key,
                "label": "CellRef",
                "citation": "Guo, M., Morley, M.P., Jiang, C. et al. Guided construction of single cell reference for human and mouse lung. Nat Commun 14, 4566 (2023).",
                "DOI": "https://doi.org/10.1038/s41467-023-40173-5",
            }
        )
    publication_cellref_vertex = vertex_collections["publication_ind"].get(_key)
    _key = "NLMP00000000001"
    if not vertex_collections["publication_cls"].has(_key):
        vertex_collections["publication_cls"].insert(
            {
                "_key": _key,
                "label": "NLM Publication",
            }
        )
    publication_cls_vertex = vertex_collections["publication_cls"].get(_key)

    # Biomarker combination class vertex
    _key = "NLMB00000000001"
    if not vertex_collections["biomarker_combination_cls"].has(_key):
        vertex_collections["biomarker_combination_cls"].insert(
            {
                "_key": _key,
                "label": "NLM Biomarker Combination",
            }
        )
    biomarker_combination_cls_vertex = vertex_collections[
        "biomarker_combination_cls"
    ].get(_key)

    return (
        anatomic_structure_cls_vertex,
        publication_hlca_vertex,
        publication_cellref_vertex,
        publication_cls_vertex,
        biomarker_combination_cls_vertex,
    )


def insert_triples(triples, edge_collections):
    """Insert triples into the ArangoDB graph.

    Parameters
    ----------
    triples : list(tuple)
    edge_collections : dict
        A dictionary with edge name keys containing
        arango.collection.EdgeCollection instance values

    Returns
    -------
    None
    """
    for triple in triples:

        from_vertex = triple[0]  # subject
        predicate = triple[1]
        to_vertex = triple[2]  # object

        # All vertices are valid or an empty dictionary
        if from_vertex == {} or to_vertex == {}:
            continue

        # Define edge, noting that all predicates are arbitrary dictionaries
        edge = {
            "_key": from_vertex["_key"] + "-" + to_vertex["_key"],
            "_from": from_vertex["_id"],
            "_to": to_vertex["_id"],
        }
        edge.update(predicate)

        # Insert triple
        from_vertex_collection = from_vertex["_id"].split("/")[0]
        to_vertex_collection = to_vertex["_id"].split("/")[0]
        edge_name = from_vertex_collection + "-" + to_vertex_collection
        if not edge_collections[edge_name].has(edge):
            edge_collections[edge_name].insert(edge)


def insert_vertices_and_edges_from_mdata_row(
    row, vertex_collections, edge_collections, nm2ids, g2t
):
    """Insert vertices and edges from each row of the manually curated
    data.

    Parameters
    ----------
    row : pd.Series
        Row Series from the DataFrame containing manually curated data
    vertex_collections : dict
        A dictionary with vertex name keys containing
        arango.collection.VertexCollection instance values
    edge_collections : dict
        A dictionary with edge name keys containing
        arango.collection.EdgeCollection instance values
    g2t : dict
        Dictionary mapping gene to transcript class

    Returns
    -------
    None
    """
    # Insert all individual vertices
    (
        anatomic_structure_cls_vertex,
        publication_hlca_vertex,
        publication_cellref_vertex,
        publication_cls_vertex,
        biomarker_combination_cls_vertex,
    ) = insert_individual_vertices(vertex_collections)

    # Insert cell set class vertices
    cell_set_cls_vertex = {}
    if not (pd.isna(row["CL_cell_type"]) or pd.isna(row["CL_PURL"])):

        # Get the number from the current CL term
        number = None
        CL_PURL = row["CL_PURL"]
        if "http" in CL_PURL:
            _, number, _, _, _ = parse_term(CL_PURL)

        elif ":" in CL_PURL:
            number = CL_PURL.split(":")[1]

        else:
            print(f"Could not get number from CL PURL: {CL_PURL}")

        # Update the CL term vertex, if it, and the proposed addition
        # exists
        _key = number
        cell_set_cls_vertex = vertex_collections["CL"].get(_key)
        if cell_set_cls_vertex and pd.isna(
            row["Proposed addition to CL definition or annotation property"]
        ):
            cell_set_cls_vertex["proposed definition"] = str(
                row["Proposed addition to CL definition or annotation property"]
            )
            vertex_collections["CL"].update(cell_set_cls_vertex)
        else:
            cell_set_cls_vertex = {}

    # Insert biomarker combination, and transcript individual vertices
    hlca_transcript_names = []
    biomarker_combination_hlca_vertex = {}
    if not pd.isna(row["HLCA_NSForestMarkers"]):
        hlca_transcript_names = ast.literal_eval(row["HLCA_NSForestMarkers"])
        _key = "hlca-" + row["uuid"]
        if not vertex_collections["biomarker_combination_ind"].has(_key):
            vertex_collections["biomarker_combination_ind"].insert(
                {"_key": _key, "label": hlca_transcript_names}
            )
        biomarker_combination_hlca_vertex = vertex_collections[
            "biomarker_combination_ind"
        ].get(_key)
    cellref_transcript_names = []
    biomarker_combination_cellref_vertex = {}
    if not pd.isna(row["CellRef_NSForestMarkers"]):
        cellref_transcript_names = ast.literal_eval(row["CellRef_NSForestMarkers"])
        _key = "cellref-" + row["uuid"]
        if not vertex_collections["biomarker_combination_ind"].has(_key):
            vertex_collections["biomarker_combination_ind"].insert(
                {
                    "_key": _key,
                    "label": cellref_transcript_names,
                }
            )
        biomarker_combination_cellref_vertex = vertex_collections[
            "biomarker_combination_ind"
        ].get(_key)

    # Insert transcript individual and class, and gene vertices
    transcript_ind_vertices = {}
    for transcript in hlca_transcript_names:
        _key = "hlca-" + transcript
        if not vertex_collections["transcript_ind"].has(_key):
            vertex_collections["transcript_ind"].insert(
                {"_key": _key, "label": transcript, "publication": "HLCA"}
            )
        transcript_ind_vertices[transcript] = vertex_collections["transcript_ind"].get(
            _key
        )
    transcript_ind_vertices = {}
    for transcript in cellref_transcript_names:
        _key = "cellref-" + transcript
        if not vertex_collections["transcript_ind"].has(_key):
            vertex_collections["transcript_ind"].insert(
                {"_key": _key, "label": transcript, "publication": "CellRef"}
            )
        transcript_ind_vertices[transcript] = vertex_collections["transcript_ind"].get(
            _key
        )
    transcript_cls_vertices = {}
    gene_cls_vertices = {}
    if "count" not in g2t:
        g2t["count"] = 0
    for transcript in hlca_transcript_names + cellref_transcript_names:
        gene = transcript
        g2t["count"] += 1

        # Transcript class
        _key = f"NLMT{g2t['count']:011d}"
        g2t[gene] = _key  # TODO: Not used, but might be useful later. Remove?
        if not vertex_collections["transcript_cls"].has(_key):
            vertex_collections["transcript_cls"].insert(
                {
                    "_key": _key,
                    "label": transcript,
                }
            )
        transcript_cls_vertices[gene] = vertex_collections["transcript_cls"].get(_key)

        # Gene class
        ids = map_gene_name_to_ids(gene, nm2ids)
        if len(ids) == 0:
            gene_cls_vertices[gene] = {}
        else:
            _key = ids[0]
            if not vertex_collections["gene_cls"].has(_key):
                vertex_collections["gene_cls"].insert(
                    {
                        "_key": _key,
                        "label": gene,
                    }
                )
            gene_cls_vertices[gene] = vertex_collections["gene_cls"].get(_key)

    # Insert cell set vertices
    cell_set_hlca_vertex = {}
    if not pd.isna(row["HLCA_cellset"]):
        _key = "hlca-" + row["uuid"]
        if not vertex_collections["cell_set_ind"].has(_key):
            vertex_collections["cell_set_ind"].insert(
                {
                    "_key": _key,
                    "label": row["HLCA_cellset"],
                }
            )
        cell_set_hlca_vertex = vertex_collections["cell_set_ind"].get(_key)
    cell_set_cellref_vertex = {}
    if not pd.isna(row["CellRef_cellset"]):
        _key = "cellref-" + row["uuid"]
        if not vertex_collections["cell_set_ind"].has(_key):
            vertex_collections["cell_set_ind"].insert(
                {
                    "_key": _key,
                    "label": row["CellRef_cellset"],
                }
            )
        cell_set_cellref_vertex = vertex_collections["cell_set_ind"].get(_key)

    # Initialize triples
    triples = []

    # Define and collect (cell_set_cls (CL), PART_OF, anatomical_structure_cls) triples
    triples.append(
        (cell_set_cls_vertex, {"label": "PART_OF"}, anatomic_structure_cls_vertex),
    )

    # Define and collect (cell_set_cls (CL), EXPRESSES, gene_cls) triples
    for transcript, gene_cls_vertex in gene_cls_vertices.items():
        triples.append(
            (
                cell_set_cls_vertex,
                {"label": "EXPRESSES", "is_shortcut": True},
                gene_cls_vertex,
            )
        )

        # Define and collect (gene_cls, PRODUCES, transcript_cls) triples
        transcript_cls_vertex = transcript_cls_vertices[transcript]
        triples.append((gene_cls_vertex, {"label": "PRODUCES"}, transcript_cls_vertex))

    # Define and collect (publication_ind, IS_INSTANCE, publication_cls) triples
    triples.extend(
        [
            (publication_hlca_vertex, {"label": "IS_INSTANCE"}, publication_cls_vertex),
            (
                publication_cellref_vertex,
                {"label": "IS_INSTANCE"},
                publication_cls_vertex,
            ),
        ]
    )

    # Define and collect (cell_set_ind, IS_INSTANCE, cell_set_cls (CL)) triples
    if row["predicate_HLCA"] in [
        "skos:broadMatch",
        "skos:exactMatch",
        "skos:relatedMatch",
    ]:
        triples.append(
            (
                cell_set_hlca_vertex,
                {"label": "IS_INSTANCE", "match": row["predicate_HLCA"]},
                cell_set_cls_vertex,
            ),
        )
    if row["predicate_CellRef"] in [
        "skos:broadMatch",
        "skos:exactMatch",
        "skos:relatedMatch",
    ]:
        triples.append(
            (
                cell_set_cellref_vertex,
                {"label": "IS_INSTANCE", "match": row["predicate_CellRef"]},
                cell_set_cls_vertex,
            ),
        )

    # Define and collect (transcript_ind, IS_INSTANCE, transcript_cls) triples
    for transcript, transcript_ind_vertex in transcript_ind_vertices.items():
        transcript_cls_vertex = transcript_cls_vertices[transcript]
        triples.append(
            (transcript_ind_vertex, {"label": "IS_INSTANCE"}, transcript_cls_vertex)
        )

    # Define and collect (biomarker_combination_ind, IS_INSTANCE, biomarker_combination_cls) triples
    triples.extend(
        [
            (
                biomarker_combination_hlca_vertex,
                {"label": "IS_INSTANCE"},
                biomarker_combination_cls_vertex,
            ),
            (
                biomarker_combination_cellref_vertex,
                {"label": "IS_INSTANCE"},
                biomarker_combination_cls_vertex,
            ),
        ]
    )

    # Define and collect (cell_set_ind, SOURCE, publication_ind) triples
    triples.extend(
        [
            (cell_set_hlca_vertex, {"label": "SOURCE"}, publication_hlca_vertex),
            (cell_set_cellref_vertex, {"label": "SOURCE"}, publication_cellref_vertex),
        ]
    )

    # Define and collect:
    #   (cell_set_ind, HAS_PART, transcript_ind), and
    #   (transcript_ind, MEMBER_OF, biomarker_combination_ind) triples
    for transcript_ind_vertex in transcript_ind_vertices.values():
        publication = transcript_ind_vertex["publication"]
        if publication == "HLCA":
            triples.extend(
                [
                    (
                        cell_set_hlca_vertex,
                        {"label": "HAS_PART"},
                        transcript_ind_vertex,
                    ),
                    (
                        transcript_ind_vertex,
                        {"label": "MEMBER_OR"},
                        biomarker_combination_hlca_vertex,
                    ),
                ]
            )

        elif publication == "CellRef":
            triples.extend(
                [
                    (
                        cell_set_cellref_vertex,
                        {"label": "HAS_PART"},
                        transcript_ind_vertex,
                    ),
                    (
                        transcript_ind_vertex,
                        {"label": "MEMBER_OF"},
                        biomarker_combination_cellref_vertex,
                    ),
                ]
            )

        else:
            print(f"Unknown publication: {publication}")

    # Define and collect (biomarker_combination_ind, IS_MARKER_FOR, cell_set_ind) triples
    triples.extend(
        [
            (
                biomarker_combination_hlca_vertex,
                {"label": "IS_MARKER_FOR", "Fbeta": row["HLCA_Fbeta"]},
                cell_set_hlca_vertex,
            ),
            (
                biomarker_combination_cellref_vertex,
                {"label": "IS_MARKER_FOR", "Fbeta": row["CellRef_Fbeta"]},
                cell_set_cellref_vertex,
            ),
        ]
    )

    # Insert all triples
    insert_triples(triples, edge_collections)


def insert_vertices_and_edges_from_gdata_for_gene_name_and_id(
    gdata, gene_name, gene_id, vertex_collections, edge_collections
):
    """Insert vertices and edges for each gene name and id from the
    gget data.

    Parameters
    ----------
    gdata : dict
        Dictionary containing cell types, marker names, and marker ids
        for each marker name all for each source publication, and all
        unique marker ids
    gene_name : str
        Gene name
    gene_id : str
        Ensemble gene stable identifier
    vertex_collections : dict
        A dictionary with vertex name keys containing
        arango.collection.VertexCollection instance values
    edge_collections : dict
        A dictionary with edge name keys containing
        arango.collection.EdgeCollection instance values

    Returns
    -------
    None
    """
    # Initialize triples
    triples = []

    # Get gene name vertex
    gene_cls_vertex = {}
    if vertex_collections["gene_cls"].has(gene_name):
        gene_cls_vertex = vertex_collections["gene_cls"].get(gene_name)

    # Get diseases and drugs
    diseases = []
    drugs = []
    if gene_id in gdata["gene_ids"]:
        diseases = gdata["gene_ids"][gene_id]["resources"]["diseases"]
        drugs = gdata["gene_ids"][gene_id]["resources"]["drugs"]

    for disease in diseases:
        # TODO: Use keyword argument with suitable default, and set in command line arguments
        if disease["score"] < 0.5:
            continue

        # Insert disease vertex
        _key = disease["id"]
        if not vertex_collections["disease_cls"].has(_key):
            disease["_key"] = _key
            disease["label"] = disease["name"]
            vertex_collections["disease_cls"].insert(disease)
        disease_cls_vertex = vertex_collections["disease_cls"].get(_key)

        # Define and collect (gene_cls, IS_GENETIC_BASIS_FOR_CONDITION, disease_cls) triples
        triples.append(
            (
                gene_cls_vertex,
                {"label": "IS_GENETIC_BASIS_FOR_CONDITION"},
                disease_cls_vertex,
            )
        )

    for drug in drugs:

        # Insert drug product vertex
        _key = drug["id"]
        if not vertex_collections["drug_product_cls"].has(_key):
            drug["_key"] = _key
            drug["label"] = drug["name"]
            vertex_collections["drug_product_cls"].insert(drug)
        drug_product_cls_vertex = vertex_collections["drug_product_cls"].get(_key)

        # Define and collect (drug_product_cls, MOLECULARLY_INTERACTS_WITH, gene_cls) triples
        triples.append(
            (
                drug_product_cls_vertex,
                {"label": "MOLECULARLY_INTERACTS_WITH"},
                gene_cls_vertex,
            )
        )

        # Get diseases vertex
        _key = drug["disease_id"]
        if vertex_collections["disease_cls"].has(_key):
            disease_cls_vertex = vertex_collections["disease_cls"].get(_key)

            # Define and collect (drug_product_cls, IS_SUBSTANCE_THAT_TREATS, disease_cls) triples
            triples.append(
                (
                    drug_product_cls_vertex,
                    {"label": "IS_SUBSTANCE_THAT_TREATS"},
                    disease_cls_vertex,
                )
            )

    # Insert all triples
    insert_triples(triples, edge_collections)


def main():

    parser = argparse.ArgumentParser(description="Load Manually Curated Data")
    parser.add_argument(
        "--mdata-dirname",
        default=Path("../data/cell-kn"),
        help="name of directory containing manually curated data",
    )
    parser.add_argument(
        "--mdata-filename",
        default="HLCA_CellRef-ver-0.4.0.xlsm",
        help="name of file containing manually curated data",
    )
    parser.add_argument(
        "--label",
        default="",
        help="label to add to database_name",
    )
    parser.add_argument(
        "--include-bnodes",
        action="store_true",
        help="assume BNodes were included when loading",
    )
    group = parser.add_argument_group(
        "Cell Ontology (CL)", "Assume this version of the CL has been loaded"
    )
    exclusive_group = group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument(
        "--slim", action="store_true", help="assume the slim ontology has been loaded"
    )
    exclusive_group.add_argument(
        "--full", action="store_true", help="assume the full ontology has been loaded"
    )

    args = parser.parse_args()

    gdata_dirname = args.mdata_dirname
    gdata_filename = args.mdata_filename.replace(".xlsm", ".json")

    if args.slim:
        db_name = "BioPortal-Slim"
        graph_name = "CL-Slim"

    if args.full:
        db_name = "BioPortal-Full"
        graph_name = "CL-Full"

    if args.include_bnodes:
        db_name += "-BNodes"
        graph_name += "-BNodes"

    if args.label:
        db_name += f"-{args.label}"

    print(f"Loading {args.mdata_dirname / args.mdata_filename}")
    mdata = load_mdata(args.mdata_dirname, args.mdata_filename)

    print("Code implements HLCA-CellRef-Triples-ver-0.6.0")

    print(f"Loading {gdata_dirname / gdata_filename}")
    nm2ids = get_gene_name_to_ids_map()
    gdata = load_gdata(mdata, gdata_dirname, gdata_filename, nm2ids)

    print(f"Getting graph {graph_name} from {db_name}")
    adb_graph = get_graph(db_name, graph_name)

    print("Defining and creating vertex and edge collections")
    vertex_collections, edge_collections = init_collections(adb_graph)

    g2t = {}
    n_row = mdata.shape[0]
    for i_row, row in mdata.iterrows():
        print(f"Inserting mdata vertices and edges from row {i_row} (of {n_row})")
        insert_vertices_and_edges_from_mdata_row(
            row, vertex_collections, edge_collections, nm2ids, g2t
        )

    for source, s in gdata["source"].items():
        for cell_type, t in s["cell_types"].items():
            for gene_name, n in t["marker_names"].items():
                for gene_id in n["marker_ids"]:
                    print(
                        f"Inserting gdata vertices and edges for source {source} gene name {gene_name} and id {gene_id}"
                    )
                    insert_vertices_and_edges_from_gdata_for_gene_name_and_id(
                        gdata, gene_name, gene_id, vertex_collections, edge_collections
                    )


if __name__ == "__main__":
    main()
