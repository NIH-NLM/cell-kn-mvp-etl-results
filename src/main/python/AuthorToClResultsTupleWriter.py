import ast
import json
from pathlib import Path
from urllib.parse import urlparse

from rdflib.term import Literal, URIRef

from E_Utilities import get_data_for_pmid
from ExternalApiResultsFetcher import CELLXGENE_PATH
from LoaderUtilities import (
    DEPRECATED_TERMS,
    MIN_CLUSTER_SIZE,
    PURLBASE,
    RDFSBASE,
    collect_results_sources_data,
    load_results,
    hyphenate,
)

TUPLES_DIRPATH = Path("../../../data/tuples")


def create_tuples_from_author_to_cl(author_to_cl_results, cellxgene_results):
    """Creates tuples from manual author cell set to CL term mapping
    consistent with schema v0.7. Exclude clusters smaller than the
    minimum size.

    Parameters
    ----------
    author_to_cl_results : pd.DataFrame
        DataFrame containing author to CL results
    cellxgene_results : dict
        Dictionary containing cellxgene results

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    """
    tuples = []

    # Nodes for these results
    csd_term = f"CSD_{author_to_cl_results['dataset_id'][0]}"
    pub_term = f"PUB_{author_to_cl_results['DOI'][0].replace('/', '-')}"
    ds_term = f"DS_{author_to_cl_results['dataset_source'][0]}"

    # Cell_set_dataset_Ind, SOURCE, Publication_Ind
    # IAO:0000100, dc:source, IAO:0000311
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}/dc#Source"),
            URIRef(f"{PURLBASE}/{pub_term}"),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#Source"),
            Literal("Manual Mapping"),
        )
    )

    # Cell_set_dataset_Ind, INSTANCE_OF, Cell_set_dataset_Class
    # IAO:0000100, rdf:type, IAO:0000100
    tuples.append(
        (
            URIRef(f"{PURLBASE}/csd_term"),
            URIRef(f"{RDFSBASE}/rdf#type"),
            URIRef(f"{PURLBASE}/ds_term"),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/csd_term"),
            URIRef(f"{PURLBASE}/ds_term"),
            URIRef(f"{RDFSBASE}#Source"),
            Literal("Manual Mapping"),
        )
    )

    # CSD node annotations
    pmid_data = get_data_for_pmid(author_to_cl_results["PMID"][0])
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}#Citation"),
            Literal(pmid_data["Citation"]),
        )
    )
    keys = [
        "Link_to_publication",
        "Link_to_CELLxGENE_collection",
        "Link_to_CELLxGENE_dataset",
        "Dataset_name",
        "Number_of_cells",
        "Organism",
        "Tissue",
        "Disease_status",
        "Collection_ID",
        "Collection_version_ID",
        "Dataset_ID",
        "Dataset_version_ID",
        "Zenodo/Nextflow_workflow/Notebook",
    ]
    for key in keys:
        value = cellxgene_results[key]
        if isinstance(value, str):
            value = value.replace("https://", "")
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{csd_term}"),
                URIRef(f"{RDFSBASE}#{key.replace(' ', '_')}"),
                Literal(value),
            )
        )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}#Cell_type"),
            Literal(str(author_to_cl_results["author_category"][0])),
        )
    )

    # PUB node annotations
    for key in pmid_data.keys():
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{pub_term}"),
                URIRef(f"{RDFSBASE}#{key.capitalize().replace(' ', '_')}"),
                Literal(pmid_data[key]),
            )
        )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#PMID"),
            Literal(str(author_to_cl_results["PMID"][0])),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#PMCID"),
            Literal(str(author_to_cl_results["PMCID"][0])),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#DOI"),
            Literal(author_to_cl_results["DOI"][0]),
        )
    )

    # Nodes for each cell type or cell set
    for _, row in author_to_cl_results.iterrows():
        uuid = row["uuid"]
        cl_term = urlparse(row["cell_ontology_id"]).path.replace("/obo/", "")
        if cl_term in DEPRECATED_TERMS:
            print(f"Warning: CL term {cl_term} deprecated")
        uberon_term = urlparse(row["uberon_entity_id"]).path.replace("/obo/", "")
        if uberon_term in DEPRECATED_TERMS:
            print(f"Warning: UBERON term {uberon_term} deprecated")
        author_cell_set = hyphenate(row["author_cell_set"])
        cluster_size = row["clusterSize"]
        if cluster_size < MIN_CLUSTER_SIZE:
            continue
        cs_term = f"CS_{author_cell_set}-{uuid}"
        bmc_term = f"BMC_{uuid}"
        bgs_term = f"BGS_{uuid}"

        # Cell_type_Class, PART_OF, Anatomical_structure_Class
        # CL:0000000, BFO:0000050, UBERON:0001062
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{PURLBASE}/BFO_0000050"),
                URIRef(f"{PURLBASE}/{uberon_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{PURLBASE}/{uberon_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Cell_type_Class, HAS_EXEMPLAR_DATA, Cell_set_dataset_Ind
        # CL:0000000, RO:0015001, IAO:0000100
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{PURLBASE}/RO_0015001"),
                URIRef(f"{PURLBASE}/{csd_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{PURLBASE}/{csd_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Cell_set_Ind, DERIVES_FROM, Anatomical_structure_Cls
        # -, RO:0001000, UBERON:0001062
        # TODO: Add Anatomical_structure_Ind annotations, remove, or replace?
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/RO_0001000"),
                URIRef(f"{PURLBASE}/{uberon_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{uberon_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Cell_set_Ind, SOURCE, Cell_set_dataset_Ind
        # -, dc:source, IAO:0000100
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}/dc#Source"),
                URIRef(f"{PURLBASE}/{csd_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{csd_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Cell_set_Ind, COMPOSED_PRIMARILY_OF, Cell_type_Class
        # -, RO:0002473, CL:0000000
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/RO_0002473"),
                URIRef(f"{PURLBASE}/{cl_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Cell_set, RO:0002292 (EXPRESSES), Binary_gene_set
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/RO_0002292"),
                URIRef(f"{PURLBASE}/{bgs_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bgs_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Biomarker_combination_Ind, IS_CHARACTERIZING_MARKER_SET_FOR, Cell_type_Class
        # TODO: Update and use RO term
        # -, RO:0015004, CL:0000000
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{PURLBASE}/IS_CHARACTERIZING_MARKER_SET_FOR"),
                URIRef(f"{PURLBASE}/{cl_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # Node annotations
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Author_cell_term"),
                Literal(row["author_cell_term"]),
            )
        )
        keys = [
            "Link_to_publication",
            "Link_to_CELLxGENE_collection",
            "Link_to_CELLxGENE_dataset",
            "Dataset_name",
        ]
        for key in keys:
            value = cellxgene_results[key]
            if isinstance(value, str):
                value = value.replace("https://", "")
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{RDFSBASE}#{key.replace(' ', '_')}"),
                    Literal(value),
                )
            )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Cell_type"),
                Literal(cl_term),
            )
        )

        # Edge annotations
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{RDFSBASE}#Match"),
                Literal(row["match"]),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{RDFSBASE}#Mapping_method"),
                Literal(row["mapping_method"]),
            )
        )

        # Nodes for each cell type and marker gene
        marker_genes = ast.literal_eval(row["NSForest_markers"])
        for gene in marker_genes:
            gs_term = f"GS_{gene}"

            # Gene_Class, PART_OF, Cell_type_Class
            # SO:0000704, BFO:0000050, CL:0000000
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/BFO_0000050"),
                    URIRef(f"{PURLBASE}/{cl_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Manual Mapping"),
                )
            )

        # Nodes for each cell type, and marker and binary gene
        binary_genes = ast.literal_eval(row["binary_genes"])
        for gene in marker_genes + binary_genes:
            gs_term = f"GS_{gene}"

            # Cell_type_Class, SELECTIVELY EXPRESS, Gene_Class
            # TODO: Update and use RO term
            # CL:0000000, RO:0002292, SO:0000704
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{RDFSBASE}#SELECTIVELY_EXPRESS"),
                    URIRef(f"{PURLBASE}/{gs_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Manual Mapping"),
                )
            )

            # Gene_Class, PART_OF, Cell_type_Class
            # SO:0000704, BFO:0000050, CL:0000000
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/BFO_0000050"),
                    URIRef(f"{PURLBASE}/{cl_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Manual Mapping"),
                )
            )

            # Gene_Class, EXPRESSED_IN, Anatomical_structure_Class
            # SO:0000704, RO:0002206, UBERON:0001062
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/RO_0002206"),
                    URIRef(f"{PURLBASE}/{uberon_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{uberon_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Manual Mapping"),
                )
            )

    return tuples


def main(summarize=False):
    """Collect paths to all NSForest results, and author cell set to
    CL term mappings identified in the results sources, and
    dataset_version_id used for creating the NSForest results in order
    to create tuples consistent with schema v0.7, and write the result
    to a JSON file. If summarizing, retain the first row only, and
    include results in output.

    Parameters
    ----------
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    None
    """
    # Collect paths to all NSForest results, and author cell set to CL
    # term mappings identified in the results sources, and load the
    # CELLxGENE results.
    (
        nsforest_paths,
        author_to_cl_paths,
        _dataset_version_ids,
        _cl_terms,
        _gene_names,
        _gene_ensembl_ids,
        _gene_entrez_ids,
    ) = collect_results_sources_data()
    with open(CELLXGENE_PATH, "r") as fp:
        cellxgene_results_sets = json.load(fp)
    for author_to_cl_path, nsforest_path, cellxgene_results in zip(
        author_to_cl_paths, nsforest_paths, cellxgene_results_sets
    ):
        if author_to_cl_path == []:
            print(
                f"No author cell set to CL term map for NSForest results {nsforest_path}"
            )
            continue

        # Load author cell set to CL term mapping, dropping "uuid"
        # column in order to merge "uuid" column from NSForest results
        author_to_cl_results = (
            load_results(author_to_cl_path)
            .sort_values("author_cell_set", ignore_index=True)
            .drop(columns=["uuid"])
        )

        # Load NSForest results
        nsforest_results = load_results(nsforest_path).sort_values(
            "clusterName", ignore_index=True
        )
        if summarize:
            nsforest_results = nsforest_results.head(1)

        # Merge NSForest results with manual author cell set to CL
        # term mapping since author cell sets may not align exactly
        author_to_cl_results = author_to_cl_results.merge(
            nsforest_results[
                ["clusterName", "clusterSize", "NSForest_markers", "binary_genes", "uuid"]
            ].copy(),
            left_on="author_cell_set",
            right_on="clusterName",
        )

        print(f"Creating tuples from {author_to_cl_path}")
        author_to_cl_tuples = create_tuples_from_author_to_cl(
            author_to_cl_results, cellxgene_results
        )
        if summarize:
            output_dirpath = TUPLES_DIRPATH / "summaries"
        else:
            output_dirpath = TUPLES_DIRPATH
        with open(
            output_dirpath / author_to_cl_path.name.replace(".csv", ".json"),
            "w",
        ) as f:
            data = {}
            if summarize:
                data["results"] = author_to_cl_results.to_dict()
            data["tuples"] = author_to_cl_tuples
            json.dump(data, f, indent=4)

        if summarize:
            break


if __name__ == "__main__":
    main(summarize=True)
    main()
