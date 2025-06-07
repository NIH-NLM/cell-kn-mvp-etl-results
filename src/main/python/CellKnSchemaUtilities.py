import json
from pathlib import Path

import numpy as np
import pandas as pd

import ArangoDbUtilities as adb
from LoaderUtilities import PURLBASE


def read_schema(schema_path):
    """Reads the schema, then resolves a number of computationally
    incompatible, hand entered, contingent anachronisms.

    Parameters
    ----------
    schema_path : Path
        Path to schema Excel file

    Returns
    -------
    schema : pd.DataFrame
        DataFrame containing the schema triples
    terms : pd.DataFrame
        DataFrame containing the mapping from terms to CURIEs
    """
    # Read schema, and mapping from term to CURIE
    schema = pd.read_excel(schema_path, 0)
    terms = pd.read_excel(schema_path, 2)

    # Drop subtype, child, parent, or pathway since not present in the
    # classes/relations
    schema["Subject Node"] = schema["Subject Node"].str.replace(" (subtype/child)", "")
    schema["Object Node"] = schema["Object Node"].str.replace(" (parent)", "")
    schema["Object Node"] = schema["Object Node"].str.replace("/pathway", "")

    # Drop class designation since redundant given the connections
    schema["Subject Node"] = schema["Subject Node"].apply(
        lambda n: n.replace("_class", "")
    )
    schema["Object Node"] = schema["Object Node"].apply(
        lambda n: n.replace("_class", "")
    )

    # Drop cellular component since not in the KG
    schema = schema[
        (schema["Subject Node"] != "Cellular_component")
        & (schema["Object Node"] != "Cellular_component")
    ]

    # Drop unknown predicate relation
    # schema = schema[schema["Predicate Relation"] != "???"]

    # TODO: Remove
    # Drop predicates since not present in the classes/relations
    # schema = schema[schema["Predicate Relation"] != "IS_MARKER_FOR_PROTEIN"]
    # schema = schema[schema["Predicate Relation"] != "IS_A"]

    # Organism and species appear in the schema as separate entries,
    # but as combined in the classes/relations
    combined = terms[terms["Schema Name"] == "Organism/Species"]
    organism = combined.copy()
    species = combined.copy()
    organism["Schema Name"] = organism["Schema Name"].str.replace("/Species", "")
    species["Schema Name"] = species["Schema Name"].str.replace("Organism/", "")
    terms[terms["Schema Name"] == "Organism/Species"] = organism
    terms = pd.concat((terms, species))

    # Replace contingent predicate with placeholder
    schema["Predicate Relation"] = schema["Predicate Relation"].str.replace(
        "??? Need looser relationship to express that the two are merely associated",
        "ASSOCIATED_WITH",
    )

    # Ensure schema subject and object nodes, and predicate terms
    # match classes/relations schema names
    missing_subjects = set(schema["Subject Node"]) - set(terms["Schema Name"])
    if missing_subjects != set():
        print(f"Unexpected missing subjects: {missing_subjects}")
    missing_objects = set(schema["Object Node"]) - set(terms["Schema Name"])
    if missing_objects != set():
        print(f"Unexpected missing objects: {missing_objects}")
    missing_predicates = set(schema["Predicate Relation"]) - set(terms["Schema Name"])
    if missing_predicates != set():
        print(f"Unexpected missing predicates: {missing_predicates}")

    # Add subject and object nodes with their type: class or
    # individual
    connections = schema["Connections"].str.split("-", expand=True)
    schema["Subject Node Type"] = schema["Subject Node"] + "_" + connections[0]
    schema["Object Node Type"] = schema["Object Node"] + "_" + connections[1]

    # Add subject and object nodes, and predicate terms with their
    # CURIE
    schema["Subject Node Curie"] = schema["Subject Node"].apply(
        lambda n: (
            terms["CURIE"][terms["Schema Name"] == n].iloc[0]
            if (terms["Schema Name"] == n).any()
            else "NA"
        )
    )
    schema["Object Node Curie"] = schema["Object Node"].apply(
        lambda n: (
            terms["CURIE"][terms["Schema Name"] == n].iloc[0]
            if (terms["Schema Name"] == n).any()
            else "NA"
        )
    )
    schema["Predicate Relation Curie"] = schema["Predicate Relation"].apply(
        lambda n: (
            terms["CURIE"][terms["Schema Name"] == n].iloc[0]
            if (terms["Schema Name"] == n).any()
            else "NA"
        )
    )

    return schema, terms


def create_tuples(schema):
    """Create tuples suitable for loading the schema into ArangoDB.

    Parameters
    ----------
    schema : pd.DataFrame
        DataFrame containing the schema triples

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    """
    tuples = []

    # Resolve even more computationally incompatible, hand entered,
    # contingent anachronisms
    df = schema[["Subject Node Curie", "Predicate Relation Curie", "Object Node Curie"]]
    df = df.map(lambda e: e.replace("MONDO:0000001 or MONDO:0021178", "MONDO:0000001"))
    df = df.map(
        lambda e: e.replace(
            "PATO:0000068, MONDO:0000001 (disease), or MOND...", "PATO:0000068"
        )
    )
    df = df.map(
        lambda e: e.replace("HsapDv:0000000 or MmusDv:0000000", "HsapDv:0000000")
    )
    df = df.map(lambda e: e.replace("EFO:0002772 or EFO:0010183", "EFO:0002772"))
    df = df.map(
        lambda e: e.replace(
            "PATO:0000068, MONDO:0000001 (disease), or MONDO:0021178 (injury)",
            "PATO:0000068",
        )
    )
    df = df.map(lambda e: e.replace(":", "_"))

    # Create the tuples
    for _, row in df.iterrows():
        s, o, p = row
        tuples.append(
            (
                f"{PURLBASE}/{s}",
                f"{PURLBASE}/{o}",
                f"{PURLBASE}/{p}",
            )
        )

    return tuples


def identify_unique_classes(schema):
    """Identify unique subject, object, and vertex classes. Vertex
    classes equal the union of subject and object classes.

    Parameters
    ----------
    schema : pd.DataFrame
        DataFrame containing the schema triples

    Returns
    -------
    subjects : np.ndarray
        Numpy array containing unique subject classes
    objects : np.ndarray
        Numpy array containing unique object classes
    vertices : np.ndarray
        Numpy array containing unique vertex classes
    """
    subjects = np.unique(schema["Subject Node Type"].values)
    objects = np.unique(schema["Object Node Type"].values)
    vertices = np.unique(np.concatenate((subjects, objects)))
    return subjects, objects, vertices


def identify_nsforest_triples(schema, subjects, objects, vertices, triples_path):
    """Identify triples which contain vertices corresponding to
    NSForest results, then write the result to an Excel spreadsheet.

    Parameters
    ----------
    schema : pd.DataFrame
        DataFrame containing the schema triples
    subjects : np.ndarray
        Numpy array containing unique subject classes
    objects : np.ndarray
        Numpy array containing unique object classes
    vertices : np.ndarray
        Numpy array containing unique vertex classes
    triples_path : Path
        Path to which to write triples

    Returns
    -------
    None
    """
    # Assign vertices corresponding to NSForest results
    selected_vertices = [
        "Biomarker_combination",
        "Binary_gene_combination",
        "Cell_set",
        "Gene",
    ]

    # Identify triples which contain selected vertices
    triples_with_names = schema.loc[
        schema["Subject Node"].isin(selected_vertices)
        | schema["Object Node"].isin(selected_vertices),
        ["Subject Node Type", "Predicate Relation", "Object Node Type"],
    ]
    triples_with_curies = schema.loc[
        schema["Subject Node"].isin(selected_vertices)
        | schema["Object Node"].isin(selected_vertices),
        ["Subject Node Curie", "Predicate Relation Curie", "Object Node Curie"],
    ]

    # Write the result to an Excel spreadsheet
    with pd.ExcelWriter(triples_path) as writer:
        pd.DataFrame(subjects, columns=["Subjects"]).to_excel(
            writer, sheet_name="Subjects"
        )
        pd.DataFrame(objects, columns=["Objects"]).to_excel(
            writer, sheet_name="Objects"
        )
        pd.DataFrame(vertices, columns=["Vertices"]).to_excel(
            writer, sheet_name="Vertices"
        )
        triples_with_names.to_excel(writer, sheet_name="Triples with Names")
        triples_with_curies.to_excel(writer, sheet_name="Triples with CURIEs")


def identify_author_to_cl_triples(schema, subjects, objects, vertices, triples_path):
    """Identify triples which contain vertices corresponding to the
    manual mapping of author cell set to CL term results, then write
    the result to an Excel spreadsheet.

    Parameters
    ----------
    schema : pd.DataFrame
        DataFrame containing the schema triples
    subjects : np.ndarray
        Numpy array containing unique subject classes
    objects : np.ndarray
        Numpy array containing unique object classes
    vertices : np.ndarray
        Numpy array containing unique vertex classes
    triples_path : Path
        Path to which to write triples

    Returns
    -------
    None
    """
    # Assign vertices corresponding to the manual mapping of author
    # cell set to CL term results
    selected_vertices = [
        "Anatomical_structure",
        "Cell_set",
        "Cell_set_dataset",
        "Cell_type",
        "Publication",
    ]

    # Identify triples which contain selected vertices
    triples_with_names = schema.loc[
        schema["Subject Node"].isin(selected_vertices)
        | schema["Object Node"].isin(selected_vertices),
        ["Subject Node Type", "Predicate Relation", "Object Node Type"],
    ]
    triples_with_curies = schema.loc[
        schema["Subject Node"].isin(selected_vertices)
        | schema["Object Node"].isin(selected_vertices),
        ["Subject Node Curie", "Predicate Relation Curie", "Object Node Curie"],
    ]

    # Write the result to an Excel spreadsheet
    with pd.ExcelWriter(triples_path) as writer:
        pd.DataFrame(subjects, columns=["Subjects"]).to_excel(
            writer, sheet_name="Subjects"
        )
        pd.DataFrame(objects, columns=["Objects"]).to_excel(
            writer, sheet_name="Objects"
        )
        pd.DataFrame(vertices, columns=["Vertices"]).to_excel(
            writer, sheet_name="Vertices"
        )
        triples_with_names.to_excel(writer, sheet_name="Triples with Names")
        triples_with_curies.to_excel(writer, sheet_name="Triples with CURIEs")


def load_graph(graph, schema):
    """Create vertices from unique schema subject and object nodes,
    then create edges from all schema triples.

    Parameters
    ----------
    adb_graph : arango.graph.Graph
        An ArangoDB graph instance
    schema : pd.DataFrame
        DataFrame containing the schema triples

    Returns
    -------
    None
    """
    # Create vertices from unique schema subject and object nodes
    for node in pd.concat([schema["Subject Node"], schema["Object Node"]]).unique():
        vertex_collection = adb.create_or_get_vertex_collection(graph, node)
        vertex = {"_key": node, "label": node.replace("_", " ").title()}
        vertex_collection.insert(vertex)

    # Create edges from all schema triples
    for _, row in schema.iterrows():
        s, p, o = row[["Subject Node", "Predicate Relation", "Object Node"]]
        edge_collection, _ = adb.create_or_get_edge_collection(graph, s, o)
        edge = {
            "_key": f"{s}-{p}-{o}",
            "_from": f"{s}/{s}",
            "_to": f"{o}/{o}",
            "label": p,
        }
        edge_collection.insert(edge)


def main():
    """Read the schema, resolve a number of computationally
    incompatible, hand entered, contingent anachronisms, then create
    and write tuples suitable for loading the schema into ArangoDB.


    Parameters
    ----------

    Returns
    -------
    """
    # Read the schema, and create the tuples
    schema_path = Path("../../../data/schema/cell-kn-schema-2024-04-16.xlsx")
    schema, _terms = read_schema(schema_path)
    schema_tuples = create_tuples(schema)

    # Write the tuples
    print(f"Creating tuples from {schema_path}")
    schema_tuples_path = Path(str(schema_path.resolve()).replace(".xlsx", ".json"))
    with open(schema_tuples_path, "w") as f:
        results = {}
        results["tuples"] = schema_tuples
        json.dump(results, f, indent=4)

    # Identify unique subject, object, and vertex classes
    subjects, objects, vertices = identify_unique_classes(schema)

    # Identify and write triples which contain vertices corresponding
    # to NSForest results
    nsforest_triples_path = Path(
        str(schema_path.resolve()).replace(".xlsx", "-nsforest.xlsx")
    )
    identify_nsforest_triples(
        schema, subjects, objects, vertices, nsforest_triples_path
    )

    # Identify and write triples which contain vertices corresponding
    # to the manual mapping of author cell set to CL term results
    author_to_cl_triples_path = Path(
        str(schema_path.resolve()).replace(".xlsx", "-author-to-cl.xlsx")
    )
    identify_author_to_cl_triples(
        schema, subjects, objects, vertices, author_to_cl_triples_path
    )

    # Load the schema graph
    db_name = "Cell-KN-Schema"
    graph_name = "KN-Schema-2024-04-16"
    print(f"Loading {db_name}/{graph_name}")
    adb.delete_database(db_name)
    db = adb.create_or_get_database(db_name)
    graph = adb.create_or_get_graph(db, graph_name)
    load_graph(graph, schema)

    # TODO: Relocate?
    # Create analyzers and views
    for database_name in ["Cell-KN-Ontologies", "Cell-KN-Phenotypes"]:

        # First delete
        try:
            adb.delete_view(database_name)
        except Exception as ex:
            print(f"No view in {database_name} to delete")
        try:
            adb.delete_analyzers(database_name)
        except Exception as ex:
            print(f"No analyzers in {database_name} to delete")

        # Then create
        adb.create_analyzers(database_name)
        adb.create_view(
            database_name,
            collection_maps_name="../../../cell-kn-mvp-etl-ontologies/data/cell-kn-mvp-collection-maps-2025-04-23.json",
        )


if __name__ == "__main__":
    main()
