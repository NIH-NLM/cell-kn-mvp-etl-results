import json
from pathlib import Path
from pprint import pprint

import numpy as np
import pandas as pd

from LoaderUtilities import load_results, hyphenate, PURLBASE, RDFSBASE


def read_schema(schema_path):

    schema = pd.read_excel(schema_path, 0)
    terms = pd.read_excel(schema_path, 2)

    # Drop Subtype, child, parent, or pathway since not present in the classes/relations
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

    # Drop unknown precicate relation
    schema = schema[schema["Predicate Relation"] != "???"]

    # Drop predicates since not present in the classes/relations
    # schema = schema[schema["Predicate Relation"] != "IS_MARKER_FOR_PROTEIN"]
    # schema = schema[schema["Predicate Relation"] != "IS_A"]

    # Organism and species appear in the schema as seperate entries, but as combined in the classes/relations
    combined = terms[terms["Schema Name"] == "Organism/Species"]
    organism = combined.copy()
    species = combined.copy()
    organism["Schema Name"] = organism["Schema Name"].str.replace("/Species", "")
    species["Schema Name"] = species["Schema Name"].str.replace("Organism/", "")
    terms[terms["Schema Name"] == "Organism/Species"] = organism
    terms = pd.concat((terms, species))

    # Ensure schema subject and object nodes, and predicate terms match classes/relations schema names
    missing_subjects = set(schema["Subject Node"]) - set(terms["Schema Name"])
    if missing_subjects != set():
        raise Exception(f"Unexpected missing subjects: {missing_subjects}")
    missing_objects = set(schema["Object Node"]) - set(terms["Schema Name"])
    if missing_objects != set():
        raise Exception(f"Unexpected missing objects: {missing_objects}")
    missing_predicates = set(schema["Predicate Relation"]) - set(terms["Schema Name"])
    if missing_predicates != set():
        raise Exception(f"Unexpected missing predicates: {missing_predicates}")

    # Add subject and object nodes with their type: class or individual
    connections = schema["Connections"].str.split("-", expand=True)
    schema["Subject Node Type"] = schema["Subject Node"] + "_" + connections[0]
    schema["Object Node Type"] = schema["Object Node"] + "_" + connections[1]

    # Add subject and object nodes, and predicate terms with their CURIE
    schema["Subject Node Curie"] = schema["Subject Node"].apply(
        lambda n: terms["CURIE"][terms["Schema Name"] == n].iloc[0]
    )
    schema["Object Node Curie"] = schema["Object Node"].apply(
        lambda n: terms["CURIE"][terms["Schema Name"] == n].iloc[0]
    )
    schema["Predicate Relation Curie"] = schema["Predicate Relation"].apply(
        lambda n: terms["CURIE"][terms["Schema Name"] == n].iloc[0]
    )

    return schema, terms


def create_tuples(schema):

    tuples = []

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

    # Identify unique subjects, objects, and vertices
    subjects = np.unique(schema["Subject Node Type"].values)
    objects = np.unique(schema["Object Node Type"].values)
    vertices = np.unique(np.concatenate((subjects, objects)))

    return subjects, objects, vertices


def identify_nsforest_triples(schema, subjects, objects, vertices, triples_path):
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


def main():
    schema_path = Path("../../../data/schema/cell-kn-schema-v0.7.0.xlsx")
    schema, terms = read_schema(schema_path)

    schema_tuples_path = Path(str(schema_path.resolve()).replace(".xlsx", ".json"))
    schema_tuples = create_tuples(schema)

    with open(schema_tuples_path, "w") as f:
        results = {}
        results["tuples"] = schema_tuples
        json.dump(results, f, indent=4)

    subjects, objects, vertices = identify_unique_classes(schema)

    nsforest_triples_path = Path(str(schema_path.resolve()).replace(".xlsx", "-nsforest.xlsx"))
    identify_nsforest_triples(schema, subjects, objects, vertices, nsforest_triples_path)

    author_to_cl_triples_path = Path(
        str(schema_path.resolve()).replace(".xlsx", "-author-to-cl.xlsx")
    )
    identify_author_to_cl_triples(schema, subjects, objects, vertices, author_to_cl_triples_path)


if __name__ == "__main__":
    main()
