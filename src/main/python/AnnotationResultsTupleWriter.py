import json
from pathlib import Path
from pprint import pprint

from CellKnSchemaUtilities import read_schema


def write_triple_components(annotation_results_path, terms):
    """Load annotation results and determine unique subject,
    predicate, and object types, and unique subject and object names
    and identifiers. Use the schema mapping of terms to CURIEs.

    Parameters
    ----------
    annotation_results_path : Path
        Path to annotation results
    terms : pd.DataFrame
        DataFrame containing schema mapping from term to CURIE

    Returns
    -------
    None
    """
    # Load annotation results
    with open(annotation_results_path, "r") as f:
        annotation_results = json.load(f)

    # Identify unique subject, predicate, and object types
    subject_types = set()
    relations = set()
    object_types = set()
    keys = annotation_results[0].keys()
    for triple in annotation_results:
        subject_types.add(triple["subject_type"])
        relations.add(triple["relation"])
        object_types.add(triple["object_type"])

        # Also ensure all triples have the same keys
        if annotation_results[0].keys() != keys:
            raise Exception("Not all triples have the same keys")

    # Identify unique subject, and object names and identifiers
    names = {}
    identifiers = {}
    for node_type in list(subject_types) + list(object_types):
        names[node_type] = set()
        identifiers[node_type] = set()
    for triple in annotation_results:
        for node in ["subject", "object"]:
            node_type = triple[node + "_type"]
            name = triple[node + "_name"]
            identifier = triple[node + "_identifier"]
            names[node_type].add(name)
            identifiers[node_type].add(identifier)

    # Write triple component descriptions
    annotation_components_path = Path(
        str(annotation_results_path).replace(".json", ".out")
    )
    with open(annotation_components_path, "w") as f:

        f.write("\n=== Subjects and their CURIE\n\n")
        for subject_type in subject_types:
            f.write(
                f"{subject_type}, {terms['CURIE'][terms['Schema Name'] == subject_type].values}\n"
            )

        f.write("\n=== Predicates and their CURIE\n\n")
        for relation in relations:
            f.write(
                f"{relation}, {terms['CURIE'][terms['Schema Name'] == relation].values}\n"
            )

        f.write("\n=== Objects and their CURIE\n\n")
        for object_type in object_types:
            f.write(
                f"{object_type}, {terms['CURIE'][terms['Schema Name'] == object_type].values}\n"
            )

        f.write("\n=== Types and their names\n\n")
        pprint(names, stream=f)

        f.write("\n=== Types and their identifiers\n\n")
        pprint(identifiers, stream=f)


def main():
    """Load schema, and annotation results to determine unique
    annotation subject, predicate, and object types, and unique
    annotation subject object names and identifiers, using the schema
    mapping of terms to CURIEs.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Load the schema, and mapping of term to CURIES
    schema_path = Path("../../../data/schema/cell-kn-schema-v0.7.0.xlsx")
    _schema, terms = read_schema(schema_path)

    # Load the annotation results, and write triple component descriptions
    annotation_results_path = Path(
        "../../../data/results/cell-kn-mvp-annotation-results-2025-03-14.json"
    )
    write_triple_components(annotation_results_path, terms)


if __name__ == "__main__":
    main()
