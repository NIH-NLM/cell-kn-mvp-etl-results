import json
from pathlib import Path
from pprint import pprint

from rdflib.term import URIRef, Literal

from CellKnSchemaUtilities import read_schema
from LoaderUtilities import (
    get_mesh_to_mondo_map,
    hyphenate,
    map_mesh_to_mondo,
    PURLBASE,
    RDFSBASE,
)

TUPLES_DIRPATH = Path("../../../data/tuples")


def write_triple_components(annotation_results, terms, components_path):
    """Determine unique subject, predicate, and object types, and
    unique subject and object names and identifiers. Use the schema
    mapping of terms to CURIEs.

    Parameters
    ----------
    annotation_results : dict
        Annotation results
    terms : pd.DataFrame
        DataFrame containing schema mapping from term to CURIE
    component_paths : Path
        Path to which to write components

    Returns
    -------
    None
    """
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
    with open(components_path, "w") as f:

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


def normalize_term(annotation, term, mesh2mondo):
    """Normalize annoation by replacing colons with underscores,
    creating unique identifiers from annotation components, using
    common vertex names, mapping MESH to MONDO terms, mapping PMID to
    DOIs, and generally cleaning up.

    Parameters
    ----------
    annotation : dict
        A single annotation
    terms : str
        The term of the annoation: "subject" or "predicate"
    mesh2mondo : dict
        Dictionary mapping MeSH term to MONDO term

    Returns
    -------
    str
        The normalized term
    """
    # Replace unicode characters
    if annotation[f"{term}_name"] is not None:
        annotation[f"{term}_name"] = annotation[f"{term}_name"].replace(
            "\u03b3\u03b4", "gamma-delta"
        )
        annotation[f"{term}_name"] = annotation[f"{term}_name"].replace("\u2212", "-")
    if annotation[f"{term}_identifier"] is not None:
        annotation[f"{term}_identifier"] = annotation[f"{term}_identifier"].replace(
            "\u2212", "-"
        )

    # Normalize by annotation type
    atype = annotation[f"{term}_type"]
    if atype == "Anatomical_structure":
        return annotation[f"{term}_identifier"].replace(":", "_")

    elif atype == "Assay":
        return annotation[f"{term}_identifier"].replace(":", "_")

    elif atype == "Biomarker_combination":  # Never a subject
        return (
            "BMC_"
            + annotation[f"{term}_name"]
            + "-"
            + annotation[f"subject_identifier"].split("-")[0]
        )

    elif atype == "Cell_set":
        return (
            "CS_"
            + hyphenate(annotation[f"{term}_name"])
            + "-"
            + annotation[f"subject_identifier"].split("-")[0]
        )

    elif atype == "Cell_set_dataset":
        return annotation[f"{term}_identifier"].replace("NLP_dataset", "CSD")

    elif atype == "Cell_type":
        return (
            annotation[f"{term}_identifier"]
            .replace("<skos:related>", "")
            .replace(":", "_")
        )

    elif atype == "Disease":
        return map_mesh_to_mondo(annotation[f"{term}_identifier"], mesh2mondo)

    elif atype == "Gene":
        gene_symbol = annotation[f"{term}_name"].replace("Myelin basic protein", "MBP")
        if gene_symbol == gene_symbol.upper():
            return "GS_" + gene_symbol
        else:
            return None

    elif atype == "Publication":
        if annotation[f"{term}_identifier"] == "37824655":
            # Jorstad: 0,37824655,PMC11687949,doi.org/10.1126/science.adf6812
            return "PUB_doi.org/10.1126/science.adf6812"

        elif annotation[f"{term}_identifier"] == "37516747":
            # Guo : 0,37516747,PMC10387117,doi.org/10.1038/s41467-023-40173-5
            return "PUB_doi.org/10.1038/s41467-023-40173-5"

        elif annotation[f"{term}_identifier"] == "37291214":
            # Sikkema: 0,37291214,PMC10287567,doi.org/10.1038/s41591-023-02327-2
            return "PUB_doi.org/10.1038/s41591-023-02327-2"

        elif annotation[f"{term}_identifier"] == "38014002":  # Never expected
            # Li: 0,38014002,PMC10680922,doi.org/10.1101/2023.11.07.566105
            return "PUB_doi.org/10.1101/2023.11.07.566105"

        else:
            return None


def create_tuples_from_annotation(annotation_results):
    """Create tuples from manual annotation of selected articles
    consistent with schema v0.7.

    Parameters
    ----------
    annotation_results : dict
        Dictionary of manual annotation results in semantic triple
        form

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples only) created
    """
    mesh2mondo = get_mesh_to_mondo_map(
        "../../../cell-kn-mvp-etl-ontologies/data/obo", "mondo-simple.owl"
    )
    tuples = []
    for annotation in annotation_results:
        subject = normalize_term(annotation, "subject", mesh2mondo)
        predicate = annotation["relation"]
        object = normalize_term(annotation, "object", mesh2mondo)
        if subject is None or object is None:
            continue

        s = URIRef(f"{PURLBASE}/{subject}")
        p = URIRef(f"{RDFSBASE}#{predicate}")
        o = URIRef(f"{PURLBASE}/{object}")
        tuples.append((s, p, o))

        p = URIRef(f"{RDFSBASE}#source")
        l = Literal("Manual Annotation")
        tuples.append((s, o, p, l))

    return tuples


def main(summarize=False):
    """Load schema, and annotation results to determine unique
    annotation subject, predicate, and object types, and unique
    annotation subject object names and identifiers, using the schema
    mapping of terms to CURIEs. Create tuples consistent with schema
    v0.7, and write the result to a JSON file. If summarizing, retain
    the first row only, and include results in output.

    Parameters
    ----------
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    None
    """
    # Load schema, and mapping of term to CURIES
    schema_path = Path("../../../data/schema/cell-kn-schema-v0.7.0.xlsx")
    _schema, terms = read_schema(schema_path)

    # Load annotation results
    annotation_path = Path(
        "../../../data/results/cell-kn-mvp-annotation-results-2025-03-14.json"
    )
    with open(annotation_path, "r") as fp:
        annotation_results = json.load(fp)

    print(f"Determining components from {annotation_path}")
    components_path = Path(str(annotation_path).replace(".json", ".out"))
    write_triple_components(annotation_results, terms, components_path)

    print(f"Creating tuples from {annotation_path}")
    if summarize:
        annotation_results = [annotation_results[0]]
    annotation_tuples = create_tuples_from_annotation(annotation_results)
    if summarize:
        output_dirpath = TUPLES_DIRPATH / "summaries"
    else:
        output_dirpath = TUPLES_DIRPATH
    data = {}
    if summarize:
        data["results"] = annotation_results
    data["tuples"] = annotation_tuples
    with open(
        output_dirpath / annotation_path.name.replace(".csv", ".json"),
        "w",
    ) as f:
        json.dump(data, f, indent=4)


if __name__ == "__main__":
    main(summarize=True)
    main()
