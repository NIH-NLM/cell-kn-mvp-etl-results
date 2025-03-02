import json
from pathlib import Path
from pprint import pprint

from CellKnSchemaParser import read_schema

schema_path = Path("../data/cell-kn-schema-v0.7.0.xlsx")
schema, terms = read_schema(schema_path)

with open(Path("../data/cell-kn-mvp-annotation-results-2025-02-18.json"), "r") as f:
    results = json.load(f)

subject_types = set()
relations = set()
object_types = set()

keys = results[0].keys()

for triple in results:

    if results[0].keys() != keys:
        raise Exception("Not all triples have the same keys")

    subject_types.add(triple["subject_type"])
    relations.add(triple["relation"])
    object_types.add(triple["object_type"])


names = {}
identifiers = {}
for type in list(subject_types) + list(object_types):
    names[type] = set()
    identifiers[type] = set()

for triple in results:
    for node in ["subject", "object"]:
        type = triple[node + "_type"]
        name = triple[node + "_name"]
        identifier = triple[node + "_identifier"]

        names[type].add(name)
        identifiers[type].add(identifier)


with open(__file__.replace(".py", ".out"), "w") as f:
    f.write("\n=== Subjects and their CURIE\n\n")
    for subject_type in subject_types:
        f.write(f"{subject_type}, {terms['CURIE'][terms['Schema Name'] == subject_type].values}\n")

    f.write("\n=== Predicates and their CURIE\n\n")
    for relation in relations:
        f.write(f"{relation}, {terms['CURIE'][terms['Schema Name'] == relation].values}\n")

    f.write("\n=== Objects and their CURIE\n\n")
    for object_type in object_types:
        f.write(f"{object_type}, {terms['CURIE'][terms['Schema Name'] == object_type].values}\n")

    f.write("\n=== Types and their names\n\n")
    pprint(names, stream=f)

    f.write("\n=== Types and their identifiers\n\n")
    pprint(identifiers, stream=f)
