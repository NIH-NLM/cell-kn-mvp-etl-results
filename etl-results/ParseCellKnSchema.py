from pathlib import Path

import numpy as np
import pandas as pd


schema = pd.read_excel(Path("../data/cell-kn-schema-v0.5.0.xlsx"), 0)
relations = pd.read_excel(Path("../data/cell-kn-schema-v0.5.0.xlsx"), 2)

# Drop Subtype, child, parent, or pathway since not present in the relations
schema["Subject Node"] = schema["Subject Node"].str.replace(" (subtype/child)", "")
schema["Object Node"] = schema["Object Node"].str.replace(" (parent)", "")
schema["Object Node"] = schema["Object Node"].str.replace("/pathway", "")

# Drop class designation since redundant given the connections
schema["Subject Node"] = schema["Subject Node"].apply(lambda n: n.replace("_class", ""))
schema["Object Node"] = schema["Object Node"].apply(lambda n: n.replace("_class", ""))

# Drop cellular component since not in the KG
schema = schema[
    (schema["Subject Node"] != "Cellular_component")
    & (schema["Object Node"] != "Cellular_component")
]

# Drop unknown precicate relation
schema = schema[schema["Predicate Relation"] != "???"]

# Organism and species appear in the schema as seperate entries, but as combined in the relations
combined = relations[relations["Schema Name"] == "Organism/Species"]
organism = combined.copy()
species = combined.copy()
organism["Schema Name"] = organism["Schema Name"].str.replace("/Species", "")
species["Schema Name"] = species["Schema Name"].str.replace("Organism/", "")
relations[relations["Schema Name"] == "Organism/Species"] = organism
relations = pd.concat((relations, species))

# Ensure schema subject and object nodes, and predicate relations match relation schema names
missing_subjects = set(schema["Subject Node"]) - set(relations["Schema Name"])
if missing_subjects != set():
    raise Exception(f"Unexpected missing subjects: {missing_subjects}")
missing_objects = set(schema["Object Node"]) - set(relations["Schema Name"])
if missing_objects != set():
    raise Exception(f"Unexpected missing objects: {missing_objects}")
missing_predicates = set(schema["Predicate Relation"]) - set(relations["Schema Name"])
if missing_predicates != set():
    raise Exception(f"Unexpected missing predicates: {missing_predicates}")

# Add subject and object nodes with their type: class or individual
connections = schema["Connections"].str.split("-", expand=True)
schema["Subject Node Type"] = schema["Subject Node"] + "_" + connections[0]
schema["Object Node Type"] = schema["Object Node"] + "_" + connections[1]

# Add subject and object nodes, and predicate relations with their CURIE
schema["Subject Node Curie"] = schema["Subject Node"].apply(
    lambda n: relations["CURIE"][relations["Schema Name"] == n].iloc[0]
)
schema["Object Node Curie"] = schema["Object Node"].apply(
    lambda n: relations["CURIE"][relations["Schema Name"] == n].iloc[0]
)
schema["Predicate Relation Curie"] = schema["Predicate Relation"].apply(
    lambda n: relations["CURIE"][relations["Schema Name"] == n].iloc[0]
)

# Identify unique subjects, objects, and vertices
subjects = np.unique(schema["Subject Node Type"].values)
objects = np.unique(schema["Object Node Type"].values)
vertices = np.unique(np.concatenate((subjects, objects)))

# Identify triples which contain selected vertices
triples_with_names = schema.loc[
    schema["Subject Node"].isin(["Biomarker_combination", "Cell_set", "Gene"])
    | schema["Object Node"].isin(["Biomarker_combination", "Cell_set", "Gene"]),
    ["Subject Node Type", "Predicate Relation", "Object Node Type"],
]
triples_with_curies = schema.loc[
    schema["Subject Node"].isin(["Biomarker_combination", "Cell_set", "Gene"])
    | schema["Object Node"].isin(["Biomarker_combination", "Cell_set", "Gene"]),
    ["Subject Node Curie", "Predicate Relation Curie", "Object Node Curie"],
]

# Write the result to an Excel spreadsheet
with pd.ExcelWriter(Path("../data/cell-kn-schema-v0.5.0-nsforest.xlsx")) as writer:
    pd.DataFrame(subjects, columns=["Subjects"]).to_excel(writer, sheet_name="Subjects")
    pd.DataFrame(objects, columns=["Objects"]).to_excel(writer, sheet_name="Objects")
    pd.DataFrame(vertices, columns=["Vertices"]).to_excel(writer, sheet_name="Vertices")
    triples_with_names.to_excel(writer, sheet_name="Triples with Names")
    triples_with_curies.to_excel(writer, sheet_name="Triples with CURIEs")
