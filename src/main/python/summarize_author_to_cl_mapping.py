import json
from pathlib import Path

import pandas as pd

SUMMARIES_DIRPATH = Path("../../../data/tuples/summaries")


with open(SUMMARIES_DIRPATH / "AuthorToClResultsLoader-guo.json", "r") as fp:
    data = json.load(fp)

df = pd.DataFrame.from_dict(data["results"])

# For Author to CL mapping
print(df.loc[:, list(df.columns[1:4])])
print(df.loc[:, list(df.columns[4:9])])
print(df.loc[:, list(df.columns[9:13])])
print(df.loc[:, list(df.columns[13:18])])
print(df.loc[:, list(df.columns[18:])])

def to_curie(e):

    if e is None:
        e = ""

    elif "http://purl.obolibrary.org/" in e:
        e = e.replace("http://purl.obolibrary.org/", "").replace("_", ":", 1)
        
    elif "http://www.w3.org/1999/02/22-rdf-syntax-ns" in e:
        e = e.replace("http://www.w3.org/1999/02/22-rdf-syntax-ns", "rdf")

    return e

df = pd.DataFrame.from_dict(data["tuples"]).map(to_curie)

print(df)
