import json
from pathlib import Path
from pprint import pprint

import pandas as pd

from ExternalApiResultsFetcher import RESOURCES

SUMMARIES_DIRPATH = Path("../../../data/tuples/summaries")


with open(SUMMARIES_DIRPATH / "ExternalApiResultsLoader-guo.json", "r") as fp:
    data = json.load(fp)

print()
print("opentargets")
for resource in RESOURCES:
    df = pd.DataFrame.from_dict(data["opentargets"]["ENSG00000108821"][resource])
    print()
    print(resource)
    if df.shape[1] > 5:
        print(df.loc[:, list(df.columns[1:6])])
        print(df.loc[:, list(df.columns[6:])])
    else:
        print(df)

print()
print("uniprot")
pprint(data["uniprot"]["ENSP00000403925"])

def to_curie(e):

    if e is None:
        e = ""

    elif "http://purl.obolibrary.org/" in e:
        e = e.replace("http://purl.obolibrary.org/", "").replace("_", ":", 1)

    elif "http://www.w3.org/1999/02/22-rdf-syntax-ns" in e:
        e = e.replace("http://www.w3.org/1999/02/22-rdf-syntax-ns", "rdf")

    return e


df = pd.DataFrame.from_dict(data["tuples"]).map(to_curie)

print()
print("tuples")
print(df)
