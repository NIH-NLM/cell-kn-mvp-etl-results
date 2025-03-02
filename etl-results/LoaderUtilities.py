import random
import string

import pandas as pd


ALPHABET = string.ascii_lowercase + string.digits
PURLBASE = "http://purl.obolibrary.org/obo"
RDFSBASE = "http://www.w3.org/1999/02/22-rdf-syntax-ns"


def get_uuid():
    """Get an eight character random string.

    Parameters
    ----------
    None

    Returns
    -------
    An eight character random string.
    """
    return "".join(random.choices(ALPHABET, k=12))


def load_results(results_path):
    """Load NSForest results CSV file and append a UUID.

    Parameters
    ----------
    results_Path : Path
        Path of CSV file containing NSForest results

    Returns
    -------
    results : pd.DataFrame
        DataFrame containing NSForest results
    """
    results = pd.read_csv(results_path)
    if not "uuid" in results.columns:
        print(f"Add UUID column to results CSV file {results_path.name}")
        results["uuid"] = [get_uuid() for idx in results.index]
        results.to_csv(results_path)
    return results


def hyphenate(name):
    for c in [" ", "_", ","]:
        name = name.replace(c, "-").replace("--", "-")
    return name
