import ast
import json
from pathlib import Path

import gget
import pandas as pd
import requests
import scanpy as sc

from LoaderUtilities import load_results
from UniprotIdMapper import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)


RESOURCES = [
    "diseases",
    "drugs",
    "interactions",
    "pharmacogenetics",
    "tractability",
    "expression",
    "depmap",
]

NSFOREST_DIRPATH = Path("../data/results")


def get_gene_name_to_ids_map():
    """Query BioMart to map gene names to ids.

    Parameters
    ----------
    None

    Returns
    -------
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene id
    """
    print("Creating gene name to ids map")
    gnm2ids = sc.queries.biomart_annotations(
        "hsapiens", ["external_gene_name", "ensembl_gene_id"], use_cache=True
    ).set_index("external_gene_name")

    return gnm2ids


def map_gene_name_to_ids(name, gnm2ids):
    """Map a gene name to a gene id list.

    Parameters
    ----------
    name : str
        Gene name
    gnm2ids : pd.DataFrame
        DataFrame indexed by name containing id

    Returns
    -------
    list
        Gene ids
    """
    if name in gnm2ids.index:
        ids = gnm2ids.loc[name, "ensembl_gene_id"]
        if isinstance(ids, pd.core.series.Series):
            ids = ids.to_list()
        else:
            ids = [ids]
        print(f"Mapped gene name {name} to ids {ids}")
    else:
        print(f"Could not find ids for name: {name}")
        ids = []
    return ids


def map_protein_id_to_accession(protein_id):

    job_id = submit_id_mapping(
        from_db="Ensembl_Protein", to_db="UniProtKB", ids=[protein_id]
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        data = get_id_mapping_results_search(link)

    if len(data["results"]) == 0:
        accession = None

    else:
        accession = data["results"][0]["to"]["primaryAccession"]

    return accession


def collect_unique_gene_symbols(nsforest_results):

    gene_symbols = set()

    for column in ["NSForest_markers", "binary_genes"]:
        for gene_list_str in nsforest_results[column]:
            gene_symbols |= set(ast.literal_eval(gene_list_str))

    return list(gene_symbols)


def collect_unique_gene_ids(gene_symbols, gnm2ids):

    gene_symbols = set(gene_symbols)

    gene_ids = set()

    for gene_symbol in gene_symbols:

        gene_ids |= set(map_gene_name_to_ids(gene_symbol, gnm2ids))

    print(
        f"Collected {len(gene_ids)} unique gene ids for {len(gene_symbols)} unique gene symbols"
    )

    return list(gene_ids)


def get_opentargets_results(nsforest_path, resources=RESOURCES):

    opentargets_path = Path(str(nsforest_path).replace(".csv", "-opentargets.json"))

    if not opentargets_path.exists():

        opentargets_results = {}

        print(f"Loading NSForest results from {nsforest_path}")
        nsforest_results = load_results(nsforest_path).sort_values(
            "clusterName", ignore_index=True
        )

        gene_symbols = collect_unique_gene_symbols(nsforest_results)
        gnm2ids = get_gene_name_to_ids_map()
        gene_ids = collect_unique_gene_ids(gene_symbols, gnm2ids)

    else:

        print(f"Loading opentargets results from {opentargets_path}")
        with open(opentargets_path, "r") as fp:
            opentargets_results = json.load(fp)

        gene_ids = opentargets_results["gene_ids"]

    n_to_dump = 25
    do_dump = False
    n_fetched = 0
    for gene_id in gene_ids:
        n_fetched += 1

        if gene_id not in opentargets_results:
            do_dump = True

            opentargets_results[gene_id] = {}

            for resource in resources:
                try:
                    opentargets_results[gene_id][resource] = gget.opentargets(
                        gene_id, resource=resource, json=True, verbose=True
                    )
                    print(
                        f"Assigned gget opentargets resource {resource} for gene id {gene_id} ({n_fetched} of {n_to_dump})"
                    )
                except Exception as exc:
                    print(
                        f"Could not assign gget opentargets resource {resource} for gene id {gene_id} ({n_fetched} of {n_to_dump})"
                    )
                    opentargets_results[gene_id][resource] = {}

        else:
            print(
                f"Already assigned gget opentargets resources for gene id {gene_id} ({n_fetched} of {n_to_dump})"
            )
            if gene_id != gene_ids[-1]:
                continue

        if do_dump and (n_fetched >= n_to_dump or gene_id == gene_ids[-1]):
            do_dump = False
            n_fetched = 0

            opentargets_results["gene_ids"] = gene_ids

            print(f"Dumping opentargets results to {opentargets_path}")
            with open(opentargets_path, "w") as fp:
                json.dump(opentargets_results, fp, indent=4)

    return opentargets_path, opentargets_results


def collect_unique_drug_names(opentargets_results):

    drug_names = set()

    for gene_id, resources in opentargets_results.items():
        if gene_id == "gene_ids":
            continue
        for drug in resources["drugs"]:
            drug_names.add(drug["name"])

    return list(drug_names)


def get_ebi_results(opentargets_path, resources=RESOURCES):

    ebi_path = Path(str(opentargets_path).replace("opentargets", "ebi"))

    if not ebi_path.exists():

        ebi_results = {}

        nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

        opentargets_path, opentargets_results = get_opentargets_results(
            nsforest_path, resources=resources
        )

        drug_names = collect_unique_drug_names(opentargets_results)

    else:

        print(f"Loading ebi results from {ebi_path}")
        with open(ebi_path, "r") as fp:
            ebi_results = json.load(fp)

        drug_names = ebi_results["drug_names"]

    n_to_dump = 25
    do_dump = False
    n_fetched = 0
    for drug_name in drug_names:
        n_fetched += 1

        if drug_name not in ebi_results:
            do_dump = True

            response = requests.get(f"https://www.ebi.ac.uk/ols/api/search?q={drug_name}&ontology=dron")
            if response.status_code == 200:
                print(
                    f"Assigned EBI results for drug name {drug_name} ({n_fetched} of {n_to_dump})"
                )
                ebi_results[drug_name] = response.json()

            else:
                print(
                    f"Could not assign EBI results for protein id {drug_name} ({n_fetched} of {n_to_dump})"
                )
                ebi_results[drug_name] = {}

        else:
            print(
                f"Already assigned EBI results for protein id {drug_name} ({n_fetched}) or {n_to_dump}"
            )
            if drug_name != drug_names[-1]:
                continue

        if do_dump and (n_fetched >= n_to_dump or drug_name == drug_names[-1]):
            do_dump = False
            n_fetched = 0

            ebi_results["drug_names"] = drug_names

            print(f"Dumping ebi results to {ebi_path}")
            with open(ebi_path, "w") as fp:
                json.dump(ebi_results, fp, indent=4)

    return ebi_path, ebi_results


def collect_unique_protein_ids(opentargets_results):

    protein_ids = set()

    for gene_id, resources in opentargets_results.items():
        if gene_id == "gene_ids":
            continue
        for interaction in resources["interactions"]:
            for key in ["protein_a_id", "protein_b_id"]:
                protein_ids |= set([interaction[key]])

    return list(protein_ids)


def get_uniprot_results(opentargets_path, resources=RESOURCES):

    uniprot_path = Path(str(opentargets_path).replace("opentargets", "uniprot"))

    if not uniprot_path.exists():

        uniprot_results = {}

        nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

        opentargets_path, opentargets_results = get_opentargets_results(
            nsforest_path, resources=resources
        )

        protein_ids = collect_unique_protein_ids(opentargets_results)

    else:

        print(f"Loading uniprot results from {uniprot_path}")
        with open(uniprot_path, "r") as fp:
            uniprot_results = json.load(fp)

        protein_ids = uniprot_results["protein_ids"]

    n_to_dump = 25
    do_dump = False
    n_fetched = 0
    for protein_id in protein_ids:
        n_fetched += 1

        if protein_id not in uniprot_results:
            do_dump = True

            if "ENSP" in protein_id:
                accession = map_protein_id_to_accession(protein_id)
                if accession is None:
                    uniprot_results[protein_id] = {}
                    continue
                print(f"Mapped accession {accession} to protein id {protein_id}")

            else:
                accession = protein_id

            response = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{accession}?fields=accession,protein_name,cc_function,ft_binding"
            )
            if response.status_code == 200:
                print(
                    f"Assigned UniProt results for protein id {protein_id} ({n_fetched} of {n_to_dump})"
                )
                uniprot_results[protein_id] = response.json()

            else:
                print(
                    f"Could not assign UniProt results for protein id {protein_id} ({n_fetched} of {n_to_dump})"
                )
                uniprot_results[protein_id] = {}

        else:
            print(
                f"Already assigned UniProt results for protein id {protein_id} ({n_fetched}) or {n_to_dump}"
            )
            if protein_id != protein_ids[-1]:
                continue

        if do_dump and (n_fetched >= n_to_dump or protein_id == protein_ids[-1]):
            do_dump = False
            n_fetched = 0

            uniprot_results["protein_ids"] = protein_ids

            print(f"Dumping uniprot results to {uniprot_path}")
            with open(uniprot_path, "w") as fp:
                json.dump(uniprot_results, fp, indent=4)

    return uniprot_path, uniprot_results


def main():

    nsforest_path = (
        NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-guo-2023-2025-02-22.csv"
    ).resolve()

    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    # uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)

    ebi_path, ebi_results = get_ebi_results(opentargets_path)


if __name__ == "__main__":
    main()
