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

# TODO: Refactor to reduce massive redundancy


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


def get_protein_id_to_accession_map(protein_ids):

    pid2acc = {}

    batch_size = 1000

    ensp_ids = []
    for protein_id in protein_ids:

        if "ENSP" in protein_id:
            ensp_ids.append(protein_id)

        if len(ensp_ids) == batch_size or (
            len(ensp_ids) > 0 and protein_id == protein_ids[-1]
        ):

            job_id = submit_id_mapping(
                from_db="Ensembl_Protein", to_db="UniProtKB", ids=ensp_ids
            )
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                data = get_id_mapping_results_search(link)

            for result in data["results"]:
                pid2acc[result["from"]] = result["to"]["primaryAccession"]

            ensp_ids = []

    return pid2acc


def map_protein_id_to_accession(ensp_id, pid2acc):

    acc = None

    if ensp_id in pid2acc:
        acc = pid2acc[ensp_id]

    return acc


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

    total_size = len(gene_ids)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for gene_id in gene_ids:
        n_in_batch += 1
        n_so_far += 1
        if gene_id not in opentargets_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            opentargets_results[gene_id] = {}

            for resource in resources:
                try:
                    opentargets_results[gene_id][resource] = gget.opentargets(
                        gene_id, resource=resource, json=True, verbose=True
                    )
                    print(
                        f"Assigned gget opentargets resource {resource} for gene id {gene_id}"
                    )
                except Exception as exc:
                    print(
                        f"Could not assign gget opentargets resource {resource} for gene id {gene_id}"
                    )
                    opentargets_results[gene_id][resource] = {}

        else:
            # print(f"Already assigned gget opentargets resources for gene id {gene_id}")
            if gene_id != gene_ids[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or gene_id == gene_ids[-1]):
            do_dump = False
            n_in_batch = 0

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

    total_size = len(drug_names)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for drug_name in drug_names:
        n_in_batch += 1
        n_so_far += 1
        if drug_name not in ebi_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            response = requests.get(
                f"https://www.ebi.ac.uk/ols/api/search?q={drug_name}&ontology=dron"
            )
            if response.status_code == 200:
                print(f"Assigned EBI results for drug name {drug_name}")
                ebi_results[drug_name] = response.json()

            else:
                print(f"Could not assign EBI results for drug name {drug_name}")
                ebi_results[drug_name] = {}

        else:
            # print(f"Already assigned EBI results for drug name {drug_name}")
            if drug_name != drug_names[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or drug_name == drug_names[-1]):
            do_dump = False
            n_in_batch = 0

            ebi_results["drug_names"] = drug_names

            print(f"Dumping ebi results to {ebi_path}")
            with open(ebi_path, "w") as fp:
                json.dump(ebi_results, fp, indent=4)

    return ebi_path, ebi_results


def get_rxnav_results(opentargets_path, resources=RESOURCES):

    rxnav_path = Path(str(opentargets_path).replace("opentargets", "rxnav"))

    if not rxnav_path.exists():

        rxnav_results = {}

        nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

        opentargets_path, opentargets_results = get_opentargets_results(
            nsforest_path, resources=resources
        )

        drug_names = collect_unique_drug_names(opentargets_results)

    else:

        print(f"Loading RxNav results from {rxnav_path}")
        with open(rxnav_path, "r") as fp:
            rxnav_results = json.load(fp)

        drug_names = rxnav_results["drug_names"]

    total_size = len(drug_names)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for drug_name in drug_names:
        n_in_batch += 1
        n_so_far += 1
        if drug_name not in rxnav_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            rxnav_results[drug_name] = {}

            urls = [
                f"https://rxnav.nlm.nih.gov/REST/rxcui.json?name={drug_name}",
                f"https://rxnav.nlm.nih.gov/REST/spellingsuggestions.json?name={drug_name}",
                f"https://rxnav.nlm.nih.gov/REST/Prescribe/drugs.json?name={drug_name}",
            ]
            for url in urls:
                response = requests.get(url)
                if response.status_code == 200:
                    content = url.split("/")[-1].split("?")[0].replace(".json", "")
                    print(f"Assigned RxNav {content} results for drug name {drug_name}")
                    rxnav_results[drug_name].update(response.json())

                else:
                    print(
                        f"Could not assign RxNav {content} results for drug name {drug_name}"
                    )

            if "rxnormId" in rxnav_results[drug_name]["idGroup"]:

                rxcui = rxnav_results[drug_name]["idGroup"]["rxnormId"][0]

                urls = [
                    f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/properties.json",
                    f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/allProperties.json?prop=names+codes",
                ]
                for url in urls:
                    response = requests.get(url)
                    if response.status_code == 200:
                        content = url.split("/")[-1].split("?")[0].replace(".json", "")
                        print(
                            f"Assigned RxNav {content} results for drug name {drug_name}"
                        )
                        rxnav_results[drug_name].update(response.json())

                    else:
                        print(
                            f"Could not assign RxNav {content} results for drug name {drug_name}"
                        )

        else:
            # print(f"Already assigned RxNav results for drug name {drug_name}")
            if drug_name != drug_names[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or drug_name == drug_names[-1]):
            do_dump = False
            n_in_batch = 0

            rxnav_results["drug_names"] = drug_names

            print(f"Dumping RxNav results to {rxnav_path}")
            with open(rxnav_path, "w") as fp:
                json.dump(rxnav_results, fp, indent=4)

    return rxnav_path, rxnav_results


def get_prop_for_drug(rxnav_results, drug_name, prop_name):

    # "DRUGBANK" or "UNII_CODE"

    prop_value = None

    if drug_name not in rxnav_results:
        print(f"No RxNav results for drug name {drug_name}")

    elif "propConceptGroup" not in rxnav_results[drug_name]:
        print(f"No property group in RxNav results for drug name {drug_name}")

    else:
        for propConcept in rxnav_results[drug_name]["propConceptGroup"]["propConcept"]:
            if propConcept["propName"] == prop_name:
                prop_value = propConcept["propValue"]
                break

    return prop_value


def get_drugbank_results(rxnav_path, resources=RESOURCES):

    drugbank_path = Path(str(rxnav_path).replace("rxnav", "drugbank"))

    if not drugbank_path.exists():

        drugbank_results = {}

        opentargets_path = Path(str(rxnav_path).replace("rxnav", "opentargets"))

        rxnav_path, rxnav_results = get_rxnav_results(
            opentargets_path, resources=resources
        )

        opentargets_path, opentargets_results = get_opentargets_results(
            opentargets_path, resources=resources
        )

        drug_names = collect_unique_drug_names(opentargets_results)

    else:

        print(f"Loading DrugBank results from {drugbank_path}")
        with open(drugbank_path, "r") as fp:
            drugbank_results = json.load(fp)

        drug_names = drugbank_results["drug_names"]

    total_size = len(drug_names)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for drug_name in drug_names:
        n_in_batch += 1
        n_so_far += 1
        if drug_name not in drugbank_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            drugbank_results[drug_name] = {}

            drugbank_id = get_prop_for_drug(rxnav_results, drug_name, "DRUGBANK")

            response = requests.get(f"https://go.drugbank.com/drugs/{drugbank_id}")
            if response.status_code == 200:
                print(f"Assigned DrugBank results for drug name {drug_name}")
                drugbank_results[drug_name].update(response.json())

            else:
                print(f"Could not assign DrugBank results for drug name {drug_name}")

        else:
            # print(f"Already assigned DrugBank results for drug name {drug_name}")
            if drug_name != drug_names[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or drug_name == drug_names[-1]):
            do_dump = False
            n_in_batch = 0

            drugbank_results["drug_names"] = drug_names

            print(f"Dumping DrugBank results to {drugbank_path}")
            with open(drugbank_path, "w") as fp:
                json.dump(drugbank_results, fp, indent=4)

    return drugbank_path, drugbank_results


def get_ncats_results(rxnav_path, resources=RESOURCES):

    ncats_path = Path(str(rxnav_path).replace("rxnav", "ncats"))

    if not ncats_path.exists():

        ncats_results = {}

        opentargets_path = Path(str(rxnav_path).replace("rxnav", "opentargets"))

        rxnav_path, rxnav_results = get_rxnav_results(
            opentargets_path, resources=resources
        )

        opentargets_path, opentargets_results = get_opentargets_results(
            opentargets_path, resources=resources
        )

        drug_names = collect_unique_drug_names(opentargets_results)

    else:

        print(f"Loading Ncats results from {ncats_path}")
        with open(ncats_path, "r") as fp:
            ncats_results = json.load(fp)

        drug_names = ncats_results["drug_names"]

    total_size = len(drug_names)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for drug_name in drug_names:
        n_in_batch += 1
        n_so_far += 1
        if drug_name not in ncats_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            ncats_results[drug_name] = {}

            unii_code = get_prop_for_drug(rxnav_results, drug_name, "UNII_CODE")

            response = requests.get(f"https://drugs.ncats.io/drug/{unii_code}")
            if response.status_code == 200:
                print(f"Assigned Ncats results for drug name {drug_name}")
                ncats_results[drug_name].update(response.json())

            else:
                print(f"Could not assign Ncats results for drug name {drug_name}")

        else:
            # print(f"Already assigned Ncats results for drug name {drug_name}")
            if drug_name != drug_names[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or drug_name == drug_names[-1]):
            do_dump = False
            n_in_batch = 0

            ncats_results["drug_names"] = drug_names

            print(f"Dumping Ncats results to {ncats_path}")
            with open(ncats_path, "w") as fp:
                json.dump(ncats_results, fp, indent=4)

    return ncats_path, ncats_results


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
        pid2acc = get_protein_id_to_accession_map(protein_ids)

    else:

        print(f"Loading uniprot results from {uniprot_path}")
        with open(uniprot_path, "r") as fp:
            uniprot_results = json.load(fp)

        protein_ids = uniprot_results["protein_ids"]
        pid2acc = uniprot_results["pid2acc"]

    total_size = len(protein_ids)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for protein_id in protein_ids:
        n_in_batch += 1
        n_so_far += 1
        if protein_id not in uniprot_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            if "ENSP" in protein_id:
                accession = map_protein_id_to_accession(protein_id, pid2acc)
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
                print(f"Assigned UniProt results for protein id {protein_id}")
                uniprot_results[protein_id] = response.json()

            else:
                print(f"Could not assign UniProt results for protein id {protein_id}")
                uniprot_results[protein_id] = {}

        else:
            # print(f"Already assigned UniProt results for protein id {protein_id}")
            if protein_id != protein_ids[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or protein_id == protein_ids[-1]):
            do_dump = False
            n_in_batch = 0

            uniprot_results["protein_ids"] = protein_ids
            uniprot_results["pid2acc"] = pid2acc

            print(f"Dumping uniprot results to {uniprot_path}")
            with open(uniprot_path, "w") as fp:
                json.dump(uniprot_results, fp, indent=4)

    return uniprot_path, uniprot_results


def main():

    nsforest_path = (
        NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-guo-2023-2025-02-22.csv"
    ).resolve()

    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    ebi_path, ebi_results = get_ebi_results(opentargets_path)

    rxnav_path, rxnav_results = get_rxnav_results(opentargets_path)

    # TODO: Restore if API becomes available
    # drugbank_path, drugbank_results = get_drugbank_results(rxnav_path)

    # TODO: Restore if API becomes available
    # ncats_path, ncats_results = get_ncats_results(rxnav_path)

    uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)

    return opentargets_results, ebi_results, rxnav_results, uniprot_results


if __name__ == "__main__":
    opentargets_results, ebi_results, rxnav_results, uniprot_results = main()
