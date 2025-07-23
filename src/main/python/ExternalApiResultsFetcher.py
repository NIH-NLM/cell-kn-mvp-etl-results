import argparse
from glob import glob
import json
import os
from pathlib import Path
import re
import shutil

import gget
import requests

from E_Utilities import get_data_for_gene_id
from OpenTargetsGGetQueries import gget_queries
from LoaderUtilities import (
    load_results,
    collect_unique_gene_ensembl_ids,
    collect_unique_gene_entrez_ids,
    collect_unique_gene_symbols,
    get_gene_name_to_and_from_entrez_id_maps,
    get_gene_name_to_ensembl_ids_map,
    get_value_or_none,
    get_values_or_none,
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
BASE_URL = "https://api.platform.opentargets.org/api/v4/graphql"

NSFOREST_DIRPATH = Path("../../../data/results")
HUBMAP_DIRPATH = Path("../../../data/hubmap")
HUBMAP_LATEST_URLS = [
    "https://lod.humanatlas.io/asct-b/allen-brain/latest/",
    "https://lod.humanatlas.io/asct-b/eye/latest/",
    "https://lod.humanatlas.io/asct-b/lung/latest/",
]

# TODO: Refactor to reduce massive redundancy


def get_cellxgene_metadata(author_to_cl_path, force=False):
    """Use the CELLxGENE curation API to fetch the metadata for the
    collection and dataset corresponding to the manual author cell set
    to CL term mapping results loaded from the specified path.

    Parameters
    ----------
    author_to_cl_path : Path
        Path to author to CL results
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    cellxgene_path : Path
        Path to cellxgene results
    cellxgene_results : list(dict)
        List of dictionaries containing cellxgene results
    """
    # Create, or load cellxgene results
    cellxgene_path = Path(str(author_to_cl_path).replace(".csv", "-cellxgene.json"))
    if not cellxgene_path.exists() or force:

        # Create results

        cellxgene_results = []

        print(f"Loading author to CL results from {author_to_cl_path}")
        author_to_cl_results = load_results(author_to_cl_path)
        collection_id = author_to_cl_results["collection_id"][0]
        collection_version_id = author_to_cl_results["collection_version_id"][0]
        datasets_ids = zip(
            author_to_cl_results["dataset_id"][0].split("--"),
            author_to_cl_results["dataset_version_id"][0].split("--"),
        )
        base_url = "https://api.cellxgene.cziscience.com/curation/v1"
        for dataset_id, dataset_version_id in datasets_ids:
            dataset_results = {}
            dataset_url = (
                f"{base_url}/collections/{collection_id}/datasets/{dataset_id}"
            )
            response = requests.get(dataset_url)
            if response.status_code == 200:
                response_json = response.json()

                print(
                    f"Assigning cellxgene metadata for collection_id {collection_id} and dataset_id {dataset_id}"
                )
                dataset_results["Link_to_publication"] = None
                dataset_results["Link_to_CELLxGENE_collection"] = None
                citation = get_value_or_none(response_json, ["citation"])
                if citation:
                    m = re.search(r"Publication:\s*(\S*)\s*Dataset Version:", citation)
                    if m:
                        dataset_results["Link_to_publication"] = m.group(1)
                    m = re.search(r"Collection:\s*(\S*)$", citation)
                    if m:
                        dataset_results["Link_to_CELLxGENE_collection"] = m.group(1)
                dataset_results["Link_to_CELLxGENE_dataset"] = response_json["assets"][
                    0
                ]["url"]
                dataset_results["Dataset_name"] = get_value_or_none(
                    response_json, ["title"]
                )
                dataset_results["Number_of_cells"] = get_value_or_none(
                    response_json, ["cell_count"]
                )
                dataset_results["Organism"] = get_values_or_none(
                    response_json, "organism", ["label"]
                )
                dataset_results["Tissue"] = get_values_or_none(
                    response_json, "tissue", ["label"]
                )
                dataset_results["Disease_status"] = get_values_or_none(
                    response_json, "disease", ["label"]
                )
                dataset_results["Collection_ID"] = collection_id
                dataset_results["Collection_version_ID"] = collection_version_id
                dataset_results["Dataset_ID"] = dataset_id
                dataset_results["Dataset_version_ID"] = dataset_version_id
                dataset_results["Zenodo/Nextflow_workflow/Notebook"] = "TBC"
                cellxgene_results.append(dataset_results)

            else:
                print(
                    f"Could not assign cellxgene metadata for collection_id {collection_id} and dataset_id {dataset_id}"
                )

            print(f"Dumping cellxgene results to {cellxgene_path}")
            with open(cellxgene_path, "w") as fp:
                json.dump(cellxgene_results, fp, indent=4)

    else:

        # Load results

        print(f"Loading cellxgene results from {cellxgene_path}")
        with open(cellxgene_path, "r") as fp:
            cellxgene_results = json.load(fp)

    return cellxgene_path, cellxgene_results


def get_opentargets_results(nsforest_path, resources=RESOURCES, force=False):
    """Use the gget opentargets command to obtain the specified
    resources for each unique Ensembl gene id mapped from each gene
    symbol in the NSForest results loaded from the specified path. The
    opentargets results are written out in batches to enable
    restarting.

    Parameters
    ----------
    nsforest_path : Path
        Path to NSForest results
    resources : list(str)
        List of resource names to use with the gget opentargets
        command
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    opentargets_path : Path
        Path to opentargets results
    opentargets_results : dict
        Dictionary containing opentargets results keyed by gene
        Ensembl id, then by resource
    """
    # Create, or load opentargets results
    opentargets_path = Path(str(nsforest_path).replace(".csv", "-opentargets.json"))
    if not opentargets_path.exists() or force:

        # Initialize results, and collect unique gene symbols and Ensembl ids

        opentargets_results = {}

        print(f"Loading NSForest results from {nsforest_path}")
        nsforest_results = load_results(nsforest_path).sort_values(
            "clusterName", ignore_index=True
        )

        gene_symbols = collect_unique_gene_symbols(nsforest_results)
        gnm2ids = get_gene_name_to_ensembl_ids_map()
        gene_ensembl_ids = collect_unique_gene_ensembl_ids(gene_symbols, gnm2ids)

    else:

        # Load results, and assign unique gene symbols and ids

        print(f"Loading opentargets results from {opentargets_path}")
        with open(opentargets_path, "r") as fp:
            opentargets_results = json.load(fp)

        gene_symbols = opentargets_results["gene_symbols"]
        gene_ensembl_ids = opentargets_results["gene_ensembl_ids"]

    # Consider each gene id, and setup to dump the results in
    # batches, and enable restarting
    total_size = len(gene_ensembl_ids)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for gene_ensembl_id in gene_ensembl_ids:
        n_in_batch += 1
        n_so_far += 1
        if gene_ensembl_id not in opentargets_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            opentargets_results[gene_ensembl_id] = {}

            try:
                query = gget_queries["target"]
                query["variables"]["ensemblId"] = gene_ensembl_id
                response = requests.post(
                    BASE_URL,
                    json={
                        "query": query["query_string"],
                        "variables": query["variables"],
                    },
                )
                response.raise_for_status()
                print(
                    f"Assigning gget opentargets resources for gene Ensembl id {gene_ensembl_id}"
                )
                data = json.loads(response.text)["data"]
                opentargets_results[gene_ensembl_id]["target"] = {}
                for key in [
                    "id",
                    "dbXrefs",
                    "proteinIds",
                    "transcriptIds",
                    "approvedSymbol",
                    "approvedName",
                ]:
                    opentargets_results[gene_ensembl_id]["target"][key] = data[
                        "target"
                    ][key]
                for resource in resources:
                    if resource == "diseases":
                        resource_data = data["target"]["associatedDiseases"]["rows"]

                    elif resource == "drugs":
                        resource_data = data["target"]["knownDrugs"]["rows"]

                    elif resource == "interactions":
                        resource_data = data["target"]["interactions"]["rows"]

                    elif resource == "pharmacogenetics":  # Not a typo
                        resource_data = data["target"]["pharmacogenomics"]

                    elif resource == "tractability":
                        resource_data = data["target"]["tractability"]

                    elif resource == "expression":
                        resource_data = data["target"]["expressions"]

                    elif resource == "depmap":
                        resource_data = data["target"]["depMapEssentiality"]

                    opentargets_results[gene_ensembl_id][resource] = resource_data

            except Exception as exc:
                print(
                    f"Could not assign gget opentargets resources for gene Ensembl id {gene_ensembl_id}"
                )
                opentargets_results[gene_ensembl_id]["target"] = {}
                for resource in resources:
                    opentargets_results[gene_ensembl_id][resource] = {}

        else:
            # print(f"Already assigned gget opentargets resources for gene Ensembl id {gene_ensembl_id}")
            if gene_ensembl_id != gene_ensembl_ids[-1]:
                continue

        if do_dump and (
            n_in_batch >= batch_size or gene_ensembl_id == gene_ensembl_ids[-1]
        ):
            do_dump = False
            n_in_batch = 0

            opentargets_results["gene_symbols"] = gene_symbols
            opentargets_results["gene_ensembl_ids"] = gene_ensembl_ids

            print(f"Dumping opentargets results to {opentargets_path}")
            with open(opentargets_path, "w") as fp:
                json.dump(opentargets_results, fp, indent=4)

    return opentargets_path, opentargets_results


def collect_unique_drug_names(opentargets_results):
    """Collect unique drug names contained in the opentargets results.

    Parameters
    ----------
    opentargets_results : dict
        Dictionary containing opentargets results keyed by gene
        Ensembl id, then by resource

    Returns
    -------
    drug_names : list(str)
        List of unique drum names
    """
    drug_names = set()

    for gene_ensembl_id, resources in opentargets_results.items():
        if gene_ensembl_id in ["gene_ensembl_ids", "gene_symbols"]:
            continue
        for drug in resources["drugs"]:
            drug_names.add(drug["approvedName"])

    return list(drug_names)


def get_ebi_results(opentargets_path, resources=RESOURCES, force=False):
    """Use an EBI API endpoint for each unique drug name in the
    opentargets results loaded from the specified path. The EBI
    results are written out in batches to enable restarting.

    Parameters
    ----------
    opentargets_path : Path
        Path to opentarget results
    resources : list(str)
        List of resource names to use with the gget opentargets
        command
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    ebi_path : Path
        Path to EBI results
    ebi_results : dict
        Dictionary containing EBI results keyed by drug name
    """
    # Create, or load EBI results
    ebi_path = Path(str(opentargets_path).replace("opentargets", "ebi"))
    if not ebi_path.exists() or force:

        # Initialize results, and collect unique drug names

        ebi_results = {}

        nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))
        opentargets_path, opentargets_results = get_opentargets_results(
            nsforest_path, resources=resources
        )

        drug_names = collect_unique_drug_names(opentargets_results)

    else:

        # Load results, and assign unique drug names

        print(f"Loading ebi results from {ebi_path}")
        with open(ebi_path, "r") as fp:
            ebi_results = json.load(fp)

        drug_names = ebi_results["drug_names"]

    # Consider each drug name, and setup to dump the results in
    # batches, and enable restarting
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
                print(f"Assigning EBI results for drug name {drug_name}")
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


def get_rxnav_results(opentargets_path, resources=RESOURCES, force=False):
    """Use an RxNav API endpoint for each unique drug name in the
    opentargets results loaded from the specified path. The RxNav
    results are written out in batches to enable restarting.

    Parameters
    ----------
    opentargets_path : Path
        Path to opentarget results
    resources : list(str)
        List of resource names to use with the gget opentargets
        command
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    rxnav_path : Path
        Path to rxnav results
    rxnav_results : dict
        Dictionary containing rxnav results keyed by drug name
    """
    # Create, or load RxNav results
    rxnav_path = Path(str(opentargets_path).replace("opentargets", "rxnav"))
    if not rxnav_path.exists() or force:

        # Initialize results, and collect unique drug names

        rxnav_results = {}

        nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

        opentargets_path, opentargets_results = get_opentargets_results(
            nsforest_path, resources=resources
        )

        drug_names = collect_unique_drug_names(opentargets_results)

    else:

        # Load results, and assign unique drug names

        print(f"Loading RxNav results from {rxnav_path}")
        with open(rxnav_path, "r") as fp:
            rxnav_results = json.load(fp)

        drug_names = rxnav_results["drug_names"]

    # Consider each drug name, and setup to dump the results in
    # batches, and enable restarting
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

            # Get mapping from drug name to RXCUI, suggested
            # spellings, and prescribable drugs information
            urls = [
                f"https://rxnav.nlm.nih.gov/REST/rxcui.json?name={drug_name}",
                f"https://rxnav.nlm.nih.gov/REST/spellingsuggestions.json?name={drug_name}",
                f"https://rxnav.nlm.nih.gov/REST/Prescribe/drugs.json?name={drug_name}",
            ]
            for url in urls:
                response = requests.get(url)
                if response.status_code == 200:
                    content = url.split("/")[-1].split("?")[0].replace(".json", "")
                    print(
                        f"Assigning RxNav {content} results for drug name {drug_name}"
                    )
                    rxnav_results[drug_name].update(response.json())

                else:
                    print(
                        f"Could not assign RxNav {content} results for drug name {drug_name}"
                    )

            # Use the RXCUI to get drug properties
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
                            f"Assigning RxNav {content} results for drug name {drug_name}"
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
    """Get the value for the specified property name contained in the
    RxNav results for the specified drug name.

    Parameters
    ----------
    rxnav_results : dict
        Dictionary containing rxnav results keyed by drug name
    drug_name : str
        Drug name key, currently only "DRUGBANK" or "UNII_CODE"
        expected
    prop_name : str
        RxNav results property name key

    Returns
    -------
    prov_value : str
        RxNav results property name value
    """
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


def get_drugbank_results(rxnav_path, resources=RESOURCES, force=False):
    """Use the DrugBank website for each unique drug name in the
    opentargets results loaded from the path corresponding to the
    specified RxNav path. The DrugBank results are written out in
    batches to enable restarting. Drug names are mapped to DrugBank
    ids using the RxNav results.

    Parameters
    ----------
    rxnav_path : Path
        Path to RxNav results
    resources : list(str)
        List of resource names to use with the gget opentargets
        command
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    durgbank_path : Path
        Path to DrugBank results
    drugbank_results : dict
        Dictionary containing RxNav results keyed by drug name

    Notes
    -----
    As currently written, this function will not work. The DrubBank
    website needs to be replaced by the DrugBank API, which requires a
    license.
    """
    # Create, or load RxNav results
    drugbank_path = Path(str(rxnav_path).replace("rxnav", "drugbank"))
    if not drugbank_path.exists() or force:

        # Initialize results, and collect unique drug names

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

        # Load results, and assign unique drug names

        print(f"Loading DrugBank results from {drugbank_path}")
        with open(drugbank_path, "r") as fp:
            drugbank_results = json.load(fp)

        drug_names = drugbank_results["drug_names"]

    # Consider each drug name, and setup to dump the results in
    # batches, and enable restarting
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

            # Map drug name to DrugBank id
            drugbank_id = get_prop_for_drug(rxnav_results, drug_name, "DRUGBANK")

            response = requests.get(f"https://go.drugbank.com/drugs/{drugbank_id}")
            if response.status_code == 200:
                print(f"Assigning DrugBank results for drug name {drug_name}")
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


def get_ncats_results(rxnav_path, resources=RESOURCES, force=False):
    """Use the NCATS website for each unique drug name in the
    opentargets results loaded from the path corresponding to the
    specified RxNav path. The NCATS results are written out in batches
    to enable restarting.

    Parameters
    ----------
    rxnav_path : Path
        Path to RxNav results
    resources : list
        List of resource names to use with the gget opentargets
        command
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    ncats_path : Path
        Path to NCATS results
    ncats_results : dict
        Dictionary containing NCATS results keyed by drug name

    Notes
    -----
    As currently written, this function will not work. The NCATS
    website needs to be replaced by an NCATS API, which is currently
    unkown.
    """
    # Create, or load NCATS results
    ncats_path = Path(str(rxnav_path).replace("rxnav", "ncats"))
    if not ncats_path.exists() or force:

        # Initialize results, and collect unique drug names

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

        # Load results, and assign unique drug names

        print(f"Loading Ncats results from {ncats_path}")
        with open(ncats_path, "r") as fp:
            ncats_results = json.load(fp)

        drug_names = ncats_results["drug_names"]

    # Consider each drug name, and setup to dump the results in
    # batches, and enable restarting
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

            # Map drug name to UNII_CODE
            unii_code = get_prop_for_drug(rxnav_results, drug_name, "UNII_CODE")

            response = requests.get(f"https://drugs.ncats.io/drug/{unii_code}")
            if response.status_code == 200:
                print(f"Assigning Ncats results for drug name {drug_name}")
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


def get_gene_results(nsforest_path, force=False):
    """Use the E-Utilities to fetch Gene data for each gene symbol in
    the NSForest results loaded from the specified path. The gene
    results are written out in batches to enable restarting.

    Parameters
    ----------
    nsforest_path : Path
        Path to NSForest results
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    gene_path : Path
        Path to gene results
    gene_results : dict
        Dictionary containing gene results keyed by gene symbol
    """
    # Create, or load gene results
    gene_path = Path(str(nsforest_path).replace(".csv", "-gene.json"))
    if not gene_path.exists() or force:

        # Initialize results, and collect unique gene symbols and Entrez ids

        gene_results = {}

        print(f"Loading NSForest results from {nsforest_path}")
        nsforest_results = load_results(nsforest_path).sort_values(
            "clusterName", ignore_index=True
        )

        gene_symbols = collect_unique_gene_symbols(nsforest_results)
        gnm2id, _gid2nm = get_gene_name_to_and_from_entrez_id_maps(gene_symbols)
        gene_entrez_ids = collect_unique_gene_entrez_ids(gene_symbols, gnm2id)

    else:

        # Load results, and assign unique gene symbols

        print(f"Loading gene results from {gene_path}")
        with open(gene_path, "r") as fp:
            gene_results = json.load(fp)

        gene_symbols = gene_results["gene_symbols"]
        gnm2id = gene_results["gnm2id"]
        gene_entrez_ids = gene_results["gene_entrez_ids"]

    # Consider each gene symbol, and setup to dump the results in
    # batches, and enable restarting
    total_size = len(gene_symbols)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for gene_symbol in gene_symbols:
        n_in_batch += 1
        n_so_far += 1
        if gene_symbol not in gene_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True
            try:
                print(f"Assigning gene data for gene symbol {gene_symbol}")
                gene_results[gene_symbol] = get_data_for_gene_id(gnm2id[gene_symbol])

            except Exception as exc:
                print(f"Could not assign gene data for gene symbol {gene_symbol}")
                gene_results[gene_symbol] = {}

        else:
            # print(f"Already assigned gene data for gene symbol {gene_symbol}")
            if gene_symbol != gene_symbols[-1]:
                continue

        if do_dump and (n_in_batch >= batch_size or gene_symbol == gene_symbols[-1]):
            do_dump = False
            n_in_batch = 0

            gene_results["gene_symbols"] = gene_symbols
            gene_results["gnm2id"] = gnm2id
            gene_results["gene_entrez_ids"] = gene_entrez_ids

            print(f"Dumping gene results to {gene_path}")
            with open(gene_path, "w") as fp:
                json.dump(gene_results, fp, indent=4)

    return gene_path, gene_results


def collect_unique_protein_accessions(gene_results):
    """Collect unique protein accessions contained in the gene
    results.

    Parameters
    ----------
    gene_results : dict
        Dictionary containing gene results keyed by gene symbol

    Returns
    -------
    protein_accessions : list
        List of unique protein accessions

    """
    protein_accessions = set()

    for gene_symbol, gene_data in gene_results.items():
        if gene_symbol in ["gene_symbols", "gnm2id", "gene_entrez_ids"] or not gene_data:
            continue
        protein_accessions |= set([gene_data["UniProt_name"]])

    return list(protein_accessions)


def get_uniprot_results(gene_path, force=False):
    """Use a UniProt API endpoint for each protein accession in the
    gene results loaded from the specified path. The UniProt results
    are written out in batches to enable restarting.

    Parameters
    ----------
    gene_path : Path
        Path to gene results
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    uniprot_path : Path
        Path to UniProt results
    uniprot_results : dict
        Dictionary containing UniProt results keyed by protein
        accession
    """
    # Create, or load UniProt results
    uniprot_path = Path(str(gene_path).replace("gene", "uniprot"))
    if not uniprot_path.exists() or force:

        # Initialize results, and collect unique protein accessions

        uniprot_results = {}

        nsforest_path = Path(str(gene_path).replace("-gene.json", ".csv"))

        gene_path, gene_results = get_gene_results(nsforest_path)

        protein_accessions = collect_unique_protein_accessions(gene_results)

    else:

        # Load results, and assign unique protein accessions

        print(f"Loading uniprot results from {uniprot_path}")
        with open(uniprot_path, "r") as fp:
            uniprot_results = json.load(fp)

        protein_accessions = uniprot_results["protein_accessions"]

    # Consider each protein accession, and setup to dump the results
    # in batches, and enable restarting
    total_size = len(protein_accessions)
    n_so_far = 0
    do_dump = False
    batch_size = 25
    n_in_batch = 0
    for protein_accession in protein_accessions:
        n_in_batch += 1
        n_so_far += 1
        if protein_accession not in uniprot_results:
            print(
                f"Fetched {n_in_batch}/{batch_size} in batch - {n_so_far}/{total_size} so far"
            )
            do_dump = True

            response = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{protein_accession}"
            )
            if response.status_code == 200:
                print(
                    f"Assigning UniProt results for protein accession {protein_accession}"
                )
                response_json = response.json()
                data = {}
                data["Protein_name"] = get_value_or_none(
                    response_json,
                    [
                        "proteinDescription",
                        "recommendedName",
                        "fullName",
                        "value",
                    ],
                )
                data["UniProt_ID"] = get_value_or_none(
                    response_json, ["primaryAccession"]
                )
                data["Gene_name"] = None
                if "genes" in response_json and len(response_json["genes"]) > 0:
                    data["Gene_name"] = get_value_or_none(
                        response_json["genes"][0],
                        [
                            "geneName",
                            "value",
                        ],
                    )
                data["Number_of_amino_acids"] = get_value_or_none(
                    response_json,
                    [
                        "sequence",
                        "length",
                    ],
                )
                data["Function"] = None
                if "comments" in response_json:
                    for comment in response_json["comments"]:
                        if (
                            "commentType" in comment
                            and comment["commentType"] == "FUNCTION"
                        ):
                            if "texts" in comment and len(comment["texts"]) > 0:
                                data["Function"] = get_value_or_none(
                                    comment["texts"][0], ["value"]
                                )
                data["Annotation_score"] = get_value_or_none(
                    response_json, ["annotationScore"]
                )
                data["Organism"] = get_value_or_none(
                    response_json,
                    [
                        "organism",
                        "scientificName",
                    ],
                )
                uniprot_results[protein_accession] = data

            else:
                print(
                    f"Could not assign UniProt results for protein accession {protein_accession}"
                )
                uniprot_results[protein_accession] = {}

        else:
            # print(f"Already assigned UniProt results for protein accession {protein_accession}")
            if protein_accession != protein_accessions[-1]:
                continue

        if do_dump and (
            n_in_batch >= batch_size or protein_accession == protein_accessions[-1]
        ):
            do_dump = False
            n_in_batch = 0

            uniprot_results["protein_accessions"] = protein_accessions

            print(f"Dumping uniprot results to {uniprot_path}")
            with open(uniprot_path, "w") as fp:
                json.dump(uniprot_results, fp, indent=4)

    return uniprot_path, uniprot_results


def get_hubmap_json_urls():
    """Get the URL to specified HuBMAP data table JSON files.

    Parameters
    ----------
    None

    Returns
    -------
    json_urls : list(tuple(str, float, str))
       List of tuples with the organ, version, and URL
    """
    json_urls = []

    # Get each HuBMAP data table latest version URL
    p_org = re.compile(r"asct-b\/(.*)\/latest")
    p_url = re.compile(r"https:\/\/.*\/v(\d\.\d)\/graph.json")
    for latest_url in HUBMAP_LATEST_URLS:
        m_org = p_org.search(latest_url)
        if m_org is not None:
            org = m_org.group(1)

        else:
            # Should never happen
            raise Exception("No organ in HuBMAP URL")
        response = requests.get(latest_url)
        if response.status_code == 200:

            # Parse the response to find version and JSON file URL
            m_url = p_url.search(response.text)
            if m_url is not None:
                json_url = m_url.group(0)
                json_ver = float(m_url.group(1))
                json_urls.append((org, json_ver, json_url))

            else:
                raise Exception("Could not find HuBMAP JSON URL or version")

        else:
            raise Exception("Could not get HuBMAP latest URL")

    return json_urls


def download_hubmap_data_tables():
    """Download specified latest HuBMAP data table JSON files,
    archiving any earlier versions.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Get the URL to all HuBMAP data table JSON files
    json_urls = get_hubmap_json_urls()
    for org, ver, url in json_urls:

        # Skip the current JSON file if it exists, otherwise, archive
        # any earlier versions
        hubmap_filepath = HUBMAP_DIRPATH / f"{org}-v{ver}.json"
        if hubmap_filepath.exists():
            print(f"HuBMAP data table {hubmap_filepath} already exists")
            continue

        else:
            for pathname in glob(str(HUBMAP_DIRPATH / f"{org}-v*.json")):
                try:
                    shutil.move(Path(pathname), HUBMAP_DIRPATH / ".archive")
                    print(f"Archived HuBMAP data table {pathname}")
                except Exception as exc:
                    # Since already archived, though should never happen
                    os.remove(pathname)
                    print(f"Removed HuBMAP data table {pathname}")

        # Download the JSON file
        response = requests.get(url)
        if response.status_code == 200:
            with open(hubmap_filepath, "w") as fp:
                fp.write(response.text)
            print(f"Downloaded HuBMAP data table {hubmap_filepath}")

        else:
            print(f"Could not download HuBMAP data table {hubmap_filepath}")


def main():
    """For datasets corresponding to the Guo et al. 2023, Jorstad et
    al. 2023, Li et al. 2023, and Sikkema, et al. 2023, publications,
    use the CELLxGENE curation API to fetch the metadata for each
    collection and dataset, and load NSForest results from processing
    each dataset, then

    - Use the gget opentargets command to obtain the diseases, drugs,
      interactions, pharmacogenetics, tractability, expression, and
      depmap resources for each unique gene id mapped from each gene
      symbol in the NSForest results

    - Use an EBI API endpoint to obtain drug ontology data for each
      unique drug name in the opentargets results

    - Use an RxNav API endpoint for each unique drug name in the
      opentargets results to obtain the mapping from drug name to
      RXCUI, suggested spellings, prescribable drugs information, and
      drug properties

    - [Use the DrugBank website for each unique drug name in the
      opentargets results loaded from the path corresponding to the
      specified RxNav path]

    - [Use the NCATS website for each unique drug name in the
      opentargets results loaded from the path corresponding to the
      specified RxNav path]

    - Use the E-Utilities to fetch Gene data for each gene symbol in
      the NSForest results

    - Use a UniProt API endpoint for each protein accession in the
      gene results to obtain protein data

    Download specified latest HuBMAP data table JSON files, archiving
    any earlier versions.

    Parameters
    ----------
    None

    Returns
    -------
    opentargets_results : dict
        Dictionary containing opentargets results keyed by gene id,
        then by resource
    ebi_results : dict
        Dictionary containing EBI results keyed by drug name
    rxnav_results : dict
        Dictionary containing rxnav results keyed by drug name
    uniprot_results : dict
        Dictionary containing UniProt results keyed by protein
        accession
    """
    # Load NSForest results and fetch external API results
    parser = argparse.ArgumentParser(description="Fetch External API Results")
    parser.add_argument(
        "-c",
        "--force-cellxgene",
        action="store_true",
        help="force fetching of cellxgene results",
    )
    parser.add_argument(
        "-o",
        "--force-opentargets",
        action="store_true",
        help="force fetching of opentargets results",
    )
    parser.add_argument(
        "-e",
        "--force-ebi",
        action="store_true",
        help="force fetching of ebi results",
    )
    parser.add_argument(
        "-r",
        "--force-rxnav",
        action="store_true",
        help="force fetching of rxnav results",
    )
    parser.add_argument(
        "-d",
        "--force-drugbank",
        action="store_true",
        help="force fetching of drugbank results",
    )
    parser.add_argument(
        "-n",
        "--force-ncats",
        action="store_true",
        help="force fetching of ncats results",
    )
    parser.add_argument(
        "-g",
        "--force-gene",
        action="store_true",
        help="force fetching of gene results",
    )
    parser.add_argument(
        "-u",
        "--force-uniprot",
        action="store_true",
        help="force fetching of uniprot results",
    )
    parser.add_argument(
        "--force-all",
        action="store_true",
        help="force fetching of all results",
    )
    args = parser.parse_args()

    nsforest_paths = [
        Path(p).resolve()
        for p in glob(str(NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-*.csv"))
    ]
    for nsforest_path in nsforest_paths:

        # Map NSForest results filename to manual author cell set to
        # CL term mapping filename
        author = re.search("results-([a-zA-Z]*)", nsforest_path.name).group(1)
        author_to_cl_path = Path(
            glob(
                str(NSFOREST_DIRPATH / f"cell-kn-mvp-map-author-to-cl-{author}-*.csv")
            )[-1]
        ).resolve()

        print(f"Fetching CELLxGENE results {author_to_cl_path}")

        get_cellxgene_metadata(
            author_to_cl_path, force=args.force_cellxgene or args.force_all
        )

        print(f"Fetching results for {nsforest_path}")

        opentargets_path, _opentargets_results = get_opentargets_results(
            nsforest_path, force=args.force_opentargets or args.force_all
        )

        _ebi_path, _ebi_results = get_ebi_results(
            opentargets_path, force=args.force_ebi or args.force_all
        )

        _rxnav_path, _rxnav_results = get_rxnav_results(
            opentargets_path, force=args.force_rxnav or args.force_all
        )

        # TODO: Restore if API becomes available
        # drugbank_path, drugbank_results = get_drugbank_results(
        #     rxnav_path, force=args.force_drugbank or args.force_all
        # )

        # TODO: Restore if API becomes available
        # ncats_path, ncats_results = get_ncats_results(
        #     rxnav_path, force=args.force_ncats or args.force_all
        # )

        gene_path, _gene_results = get_gene_results(
            nsforest_path, force=args.force_gene or args.force_all
        )

        _uniprot_path, _uniprot_results = get_uniprot_results(
            gene_path, force=args.force_uniprot or args.force_all
        )

    # Download specified latest HuBMAP data table JSON files
    download_hubmap_data_tables()


if __name__ == "__main__":
    main()
