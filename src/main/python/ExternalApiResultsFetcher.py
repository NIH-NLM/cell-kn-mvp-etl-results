import argparse
from glob import glob
import json
import os
from pathlib import Path
import re
import shutil

import gget
import requests

from LoaderUtilities import (
    load_results,
    collect_unique_gene_ids,
    collect_unique_gene_symbols,
    get_accession_to_protein_id_map,
    get_gene_name_to_ids_map,
    get_protein_id_to_accession_map,
    map_protein_id_to_accession,
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

NSFOREST_DIRPATH = Path("../../../data/results")
HUBMAP_DIRPATH = Path("../../../data/hubmap")
HUBMAP_LATEST_URLS = [
    "https://lod.humanatlas.io/asct-b/allen-brain/latest/",
    "https://lod.humanatlas.io/asct-b/eye/latest/",
    "https://lod.humanatlas.io/asct-b/lung/latest/",
]

# TODO: Refactor to reduce massive redundancy


def get_opentargets_results(nsforest_path, resources=RESOURCES, force=False):
    """Use the gget opentargets command to obtain the specified
    resources for each unique gene id mapped from each gene symbol in
    the NSForest results loaded from the specified path. The
    opentarget results are written out in batches to enable
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
        Dictionary containg opentargets results keyed by gene id, then
        by resource
    """
    # Create, or load opentargets results
    opentargets_path = Path(str(nsforest_path).replace(".csv", "-opentargets.json"))
    if not opentargets_path.exists() or force:

        # Initialize results, and collect unique gene symbols and ids

        opentargets_results = {}

        print(f"Loading NSForest results from {nsforest_path}")
        nsforest_results = load_results(nsforest_path).sort_values(
            "clusterName", ignore_index=True
        )

        gene_symbols = collect_unique_gene_symbols(nsforest_results)
        gnm2ids = get_gene_name_to_ids_map()
        gene_ids = collect_unique_gene_ids(gene_symbols, gnm2ids)

    else:

        # Load results, and assign unique gene symbols and ids

        print(f"Loading opentargets results from {opentargets_path}")
        with open(opentargets_path, "r") as fp:
            opentargets_results = json.load(fp)

        gene_symbols = opentargets_results["gene_symbols"]
        gene_ids = opentargets_results["gene_ids"]

    # Consider each gene id, and setup to dump the results in
    # batches, and enable restarting
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

            opentargets_results["gene_symbols"] = gene_symbols
            opentargets_results["gene_ids"] = gene_ids

            print(f"Dumping opentargets results to {opentargets_path}")
            with open(opentargets_path, "w") as fp:
                json.dump(opentargets_results, fp, indent=4)

    return opentargets_path, opentargets_results


def collect_unique_drug_names(opentargets_results):
    """Collect unique drug names contained in the opentargets results.

    Parameters
    ----------
    opentargets_results : dict
        Dictionary containg opentargets results keyed by gene id, then
        by resource

    Returns
    -------
    drug_names : list(str)
        List of unique drum names
    """
    drug_names = set()

    for gene_id, resources in opentargets_results.items():
        if gene_id in ["gene_ids", "gene_symbols"]:
            continue
        for drug in resources["drugs"]:
            drug_names.add(drug["name"])

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
        Dictionary containg EBI results keyed by drug name
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
        Dictionary containg rxnav results keyed by drug name
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
                    print(f"Assigned RxNav {content} results for drug name {drug_name}")
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
    """Get the value for the specified property name contained in the
    RxNav results for the specified drug name.

    Parameters
    ----------
    rxnav_results : dict
        Dictionary containg rxnav results keyed by drug name
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
        Dictionary containg RxNav results keyed by drug name

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
        Dictionary containg NCATS results keyed by drug name

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
    """Collect unique protein ids contained in the opentargets results.

    Parameters
    ----------
    opentargets_results : dict
        Dictionary containg opentargets results keyed by gene id, then
        by resource

    Returns
    -------
    protein_ids : list
        List of unique protein ids

    Notes
    -----
    Protein ids may be either Ensembl ids or UniProt accessions
    """
    protein_ids = set()

    for gene_id, resources in opentargets_results.items():
        if gene_id in ["gene_ids", "gene_symbols"]:
            continue
        for interaction in resources["interactions"]:
            for key in ["protein_a_id", "protein_b_id"]:
                protein_ids |= set([interaction[key]])

    return list(protein_ids)


def get_uniprot_results(opentargets_path, resources=RESOURCES, force=False):
    """Use a UniProt API endpoint for each protein id in the
    opentargets results loaded from the specified path. The UniProt
    results are written out in batches to enable restarting.

    Parameters
    ----------
    opentargets_path : Path
        Path to opentarget results
    resources : list
        List of resource names to use with the gget opentargets
        command
    force : bool
        Flag to force fetching, or not

    Returns
    -------
    uniprot_path : Path
        Path to UniProt results
    uniprot_results : dict
        Dictionary containg UniProt results keyed by protein id
    """
    # Create, or load UniProt results
    uniprot_path = Path(str(opentargets_path).replace("opentargets", "uniprot"))
    if not uniprot_path.exists() or force:

        # Initialize results, and collect unique protein ids, and
        # their mapping to and from Ensembl protein ids and UniProg
        # accessions

        uniprot_results = {}

        nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

        opentargets_path, opentargets_results = get_opentargets_results(
            nsforest_path, resources=resources
        )

        protein_ids = collect_unique_protein_ids(opentargets_results)
        ensp2accn = get_protein_id_to_accession_map(protein_ids)
        accn2ensp = get_accession_to_protein_id_map(protein_ids)

    else:

        # Load results, and assign unique protein ids, and their
        # mapping to and from Ensembl protein ids and UniProg
        # accessions

        print(f"Loading uniprot results from {uniprot_path}")
        with open(uniprot_path, "r") as fp:
            uniprot_results = json.load(fp)

        protein_ids = uniprot_results["protein_ids"]
        ensp2accn = uniprot_results["ensp2accn"]
        accn2ensp = uniprot_results["accn2ensp"]

    # Consider each protein id, and setup to dump the results in
    # batches, and enable restarting
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

                # Map Ensembl id to UniProt accession
                accession = map_protein_id_to_accession(protein_id, ensp2accn)
                if accession is None:
                    uniprot_results[protein_id] = {}
                    continue
                print(f"Mapped accession {accession} to protein id {protein_id}")

            else:

                # Assume protein id is a UniProt accession
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
            uniprot_results["ensp2accn"] = ensp2accn
            uniprot_results["accn2ensp"] = accn2ensp

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
    """Load NSForest results from processing datasets corresponding to
    the Guo et al. 2023, Jorstad et al. 2023, Li et al. 2023, and
    Sikkema, et al. 2023, publications, then

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

    - Use a UniProt API endpoint for each protein id in the
      opentargets results to obtain other protein ids, descriptions,
      and comments

    Download specified latest HuBMAP data table JSON files, archiving
    any earlier versions.

    Parameters
    ----------
    None

    Returns
    -------
    opentargets_results : dict
        Dictionary containg opentargets results keyed by gene id, then
        by resource
    ebi_results : dict
        Dictionary containg EBI results keyed by drug name
    rxnav_results : dict
        Dictionary containg rxnav results keyed by drug name
    uniprot_results : dict
        Dictionary containg UniProt results keyed by protein id
    """
    # Load NSForest results and fetch external API results
    parser = argparse.ArgumentParser(description="Fetch External API Results")
    parser.add_argument(
        "--force-opentargets",
        action="store_true",
        help="force fetching of opentargets results",
    )
    parser.add_argument(
        "--force-ebi",
        action="store_true",
        help="force fetching of ebi results",
    )
    parser.add_argument(
        "--force-rxnav",
        action="store_true",
        help="force fetching of rxnav results",
    )
    parser.add_argument(
        "--force-drugbank",
        action="store_true",
        help="force fetching of drugbank results",
    )
    parser.add_argument(
        "--force-ncats",
        action="store_true",
        help="force fetching of ncats results",
    )
    parser.add_argument(
        "--force-uniprot",
        action="store_true",
        help="force fetching of uniprot results",
    )
    args = parser.parse_args()

    nsforest_paths = [
        Path(p).resolve()
        for p in glob(str(NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-*.csv"))
    ]
    for nsforest_path in nsforest_paths:

        print(f"Fetching results for {nsforest_path}")

        opentargets_path, _opentargets_results = get_opentargets_results(
            nsforest_path, force=args.force_opentargets
        )

        _ebi_path, _ebi_results = get_ebi_results(
            opentargets_path, force=args.force_ebi
        )

        _rxnav_path, _rxnav_results = get_rxnav_results(
            opentargets_path, force=args.force_rxnav
        )

        # TODO: Restore if API becomes available
        # drugbank_path, drugbank_results = get_drugbank_results(rxnav_path, force=args.force_drugbank)

        # TODO: Restore if API becomes available
        # ncats_path, ncats_results = get_ncats_results(rxnav_path, force=args.force_ncats)

        _uniprot_path, _uniprot_results = get_uniprot_results(
            opentargets_path, force=args.force_uniprot
        )

    # Download specified latest HuBMAP data table JSON files
    download_hubmap_data_tables()


if __name__ == "__main__":
    main()
