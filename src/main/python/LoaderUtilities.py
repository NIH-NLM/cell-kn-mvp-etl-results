import ast
import random
import string

import pandas as pd
import scanpy as sc

from UniProtIdMapper import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)

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
    """Load results CSV file and append a UUID.

    Parameters
    ----------
    results_Path : Path
        Path of results CSV file

    Returns
    -------
    results : pd.DataFrame
        DataFrame containing results
    """
    results = pd.read_csv(results_path)
    if "uuid" not in results.columns:
        print(f"Add UUID column to results CSV file {results_path.name}")
        results["uuid"] = [get_uuid() for idx in results.index]
        results.to_csv(results_path)
    return results


def hyphenate(iname):
    """Replace spaces, underscores, commas and forward slashes with
    hyphens, but only one.

    Parameters
    ----------
    iname : str
        Input name

    Returns
    -------
    oname : str
        Output name
    """
    cname = iname
    for c in [" ", "_", ",", "/"]:
        cname = cname.replace(c, "-")
        oname = cname.replace("--", "-")
        while cname != oname:
            cname = oname
            oname = oname.replace("--", "-")
    return oname


def get_gene_names_and_ids():
    """Query BioMart to get gene names and ids.

    Parameters
    ----------
    None

    Returns
    -------
    gnmsids : pd.DataFrame
        DataFrame with columns containing gene names and ids
    """
    print("Getting gene names and ids from BioMart")
    gnmsids = sc.queries.biomart_annotations(
        "hsapiens", ["external_gene_name", "ensembl_gene_id"], use_cache=True
    )
    return gnmsids


def get_gene_name_to_ids_map():
    """Get gene name to ids map.

    Parameters
    ----------
    None

    Returns
    -------
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene id
    """
    print("Creating gene name to ids map")
    gnmsids = get_gene_names_and_ids()
    gnm2ids = gnmsids.set_index("external_gene_name")
    return gnm2ids


def map_gene_name_to_ids(name, gnm2ids):
    """Map a gene name to a gene id list.

    Parameters
    ----------
    name : str
        Gene name
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene id

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
        # print(f"Mapped gene name {name} to ids {ids}")
    else:
        print(f"Could not find gene ids for gene name: {name}")
        ids = []
    return ids


def get_gene_id_to_names_map():
    """Map gene id to names.

    Parameters
    ----------
    None

    Returns
    -------
    gid2nms : pd.DataFrame
        DataFrame indexed by gene ids containing gene names
    """
    print("Creating gene id to names map")
    gnmsids = get_gene_names_and_ids()
    gid2nms = gnmsids.set_index("ensembl_gene_id")
    return gid2nms


def map_gene_id_to_names(gid, gid2nms):
    """Map a gene id to a gene name list.

    Parameters
    ----------
    gid : str
        Gene id
    gid2nms : pd.DataFrame
        DataFrame indexed by gene id containing gene name

    Returns
    -------
    list
        Gene names
    """
    if gid in gid2nms.index:
        names = gid2nms.loc[gid, "external_gene_name"]
        if isinstance(names, pd.core.series.Series):
            names = names.to_list()
        else:
            names = [names]
        # print(f"Mapped gene id {gid} to names {names}")
    else:
        print(f"Could not find gene names for gene id: {gid}")
        names = []
    return names


def get_protein_id_to_accession_map(protein_ids):
    """Map Ensembl protein ids to UniProt accession lists.

    Parameters
    ----------
    protein_ids : list(str)
        Protein ids returned by gget opentargets command

    Returns
    -------
    ensp2accn : dict
        Dictionary mapping Ensembl protein ids to UniProt accession
        lists
    """
    ensp2accn = {}

    # Submit Ensembl ids in batches to the UniProt id mapping service
    batch_size = 1000
    ensps = []
    for protein_id in protein_ids:
        if "ENSP" in protein_id:
            ensps.append(protein_id)

        if len(ensps) == batch_size or (
            len(ensps) > 0 and protein_id == protein_ids[-1]
        ):
            # Submit full, or the last batch
            job_id = submit_id_mapping(
                from_db="Ensembl_Protein", to_db="UniProtKB", ids=ensps
            )
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                data = get_id_mapping_results_search(link)

            # Collect the mapping results
            for result in data["results"]:
                ensp = result["from"]
                accn = result["to"]["primaryAccession"]
                if ensp not in ensp2accn:
                    ensp2accn[ensp] = accn
                else:
                    if not isinstance(ensp2accn[ensp], list):
                        ensp2accn[ensp] = [ensp2accn[ensp]]
                    ensp2accn[ensp].append(accn)

            # Initialize for the next batch
            ensps = []

    return ensp2accn


def map_protein_id_to_accession(ensp, ensp2accn):
    """Map Ensembl protein id to UniProt accession, selecting the
    first if more than one found.

    Parameters
    ----------
    ensp : str
        Ensembl protein id
    ensp2accn : dict
        Dictionary mapping Ensembl protein ids to UniProt accession
        lists

    Returns
    -------
    accn : str
        UniProt accession
    """
    accn = None

    if ensp in ensp2accn:
        accn = ensp2accn[ensp]
        if isinstance(accn, list):
            accn = accn[0]

    return accn


def get_accession_to_protein_id_map(protein_ids):
    """Map UniProt accession to Ensembl protein ids lists.

    Parameters
    ----------
    protein_ids : list(str)
        Protein ids returned by gget opentargets command

    Returns
    -------
    accn2esnp : dict
        Dictionary mapping UniProt accession to Ensembl protein ids
        lists
    """
    accn2esnp = {}

    # Submit UniProt accessions in batches to the UniProt id mapping
    # service
    batch_size = 1000
    accns = []
    for protein_id in protein_ids:
        if "ENSP" not in protein_id:
            accns.append(protein_id)

        if len(accns) == batch_size or (
            len(accns) > 0 and protein_id == protein_ids[-1]
        ):
            # Submit full, or the last batch
            job_id = submit_id_mapping(
                from_db="UniProtKB_AC-ID", to_db="Ensembl_Protein", ids=accns
            )
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                data = get_id_mapping_results_search(link)

            # Collect the mapping results
            for result in data["results"]:
                accn = result["from"]
                ensp = result["to"]
                if accn not in accn2esnp:
                    accn2esnp[accn] = ensp
                else:
                    if not isinstance(accn2esnp[accn], list):
                        accn2esnp[accn] = [accn2esnp[accn]]
                    accn2esnp[accn].append(ensp)

            # Initialize for the next batch
            accns = []

    return accn2esnp


def map_accession_to_protein_id(accn, accn2ensp):
    """Map UniProt accession to Ensembl protein id, selecting the
    first if more than one found.

    Parameters
    ----------
    accn : str
        UniProt accession
    accn2esnp : dict
        Dictionary mapping UniProt accession to Ensembl protein ids
        lists

    Returns
    -------
    ensp : str
        Ensembl protein id
    """
    ensp = None

    if accn in accn2ensp:
        ensp = accn2ensp[accn]
        if isinstance(ensp, list):
            ensp = ensp[0]

    return ensp


def collect_unique_gene_symbols(nsforest_results):
    """Collect unique gene symbols found in the NSForest results
    marker or binary genes.

    Parameters
    ----------
    nsforest_results : pd.DataFrame
        DataFrame containing NSForest results

    Returns
    -------
    gene_symbols : list(str)
        List of unique gene symbols
    """
    gene_symbols = set()

    for column in ["NSForest_markers", "binary_genes"]:
        for gene_list_str in nsforest_results[column]:
            gene_symbols |= set(ast.literal_eval(gene_list_str))

    return list(gene_symbols)


def collect_unique_gene_ids(gene_symbols, gnm2ids):
    """Collect unique Ensembl gene ids corresponding to the specified list of gene symbols.

    Parameters
    ----------
    gene_symbols : list(str)
        List of unique gene symbols
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene id

    Returns
    -------
    gene_ids : list(str)
        List of unique gene ids
    """
    gene_ids = set()

    gene_symbols = set(gene_symbols)
    for gene_symbol in gene_symbols:
        gene_ids |= set(map_gene_name_to_ids(gene_symbol, gnm2ids))
    print(
        f"Collected {len(gene_ids)} unique gene ids for {len(gene_symbols)} unique gene symbols"
    )

    return list(gene_ids)


def get_efo_to_mondo_map():
    """Get EFO to MONDO term map.

    Parameters
    ----------
    None

    Returns
    -------
    efo2mondo : pd.DataFrame
        DataFrame indexed by EFO containing MONDO term
    """
    print("Creating EFO to MONDO term map")
    mondo_efo_mappings_name = (
        "../../../cell-kn-etl-ontologies/data/mondo_efo_mappings.tsv"
    )
    efo2mondo = pd.read_csv(mondo_efo_mappings_name)
    efo2mondo = efo2mondo.set_index("EFO")
    return efo2mondo


def map_efo_to_mondo(efo, efo2mondo):
    """Map EFO to MONDO term.

    Parameters
    ----------
    efo : str
        EFO term
    efo2mondo : pd.DataFrame
        DataFrame indexed by EFO containing MONDO term

    Returns
    -------
    str
        MONDO term
    """
    if efo in efo2mondo.index:
        mondo = efo2mondo.loc[efo, "MONDO"]
    else:
        print(f"Could not find MONDO for EFO term: {efo}")
        return None
    return mondo
