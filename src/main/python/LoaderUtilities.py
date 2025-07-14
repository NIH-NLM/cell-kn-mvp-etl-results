import ast
import hashlib
import json
from pathlib import Path
import random
import string

from lxml import etree
import pandas as pd
import scanpy as sc

from E_Utilities import find_gene_id_for_gene_name
from OntologyParserLoader import parse_term
from UniProtIdMapper import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)

ALPHABET = string.ascii_lowercase + string.digits
PURLBASE = "http://purl.obolibrary.org/obo"
RDFSBASE = "http://www.w3.org/1999/02/22-rdf-syntax-ns"

OWL_NS = "{http://www.w3.org/2002/07/owl#}"
OBO_IN_OWL_NS = "{http://www.geneontology.org/formats/oboInOwl#}"
RDF_NS = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"


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


def get_gene_names_and_ensembl_ids():
    """Query BioMart to get gene names and Ensembl ids.

    Parameters
    ----------
    None

    Returns
    -------
    gnmsids : pd.DataFrame
        DataFrame with columns containing gene names and Ensembl ids
    """
    print("Getting gene names and Ensembl ids from BioMart")
    gnmsids = sc.queries.biomart_annotations(
        "hsapiens", ["external_gene_name", "ensembl_gene_id"], use_cache=True
    )
    return gnmsids


def get_gene_name_to_ensembl_ids_map():
    """Get gene name to Ensembl ids map.

    Parameters
    ----------
    None

    Returns
    -------
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Ensembl id
    """
    print("Creating gene name to Ensembl ids map")
    gnmsids = get_gene_names_and_ensembl_ids()
    gnm2ids = gnmsids.set_index("external_gene_name")
    return gnm2ids


def map_gene_name_to_ensembl_ids(name, gnm2ids):
    """Map a gene name to a gene Ensembl id list.

    Parameters
    ----------
    name : str
        Gene name
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Ensembl id

    Returns
    -------
    list
        Gene Ensembl ids
    """
    if name in gnm2ids.index:
        ids = gnm2ids.loc[name, "ensembl_gene_id"]
        if isinstance(ids, pd.core.series.Series):
            ids = ids.to_list()
        else:
            ids = [ids]
        # print(f"Mapped gene name {name} to Ensembl ids {ids}")
    else:
        print(f"Could not find gene Ensembl ids for gene name: {name}")
        ids = []
    return ids


def get_gene_ensembl_id_to_names_map():
    """Map gene Ensembl id to names.

    Parameters
    ----------
    None

    Returns
    -------
    gid2nms : pd.DataFrame
        DataFrame indexed by gene Ensembl ids containing gene names
    """
    print("Creating gene Ensembl id to names map")
    gnmsids = get_gene_names_and_ensembl_ids()
    gid2nms = gnmsids.set_index("ensembl_gene_id")
    return gid2nms


def map_gene_ensembl_id_to_names(gid, gid2nms):
    """Map a gene Ensembl id to a gene name list.

    Parameters
    ----------
    gid : str
        Gene Ensembl id
    gid2nms : pd.DataFrame
        DataFrame indexed by gene Ensembl id containing gene name

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
        # print(f"Mapped gene Ensembl id {gid} to names {names}")
    else:
        print(f"Could not find gene names for gene Ensembl id: {gid}")
        names = []
    return names


def get_protein_ensembl_id_to_accession_map(protein_ids):
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


def map_protein_ensembl_id_to_accession(ensp, ensp2accn):
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


def get_protein_accession_to_ensembl_id_map(protein_ids):
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


def map_accession_to_protein_ensembl_id(accn, accn2ensp):
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


def get_gene_name_to_and_from_entrez_id_maps(gene_symbols):
    """Get gene name to Entrez id map, and its reverse. Cache results
    after each success to prevent duplicate requests on restart.

    Parameters
    ----------
    gene_symbols : list(str)
        List of gene symbols

    Returns
    -------
    gnm2id : dict
        Dictionary with name keys and id values
    gid2nm : dict
        Dictionary with id keys and name values
    """
    gnm2id = {}
    gid2nm = {}

    # Initialize a cache to prevent duplicate requests on restart
    hasher = hashlib.sha256()
    gene_symbols = sorted(set(gene_symbols))
    hasher.update(":".join(gene_symbols).encode("utf-8"))
    cache_path = Path(f".cache/{hasher.hexdigest()}.json")
    if cache_path.exists():
        with open(cache_path, "r") as fp:
            cache = json.load(fp)
        gnm2id = cache["gnm2id"]
        gid2nm = cache["gid2nm"]

    else:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache = {}
        cache["gnm2id"] = gnm2id
        cache["gid2nm"] = gid2nm
        with open(cache_path, "w") as fp:
            json.dump(cache, fp, indent=4)

    print("Creating gene name to and from Entrez id maps")
    for gene_symbol in gene_symbols:

        # Skip cached results
        if gene_symbol in gnm2id:
            continue

        gene_id = find_gene_id_for_gene_name(gene_symbol)
        gnm2id[gene_symbol] = gene_id
        gid2nm[gene_id] = gene_symbol

        # Cache current results
        cache["gnm2id"] = gnm2id
        cache["gid2nm"] = gid2nm
        with open(cache_path, "w") as fp:
            json.dump(cache, fp, indent=4)

    return gnm2id, gid2nm


def collect_unique_gene_ensembl_ids(gene_symbols, gnm2ids):
    """Collect unique Ensembl gene ids corresponding to the specified
    list of gene symbols.

    Parameters
    ----------
    gene_symbols : list(str)
        List of gene symbols
    gnm2ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Ensembl id

    Returns
    -------
    gene_ids : list(str)
        List of unique gene Ensembl ids
    """
    gene_ids = set()

    gene_symbols = set(gene_symbols)
    for gene_symbol in gene_symbols:
        gene_ids |= set(map_gene_name_to_ensembl_ids(gene_symbol, gnm2ids))
    print(
        f"Collected {len(gene_ids)} unique Ensembl gene ids for {len(gene_symbols)} unique gene symbols"
    )

    return list(gene_ids)


def collect_unique_gene_entrez_ids(gene_symbols, gnm2id):
    """Collect unique Entrez gene ids corresponding to the specified
    list of gene symbols.

    Parameters
    ----------
    gene_symbols : list(str)
        List of gene symbols
    gnm2id : dict
        Dictionary with name keys and id values

    Returns
    -------
    gene_ids : list(str)
        List of unique gene Entrez ids
    """
    gene_ids = set()

    gene_symbols = set(gene_symbols)
    for gene_symbol in gene_symbols:
        gene_id = gnm2id[gene_symbol]
        if gene_id:
            gene_ids.add(gene_id)
    print(
        f"Collected {len(gene_ids)} unique Entrez gene ids for {len(gene_symbols)} unique gene symbols"
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
        "../../../cell-kn-mvp-etl-ontologies/data/mondo_efo_mappings.tsv"
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
        # print(f"Could not find MONDO for EFO term: {efo}")
        return None
    return mondo


def get_mesh_to_mondo_map(obo_dir, obo_fnm):
    """Parse MONDO ontology XML downloaded from the OBO Foundry to
    create a mapping from MeSH term to MONDO term.

    Parameters
    ----------
    obo_dir : str | Path
        Name of directory containing downloaded MONDO ontology XML
    obo_fnm : str
        Name of downloaded MONDO ontology XML file

    Returns
    -------
    mesh2mondo : dict
        Dictionary mapping MeSH term to MONDO term
    """
    mesh2mondo = {}
    root = etree.parse(Path(obo_dir) / obo_fnm)
    for class_element in root.iter(f"{OWL_NS}Class"):

        # Look for an about attribute
        uriref = class_element.get(f"{RDF_NS}about")
        if uriref is None:
            continue

        id, number, mondo_term, _, _ = parse_term(uriref)
        if id is None:
            continue

        for hasDbXref_element in class_element.iter(f"{OBO_IN_OWL_NS}hasDbXref"):
            if hasDbXref_element is None:
                continue
            mesh_term = hasDbXref_element.text
            if "MESH" in mesh_term:
                mesh2mondo[mesh_term] = mondo_term
                break

    # https://meshb.nlm.nih.gov/record/ui?ui=D000077192
    # http://purl.obolibrary.org/obo/MONDO_0004991
    mesh2mondo["MESH:D000077192"] = "MONDO_0004991"

    # https://meshb.nlm.nih.gov/record/ui?ui=D000086382
    # http://purl.obolibrary.org/obo/MONDO_0100096
    mesh2mondo["MESH:D000086382"] = "MONDO_0100096"

    # https://meshb.nlm.nih.gov/record/ui?ui=D003643
    # http://purl.obolibrary.org/obo/UBERON_0000071
    mesh2mondo["MESH:D003643"] = "UBERON_0000071"

    # https://meshb.nlm.nih.gov/record/ui?ui=D005355
    # http://purl.obolibrary.org/obo/MONDO_0002771
    mesh2mondo["MESH:D005355"] = "MONDO_0002771"

    return mesh2mondo


def map_mesh_to_mondo(mesh, mesh2mondo):
    """Map MeSH term to MONDO term.

    Parameters
    ----------
    mesh : str
        MeSH term
    mesh2mondo : dict
        Dictionary mapping MeSH term to MONDO term

    Returns
    -------
    mondo : str
        MONDO term
    """
    mondo = None

    if mesh in mesh2mondo:
        mondo = mesh2mondo[mesh]

    return mondo


def get_value_or_none(data, keys):
    """Return the value in the data corresponding to the last key, or
    None, if any key is not in the data.

    Parameters
    ----------
    data : dict
        Dictionary which may or may not contain the keys
    keys : list(str)
        List of keys to access the dictionary in order
    """
    value = None
    for key in keys:
        try:
            if value is None:
                value = data[key]
            else:
                value = value[key]
        except:
            return None
    return value


def get_values_or_none(data, list_key, value_keys):
    """Collect and return the values for each list item in the data
    corresponding to the list key, and last value key.

    Parameters
    ----------
    data : dict
        Dictionary which may or may not contain the keys
    list_key : str
        Key of the list of items
    value_keys : list(str)
        List of keys to access each item in order
    """
    values = ""
    if list_key in data:
        for item in data[list_key]:
            value = get_value_or_none(item, value_keys)
            if values == "":
                values = value
            else:
                values += ", " + value
    return values
