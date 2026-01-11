import ast
from glob import glob
import json
from pathlib import Path
import random
import re
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

with open(Path("../../../data/obo/deprecated_terms.txt"), "r") as fp:
    DEPRECATED_TERMS = fp.read().splitlines()

MIN_CLUSTER_SIZE = 10

RESULTS_SOURCES_PATH = Path("../../../data/results-sources-2026-01-06.json")


def get_cl_terms(author_to_cl_results):
    """Create a set of clean CL terms from the given author to CL results.

    Parameters
    ----------
    author_cl_results : pd.DataFrame
        DataFrame containing author to CL results

    Returns
    -------
    set(str)
        Set of clean CL terms
    """
    return set(
        author_to_cl_results.loc[
            author_to_cl_results["cell_ontology_id"].str.contains("CL"),
            "cell_ontology_id",
        ]
        .str.replace("http://purl.obolibrary.org/obo/", "")
        .str.replace("https://purl.obolibrary.org/obo/", "")
    )


def collect_results_sources_data():
    """Collect paths to all NSForest results, silhouette scores, and
    author cell set to CL term mappings identified in the results
    sources. Collect the dataset_version_id used for creating the
    NSForest results. Collect the unique gene names, Ensembl
    identifiers, and Entrez identifiers corresponding to all NSForet
    results.

    Parameters
    ----------
    None

    Returns
    -------
    nsforest_paths: list(Path)
        List of paths to all NSForest results files
    silhouette_paths: list(Path)
        List of paths to all silhouette scores files
    author_to_cl_paths: list(Path | None)
        List of paths to all author cell set to CL term mapping files
    dataset_version_ids: list(str)
        List of the dataset version identifier corresponding to the
        dataset used to generate the NSForest results, or the first
        such identifier if more than one dataset was combined
    gene_names: list(str)
        List of gene names found in all NSForest results
    gene_ensembl_ids: list(str)
        List of Ensembl identifiers coooresponding to the gene names
        found in all NSForest results
    gene_entrez_ids: list(str)
        List of Entrez identifiers coooresponding to the gene names
        found in all NSForest results
    """
    nsforest_paths = []
    silhouette_paths = []
    author_to_cl_paths = []
    dataset_version_ids = []
    cl_terms = set()
    gene_names = set()
    gene_ensembl_ids = set()
    gene_entrez_ids = set()

    print(f"Loading results sources from {RESULTS_SOURCES_PATH}")
    with open(RESULTS_SOURCES_PATH, "r") as fp:
        results_sources = json.load(fp)

    for results_source in results_sources:

        print(f'Finding NSForest results in {results_source["nsforest_dirpath"]}')
        _nsforest_paths = [
            Path(p).resolve()
            for p in glob(
                str(
                    Path(results_source["nsforest_dirpath"])
                    / results_source["nsforest_pattern"]
                )
            )
        ]
        nsforest_paths.extend(_nsforest_paths)

        # Collect paths to silhouette scores and author cell set to CL
        # term mappings, dataset version identifiers and unique gene
        # names considering all NSForest results
        for _nsforest_path in _nsforest_paths:
            print(
                f"Finding silhouette scores and author cell set to CL term mapping paths, and dataset version id for {_nsforest_path}"
            )
            silhouette_path = None
            author_to_cl_path = None
            dataset_version_id = None
            match = re.search(
                results_source["identity_pattern"], str(_nsforest_path.name)
            )
            if match:
                # TODO: Put this information in the results_source.json file?
                if len(match.groups()) == 2:
                    tissue = ""
                    author = match.group(1)
                    year = match.group(2)

                elif len(match.groups()) == 3:
                    tissue = match.group(1)
                    author = match.group(2)
                    year = match.group(3)

                else:
                    raise Exception("Unexpected number of groups")

                # Find path to silhouette scores
                silhouette_path = glob(
                    str(
                        Path(results_source["silhouette_dirpath"])
                        / results_source["silhouette_pattern"].format(
                            tissue=tissue, author=author, year=year
                        )
                    )
                )
                if len(silhouette_path) != 0:
                    silhouette_path = Path(silhouette_path[0]).resolve()

                # Find path to author cell set to CL term mappings
                author_to_cl_path = glob(
                    str(
                        _nsforest_path.parent
                        / results_source["mapping_pattern"].format(
                            tissue=tissue, author=author, year=year
                        )
                    )
                )
                if len(author_to_cl_path) != 0:
                    author_to_cl_path = Path(author_to_cl_path[0]).resolve()

                    # Load author cell set to CL term mapping, and assign dataset
                    # version identifier
                    print(
                        f"Loading author cell set to CL term mapping from {author_to_cl_path}"
                    )
                    author_to_cl_results = load_results(author_to_cl_path)
                    dataset_version_id = author_to_cl_results["dataset_version_id"][
                        0
                    ].split("--")[
                        0
                    ]  # Some dataset_version_ids are a "--" separated list of dataset_version_ids
                    cl_terms = cl_terms.union(get_cl_terms(author_to_cl_results))

                else:

                    # Parse dataset identifier within NSForest results
                    # path, and assign
                    match = re.search(
                        "([0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12})",
                        str(_nsforest_path.name),
                    )
                    dataset_version_id = match.group(1)

            silhouette_paths.append(silhouette_path)
            author_to_cl_paths.append(author_to_cl_path)
            dataset_version_ids.append(dataset_version_id)

            # Collect unique gene names
            print(f"Loading NSForest results from {_nsforest_path}")
            nsforest_results = load_results(_nsforest_path).sort_values(
                "clusterName", ignore_index=True
            )
            gene_names |= set(collect_unique_gene_names(nsforest_results))

    # Collect unique gene identifiers
    gene_ensembl_ids = collect_unique_gene_ensembl_ids(gene_names)
    gene_entrez_ids = collect_unique_gene_entrez_ids(gene_names)

    return (
        nsforest_paths,
        silhouette_paths,
        author_to_cl_paths,
        dataset_version_ids,
        cl_terms,
        gene_names,
        gene_ensembl_ids,
        gene_entrez_ids,
    )


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


def get_gene_names_and_ensembl_and_entrez_ids():
    """Query BioMart to get gene names, and Ensembl and Entrez ids.

    Parameters
    ----------
    None

    Returns
    -------
    gene_names_and_ids : pd.DataFrame
        DataFrame with columns containing gene names, and Ensembl and
        Entrez ids
    """
    print("Getting gene names, and Ensembl and Entrez ids from BioMart")
    gene_names_and_ids = (
        sc.queries.biomart_annotations(
            "hsapiens",
            ["external_gene_name", "ensembl_gene_id", "entrezgene_id"],
            use_cache=True,
        )
        .dropna()
        .drop_duplicates()
    )
    gene_names_and_ids["entrezgene_id"] = (
        gene_names_and_ids["entrezgene_id"].astype(int).astype(str)
    )
    return gene_names_and_ids


def get_gene_name_to_ensembl_ids_map():
    """Get gene name to Ensembl ids map.

    Parameters
    ----------
    None

    Returns
    -------
    gene_name_to_ensembl_ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Ensembl id
    """
    print("Creating gene name to Ensembl ids map")
    gene_names_and_ids = get_gene_names_and_ensembl_and_entrez_ids()
    gene_name_to_ensembl_ids = gene_names_and_ids.set_index("external_gene_name")
    return gene_name_to_ensembl_ids


def map_gene_name_to_ensembl_ids(name, gene_name_to_ensembl_ids):
    """Map a gene name to a gene Ensembl id list.

    Parameters
    ----------
    name : str
        Gene name
    gene_name_to_ensembl_ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Ensembl id

    Returns
    -------
    list
        Gene Ensembl ids
    """
    if name in gene_name_to_ensembl_ids.index:
        ids = gene_name_to_ensembl_ids.loc[name, "ensembl_gene_id"]
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
    gene_ensembl_id_to_names : pd.DataFrame
        DataFrame indexed by gene Ensembl ids containing gene names
    """
    print("Creating gene Ensembl id to names map")
    gene_names_and_ids = get_gene_names_and_ensembl_and_entrez_ids()
    gene_ensembl_id_to_names = gene_names_and_ids.set_index("ensembl_gene_id")
    return gene_ensembl_id_to_names


def map_gene_ensembl_id_to_names(gid, gene_ensembl_id_to_names):
    """Map a gene Ensembl id to a gene name list.

    Parameters
    ----------
    gid : str
        Gene Ensembl id
    gene_ensembl_id_to_names : pd.DataFrame
        DataFrame indexed by gene Ensembl id containing gene name

    Returns
    -------
    list
        Gene names
    """
    if gid in gene_ensembl_id_to_names.index:
        names = gene_ensembl_id_to_names.loc[gid, "external_gene_name"]
        if isinstance(names, pd.core.series.Series):
            names = names.to_list()
        else:
            names = [names]
        # print(f"Mapped gene Ensembl id {gid} to names {names}")
    else:
        print(f"Could not find gene names for gene Ensembl id: {gid}")
        names = []
    return names


def get_gene_name_to_entrez_ids_map():
    """Get gene name to Entrez ids map.

    Parameters
    ----------
    None

    Returns
    -------
    gene_name_to_entrez_ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Entrez id
    """
    print("Creating gene name to Entrez ids map")
    gene_names_and_ids = get_gene_names_and_ensembl_and_entrez_ids()
    gene_name_to_entrez_ids = gene_names_and_ids.set_index("external_gene_name")
    return gene_name_to_entrez_ids


def map_gene_name_to_entrez_ids(name, gene_name_to_entrez_ids):
    """Map a gene name to a gene Entrez id list.

    Parameters
    ----------
    name : str
        Gene name
    gene_name_to_entrez_ids : pd.DataFrame
        DataFrame indexed by gene name containing gene Entrez id

    Returns
    -------
    list
        Gene Entrez ids
    """
    if name in gene_name_to_entrez_ids.index:
        ids = gene_name_to_entrez_ids.loc[name, "entrezgene_id"]
        if isinstance(ids, pd.core.series.Series):
            ids = ids.to_list()
        else:
            ids = [ids]
        # print(f"Mapped gene name {name} to Entrez ids {ids}")
    else:
        print(f"Could not find gene Entrez ids for gene name: {name}")
        ids = []
    return ids


def get_gene_entrez_id_to_names_map():
    """Map gene Entrez id to names.

    Parameters
    ----------
    None

    Returns
    -------
    gene_entrez_id_to_names : pd.DataFrame
        DataFrame indexed by gene Entrez ids containing gene names
    """
    print("Creating gene Entrez id to names map")
    gene_names_and_ids = get_gene_names_and_ensembl_and_entrez_ids()
    gene_entrez_id_to_names = gene_names_and_ids.set_index("entrezgene_id")
    return gene_entrez_id_to_names


def map_gene_entrez_id_to_names(gid, gene_entrez_id_to_names):
    """Map a gene Entrez id to a gene name list.

    Parameters
    ----------
    gid : str
        Gene Entrez id
    gene_entrez_id_to_names : pd.DataFrame
        DataFrame indexed by gene Entrez id containing gene name

    Returns
    -------
    list
        Gene names
    """
    if gid in gene_entrez_id_to_names.index:
        names = gene_entrez_id_to_names.loc[gid, "external_gene_name"]
        if isinstance(names, pd.core.series.Series):
            names = names.to_list()
        else:
            names = [names]
        # print(f"Mapped gene Entrez id {gid} to names {names}")
    else:
        print(f"Could not find gene names for gene Entrez id: {gid}")
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


def collect_unique_gene_names(nsforest_results):
    """Collect unique gene names found in the NSForest results marker
    or binary genes. Exclude clusters smaller than the minimum
    size. Return these values as a sorted list for restarting.

    Parameters
    ----------
    nsforest_results : pd.DataFrame
        DataFrame containing NSForest results

    Returns
    -------
    gene_names : list(str)
        List of unique gene names
    """
    gene_names = set()

    for column in ["NSForest_markers", "binary_genes"]:
        for gene_list_str in nsforest_results.loc[
            nsforest_results["clusterSize"] >= MIN_CLUSTER_SIZE, column
        ]:
            gene_names |= set(ast.literal_eval(gene_list_str))

    return sorted(gene_names)


def collect_unique_gene_ensembl_ids(gene_names):
    """Collect unique Ensembl gene ids corresponding to the specified
    list of gene names. Return these values as a sorted list for
    restarting.

    Note that if gene names are taken from NSForest results, gene
    names might actually be Ensembl ids.

    Parameters
    ----------
    gene_names : list(str)
        List of gene names

    Returns
    -------
    gene_ensembl_ids : list(str)
        List of unique gene Ensembl ids
    """
    gene_ensembl_ids = set()

    gene_names = set(gene_names)
    gene_name_to_ensembl_ids = get_gene_name_to_ensembl_ids_map()
    for gene_name in gene_names:
        if "ENSG" in gene_name:
            gene_ensembl_id = gene_name.split(".")[0]
        else:
            _gene_ensembl_ids = map_gene_name_to_ensembl_ids(
                gene_name, gene_name_to_ensembl_ids
            )
            if len(_gene_ensembl_ids) == 0:
                gene_ensembl_id = None
            else:
                gene_ensembl_id = _gene_ensembl_ids[0]
        if gene_ensembl_id:
            gene_ensembl_ids.add(gene_ensembl_id)
    print(
        f"Collected {len(gene_ensembl_ids)} unique Ensembl gene ids for {len(gene_names)} unique gene names"
    )

    return sorted(gene_ensembl_ids)


def collect_unique_gene_entrez_ids(gene_names):
    """Collect unique Entrez gene ids corresponding to the specified
    list of gene names. Return these values as a sorted list for
    restarting.

    Note that if gene names are taken from NSForest results, gene
    names might actually be Ensembl ids.

    Parameters
    ----------
    gene_names : list(str)
        List of gene names

    Returns
    -------
    gene_ids : list(str)
        List of unique gene Entrez ids
    """
    gene_entrez_ids = set()

    gene_names = set(gene_names)
    gene_ensembl_id_to_names = get_gene_ensembl_id_to_names_map()
    gene_name_to_entrez_ids = get_gene_name_to_entrez_ids_map()
    for gene_name in gene_names:
        if "ENSG" in gene_name:
            gene_ensembl_id = gene_name.split(".")[0]
            _gene_names = map_gene_ensembl_id_to_names(
                gene_ensembl_id, gene_ensembl_id_to_names
            )
            if len(_gene_names) == 0:
                gene_name = None
            else:
                gene_name = _gene_names[0]
        _gene_entrez_ids = map_gene_name_to_entrez_ids(
            gene_name, gene_name_to_entrez_ids
        )
        if len(_gene_entrez_ids) == 0:
            gene_entrez_id = None
        else:
            gene_entrez_id = _gene_entrez_ids[0]
        if gene_entrez_id:
            gene_entrez_ids.add(gene_entrez_id)
    print(
        f"Collected {len(gene_entrez_ids)} unique Entrez gene ids for {len(gene_names)} unique gene names"
    )

    return sorted(gene_entrez_ids)


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
        "../../../cell-kn-mvp-etl-ontologies/data/mondo_efo_mappings.csv"
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


def get_chembl_to_pubchem_map():
    """Get ChEMBL to PubChem id map.

    Parameters
    ----------
    None

    Returns
    -------
    chembl2pubchem : pd.DataFrame
        DataFrame indexed by ChEMBL id containing PubChem id
    """
    print("Creating ChEMBL to PubChem id map")
    src1src22_name = "../../../cell-kn-mvp-etl-ontologies/data/src1src22.csv"
    chembl2pubchem = pd.read_csv(src1src22_name)
    chembl2pubchem = chembl2pubchem.set_index("ChEMBL")
    return chembl2pubchem


def map_chembl_to_pubchem(chembl, chembl2pubchem):
    """Map ChEMBL to PubChem id.

    Parameters
    ----------
    chembl : str
        ChEMLB id
    chembl2pubchem : pd.DataFrame
        DataFrame indexed by ChEMBL containing PubChem id

    Returns
    -------
    str
        PubChem id
    """
    pubchem = None

    if chembl in chembl2pubchem.index:
        pubchem = chembl2pubchem.loc[chembl, "PubChem"]
        if isinstance(pubchem, pd.core.series.Series):
            pubchem = pubchem.iloc[0]

    return pubchem


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
