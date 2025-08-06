import json
import os
from pathlib import Path
import re
from time import sleep
from urllib import parse

import bs4
import requests

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1


def find_names_or_none(soup, names, attribute=None):
    """Find the text, or specified attribute, in the last named tag,
    if all previously named tags are found.

    Parameters
    ----------
    soup : bs4.element.Tag
        Any soup returned by BeautifulSoup
    names : list(str)
        List of tag names to find in order
    attribute : str
        Attribute of the last named tag

    Returns
    -------
    str
        text, or attribute, in the last named tag, or None
    """
    soup = soup.find(names[0])
    for name in names[1:]:
        if soup:
            soup = soup.find(name)
    if soup:
        if attribute:
            return soup.get(attribute)
        else:
            return soup.text
    else:
        return soup


def get_data_for_pmid(pmid, do_write=False):
    """Fetch from PubMed using a PMID to find the last name of the
    first author, journal title, article title, and article year of
    publication.

    Parameters
    ----------
    pmid : str
        The PubMed identifier to use in the fetch
    do_write : bool
        Flag to write fetched results, or not (default: False)

    Returns
    -------
    data : dict
       Dictionary containing the last name of the first author,
       journal title, article title, and article year of publication
    """
    # Need a default return value
    data = {}

    # Fetch from PubMed
    print(f"Getting title for PMID: '{pmid}'")
    fetch_url = EUTILS_URL + "efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": pmid,
        "rettype": "xml",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    sleep(NCBI_API_SLEEP)
    response = requests.get(fetch_url, params=parse.urlencode(params, safe=","))
    if response.status_code == 200:
        xml_data = response.text
        if do_write:
            with open(f"{pmid}.xml", "w") as fp:
                fp.write(bs4.BeautifulSoup(xml_data, "xml").prettify())

        # Got the page, so parse it, and search for the title
        root = bs4.BeautifulSoup(xml_data, "xml").find("Article")
        if root:
            data["Author"] = find_names_or_none(
                root, ["AuthorList", "Author", "LastName"]
            )  # First author
            if len(find_names_or_none(root, ["AuthorList"])) > 1:
                data["Author"] += " et al."
            data["Journal"] = find_names_or_none(root, ["Journal", "ISOAbbreviation"])
            data["Title"] = find_names_or_none(root, ["ArticleTitle"])
            data["Year"] = find_names_or_none(root, ["ArticleDate", "Year"])
            data["Citation"] = f"{data['Author']} ({data['Year']}) {data['Journal']}"
    else:
        print(f"Encountered error in fetching from PubMed: {response.status_code}")

    return data


def find_gene_id_for_gene_name(name, do_write=False):
    """Search Gene using a gene name to find the corresponding gene
    id.

    Parameters
    ----------
    name : str
       The gene name for which to search
    do_write : bool
        Flag to write fetched results, or not (default: False)

    Returns
    -------
    str
       The gene id
    """
    # Need a default return value
    gene_id = None

    # Search Gene
    print(f"Searching Gene for name: '{name}'")
    search_url = EUTILS_URL + "esearch.fcgi"
    params = {
        "db": "gene",
        "term": name,
        "sort": "relevance",
        "retmax": 1,
        "retmode": "json",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    sleep(NCBI_API_SLEEP)
    response = requests.get(search_url, params=parse.urlencode(params, safe=","))
    if response.status_code == 200:
        json_data = response.json()
        if do_write:
            with open(f"{name}.json", "w") as fp:
                json.dump(json_data, fp, indent=4)

        # Got the response, so assign the gene id
        if len(json_data["esearchresult"]["idlist"]) > 0:
            gene_id = json_data["esearchresult"]["idlist"][0]
            print(f"Found gene id {gene_id} while searching Gene for name {name}")

        else:
            print(f"No gene id found while searching Gene for name {name}")

    else:
        print(
            f"Encountered error in searching Gene for name {name}: {response.status_code}"
        )

    return gene_id


def get_data_for_gene_id(gene_id, do_write=False):
    """Fetch from Gene using a gene id to get the full record and find
    required values.

    Parameters
    ----------
    gene_id : str
        The Gene identifier to use in the fetch
    do_write : bool
        Flag to write fetched results, or not (default: False)

    Returns
    -------
    data : dict
       Dictionary containing the required values of the full record
    """
    # Need a default return value
    data = {}

    # Fetch from Gene
    print(f"Getting data for gene id: '{gene_id}'")
    fetch_url = EUTILS_URL + "efetch.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "xml",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    sleep(NCBI_API_SLEEP)
    response = requests.get(fetch_url, params=parse.urlencode(params, safe=","))
    if response.status_code == 200:
        xml_data = response.text
        if do_write:
            with open(f"{gene_id}.xml", "w") as fp:
                fp.write(bs4.BeautifulSoup(xml_data, "xml").prettify())

        # Got the page, so parse it, and search for the required values
        tags = bs4.BeautifulSoup(xml_data, "xml").find_all("Entrezgene")
        if len(tags) > 1:
            raise Exception("Expect a single Entrezgene element")
        root = tags[0]
        data["Official_symbol"] = find_names_or_none(
            root,
            [
                "Entrezgene_gene",
                "Gene-ref",
                "Gene-ref_formal-name",
                "Gene-nomenclature_symbol",
            ],
        )
        data["Official_full_name"] = find_names_or_none(
            root,
            [
                "Entrezgene_gene",
                "Gene-ref",
                "Gene-ref_formal-name",
                "Gene-nomenclature_name",
            ],
        )
        data["Gene_type"] = find_names_or_none(
            root, ["Entrezgene_type"], attribute="value"
        )
        for child in root.find_all("Other-source_url"):
            if "www.uniprot.org" in child.text:
                data["Link_to_UniProt_ID"] = child.text
        data["Organism"] = find_names_or_none(
            root,
            [
                "Entrezgene_source",
                "BioSource",
                "BioSource_org",
                "Org-ref",
                "Org-ref_taxname",
            ],
        )
        data["RefSeq_gene_ID"] = None
        for child in root.find_all("Gene-commentary_heading"):
            if "GCF_" in child.text:
                m = re.search(r":\s*(GCF_.*)", child.text)
                if m:
                    data["RefSeq_gene_ID"] = m.group(1)
        data["Also_known_as"] = []
        for child in root.find_all("Gene-ref_syn_E"):
            data["Also_known_as"].append(child.text)
        data["Summary"] = find_names_or_none(root, ["Entrezgene_summary"])
        pr_desc = find_names_or_none(root, ["Entrezgene_prot", "Prot-ref_desc"])
        data["UniProt_name"] = Path(
            parse.urlparse(data["Link_to_UniProt_ID"]).path
        ).stem
        for product in root.find_all("Gene-commentary_products"):
            if find_names_or_none(product, ["Gene-commentary_type"], "value") == "mRNA":
                nm_id = None
                np_id = None
                for accession in product.find_all("Gene-commentary_accession"):
                    if "NM_" in accession.text:
                        nm_id = accession.text
                    elif "NP_" in accession.text:
                        np_id = accession.text
                if nm_id and np_id and pr_desc:
                    data["mRNA_(NM)_and_protein_(NP)_sequences"] = (
                        f"{nm_id} -> {np_id}, {pr_desc}"
                    )
                break

        return data

    else:
        print(f"Encountered error in fetching from Gene: {response.status_code}")

    return data


def main():
    print(get_data_for_pmid("37291214"))
    print(get_data_for_gene_id("1080"))


if __name__ == "__main__":
    main()
