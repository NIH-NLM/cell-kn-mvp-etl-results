import json
import os
from time import sleep
from urllib import parse

from bs4 import BeautifulSoup
import requests

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1


def find_names_or_none(soup, names):
    """Find the text in the last named tag, if all previously named
    tags are found.

    Parameters
    ----------
    soup : bs4.element.Tag
        Any soup returned by BeautifulSoup
    names : list(str)
        List of tag names to find in order

    Returns
    -------
    str
        text in the last named tag, or None
    """
    soup = soup.find(names[0])
    for name in names[1:]:
        if soup:
            soup = soup.find(name)
    if soup:
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
                fp.write(BeautifulSoup(xml_data, "xml").prettify())

        # Got the page, so parse it, and search for the title
        soup = BeautifulSoup(xml_data, "xml").find("Article")
        if soup:
            data["author"] = find_names_or_none(
                soup, ["AuthorList", "Author", "LastName"]
            )  # First author
            if len(find_names_or_none(soup, ["AuthorList"])) > 1:
                data["author"] += " et al."
            data["journal"] = find_names_or_none(soup, ["Journal", "ISOAbbreviation"])
            data["title"] = find_names_or_none(soup, ["ArticleTitle"])
            data["year"] = find_names_or_none(soup, ["ArticleDate", "Year"])
            data["citation"] = f"{data['author']} ({data['year']}) {data['journal']}"
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
        gene_id = json_data["esearchresult"]["idlist"][0]

    else:
        print(f"Encountered error in searching Gene: {response.status_code}")

    return gene_id


def get_data_for_gene_id(gene_id, is_summary=False, do_write=False):
    """Fetch from Gene using a gene id to get the full, or summary,
    record.

    Parameters
    ----------
    gene_id : str
        The Gene identifier to use in the fetch
    is_summary : bool
        Flag to fetch summary, rather than full, data (default: False)
    do_write : bool
        Flag to write fetched results, or not (default: False)

    Returns
    -------
    data : dict
       Dictionary containing the full record
    """
    # Need a default return value
    data = {}

    # Fetch from Gene
    print(f"Getting data for gene id: '{gene_id}'")
    if is_summary:
        fetch_url = EUTILS_URL + "esummary.fcgi"
    else:
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
                fp.write(BeautifulSoup(xml_data, "xml").prettify())

        # Got the page, so parse it, and search for the title
        breakpoint()
        soup = BeautifulSoup(xml_data, "xml").find("Article")
        if soup:
            data["author"] = find_names_or_none(
                soup, ["AuthorList", "Author", "LastName"]
            )  # First author
            if len(find_names_or_none(soup, ["AuthorList"])) > 1:
                data["author"] += " et al."
            data["journal"] = find_names_or_none(soup, ["Journal", "ISOAbbreviation"])
            data["title"] = find_names_or_none(soup, ["ArticleTitle"])
            data["year"] = find_names_or_none(soup, ["ArticleDate", "Year"])
            data["citation"] = f"{data['author']} ({data['year']}) {data['journal']}"
    else:
        print(f"Encountered error in fetching from Gene: {response.status_code}")

    return data


def main():
    print(get_data_for_pmid("37291214"))


if __name__ == "__main__":
    main()
