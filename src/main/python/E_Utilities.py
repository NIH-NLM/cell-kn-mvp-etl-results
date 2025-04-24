import os
from time import sleep
from urllib import parse

from bs4 import BeautifulSoup
import requests


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1
PUBMED = "pubmed"


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


def get_data_for_pmid(pmid):
    """Fetch from PubMed using a PMID to find the last name of the
    first author, journal title, article title, and article year of
    publication.

    Parameters
    ----------
    pmid : str
       The PubMed identifier to use in the fetch

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
        "db": PUBMED,
        "id": pmid,
        "rettype": "xml",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    sleep(NCBI_API_SLEEP)
    response = requests.get(fetch_url, params=parse.urlencode(params, safe=","))
    if response.status_code == 200:
        xml_data = response.text
        with open(f"{pmid}.xml", "w") as fp:
            fp.write(BeautifulSoup(xml_data, "xml").prettify())

        # Got the page, so parse it, and search for the title
        soup = BeautifulSoup(xml_data, "xml").find("Article")
        if soup:
            data["author"] = find_names_or_none(
                soup, ["AuthorList", "Author", "LastName"]
            )
            data["journal"] = find_names_or_none(soup, ["Journal", "Title"])
            data["title"] = find_names_or_none(soup, ["ArticleTitle"])
            data["year"] = find_names_or_none(soup, ["ArticleDate", "Year"])

    else:
        print(f"Encountered error in fetching from PubMed: {response.status_code}")

    return data


def main():
    print(get_data_for_pmid("37291214"))


if __name__ == "__main__":
    main()
