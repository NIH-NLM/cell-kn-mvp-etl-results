"""Source from https://www.uniprot.org/help/id_mapping, comments added"""

import json
import re
import time
from urllib.parse import urlparse, parse_qs, urlencode
from xml.etree import ElementTree
import zlib

import requests
from requests.adapters import HTTPAdapter, Retry


POLLING_INTERVAL = 1.5
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    """Check requests response for HTTP errors.

    Parameters
    ----------
    response : requests.Response
        Response object from request

    Returns
    -------
    None

    Raises
    ------
    HTTPError
        If the request raises an exception
    """
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    """Submit and id mapping job.

    Parameters
    ----------
    from_db : str
        Database from which to map
    to_db : str
        Database to which to map
    ids : list(str)
        Frob database ids

    Returns
    -------
    str
        Job id of the submitted job

    Raises
    ------
    HTTPError
        If the request raises an exception
    """
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    """Get the next link in the response headers.

    Parameters
    ----------
    headers : dict
        Case-insensitive Dictionary of Response Headers

    Returns
    -------
    str
        The next link
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    """Periodically check the status of a submitted id mapping job
    until completed.

    Parameters
    ----------
    job_id : str
        Job id of a submitted job

    Returns
    -------
    bool
       Flag indicating id mapping job has completed

    Raises
    ------
    HTTPError
        If the request raises an exception
    Exception
        If the the job status is not "NEW" or "RUNNING"
    """
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    """Get and decode batch response, yielding after each batch.

    Parameters
    ----------
    batch_response : requests.Response
        Response object from request
    file_format : str
        One of "json", "tsv", "xlsx", or "xml"
    compressed : bool
        Flag indicating the response is compressed, or not

    Returns
    -------
    dict | list(str) | str
        The decoded response, dict if file_format = "json", list(str)
        if file_format = "tsv", str otherwise
    Raises
    ------
    HTTPError
        If the request raises an exception
    """
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    """Combine batch responses consistent with file_format.

    Parameters
    ----------
    all_results : dict | list(str) | str
        Batch response combined consistent with file_format
    batch_results : dict | list(str) | str
        Another batch response
    file_format : str
        One of "json", "tsv", "xlsx", or "xml"

    Returns
    -------
    all_results : dict | list(str) | str
        Batch response combined consistent with file_format
    """
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    """Get redirect link for an id mapping job.

    Parameters
    ----------
    job_id : str
        Job id of a submitted job

    Returns
    -------
    str
        Redirect link

    Raises
    ------
    HTTPError
        If the request raises an exception
    """
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    """Decompress response, if needed, then decode consisten with
    file_format.

    Parameters
    ----------
    response : requests.Response
        The request response to decode
    file_format : str
        One of "json", "tsv", "xlsx", or "xml"
    compressed : bool
        Flag indicating the response is compressed, or not

    Returns
    -------
    dict | list(str) | str
        The decoded response, dict if file_format = "json", list(str)
        if file_format = "tsv", str otherwise
    """
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    """Get the XML namespace of an element.

    Parameters
    ----------
    element : ElementTree.Element
       The element for which to get the namespace

    Returns
    -------
    str
       The XML namespace
    """
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    """Merge all entry elements.

    Parameters
    ----------
    xml_results : list(str)
        List containing XML result strings

    Returns
    -------
    str
        The merged XML string
    """
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    """Print the progress of the batch request.

    Parameters
    ----------
    batch_index : int
        The number of the batch
    size : int
        The size of the batch
    total : int
        The number of results in the batch response

    Returns
    -------
    None
    """
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    """Get id mapping request search results.

    Parameters
    ----------
    url : str
        The id mapping request URL

    Returns
    -------
    results : dict | list(str) | str
        Batch response combined consistent with file_format

    Raises
    ------
    HTTPError
        If the request raises an exception
    """
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    """Get id mapping request stream results.

    Parameters
    ----------
    url : str
        The id mapping request URL

    Returns
    -------
    dict | list(str) | str
        The decoded response, dict if file_format = "json", list(str)
        if file_format = "tsv", str otherwise

    Raises
    ------
    HTTPError
        If the request raises an exception
    """
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)


def main():
    """Provide a simple example of submitting and id mapping job, then
    waiting for the result, and printing it.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    job_id = submit_id_mapping(
        from_db="Ensembl_Protein", to_db="UniProtKB", ids=["ENSP00000484862"]
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        # Equivalently using the stream endpoint which is more demanding
        # on the API and so is less stable:
        # results = get_id_mapping_results_stream(link)

    print(results)
    # {'results': [{'from': 'P05067', 'to': 'CHEMBL2487'}], 'failedIds': ['P12345']}


if __name__ == "__main__":
    main()
