import sys
import io
import os
from typing import List
from zipfile import ZipFile

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.openapi import GeneApi as DatasetsGeneApi

zipfile_name = "gene_ds.zip"

# download the data package using the DatasetsGeneApi, and then print out protein sequences for A2M and GNAS
with DatasetsApiClient() as api_client:
    gene_ids: List[int] = [2, 2778]
    gene_api = DatasetsGeneApi(api_client)
    try:
        gene_dataset_download = gene_api.download_gene_package_without_preload_content(
            gene_ids,
            include_annotation_type=["FASTA_GENE", "FASTA_PROTEIN"],
        )

        with open(zipfile_name, "wb") as f:
            f.write(gene_dataset_download.data)
    except DatasetsApiException as e:
        sys.exit(f"Exception when calling GeneApi: {e}\n")


try:
    dataset_zip = ZipFile(zipfile_name)
    zinfo = dataset_zip.getinfo(os.path.join("ncbi_dataset/data", "protein.faa"))
    with io.TextIOWrapper(dataset_zip.open(zinfo), encoding="utf8") as fh:
        print(fh.read())
except KeyError as e:
    logger.error("File %s not found in zipfile: %s", file_name, e)
