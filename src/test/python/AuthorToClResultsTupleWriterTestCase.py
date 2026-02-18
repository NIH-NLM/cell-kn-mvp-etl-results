import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path("../../main/python").resolve()))

import pandas as pd

from AuthorToClResultsTupleWriter import create_tuples_from_author_to_cl

SUMMARIES_DIRPATH = Path("../../test/data/summaries")


def to_string_tuples(tuples):
    """Convert list of URIRef/Literal tuples to list of string lists."""
    return [list(str(x) for x in t) for t in tuples]


class CreateTuplesFromAuthorToClTestCase(unittest.TestCase):
    """Tests for create_tuples_from_author_to_cl using summary fixtures."""

    def setUp(self):
        # Load author-to-cl summary
        author_to_cl_path = (
            SUMMARIES_DIRPATH / "cell-kn-mvp-map-author-to-cl-li-2023.json"
        )
        with open(author_to_cl_path, "r") as fp:
            self.author_to_cl_summary = json.load(fp)
        self.results_df = pd.DataFrame(self.author_to_cl_summary["results"])
        self.expected_tuples = self.author_to_cl_summary["tuples"]

        # Load cellxgene results from external-api summary
        external_api_path = SUMMARIES_DIRPATH / "cell-kn-mvp-external-api-results.json"
        with open(external_api_path, "r") as fp:
            external_api_summary = json.load(fp)
        self.cellxgene_results = external_api_summary["results"]["cellxgene"]

        # PMID metadata extracted from expected tuples
        self.pmid_data = {
            "Author": "Li et al.",
            "Journal": "Res Sq",
            "Title": "Integrated multi-omics single cell atlas of the human retina.",
            "Year": "2023",
            "Citation": "Li et al. (2023) Res Sq",
        }

    @patch("AuthorToClResultsTupleWriter.get_data_for_pmid")
    def test_create_tuples_from_author_to_cl(self, mock_get_data):
        """Tuples created from summary results match expected tuples."""
        mock_get_data.return_value = self.pmid_data
        actual_tuples = create_tuples_from_author_to_cl(
            self.results_df, self.cellxgene_results
        )
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    @patch("AuthorToClResultsTupleWriter.get_data_for_pmid")
    def test_tuple_count(self, mock_get_data):
        """Number of tuples matches expected count."""
        mock_get_data.return_value = self.pmid_data
        actual_tuples = create_tuples_from_author_to_cl(
            self.results_df, self.cellxgene_results
        )
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    @patch("AuthorToClResultsTupleWriter.get_data_for_pmid")
    def test_first_tuple_is_csd_citation(self, mock_get_data):
        """First tuple is a CSD Citation annotation."""
        mock_get_data.return_value = self.pmid_data
        actual_tuples = create_tuples_from_author_to_cl(
            self.results_df, self.cellxgene_results
        )
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("CSD_", first[0])
        self.assertIn("Citation", first[1])
        self.assertEqual(first[2], "Li et al. (2023) Res Sq")

    @patch("AuthorToClResultsTupleWriter.get_data_for_pmid")
    def test_last_tuple_is_gene_source(self, mock_get_data):
        """Last tuple is a Gene PART_OF Cell type Source edge annotation."""
        mock_get_data.return_value = self.pmid_data
        actual_tuples = create_tuples_from_author_to_cl(
            self.results_df, self.cellxgene_results
        )
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("GS_ITIH5", last[0])
        self.assertIn("CL_4030027", last[1])
        self.assertIn("Source", last[2])
        self.assertEqual(last[3], "NSForest")
