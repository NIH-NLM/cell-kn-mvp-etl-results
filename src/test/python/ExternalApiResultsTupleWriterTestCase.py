from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path("../../main/python").resolve()))

import pandas as pd

from ExternalApiResultsTupleWriter import get_mondo_term, get_protein_term, remove_protocols


class RemoveProtocolsTestCase(unittest.TestCase):
    """Tests for remove_protocols."""

    def test_removes_https(self):
        self.assertEqual(remove_protocols("https://example.com"), "example.com")

    def test_removes_http(self):
        self.assertEqual(remove_protocols("http://example.com"), "example.com")

    def test_non_string_unchanged(self):
        """Non-string values pass through unchanged."""
        self.assertEqual(remove_protocols(42), 42)
        self.assertIsNone(remove_protocols(None))
        self.assertEqual(remove_protocols([1, 2]), [1, 2])

    def test_string_without_protocol(self):
        """String with no protocol is returned as-is."""
        self.assertEqual(remove_protocols("example.com"), "example.com")


class GetMondoTermTestCase(unittest.TestCase):
    """Tests for get_mondo_term."""

    def setUp(self):
        self.efo2mondo = pd.DataFrame(
            {"EFO": ["EFO_0000270"], "MONDO": ["MONDO_0004992"]}
        ).set_index("EFO")

    def test_mondo_term_passes_through(self):
        """MONDO term is returned directly."""
        result = get_mondo_term("MONDO_0004992", self.efo2mondo)
        self.assertEqual(result, "MONDO_0004992")

    def test_efo_term_maps_to_mondo(self):
        """EFO term is mapped to MONDO via efo2mondo."""
        result = get_mondo_term("EFO_0000270", self.efo2mondo)
        self.assertEqual(result, "MONDO_0004992")

    def test_efo_term_not_found(self):
        """Unknown EFO term returns None."""
        result = get_mondo_term("EFO_9999999", self.efo2mondo)
        self.assertIsNone(result)

    def test_non_mondo_non_efo_returns_none(self):
        """Term that is neither MONDO nor EFO returns None."""
        result = get_mondo_term("DOID_1234", self.efo2mondo)
        self.assertIsNone(result)

    @patch("ExternalApiResultsTupleWriter.DEPRECATED_TERMS", ["MONDO_0004992"])
    def test_deprecated_mondo_returns_none(self):
        """Deprecated MONDO term returns None."""
        result = get_mondo_term("MONDO_0004992", self.efo2mondo)
        self.assertIsNone(result)

    @patch("ExternalApiResultsTupleWriter.DEPRECATED_TERMS", ["MONDO_0004992"])
    def test_efo_mapping_to_deprecated_mondo_returns_none(self):
        """EFO mapped to deprecated MONDO term returns None."""
        result = get_mondo_term("EFO_0000270", self.efo2mondo)
        self.assertIsNone(result)


class GetProteinTermTestCase(unittest.TestCase):
    """Tests for get_protein_term."""

    def setUp(self):
        self.ensp2accn = {
            "ENSP00000269305": "P04637",
            "ENSP00000000001": ["P11111", "P22222"],
        }

    def test_ensp_maps_to_pr_term(self):
        """ENSP id maps through ensp2accn to PR_accession."""
        result = get_protein_term("ENSP00000269305", self.ensp2accn)
        self.assertEqual(result, "PR_P04637")

    def test_ensp_with_list_uses_first(self):
        """ENSP with multiple accessions uses first."""
        result = get_protein_term("ENSP00000000001", self.ensp2accn)
        self.assertEqual(result, "PR_P11111")

    def test_non_ensp_passes_through_as_accession(self):
        """Non-ENSP id is used directly as accession."""
        result = get_protein_term("P04637", self.ensp2accn)
        self.assertEqual(result, "PR_P04637")

    def test_ensp_not_in_map_returns_none(self):
        """Unknown ENSP id returns None."""
        result = get_protein_term("ENSP99999999999", self.ensp2accn)
        self.assertIsNone(result)
