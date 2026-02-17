import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path("../../main/python").resolve()))

import pandas as pd

from NSForestResultsTupleWriter import create_tuples_from_nsforest

SUMMARIES_DIRPATH = Path("../../test/data/summaries")


def to_string_tuples(tuples):
    """Convert list of URIRef/Literal tuples to list of string lists."""
    return [list(str(x) for x in t) for t in tuples]


class CreateTuplesFromNSForestTestCase(unittest.TestCase):
    """Tests for create_tuples_from_nsforest using summary fixture."""

    def setUp(self):
        summary_path = SUMMARIES_DIRPATH / "cell-kn-mvp-nsforest-results-li-2023.json"
        with open(summary_path, "r") as fp:
            self.summary = json.load(fp)
        self.results_df = pd.DataFrame(self.summary["results"])
        self.expected_tuples = self.summary["tuples"]

    def test_create_tuples_from_nsforest(self):
        """Tuples created from summary results match expected tuples."""
        actual_tuples = create_tuples_from_nsforest(self.results_df)
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    def test_tuple_count(self):
        """Number of tuples matches expected count."""
        actual_tuples = create_tuples_from_nsforest(self.results_df)
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    def test_contains_bmc_type_tuple(self):
        """First tuple is the BMC INSTANCE_OF triple."""
        actual_tuples = create_tuples_from_nsforest(self.results_df)
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("BMC_", first[0])
        self.assertIn("rdf#type", first[1])
        self.assertIn("SO_0001260", first[2])

    def test_contains_marker_count_tuple(self):
        """Last tuple is the Marker_count edge annotation."""
        actual_tuples = create_tuples_from_nsforest(self.results_df)
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("CS_", last[0])
        self.assertIn("BMC_", last[1])
        self.assertIn("Marker_count", last[2])
        self.assertEqual(last[3], "2")
