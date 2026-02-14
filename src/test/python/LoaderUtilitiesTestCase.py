import os
from pathlib import Path
import string
import sys
import unittest

sys.path.insert(0, str(Path("../../main/python").resolve()))

import pandas as pd

import LoaderUtilities as lu


class LoaderUtilitiesTestCase(unittest.TestCase):
    """Pure unit tests for LoaderUtilities functions."""

    # -- get_uuid tests --

    def test_get_uuid_length(self):
        """Returns a 12-character string."""
        uuid = lu.get_uuid()
        self.assertEqual(len(uuid), 12)

    def test_get_uuid_characters(self):
        """All characters are lowercase alphanumeric."""
        uuid = lu.get_uuid()
        allowed = set(string.ascii_lowercase + string.digits)
        for c in uuid:
            self.assertIn(c, allowed)

    def test_get_uuid_uniqueness(self):
        """Two calls return different values."""
        self.assertNotEqual(lu.get_uuid(), lu.get_uuid())

    # -- hyphenate tests --

    def test_hyphenate_space(self):
        """Replaces space with hyphen."""
        self.assertEqual(lu.hyphenate("hello world"), "hello-world")

    def test_hyphenate_space(self):
        """Replaces spaces with hyphen."""
        self.assertEqual(lu.hyphenate("hello   world"), "hello-world")

    def test_hyphenate_underscore(self):
        """Replaces underscore with hyphen."""
        self.assertEqual(lu.hyphenate("hello_world"), "hello-world")

    def test_hyphenate_underscore(self):
        """Replaces underscores with hyphen."""
        self.assertEqual(lu.hyphenate("hello___world"), "hello-world")

    def test_hyphenate_comma(self):
        """Replaces commas with hyphen."""
        self.assertEqual(lu.hyphenate("hello,world"), "hello-world")

    def test_hyphenate_comma(self):
        """Replaces commas with hyphens."""
        self.assertEqual(lu.hyphenate("hello,,,world"), "hello-world")

    def test_hyphenate_slash(self):
        """Replaces forward slashe with hyphen."""
        self.assertEqual(lu.hyphenate("hello/world"), "hello-world")

    def test_hyphenate_slash(self):
        """Replaces forward slashes with hyphen."""
        self.assertEqual(lu.hyphenate("hello///world"), "hello-world")

    def test_hyphenate_multiple_separators(self):
        """Handles multiple different separators."""
        self.assertEqual(lu.hyphenate("a b_c/d"), "a-b-c-d")

    # -- get_value_or_none tests --

    def test_get_value_or_none_nested(self):
        """Accesses nested dict value."""
        data = {"a": {"b": {"c": 42}}}
        self.assertEqual(lu.get_value_or_none(data, ["a", "b", "c"]), 42)

    def test_get_value_or_none_partial(self):
        """Returns intermediate dict."""
        data = {"a": {"b": {"c": 42}}}
        self.assertEqual(lu.get_value_or_none(data, ["a", "b"]), {"c": 42})

    def test_get_value_or_none_missing_key(self):
        """Returns None for missing key."""
        data = {"a": {"b": {"c": 42}}}
        self.assertIsNone(lu.get_value_or_none(data, ["a", "x"]))

    def test_get_value_or_none_missing_first_key(self):
        """Returns None for missing first key."""
        data = {"a": 1}
        self.assertIsNone(lu.get_value_or_none(data, ["x"]))

    def test_get_value_or_none_empty_dict(self):
        """Returns None for empty dict."""
        self.assertIsNone(lu.get_value_or_none({}, ["a"]))

    # -- get_values_or_none tests --

    def test_get_values_or_none_collects(self):
        """Collects comma-separated values from list items."""
        data = {"items": [{"name": "Alice"}, {"name": "Bob"}]}
        self.assertEqual(lu.get_values_or_none(data, "items", ["name"]), "Alice, Bob")

    def test_get_values_or_none_missing_key(self):
        """Returns empty string for missing list key."""
        data = {"items": [{"name": "Alice"}]}
        self.assertEqual(lu.get_values_or_none(data, "missing", ["name"]), "")

    def test_get_values_or_none_single_item(self):
        """Single item returns just that value."""
        data = {"items": [{"name": "Alice"}]}
        self.assertEqual(lu.get_values_or_none(data, "items", ["name"]), "Alice")

    # -- map_gene_name_to_ensembl_ids tests --

    def test_map_gene_name_to_ensembl_ids_single(self):
        """Maps gene name to single Ensembl id."""

        df = pd.DataFrame(
            {
                "external_gene_name": ["BRCA1"],
                "ensembl_gene_id": ["ENSG00000012048"],
            }
        ).set_index("external_gene_name")
        ids = lu.map_gene_name_to_ensembl_ids("BRCA1", df)
        self.assertEqual(ids, ["ENSG00000012048"])

    def test_map_gene_name_to_ensembl_ids_multiple(self):
        """Maps gene name to multiple Ensembl ids."""

        df = pd.DataFrame(
            {
                "external_gene_name": ["TP53", "TP53"],
                "ensembl_gene_id": ["ENSG00000141510", "ENSG00000999999"],
            }
        ).set_index("external_gene_name")
        ids = lu.map_gene_name_to_ensembl_ids("TP53", df)
        self.assertIsInstance(ids, list)
        self.assertIn("ENSG00000141510", ids)
        self.assertIn("ENSG00000999999", ids)

    def test_map_gene_name_to_ensembl_ids_missing(self):
        """Returns empty list for unknown gene name."""

        df = pd.DataFrame(
            {
                "external_gene_name": ["BRCA1"],
                "ensembl_gene_id": ["ENSG00000012048"],
            }
        ).set_index("external_gene_name")
        ids = lu.map_gene_name_to_ensembl_ids("NONEXISTENT", df)
        self.assertEqual(ids, [])

    # -- map_gene_ensembl_id_to_names tests --

    def test_map_gene_ensembl_id_to_names(self):
        """Maps Ensembl id to single gene name."""

        df = pd.DataFrame(
            {
                "ensembl_gene_id": ["ENSG00000012048"],
                "external_gene_name": ["BRCA1"],
            }
        ).set_index("ensembl_gene_id")
        names = lu.map_gene_ensembl_id_to_names("ENSG00000012048", df)
        self.assertEqual(names, ["BRCA1"])

    def test_map_gene_ensembl_id_to_names(self):
        """Maps Ensembl id to multiple gene names."""

        df = pd.DataFrame(
            {
                "ensembl_gene_id": ["ENSG00000012048", "ENSG00000012048"],
                "external_gene_name": ["BRCA1", "BRCA9"],
            }
        ).set_index("ensembl_gene_id")
        names = lu.map_gene_ensembl_id_to_names("ENSG00000012048", df)
        self.assertIsInstance(names, list)
        self.assertIn("BRCA1", names)
        self.assertIn("BRCA9", names)

    def test_map_gene_ensembl_id_to_names_missing(self):
        """Returns empty list for unknown Ensembl id."""

        df = pd.DataFrame(
            {
                "ensembl_gene_id": ["ENSG00000012048"],
                "external_gene_name": ["BRCA1"],
            }
        ).set_index("ensembl_gene_id")
        names = lu.map_gene_ensembl_id_to_names("ENSG99999999999", df)
        self.assertEqual(names, [])

    # -- map_gene_name_to_entrez_ids tests --

    def test_map_gene_name_to_entrez_ids(self):
        """Maps gene name to single Entrez id."""

        df = pd.DataFrame(
            {
                "external_gene_name": ["BRCA1"],
                "entrezgene_id": ["672"],
            }
        ).set_index("external_gene_name")
        ids = lu.map_gene_name_to_entrez_ids("BRCA1", df)
        self.assertEqual(ids, ["672"])

    def test_map_gene_name_to_entrez_ids(self):
        """Maps gene name to multiple Entrez ids."""

        df = pd.DataFrame(
            {
                "external_gene_name": ["BRCA1", "BRCA1"],
                "entrezgene_id": ["672", "999"],
            }
        ).set_index("external_gene_name")
        ids = lu.map_gene_name_to_entrez_ids("BRCA1", df)
        self.assertIsInstance(ids, list)
        self.assertIn("672", ids)
        self.assertIn("999", ids)

    def test_map_gene_name_to_entrez_ids_missing(self):
        """Returns empty list for unknown gene name."""

        df = pd.DataFrame(
            {
                "external_gene_name": ["BRCA1"],
                "entrezgene_id": ["672"],
            }
        ).set_index("external_gene_name")
        ids = lu.map_gene_name_to_entrez_ids("NONEXISTENT", df)
        self.assertEqual(ids, [])

    # -- map_gene_entrez_id_to_names tests --

    def test_map_gene_entrez_id_to_names(self):
        """Maps Entrez id to single gene name."""

        df = pd.DataFrame(
            {
                "entrezgene_id": ["672"],
                "external_gene_name": ["BRCA1"],
            }
        ).set_index("entrezgene_id")
        names = lu.map_gene_entrez_id_to_names("672", df)
        self.assertEqual(names, ["BRCA1"])

    def test_map_gene_entrez_id_to_names(self):
        """Maps Entrez id to multiple gene names."""

        df = pd.DataFrame(
            {
                "entrezgene_id": ["672", "672"],
                "external_gene_name": ["BRCA1", "BRCA9"],
            }
        ).set_index("entrezgene_id")
        names = lu.map_gene_entrez_id_to_names("672", df)
        self.assertIsInstance(names, list)
        self.assertIn("BRCA1", names)
        self.assertIn("BRCA9", names)

    def test_map_gene_entrez_id_to_names_missing(self):
        """Returns empty list for unknown Entrez id."""

        df = pd.DataFrame(
            {
                "entrezgene_id": ["672"],
                "external_gene_name": ["BRCA1"],
            }
        ).set_index("entrezgene_id")
        names = lu.map_gene_entrez_id_to_names("99999", df)
        self.assertEqual(names, [])

    # -- map_protein_ensembl_id_to_accession tests --

    def test_map_protein_ensembl_id_to_accession_single(self):
        """Maps Ensembl protein id to single accession."""
        ensp2accn = {"ENSP001": "P12345"}
        self.assertEqual(
            lu.map_protein_ensembl_id_to_accession("ENSP001", ensp2accn), "P12345"
        )

    def test_map_protein_ensembl_id_to_accession_list(self):
        """Maps Ensembl protein id to first of multiple accessions."""
        ensp2accn = {"ENSP001": ["P11111", "P22222"]}
        self.assertEqual(
            lu.map_protein_ensembl_id_to_accession("ENSP001", ensp2accn), "P11111"
        )

    def test_map_protein_ensembl_id_to_accession_missing(self):
        """Returns None for unknown protein id."""
        ensp2accn = {"ENSP001": "P12345"}
        self.assertIsNone(lu.map_protein_ensembl_id_to_accession("ENSP999", ensp2accn))

    # -- map_accession_to_protein_ensembl_id tests --

    def test_map_accession_to_protein_ensembl_id_single(self):
        """Maps accession to single Ensembl protein id."""
        accn2ensp = {"P12345": "ENSP001"}
        self.assertEqual(
            lu.map_accession_to_protein_ensembl_id("P12345", accn2ensp), "ENSP001"
        )

    def test_map_accession_to_protein_ensembl_id_list(self):
        """Maps accession to first of multiple Ensembl protein ids."""
        accn2ensp = {"P12345": ["ENSP001", "ENSP002"]}
        self.assertEqual(
            lu.map_accession_to_protein_ensembl_id("P12345", accn2ensp), "ENSP001"
        )

    def test_map_accession_to_protein_ensembl_id_missing(self):
        """Returns None for unknown accession."""
        accn2ensp = {"P12345": "ENSP001"}
        self.assertIsNone(lu.map_accession_to_protein_ensembl_id("P99999", accn2ensp))

    # -- map_efo_to_mondo tests --

    def test_map_efo_to_mondo(self):
        """Maps EFO term to MONDO term."""

        efo2mondo = pd.DataFrame(
            {"EFO": ["EFO_0000270"], "MONDO": ["MONDO_0004992"]}
        ).set_index("EFO")
        self.assertEqual(lu.map_efo_to_mondo("EFO_0000270", efo2mondo), "MONDO_0004992")

    def test_map_efo_to_mondo_missing(self):
        """Returns None for unknown EFO term."""

        efo2mondo = pd.DataFrame(
            {"EFO": ["EFO_0000270"], "MONDO": ["MONDO_0004992"]}
        ).set_index("EFO")
        self.assertIsNone(lu.map_efo_to_mondo("EFO_9999999", efo2mondo))

    # -- map_mesh_to_mondo tests --

    def test_map_mesh_to_mondo(self):
        """Maps MeSH term to MONDO term."""
        mesh2mondo = {"MESH:D008264": "MONDO_0004992"}
        self.assertEqual(
            lu.map_mesh_to_mondo("MESH:D008264", mesh2mondo), "MONDO_0004992"
        )

    def test_map_mesh_to_mondo_missing(self):
        """Returns None for unknown MeSH term."""
        mesh2mondo = {"MESH:D008264": "MONDO_0004992"}
        self.assertIsNone(lu.map_mesh_to_mondo("MESH:D999999", mesh2mondo))

    # -- map_chembl_to_pubchem tests --

    def test_map_chembl_to_pubchem_single(self):
        """Maps ChEMBL id to single PubChem id."""

        chembl2pubchem = pd.DataFrame(
            {"ChEMBL": ["CHEMBL25"], "PubChem": ["2244"]}
        ).set_index("ChEMBL")
        self.assertEqual(lu.map_chembl_to_pubchem("CHEMBL25", chembl2pubchem), "2244")

    def test_map_chembl_to_pubchem_list(self):
        """Maps ChEMBL id to first of multiple PubChem id."""

        chembl2pubchem = pd.DataFrame(
            {"ChEMBL": ["CHEMBL25", "CHEMBL25"], "PubChem": ["2244", "3344"]}
        ).set_index("ChEMBL")
        self.assertEqual(lu.map_chembl_to_pubchem("CHEMBL25", chembl2pubchem), "2244")

    def test_map_chembl_to_pubchem_missing(self):
        """Returns None for unknown ChEMBL id."""

        chembl2pubchem = pd.DataFrame(
            {"ChEMBL": ["CHEMBL25"], "PubChem": ["2244"]}
        ).set_index("ChEMBL")
        self.assertIsNone(lu.map_chembl_to_pubchem("CHEMBL99999", chembl2pubchem))

    # -- collect_unique_gene_names tests --

    def test_collect_unique_gene_names(self):
        """Extracts unique gene names from NSForest markers and binary genes."""

        df = pd.DataFrame(
            {
                "clusterName": ["cluster1", "cluster2", "small_cluster"],
                "clusterSize": [100, 200, 5],
                "NSForest_markers": ["['TP53', 'BRCA1']", "['EGFR']", "['MYC']"],
                "binary_genes": ["['TP53', 'EGFR']", "['BRCA2']", "['KRAS']"],
            }
        )
        genes = lu.collect_unique_gene_names(df)
        self.assertIn("TP53", genes)
        self.assertIn("BRCA1", genes)
        self.assertIn("EGFR", genes)
        self.assertIn("BRCA2", genes)
        # Small cluster (size 5 < MIN_CLUSTER_SIZE=10) should be excluded
        self.assertNotIn("MYC", genes)
        self.assertNotIn("KRAS", genes)
        # Should be sorted
        self.assertEqual(genes, sorted(genes))
