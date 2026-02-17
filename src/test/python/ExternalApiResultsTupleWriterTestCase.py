import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path("../../main/python").resolve()))

import pandas as pd

from ExternalApiResultsTupleWriter import (
    create_tuples_from_cellxgene,
    create_tuples_from_gene,
    create_tuples_from_hubmap,
    create_tuples_from_opentargets,
    create_tuples_from_uniprot,
    get_mondo_term,
    get_protein_term,
    remove_protocols,
)

SUMMARIES_DIRPATH = Path("../../test/data/summaries")


def to_string_tuples(tuples):
    """Convert list of URIRef/Literal tuples to list of string lists."""
    return [list(str(x) for x in t) for t in tuples]


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
            {"EFO": ["EFO_0008524"], "MONDO": ["MONDO_0002120"]}
        ).set_index("EFO")

    def test_mondo_term_passes_through(self):
        """MONDO term is returned directly."""
        result = get_mondo_term("MONDO_0002120", self.efo2mondo)
        self.assertEqual(result, "MONDO_0002120")

    def test_efo_term_maps_to_mondo(self):
        """EFO term is mapped to MONDO via efo2mondo."""
        result = get_mondo_term("EFO_0008524", self.efo2mondo)
        self.assertEqual(result, "MONDO_0002120")

    def test_efo_term_not_found(self):
        """Unknown EFO term returns None."""
        result = get_mondo_term("EFO_9999999", self.efo2mondo)
        self.assertIsNone(result)

    def test_non_mondo_non_efo_returns_none(self):
        """Term that is neither MONDO nor EFO returns None."""
        result = get_mondo_term("DOID_1234", self.efo2mondo)
        self.assertIsNone(result)

    @patch("ExternalApiResultsTupleWriter.DEPRECATED_TERMS", ["MONDO_0002120"])
    def test_deprecated_mondo_returns_none(self):
        """Deprecated MONDO term returns None."""
        result = get_mondo_term("MONDO_0002120", self.efo2mondo)
        self.assertIsNone(result)

    @patch("ExternalApiResultsTupleWriter.DEPRECATED_TERMS", ["MONDO_0002120"])
    def test_efo_mapping_to_deprecated_mondo_returns_none(self):
        """EFO mapped to deprecated MONDO term returns None."""
        result = get_mondo_term("EFO_0008524", self.efo2mondo)
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


class CreateTuplesFromCellxgeneTestCase(unittest.TestCase):
    """Tests for create_tuples_from_cellxgene using summary fixture."""

    def setUp(self):
        summary_path = SUMMARIES_DIRPATH / "cell-kn-mvp-external-api-results.json"
        with open(summary_path, "r") as fp:
            self.summary = json.load(fp)
        self.cellxgene_results = self.summary["cellxgene"]

        # Extract expected cellxgene tuples (CSD_ or PUB_ subjects) from combined tuples
        self.expected_tuples = [
            t for t in self.summary["tuples"] if "CSD_" in t[0] or "PUB_" in t[0]
        ]

    def test_create_tuples_from_cellxgene(self):
        """Tuples created from summary cellxgene data match expected."""
        actual_tuples, _ = create_tuples_from_cellxgene(
            self.cellxgene_results
        )
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    def test_tuple_count(self):
        """Number of tuples matches expected count."""
        actual_tuples, _ = create_tuples_from_cellxgene(
            self.cellxgene_results
        )
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    def test_first_tuple_is_csd_source(self):
        """First tuple is a CSD dc:Source PUB relation."""
        actual_tuples, _ = create_tuples_from_cellxgene(
            self.cellxgene_results
        )
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("CSD_", first[0])
        self.assertIn("dc#Source", first[1])
        self.assertIn("PUB_", first[2])

    def test_last_tuple_is_zenodo_annotation(self):
        """Last tuple is a CSD Zenodo/Nextflow_workflow/Notebook annotation."""
        actual_tuples, _ = create_tuples_from_cellxgene(
            self.cellxgene_results
        )
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("CSD_", last[0])
        self.assertIn("Zenodo/Nextflow_workflow/Notebook", last[1])
        self.assertEqual(last[2], "TBC")


class CreateTuplesFromOpenTargetsTestCase(unittest.TestCase):
    """Tests for create_tuples_from_opentargets using summary fixture."""

    def setUp(self):
        summary_path = SUMMARIES_DIRPATH / "cell-kn-mvp-external-api-results.json"
        with open(summary_path, "r") as fp:
            self.summary = json.load(fp)
        self.opentargets_data = self.summary["opentargets"]

        # Reconstruct the input format expected by create_tuples_from_opentargets
        gene_ensembl_ids = [
            k for k in self.opentargets_data.keys() if k.startswith("ENSG")
        ]
        self.opentargets_results = dict(self.opentargets_data)
        self.opentargets_results["gene_ensembl_ids"] = gene_ensembl_ids

        # Gene results for UniProt lookups
        self.gene_results = {
            "gene_entrez_ids": ["1080"],
            "1080": {
                "UniProt_name": "P13569",
                "Link_to_UniProt_ID": "https://www.uniprot.org/uniprot/P13569",
            },
        }

        # Extract expected opentargets tuples from combined tuples
        opentargets_subjects = (
            "GS_CFTR",
            "MONDO_",
            "CHEMBL_",
            "RS_",
            "SO_",
        )
        self.expected_tuples = [
            t
            for t in self.summary["tuples"]
            if any(s in t[0] for s in opentargets_subjects)
        ]

        # Mapping DataFrames
        self.gene_ensembl_id_to_names = pd.DataFrame(
            {"external_gene_name": ["CFTR"]},
            index=pd.Index(["ENSG00000001626"], name="ensembl_gene_id"),
        )
        self.gene_name_to_entrez_ids = pd.DataFrame(
            {"entrezgene_id": ["1080"]},
            index=pd.Index(["CFTR"], name="external_gene_name"),
        )
        self.efo2mondo = pd.DataFrame(
            {"MONDO": ["MONDO_0002120"]},
            index=pd.Index(["EFO_0008524"], name="EFO"),
        )
        self.chembl2pubchem = pd.DataFrame(
            {"PubChem": [16220172]},
            index=pd.Index(["CHEMBL2010601"], name="ChEMBL"),
        )

    def _create_tuples(self):
        """Call create_tuples_from_opentargets with mocked mappings."""
        with (
            patch(
                "ExternalApiResultsTupleWriter.DEPRECATED_TERMS",
                [],
            ),
            patch(
                "ExternalApiResultsTupleWriter.get_gene_ensembl_id_to_names_map",
                return_value=self.gene_ensembl_id_to_names,
            ),
            patch(
                "ExternalApiResultsTupleWriter.get_gene_name_to_entrez_ids_map",
                return_value=self.gene_name_to_entrez_ids,
            ),
            patch(
                "ExternalApiResultsTupleWriter.get_efo_to_mondo_map",
                return_value=self.efo2mondo,
            ),
            patch(
                "ExternalApiResultsTupleWriter.get_chembl_to_pubchem_map",
                return_value=self.chembl2pubchem,
            ),
        ):
            actual_tuples, _ = create_tuples_from_opentargets(
                self.opentargets_results, self.gene_results
            )
        return actual_tuples

    def test_create_tuples_from_opentargets(self):
        """Tuples created from summary opentargets data match expected."""
        actual_tuples = self._create_tuples()
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    def test_tuple_count(self):
        """Number of tuples matches expected count."""
        actual_tuples = self._create_tuples()
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    def test_first_tuple_is_genetic_basis_for(self):
        """First tuple is a Gene GENETIC_BASIS_FOR Disease relation."""
        actual_tuples = self._create_tuples()
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("GS_CFTR", first[0])
        self.assertIn("GENETIC_BASIS_FOR", first[1])
        self.assertIn("MONDO_0009061", first[2])

    def test_last_tuple_is_variant_consequence(self):
        """Last tuple is a Variant_consequence_label annotation."""
        actual_tuples = self._create_tuples()
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("SO_0001583", last[0])
        self.assertIn("Variant_consequence_label", last[1])
        self.assertIn("missense_variant", last[2])


class CreateTuplesFromGeneTestCase(unittest.TestCase):
    """Tests for create_tuples_from_gene using summary fixture."""

    def setUp(self):
        summary_path = SUMMARIES_DIRPATH / "cell-kn-mvp-external-api-results.json"
        with open(summary_path, "r") as fp:
            self.summary = json.load(fp)
        self.gene_data = self.summary["gene"]

        # Reconstruct the input format expected by create_tuples_from_gene
        gene_entrez_ids = list(self.gene_data.keys())
        self.gene_results = dict(self.gene_data)
        self.gene_results["gene_entrez_ids"] = gene_entrez_ids

        # Extract expected gene tuples (GS_CDH2 subjects) from combined tuples
        self.expected_tuples = [t for t in self.summary["tuples"] if "GS_CDH2" in t[0]]

        # Mapping: gene entrez id "1000" -> gene name "CDH2"
        self.gene_entrez_id_to_names = pd.DataFrame(
            {"external_gene_name": ["CDH2"]},
            index=pd.Index(["1000"], name="entrezgene_id"),
        )

    def _create_tuples(self):
        """Call create_tuples_from_gene with mocked mapping."""
        with patch(
            "ExternalApiResultsTupleWriter.get_gene_entrez_id_to_names_map",
            return_value=self.gene_entrez_id_to_names,
        ):
            actual_tuples, _ = create_tuples_from_gene(
                self.gene_results
            )
        return actual_tuples

    def test_create_tuples_from_gene(self):
        """Tuples created from summary gene data match expected."""
        actual_tuples = self._create_tuples()
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    def test_tuple_count(self):
        """Number of tuples matches expected count."""
        actual_tuples = self._create_tuples()
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    def test_first_tuple_is_produces_relation(self):
        """First tuple is a Gene PRODUCES Protein relation."""
        actual_tuples = self._create_tuples()
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("GS_CDH2", first[0])
        self.assertIn("PRODUCES", first[1])
        self.assertIn("PR_P19022", first[2])

    def test_last_tuple_is_mrna_sequences(self):
        """Last tuple is a mRNA/protein sequences annotation."""
        actual_tuples = self._create_tuples()
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("GS_CDH2", last[0])
        self.assertIn("mRNA_(NM)_and_protein_(NP)_sequences", last[1])
        self.assertIn("NM_001308176", last[2])


class CreateTuplesFromUniprotTestCase(unittest.TestCase):
    """Tests for create_tuples_from_uniprot using summary fixture."""

    def setUp(self):
        summary_path = SUMMARIES_DIRPATH / "cell-kn-mvp-external-api-results.json"
        with open(summary_path, "r") as fp:
            self.summary = json.load(fp)
        self.uniprot_data = self.summary["uniprot"]

        # Reconstruct the input format expected by create_tuples_from_uniprot
        protein_accessions = list(self.uniprot_data.keys())
        self.uniprot_results = dict(self.uniprot_data)
        self.uniprot_results["protein_accessions"] = protein_accessions

        # Extract expected uniprot tuples (PR_ subjects) from combined tuples
        self.expected_tuples = [t for t in self.summary["tuples"] if "PR_" in t[0]]

    def test_create_tuples_from_uniprot(self):
        """Tuples created from summary uniprot data match expected."""
        actual_tuples, _ = create_tuples_from_uniprot(
            self.uniprot_results
        )
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    def test_tuple_count(self):
        """Number of tuples matches expected count."""
        actual_tuples, _ = create_tuples_from_uniprot(
            self.uniprot_results
        )
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    def test_first_tuple_is_protein_name(self):
        """First tuple is a Protein_name annotation."""
        actual_tuples, _ = create_tuples_from_uniprot(
            self.uniprot_results
        )
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("PR_P55017", first[0])
        self.assertIn("Protein_name", first[1])
        self.assertEqual(first[2], "Solute carrier family 12 member 3")

    def test_last_tuple_is_organism(self):
        """Last tuple is an Organism annotation."""
        actual_tuples, _ = create_tuples_from_uniprot(
            self.uniprot_results
        )
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("PR_P55017", last[0])
        self.assertIn("Organism", last[1])
        self.assertEqual(last[2], "Homo sapiens")


class CreateTuplesFromHubmapTestCase(unittest.TestCase):
    """Tests for create_tuples_from_hubmap using summary fixture."""

    def setUp(self):
        summary_path = SUMMARIES_DIRPATH / "hubmap-allen-brain-v1.7.json"
        with open(summary_path, "r") as fp:
            self.summary = json.load(fp)
        self.expected_tuples = self.summary["tuples"]

        # Reconstruct the original HuBMAP input format from summary data
        self.hubmap_data = {
            "data": {
                "anatomical_structures": [
                    self.summary["hubmap"]["anatomical_structure"]
                ],
                "cell_types": [self.summary["hubmap"]["cell_type"]],
            }
        }

        # Extract cl_terms from the summary cell_type id
        cl_id = self.summary["hubmap"]["cell_type"]["id"]
        self.cl_terms = {cl_id.replace(":", "_")}

    @patch("ExternalApiResultsTupleWriter.DEPRECATED_TERMS", [])
    def test_create_tuples_from_hubmap(self):
        """Tuples created from reconstructed HuBMAP data match expected."""
        actual_tuples, _ = create_tuples_from_hubmap(
            self.hubmap_data, self.cl_terms
        )
        actual_as_strings = to_string_tuples(actual_tuples)
        self.assertEqual(actual_as_strings, self.expected_tuples)

    @patch("ExternalApiResultsTupleWriter.DEPRECATED_TERMS", [])
    def test_tuple_count(self):
        """Number of tuples matches expected count."""
        actual_tuples, _ = create_tuples_from_hubmap(
            self.hubmap_data, self.cl_terms
        )
        self.assertEqual(len(actual_tuples), len(self.expected_tuples))

    def test_first_tuple_is_brain_part_of_body_proper(self):
        """First tuple is Brain PART_OF Body_proper relation."""
        actual_tuples, _ = create_tuples_from_uniprot(self.hubmap_data)
        first = list(str(x) for x in actual_tuples[0])
        self.assertIn("UBERON_0000955", first[0])
        self.assertIn("PART_OF", first[1])
        self.assertEqual(first[2], "UBERON_0013702")

    def test_last_tuple_is_b_cell_part_of_brain_source(self):
        """Last tuple is B cell PART_OF Brain Source is HuBMAP edge annotation."""
        actual_tuples, _ = create_tuples_from_uniprot(self.hubmap_data)
        last = list(str(x) for x in actual_tuples[-1])
        self.assertIn("CL_0000236", last[0])
        self.assertIn("UBERON_0000955", last[1])
        self.assertIn("Source", last[2])
        self.assertEqual(last[3], "HuBMAP")
