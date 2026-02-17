from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path("../../main/python").resolve()))

from AnnotationResultsTupleWriter import normalize_term


class NormalizeTermTestCase(unittest.TestCase):
    """Tests for normalize_term covering all annotation type branches."""

    def setUp(self):
        self.mesh2mondo = {
            "MESH:D000077192": "MONDO_0004991",
        }

    def _make_annotation(self, term, atype, name, identifier):
        """Helper to build an annotation dict for a given term role."""
        other = "object" if term == "subject" else "subject"
        return {
            f"{term}_type": atype,
            f"{term}_name": name,
            f"{term}_identifier": identifier,
            f"{other}_type": "Gene",
            f"{other}_name": "TP53",
            f"{other}_identifier": "some-uuid-001",
        }

    # -- Anatomical_structure --

    def test_anatomical_structure(self):
        """Replaces colon with underscore."""
        ann = self._make_annotation(
            "subject", "Anatomical_structure", "lung", "UBERON:0002048"
        )
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertEqual(result, "UBERON_0002048")

    # -- Assay --

    def test_assay(self):
        """Replaces colon with underscore."""
        ann = self._make_annotation("subject", "Assay", "some assay", "EFO:0002772")
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertEqual(result, "EFO_0002772")

    # -- Biomarker_combination --

    def test_biomarker_combination(self):
        """Builds BMC_ term from name and subject identifier prefix."""
        ann = self._make_annotation(
            "object", "Biomarker_combination", "marker_set", "some-id"
        )
        ann["subject_identifier"] = "abc123-rest"
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "BMC_marker_set-abc123")

    # -- Cell_set --

    def test_cell_set(self):
        """Builds CS_ term with hyphenated name and subject identifier prefix."""
        ann = self._make_annotation("subject", "Cell_set", "T cells alpha", "some-id")
        ann["subject_identifier"] = "uuid123-rest"
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertEqual(result, "CS_T-cells-alpha-uuid123")

    # -- Cell_set_dataset --

    def test_cell_set_dataset(self):
        """Replaces NLP_dataset with CSD."""
        ann = self._make_annotation(
            "subject", "Cell_set_dataset", "dataset1", "NLP_dataset_v1"
        )
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertEqual(result, "CSD_v1")

    # -- Cell_type --

    def test_cell_type(self):
        """Removes skos:related tag and replaces colon with underscore."""
        ann = self._make_annotation(
            "object", "Cell_type", "macrophage", "<skos:related>CL:0000235"
        )
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "CL_0000235")

    def test_cell_type_no_skos(self):
        """Works without skos:related tag."""
        ann = self._make_annotation("object", "Cell_type", "macrophage", "CL:0000235")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "CL_0000235")

    # -- Disease --

    def test_disease_maps_mesh_to_mondo(self):
        """Maps MESH identifier to MONDO term."""
        ann = self._make_annotation("object", "Disease", "some disease", "MESH:D000077192")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "MONDO_0004991")

    def test_disease_unknown_mesh_returns_none(self):
        """Unknown MESH term returns None."""
        ann = self._make_annotation(
            "object", "Disease", "unknown disease", "MESH:D999999"
        )
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertIsNone(result)

    # -- Gene --

    def test_gene_uppercase(self):
        """Uppercase gene name returns GS_ prefixed term."""
        ann = self._make_annotation("subject", "Gene", "TP53", "some-id")
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertEqual(result, "GS_TP53")

    def test_gene_mixed_case_returns_none(self):
        """Mixed case gene name returns None."""
        ann = self._make_annotation("subject", "Gene", "Tp53", "some-id")
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertIsNone(result)

    def test_gene_mbp_substitution(self):
        """Myelin basic protein is replaced with MBP."""
        ann = self._make_annotation(
            "subject", "Gene", "Myelin basic protein", "some-id"
        )
        result = normalize_term(ann, "subject", self.mesh2mondo)
        self.assertEqual(result, "GS_MBP")

    # -- Publication --

    def test_publication_known_pmid_jorstad(self):
        """Known PMID 37824655 returns correct DOI-based term."""
        ann = self._make_annotation("object", "Publication", "Jorstad", "37824655")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "PUB_doi.org/10.1126/science.adf6812")

    def test_publication_known_pmid_guo(self):
        """Known PMID 37516747 returns correct DOI-based term."""
        ann = self._make_annotation("object", "Publication", "Guo", "37516747")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "PUB_doi.org/10.1038/s41467-023-40173-5")

    def test_publication_known_pmid_sikkema(self):
        """Known PMID 37291214 returns correct DOI-based term."""
        ann = self._make_annotation("object", "Publication", "Sikkema", "37291214")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "PUB_doi.org/10.1038/s41591-023-02327-2")

    def test_publication_known_pmid_li(self):
        """Known PMID 38014002 returns correct DOI-based term."""
        ann = self._make_annotation("object", "Publication", "Li", "38014002")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertEqual(result, "PUB_doi.org/10.1101/2023.11.07.566105")

    def test_publication_unknown_pmid_returns_none(self):
        """Unknown PMID returns None."""
        ann = self._make_annotation("object", "Publication", "Unknown", "99999999")
        result = normalize_term(ann, "object", self.mesh2mondo)
        self.assertIsNone(result)

    # -- Unicode handling --

    def test_unicode_gamma_delta_replacement(self):
        """Unicode gamma-delta characters are replaced."""
        ann = self._make_annotation("subject", "Gene", "\u03b3\u03b4", "some-id")
        result = normalize_term(ann, "subject", self.mesh2mondo)
        # "gamma-delta" is mixed case so Gene branch returns None
        self.assertIsNone(result)
        # But verify the replacement happened in the annotation
        self.assertEqual(ann["subject_name"], "gamma-delta")

    def test_unicode_minus_sign_replacement(self):
        """Unicode minus sign is replaced with hyphen."""
        ann = self._make_annotation("subject", "Gene", "CD4\u2212", "some\u2212id")
        # Verify replacements happened
        normalize_term(ann, "subject", self.mesh2mondo)
        self.assertNotIn("\u2212", ann["subject_name"])
        self.assertNotIn("\u2212", ann["subject_identifier"])
