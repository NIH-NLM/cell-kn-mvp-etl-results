import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest

from arango import ArangoClient
from rdflib import Graph
from rdflib.term import BNode, Literal, URIRef

sys.path.insert(0, str(Path("../../main/python").resolve()))

import OntologyParserLoader as opl

SH_DIR = Path(__file__).parents[2] / "main" / "shell"
ARANGODB_DIR = Path(__file__).parent / "arangodb"
ARANGO_URL = "http://localhost:8529"

# Test data directory containing macrophage.owl and ro.owl
TEST_OBO_DIR = Path(__file__).parents[2].parent / "src" / "test" / "data" / "obo"


class OntologyParserLoaderUnitTestCase(unittest.TestCase):
    """Pure unit tests that do NOT require ArangoDB."""

    def setUp(self):
        # Reset module-level mutable state between tests
        opl.SKIPPED_VERTICES.clear()

    # -- parse_term tests --

    def test_parse_term_obo_class(self):
        """Standard OBO class URIRef parses correctly."""
        term = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        oid, number, t, label, term_type = opl.parse_term(term)
        self.assertEqual(oid, "CL")
        self.assertEqual(number, "0000235")
        self.assertEqual(t, "CL_0000235")
        self.assertIsNone(label)
        self.assertEqual(term_type, "class")

    def test_parse_term_with_ro_lookup(self):
        """OBO class with RO dict lookup returns label."""
        ro = {"RO_0002175": "present in taxon"}
        term = URIRef("http://purl.obolibrary.org/obo/RO_0002175")
        oid, number, t, label, term_type = opl.parse_term(term, ro=ro)
        self.assertEqual(oid, "RO")
        self.assertEqual(number, "0002175")
        self.assertEqual(t, "RO_0002175")
        self.assertEqual(label, "present in taxon")
        self.assertEqual(term_type, "class")

    def test_parse_term_predicate_fragment(self):
        """URIRef with fragment parses as predicate."""
        term = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        oid, number, t, fragment, term_type = opl.parse_term(term)
        self.assertIsNone(oid)
        self.assertEqual(fragment, "subClassOf")
        self.assertEqual(term_type, "predicate")

    def test_parse_term_bnode(self):
        """BNode parses as class with pseudo identifiers."""
        term = BNode("N1234567890ab")
        oid, number, t, label, term_type = opl.parse_term(term)
        self.assertEqual(oid, "BNode")
        self.assertEqual(term_type, "class")

    def test_parse_term_gorel_rejection(self):
        """GOREL identifier returns all None."""
        term = URIRef("http://purl.obolibrary.org/obo/GOREL_0001006")
        result = opl.parse_term(term)
        self.assertEqual(result, (None, None, None, None, None))

    def test_parse_term_empty_oid_or_number(self):
        """Empty oid or number returns all None."""
        term = URIRef("http://purl.obolibrary.org/obo/_0000235")
        result = opl.parse_term(term)
        self.assertEqual(result, (None, None, None, None, None))

    def test_parse_term_literal_fallback(self):
        """URIRef with no fragment and no OBO match parses as literal."""
        term = URIRef("http://example.org/somepath/somevalue")
        _, _, _, value, term_type = opl.parse_term(term)
        self.assertEqual(term_type, "literal")
        self.assertEqual(value, "somevalue")

    # -- parse_obo tests --

    def test_parse_obo(self):
        """parse_obo extracts term-to-label mapping from OWL file."""
        if not TEST_OBO_DIR.exists():
            self.skipTest(f"Test OBO dir not found: {TEST_OBO_DIR}")
        t2l, l2t, ids = opl.parse_obo(TEST_OBO_DIR, "macrophage.owl")
        self.assertIn("CL_0000235", t2l)
        self.assertEqual(t2l["CL_0000235"], "macrophage")
        self.assertIn("macrophage", l2t)
        self.assertEqual(l2t["macrophage"], "CL_0000235")
        self.assertIn("CL", ids)

    # -- get_fnode tests --

    def test_get_fnode_bnode_subject(self):
        """When subject is BNode, returns object."""
        bnode = BNode("abc123")
        uriref = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        self.assertEqual(opl.get_fnode(bnode, uriref), uriref)

    def test_get_fnode_bnode_object(self):
        """When object is BNode, returns subject."""
        bnode = BNode("abc123")
        uriref = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        self.assertEqual(opl.get_fnode(uriref, bnode), uriref)

    def test_get_fnode_both_blank_raises(self):
        """Both BNode raises Exception."""
        with self.assertRaises(Exception):
            opl.get_fnode(BNode("a"), BNode("b"))

    def test_get_fnode_both_filled_raises(self):
        """Both filled raises Exception."""
        u1 = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        u2 = URIRef("http://purl.obolibrary.org/obo/GO_0031268")
        with self.assertRaises(Exception):
            opl.get_fnode(u1, u2)

    # -- count_triple_types tests --

    def test_count_triple_types(self):
        """Counts triple types correctly."""
        g = Graph()
        s1 = URIRef("http://example.org/s1")
        p1 = URIRef("http://example.org/p1")
        o1 = URIRef("http://example.org/o1")
        o2 = Literal("value")
        b = BNode()
        g.add((s1, p1, o1))
        g.add((s1, p1, o2))
        g.add((b, p1, o1))

        types = opl.count_triple_types(g)
        self.assertEqual(types[(URIRef, URIRef, URIRef)], 1)
        self.assertEqual(types[(URIRef, URIRef, Literal)], 1)
        self.assertEqual(types[(BNode, URIRef, URIRef)], 1)

    # -- collect_fnode_triples tests --

    def test_collect_fnode_triples(self):
        """Excludes triples with BNode subject or object."""
        g = Graph()
        s1 = URIRef("http://example.org/s1")
        p1 = URIRef("http://example.org/p1")
        o1 = URIRef("http://example.org/o1")
        b = BNode()
        g.add((s1, p1, o1))  # fnode triple
        g.add((b, p1, o1))  # has BNode subject
        g.add((s1, p1, b))  # has BNode object

        fnode_triples = opl.collect_fnode_triples(g)
        self.assertEqual(len(fnode_triples), 1)
        self.assertEqual(fnode_triples[0], (s1, p1, o1))

    # -- collect_bnode_triple_sets tests --

    def test_collect_bnode_triple_sets(self):
        """Categorizes BNode triples into relation/annotation/other."""
        g = Graph()
        b = BNode()
        s_class = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        o_class = URIRef("http://purl.obolibrary.org/obo/CL_0000113")

        # Relation predicates
        p_subclass = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        p_on_prop = URIRef("http://www.w3.org/2002/07/owl#onProperty")
        p_some = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")

        # BNode as subject
        g.add((b, p_subclass, s_class))
        g.add((b, p_on_prop, o_class))
        g.add((b, p_some, o_class))

        triple_sets = {}
        opl.collect_bnode_triple_sets(g, triple_sets, use="subject")
        self.assertIn(b, triple_sets)
        self.assertEqual(len(triple_sets[b]["relation"]), 3)

    # -- create_bnode_triples_from_bnode_triple_set tests --

    def test_create_bnode_triples_from_bnode_triple_set_relation(self):
        """Creates triple from relation set with exactly 3 triples."""
        b = BNode()
        s_class = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p_rel = URIRef("http://purl.obolibrary.org/obo/RO_0002175")
        o_class = URIRef("http://purl.obolibrary.org/obo/CL_0000113")

        p_subclass = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        p_on_prop = URIRef("http://www.w3.org/2002/07/owl#onProperty")
        p_some = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")

        triple_set = {
            "relation": [
                (s_class, p_subclass, b),  # s_class is the filled node → subject
                (b, p_on_prop, p_rel),  # p_rel is the filled node → predicate
                (b, p_some, o_class),  # o_class is the filled node → object
            ],
            "annotation": [],
            "literal": [],
            "class": [],
            "other": [],
        }

        triples, ignored = opl.create_bnode_triples_from_bnode_triple_set(
            triple_set, "relation"
        )
        self.assertEqual(len(triples), 1)
        self.assertEqual(len(ignored), 0)
        self.assertEqual(triples[0], (s_class, p_rel, o_class))

    def test_create_bnode_triples_from_bnode_triple_set_invalid_type(self):
        """Invalid set_type raises Exception."""
        with self.assertRaises(Exception):
            opl.create_bnode_triples_from_bnode_triple_set({}, "invalid")

    def test_create_bnode_triples_from_bnode_triple_set_empty(self):
        """Empty relation set returns empty lists."""
        triple_set = {
            "relation": [],
            "annotation": [],
            "literal": [],
            "class": [],
            "other": [],
        }
        triples, ignored = opl.create_bnode_triples_from_bnode_triple_set(
            triple_set, "relation"
        )
        self.assertEqual(triples, [])
        self.assertEqual(ignored, [])

    # -- create_bnode_triples_from_bnode_triple_sets tests --

    def test_create_bnode_triples_from_bnode_triple_sets(self):
        """Processes multiple BNode triple sets."""
        b1 = BNode()
        b2 = BNode()
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://purl.obolibrary.org/obo/RO_0002175")
        o = URIRef("http://purl.obolibrary.org/obo/CL_0000113")

        p_sub = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        p_on = URIRef("http://www.w3.org/2002/07/owl#onProperty")
        p_some = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")

        triple_sets = {
            b1: {
                "relation": [
                    (s, p_sub, b1),
                    (b1, p_on, p),
                    (b1, p_some, o),
                ],
                "annotation": [],
                "literal": [],
                "class": [],
                "other": [],
            },
            b2: {
                "relation": [],
                "annotation": [],
                "literal": [],
                "class": [(s, p, o)],
                "other": [],
            },
        }

        bnode_triples, ignored = opl.create_bnode_triples_from_bnode_triple_sets(
            triple_sets
        )
        self.assertEqual(len(bnode_triples), 1)
        # b2's class triples go to ignored
        self.assertEqual(len(ignored), 1)

    # -- create_or_get_vertex tests --

    def test_create_or_get_vertex_creates(self):
        """Creates new vertex in empty collections."""
        vc = {}
        v = opl.create_or_get_vertex(vc, "CL", "0000235", "CL_0000235")
        self.assertEqual(v, {"_key": "0000235", "term": "CL_0000235"})
        self.assertIn("CL", vc)
        self.assertIn("0000235", vc["CL"])

    def test_create_or_get_vertex_gets_existing(self):
        """Returns existing vertex on second call."""
        vc = {}
        v1 = opl.create_or_get_vertex(vc, "CL", "0000235", "CL_0000235")
        v2 = opl.create_or_get_vertex(vc, "CL", "0000235", "CL_0000235")
        self.assertIs(v1, v2)

    def test_create_or_get_vertex_invalid_name(self):
        """Invalid vertex name returns None."""
        vc = {}
        v = opl.create_or_get_vertex(vc, "INVALID", "001", "INVALID_001")
        self.assertIsNone(v)

    # -- create_or_get_vertices_from_triple tests --

    def test_create_or_get_vertices_from_triple(self):
        """Creates vertices for subject and object."""
        vc = {}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        o = URIRef("http://purl.obolibrary.org/obo/CL_0000113")
        vertices = opl.create_or_get_vertices_from_triple(vc, s, p, o)
        self.assertEqual(len(vertices), 2)
        self.assertIn("CL", vc)
        self.assertIn("0000235", vc["CL"])
        self.assertIn("0000113", vc["CL"])

    def test_create_or_get_vertices_from_triple_literal_object(self):
        """Returns None for Literal object."""
        vc = {}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        o = Literal("macrophage")
        result = opl.create_or_get_vertices_from_triple(vc, s, p, o)
        self.assertIsNone(result)

    # -- create_or_get_edge_from_triple tests --

    def test_create_or_get_edge_from_triple(self):
        """Creates edge for valid triple."""
        vc = {}
        ec = {}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        o = URIRef("http://purl.obolibrary.org/obo/CL_0000113")
        edge = opl.create_or_get_edge_from_triple(vc, ec, s, p, o)
        self.assertIsNotNone(edge)
        self.assertEqual(edge["_from"], "CL/0000235")
        self.assertEqual(edge["_to"], "CL/0000113")
        self.assertEqual(edge["label"], "subClassOf")

    def test_create_or_get_edge_from_triple_literal_object(self):
        """Returns None for Literal object."""
        vc = {}
        ec = {}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        o = Literal("macrophage")
        result = opl.create_or_get_edge_from_triple(vc, ec, s, p, o)
        self.assertIsNone(result)

    # -- create_or_get_edge tests --

    def test_create_or_get_edge(self):
        """Creates edge document with correct keys."""
        vc = {}
        ec = {}
        edge = opl.create_or_get_edge(
            vc,
            ec,
            "CL",
            "0000235",
            "CL_0000235",
            "CL",
            "0000113",
            "CL_0000113",
            "subClassOf",
        )
        self.assertIsNotNone(edge)
        self.assertEqual(edge["_key"], "0000235-0000113")
        self.assertEqual(edge["_from"], "CL/0000235")
        self.assertEqual(edge["_to"], "CL/0000113")
        self.assertEqual(edge["label"], "subClassOf")
        self.assertIn("CL-CL", ec)

    def test_create_or_get_edge_returns_existing(self):
        """Second call returns existing edge."""
        vc = {}
        ec = {}
        e1 = opl.create_or_get_edge(
            vc,
            ec,
            "CL",
            "0000235",
            "CL_0000235",
            "CL",
            "0000113",
            "CL_0000113",
            "subClassOf",
        )
        e2 = opl.create_or_get_edge(
            vc,
            ec,
            "CL",
            "0000235",
            "CL_0000235",
            "CL",
            "0000113",
            "CL_0000113",
            "subClassOf",
        )
        self.assertIs(e1, e2)

    # -- update_vertex_from_triple tests --

    def test_update_vertex_from_triple_adds_annotation(self):
        """Adds literal value as vertex annotation."""
        vc = {"CL": {"0000235": {"_key": "0000235", "term": "CL_0000235"}}}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        o = Literal("macrophage")
        v = opl.update_vertex_from_triple(vc, s, p, o)
        self.assertIsNotNone(v)
        self.assertEqual(v["label"], "macrophage")

    def test_update_vertex_from_triple_creates_list_on_duplicate(self):
        """Second annotation with same predicate creates list."""
        vc = {"CL": {"0000235": {"_key": "0000235", "term": "CL_0000235"}}}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        opl.update_vertex_from_triple(vc, s, p, Literal("macrophage"))
        v = opl.update_vertex_from_triple(vc, s, p, Literal("histiocyte"))
        self.assertIsInstance(v["label"], list)
        self.assertIn("macrophage", v["label"])
        self.assertIn("histiocyte", v["label"])

    def test_update_vertex_from_triple_non_literal_returns_none(self):
        """Non-Literal object returns None."""
        vc = {"CL": {"0000235": {"_key": "0000235", "term": "CL_0000235"}}}
        s = URIRef("http://purl.obolibrary.org/obo/CL_0000235")
        p = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        o = URIRef("http://purl.obolibrary.org/obo/CL_0000113")
        result = opl.update_vertex_from_triple(vc, s, p, o)
        self.assertIsNone(result)

    # -- find_obo_version tests --

    def test_find_obo_version(self):
        """Extracts version from macrophage.owl."""
        if not TEST_OBO_DIR.exists():
            self.skipTest(f"Test OBO dir not found: {TEST_OBO_DIR}")
        version = opl.find_obo_version(TEST_OBO_DIR / "macrophage.owl")
        self.assertIsNotNone(version)
        # Version should be a date string
        self.assertRegex(version, r"\d{4}-\d{2}-\d{2}")


class OntologyParserLoaderTestCase(unittest.TestCase):
    """Integration tests that require a running ArangoDB instance."""

    @classmethod
    def setUpClass(cls):

        # Stop any ArangoDB instance
        subprocess.run(["./stop-arangodb.sh"], cwd=SH_DIR)

        # Start an ArangoDB instance using the test data directory
        os.environ["ARANGO_DB_HOME"] = str(ARANGODB_DIR)
        subprocess.run(["./start-arangodb.sh"], cwd=SH_DIR)

    def test_main(self):
        """Compare actual and expected macrophage vertex and edges,
        obtaining expected values by manual inspection of the
        macrophage OWL file.
        """
        arango_client = ArangoClient(hosts=ARANGO_URL)
        arango_root_password = os.environ["ARANGO_DB_PASSWORD"]

        # Parse macrophage OWL file and load the result into ArangoDB
        opl.main(parameters=["--test"])

        # Connect to ArangoDB
        db = arango_client.db(
            "Cell-KN-v1.5", username="root", password=arango_root_password
        )
        graph = db.graph("CL-Test")

        # Get the actual macrophage vertex
        vertex_collection = graph.vertex_collection("CL")
        key = "0000235"
        self.assertTrue(vertex_collection.has(key))
        a_vertex = vertex_collection.get(key)

        # Define the expected macrophage vertex
        e_vertex = {
            "_key": "0000235",
            "_id": "CL/0000235",
            "_rev": "_jPkjjC6---",
            "term": "CL_0000235",
            "hasExactSynonym": "histiocyte",
            "hasDbXref": [
                "MESH:D008264",
                "ZFA:0009141",
                "FMA:63261",
                "BTO:0000801",
                "CALOHA:TS-0587",
                "FMA:83585",
                "PMID:16213494",
                "GOC:tfm",
                "GO_REF:0000031",
                "PMID:1919437",
                "GOC:add",
            ],
            "comment": "Morphology: Diameter 30_M-80 _M, abundant cytoplasm, low N/C ratio, eccentric nucleus. Irregular shape with pseudopods, highly adhesive. Contain vacuoles and phagosomes, may contain azurophilic granules; markers: Mouse & Human: CD68, in most cases CD11b. Mouse: in most cases F4/80+; role or process: immune, antigen presentation, & tissue remodelling; lineage: hematopoietic, myeloid.",
            "label": "macrophage",
            "id": "CL:0000235",
            "definition": [
                "A mononuclear phagocyte present in variety of tissues, typically differentiated from monocytes, capable of phagocytosing a variety of extracellular particulate material, including immune complexes, microorganisms, and dead cells."
            ],
        }

        # Assert equal vertex keys
        self.assertTrue(sorted(a_vertex.keys()) == sorted(e_vertex.keys()))

        # Assert equal vertex values, ignoring the revision
        for e_key, e_value in e_vertex.items():
            if e_key == "_rev":
                continue
            a_value = a_vertex[e_key]
            if type(e_value) == list:
                self.assertTrue(type(a_value) == list)
                self.assertTrue(sorted(a_value) == sorted(e_value))
            else:
                self.assertTrue(a_value == e_value)

        # Get macrophage edges to CL terms, then assert equal labels
        edge_collection = graph.edge_collection("CL-CL")
        keys = ["0000235-0000113", "0000235-0000145", "0000235-0000766"]
        for key in keys:
            self.assertTrue(edge_collection.has(key))
            edge = edge_collection.get(key)
            self.assertTrue(edge["label"] == "subClassOf")
        key = "0000235-0000576"
        self.assertTrue(edge_collection.has(key))
        edge = edge_collection.get(key)
        self.assertTrue(edge["label"] == "develops from")

        # Get macrophage edges to GO terms, then assert equal labels
        edge_collection = graph.edge_collection("CL-GO")
        key = "0000235-0031268"
        self.assertTrue(edge_collection.has(key))
        edge = edge_collection.get(key)
        self.assertTrue(edge["label"] == "capable of")

        # Get macrophage edges to NCBITaxon terms, then assert equal labels
        edge_collection = graph.edge_collection("CL-NCBITaxon")
        key = "0000235-9606"
        self.assertTrue(edge_collection.has(key))
        edge = edge_collection.get(key)
        self.assertTrue(edge["label"] == "present in taxon")

    @classmethod
    def tearDownClass(cls):

        # Stop the ArangoDB instance using the test data directory
        subprocess.run(["./stop-arangodb.sh"], cwd=SH_DIR)

        # Remove ArangoDB test data directory
        if ARANGODB_DIR.exists():
            shutil.rmtree(ARANGODB_DIR)
