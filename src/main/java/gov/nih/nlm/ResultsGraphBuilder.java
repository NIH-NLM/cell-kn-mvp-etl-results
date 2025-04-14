package gov.nih.nlm;

import com.arangodb.ArangoDatabase;
import com.arangodb.ArangoEdgeCollection;
import com.arangodb.ArangoGraph;
import com.arangodb.ArangoVertexCollection;
import com.arangodb.entity.BaseDocument;
import com.arangodb.entity.BaseEdgeDocument;
import org.apache.jena.graph.Node;
import org.apache.jena.graph.NodeFactory;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static gov.nih.nlm.OntologyElementParser.parseOntologyElements;
import static gov.nih.nlm.OntologyTripleLoader.*;
import static gov.nih.nlm.PathUtilities.listFilesMatchingPattern;

public class ResultsGraphBuilder {

    // Assign location of ontology files
    private static final Path usrDir = Paths.get(System.getProperty("user.dir"));
    public static final Path tuplesDir = usrDir.resolve("data/tuples");
    public static final Path schemaDir = usrDir.resolve("data/schema");

    // Connect to a local ArangoDB server instance
    private static final ArangoDbUtilities arangoDbUtilities = new ArangoDbUtilities();

    // Assign triple indices
    private static final int tripleSubjectIdx = 0;
    private static final int triplePredicateIdx = 1;
    private static final int tripleObjectIdx = 2;

    // Assign quadruple indices
    private static final int quadrupleSubjectIdx = 0;
    private static final int quadrupleObjectIdx = 1;
    private static final int quadruplePredicateIdx = 2;
    private static final int quadrupleLiteralIdx = 3;

    public static ArrayList<ArrayList<Node>> readJsonFile(String jsonFilePath) throws IOException {
        ArrayList<ArrayList<Node>> tuplesArrayList = new ArrayList<>();
        try {
            String content = new String(Files.readAllBytes(Paths.get(jsonFilePath)));
            JSONObject jsonObject = new JSONObject(content);
            JSONArray tuplesJsonArray = (JSONArray) jsonObject.get("tuples");
            for (int iTuple = 0; iTuple < tuplesJsonArray.length(); iTuple++) {
                ArrayList<Node> tupleArrayList = new ArrayList<>();
                JSONArray tupleJsonArray = (JSONArray) tuplesJsonArray.get(iTuple);
                for (int iElement = 0; iElement < tupleJsonArray.length(); iElement++) {
                    String value = tupleJsonArray.get(iElement).toString();
                    Node node;
                    if (value.contains("http")) {
                        node = NodeFactory.createURI(value);
                    } else {
                        node = NodeFactory.createLiteral(value);
                    }
                    tupleArrayList.add(node);
                }
                if (tupleArrayList.size() == 3 && !(
                        tupleArrayList.get(tripleSubjectIdx).isURI() && tupleArrayList.get(triplePredicateIdx).isURI() &&
                                (tupleArrayList.get(tripleObjectIdx).isURI() || tupleArrayList.get(tripleObjectIdx).isLiteral()))) {
                    throw new IOException("Invalid triple " + tupleArrayList);
                }
                if (tupleArrayList.size() == 4 && !(
                        tupleArrayList.get(quadrupleSubjectIdx).isURI() && tupleArrayList.get(quadrupleObjectIdx).isURI() &&
                                tupleArrayList.get(quadruplePredicateIdx).isURI() && tupleArrayList.get(quadrupleLiteralIdx).isLiteral())) {
                    throw new IOException("Invalid quadruple " + tupleArrayList);
                }
                tuplesArrayList.add(tupleArrayList);
            }
        } catch (IOException e) {
            System.err.println("Error reading the file: " + e.getMessage());
        } catch (org.json.JSONException e) {
            System.err.println("Error parsing JSON: " + e.getMessage());
        }
        return tuplesArrayList;
    }

    /**
     * Construct vertices using tuples parsed from a results file that
     * contain a filled subject and object which contain an ontology ID contained in
     * the valid vertices collection.
     *
     * @param tuplesArrayList   list of tuples parsed from a results file
     * @param graph             ArangoDB graph in which to create vertex collections
     * @param vertexCollections ArangoDB vertex collections
     * @param vertexDocuments   ArangoDB vertex documents
     */
    public static void constructVertices(ArrayList<ArrayList<Node>> tuplesArrayList, ArangoGraph graph,
                                         Map<String, Set<String>> vertexKeys, Map<String, ArangoVertexCollection> vertexCollections,
                                         Map<String, Map<String, BaseDocument>> vertexDocuments) {

        int nVertices = 0;
        System.out.println("Constructing vertices using " + tuplesArrayList.size() + " tuples");
        long startTime = System.nanoTime();
        for (ArrayList<Node> tupleArrayList : tuplesArrayList) {

            // Only construct vertices using triples
            if (tupleArrayList.size() != 3) continue;

            for (Node n : tupleArrayList) {

                // Only construct valid vertices
                OntologyTripleLoader.VTuple vtuple = createVTuple(n);
                if (!vtuple.isValidVertex()) continue;

                // Create a vertex collection, if needed
                if (!vertexCollections.containsKey(vtuple.id())) {
                    vertexCollections.put(vtuple.id(), arangoDbUtilities.createOrGetVertexCollection(graph, vtuple.id()));
                    vertexDocuments.put(vtuple.id(), new HashMap<>());
                    vertexKeys.put(vtuple.id(), new HashSet<>());
                }

                // Construct the vertex, if needed
                if (!vertexKeys.get(vtuple.id()).contains(vtuple.number())) {
                    nVertices++;
                    BaseDocument doc = new BaseDocument(vtuple.number());
                    vertexDocuments.get(vtuple.id()).put(vtuple.number(), doc);
                    vertexKeys.get(vtuple.id()).add(vtuple.number());
                }
            }
        }
        long stopTime = System.nanoTime();
        System.out.println("Constructed " + nVertices + " vertices using " + tuplesArrayList.size() + " tuples in " + (stopTime - startTime) / 1e9 + " s");
    }

    /**
     * Update vertices using tuples parsed from a results file that
     * contain a filled subject which contains an ontology ID contained in the valid
     * vertices collection, and a filled object literal.
     *
     * @param tuplesArrayList     list of tuples parsed from a results file
     * @param ontologyElementMaps Maps terms and labels
     * @param vertexDocuments     ArangoDB vertex documents
     */
    public static void updateVertices(ArrayList<ArrayList<Node>> tuplesArrayList, Map<String, OntologyElementMap> ontologyElementMaps,
                                      Map<String, Map<String, BaseDocument>> vertexDocuments) throws RuntimeException {

        Set<String> updatedVertices = new HashSet<>(); // For counting only
        System.out.println("Updating vertices using " + tuplesArrayList.size() + " tuples");
        long startTime = System.nanoTime();
        for (ArrayList<Node> tupleArrayList : tuplesArrayList) {

            // Only update vertices using triples
            if (tupleArrayList.size() != 3) continue;

            // Ensure the object contains a literal
            Node o = tupleArrayList.get(tripleObjectIdx);
            if (!o.isLiteral()) {
                continue;
            }

            // Parse the object
            String literal = o.getLiteralValue().toString();

            // Ensure the subject contains a valid ontology ID
            OntologyTripleLoader.VTuple vtuple = createVTuple(tupleArrayList.get(tripleSubjectIdx));
            if (!vtuple.isValidVertex()) continue;

            // Parse the predicate
            String attribute = parsePredicate(ontologyElementMaps, tupleArrayList.get(triplePredicateIdx));

            // Update the corresponding vertex
            if (!vertexDocuments.get(vtuple.id()).containsKey(vtuple.number()))
                throw new RuntimeException("No vertex for VTuple " + vtuple);
            updatedVertices.add(vtuple.id() + "-" + vtuple.number() + "-" + attribute + "-" + literal);
            BaseDocument doc = vertexDocuments.get(vtuple.id()).get(vtuple.number());
            Set<String> literals;
            if (doc.getAttribute(attribute) == null) {
                literals = new HashSet<>();
                doc.addAttribute(attribute, literals);
            } else {
                literals = (HashSet<String>) doc.getAttribute(attribute);
            }
            literals.add(literal);
        }
        long stopTime = System.nanoTime();
        System.out.println("Updated " + updatedVertices.size() + " vertices using " + tuplesArrayList.size() + " tuples in " + (stopTime - startTime) / 1e9 + " s");
    }

    /**
     * @param tuplesArrayList     list of tuples parsed from a results file
     * @param ontologyElementMaps Maps terms and labels
     * @param graph               ArangoDB graph
     * @param edgeCollections     ArangoDB edge collections
     * @param edgeDocuments       ArangoDB edge documents
     */
    public static void constructEdges(ArrayList<ArrayList<Node>> tuplesArrayList, Map<String, OntologyElementMap> ontologyElementMaps, ArangoGraph graph,
                                      Map<String, Set<String>> edgeKeys, Map<String, ArangoEdgeCollection> edgeCollections, Map<String, Map<String, BaseEdgeDocument>> edgeDocuments) throws RuntimeException {

        int nEdges = 0;
        System.out.println("Constructing edges using " + tuplesArrayList.size() + " tuples");
        long startTime = System.nanoTime();
        for (ArrayList<Node> tupleArrayList : tuplesArrayList) {

            // Only construct edges using triples
            if (tupleArrayList.size() != 3) continue;

            // Ensure the subject contains a valid ontology ID
            OntologyTripleLoader.VTuple s_vtuple = createVTuple(tupleArrayList.get(tripleSubjectIdx));
            if (!s_vtuple.isValidVertex()) continue;

            // Ensure the object contains a valid ontology ID
            OntologyTripleLoader.VTuple o_vtuple = createVTuple(tupleArrayList.get(tripleObjectIdx));
            if (!o_vtuple.isValidVertex()) continue;

            // Parse the predicate
            String label = parsePredicate(ontologyElementMaps, tupleArrayList.get(triplePredicateIdx));

            // Create an edge collection, if needed
            String idPair = s_vtuple.id() + "-" + o_vtuple.id();
            if (!edgeCollections.containsKey(idPair)) {
                edgeCollections.put(idPair, arangoDbUtilities.createOrGetEdgeCollection(graph, s_vtuple.id(), o_vtuple.id()));
                edgeDocuments.put(idPair, new HashMap<>());
                edgeKeys.put(idPair, new HashSet<>());
            }

            // Construct the edge, if needed
            String key = s_vtuple.number() + "-" + o_vtuple.number();
            if (!edgeKeys.get(idPair).contains(key)) {
                nEdges++;
                BaseEdgeDocument doc = new BaseEdgeDocument(key, s_vtuple.id() + "/" + s_vtuple.number(), o_vtuple.id() + "/" + o_vtuple.number());
                doc.addAttribute("label", label);
                edgeDocuments.get(idPair).put(key, doc);
                edgeKeys.get(idPair).add(key);
            }
        }
        long stopTime = System.nanoTime();
        System.out.println("Constructed " + nEdges + " edges using " + tuplesArrayList.size() + " tuples in " + (stopTime - startTime) / 1e9 + " s");
    }

    /**
     * Update edges using tuples parsed from a results file that
     * contain a filled subject and object each which contains an ontology ID contained in the valid
     * vertices collection, and two filled predicate literals.
     *
     * @param tuplesArrayList     list of tuples parsed from a results file
     * @param ontologyElementMaps Maps terms and labels
     * @param edgeDocuments       ArangoDB edge documents
     */
    public static void updateEdges(ArrayList<ArrayList<Node>> tuplesArrayList, Map<String, OntologyElementMap> ontologyElementMaps,
                                   Map<String, Map<String, BaseEdgeDocument>> edgeDocuments) throws RuntimeException {

        Set<String> updatedEdges = new HashSet<>(); // For counting only
        System.out.println("Updating edges using " + tuplesArrayList.size() + " tuples");
        long startTime = System.nanoTime();
        for (ArrayList<Node> tupleArrayList : tuplesArrayList) {

            // Only update edges using quadruples
            if (tupleArrayList.size() != 4) continue;

            // Ensure the subject contains a valid ontology ID
            OntologyTripleLoader.VTuple s_vtuple = createVTuple(tupleArrayList.get(quadrupleSubjectIdx));
            if (!s_vtuple.isValidVertex()) continue;

            // Ensure the object contains a valid ontology ID
            OntologyTripleLoader.VTuple o_vtuple = createVTuple(tupleArrayList.get(quadrupleObjectIdx));
            if (!o_vtuple.isValidVertex()) continue;

            // Parse the predicate
            String attribute = parsePredicate(ontologyElementMaps, tupleArrayList.get(quadruplePredicateIdx));

            // Parse the literal
            String literal = tupleArrayList.get(quadrupleLiteralIdx).getLiteralValue().toString();

            // Update the corresponding edge
            String idPair = s_vtuple.id() + "-" + o_vtuple.id();
            String key = s_vtuple.number() + "-" + o_vtuple.number();
            if (!edgeDocuments.get(idPair).containsKey(key))
                throw new RuntimeException("Invalid edge in collection " + idPair + " with key " + key);
            updatedEdges.add(s_vtuple.id() + "/" + s_vtuple.number() + "-" + o_vtuple.id() + "/" + o_vtuple.number() + "-" + attribute + "-" + literal);
            BaseEdgeDocument doc = edgeDocuments.get(idPair).get(key);
            Set<String> literals;
            if (doc.getAttribute(attribute) == null) {
                literals = new HashSet<>();
                doc.addAttribute(attribute, literals);
            } else {
                literals = (HashSet<String>) doc.getAttribute(attribute);
            }
            literals.add(literal);
        }
        long stopTime = System.nanoTime();
        System.out.println("Updated " + updatedEdges.size() + " edges using " + tuplesArrayList.size() + " tuples in " + (stopTime - startTime) / 1e9 + " s");
    }

    /**
     * Load tuples parsed from a schema file into a local ArangoDB server instance.
     *
     * @param cellKnGraph         Cell-KN ArangoDB graph
     * @param ontologyElementMaps Maps terms and labels
     */
    public static void loadSchemaTuples(ArangoGraph cellKnGraph, Map<String, OntologyElementMap> ontologyElementMaps) {

        // Create the database and graph
        String cellKnSchemaDbName = "Cell-KN-Schema";
        arangoDbUtilities.deleteDatabase(cellKnSchemaDbName);
        ArangoDatabase cellKnSchemaDb = arangoDbUtilities.createOrGetDatabase(cellKnSchemaDbName);
        String cellKnSchemaGraphName = "KN-Schema-v0.7";
        arangoDbUtilities.deleteGraph(cellKnSchemaDb, cellKnSchemaGraphName);
        ArangoGraph cellKnSchemaGraph = arangoDbUtilities.createOrGetGraph(cellKnSchemaDb, cellKnSchemaGraphName);

        // Collect vertex keys for each vertex collection to prevent constructing
        // duplicate vertices in the vertex collection
        Map<String, Set<String>> vertexKeys = new HashMap<>();

        // Collect edge keys in each edge collection to prevent constructing duplicate
        // edges in the edge collection
        Map<String, Set<String>> edgeKeys = new HashMap<>();

        // Collect all vertices and edges before inserting them into the graph for improved performace
        Map<String, ArangoVertexCollection> vertexCollections = new HashMap<>();
        Map<String, Map<String, BaseDocument>> vertexDocuments = new HashMap<>();
        Map<String, ArangoEdgeCollection> edgeCollections = new HashMap<>();
        Map<String, Map<String, BaseEdgeDocument>> edgeDocuments = new HashMap<>();

        Path tuplesFile = schemaDir.resolve("cell-kn-schema-v0.7.0.json");
        System.out.println("Processing tuples file " + tuplesFile);

        // Read the tuples file
        ArrayList<ArrayList<Node>> tuplesArrayList;
        try {
            tuplesArrayList = readJsonFile(tuplesFile.toString());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // Construct vertices
        constructVertices(tuplesArrayList, cellKnSchemaGraph, vertexKeys, vertexCollections, vertexDocuments);

        // Construct edges
        constructEdges(tuplesArrayList, ontologyElementMaps, cellKnSchemaGraph, edgeKeys, edgeCollections, edgeDocuments);

        // Insert vertices, and edges
        insertVertices(vertexCollections, vertexDocuments);
        insertEdges(edgeCollections, edgeDocuments);

        // Replace Cell-KN schema vertices which have no attributes with Cell-KN vertices which do
        for (String id : vertexDocuments.keySet()) {
            for (String number : vertexDocuments.get(id).keySet()) {
                BaseDocument cellKnDoc = cellKnGraph.vertexCollection(id).getVertex(number, BaseDocument.class);
                if (cellKnDoc != null) {
                    vertexCollections.get(id).replaceVertex(number, cellKnDoc);
                }
            }
        }
    }

    /**
     * Load tuples parsed from a results file into a local ArangoDB server instance.
     *
     * @param args (None expected)
     */
    public static void main(String[] args) {

        // Identify the results tuples files
        String tuplesPath;
        if (args.length > 0) {
            tuplesPath = args[0];
        } else {
            tuplesPath = tuplesDir.toString();
        }
        String tuplesPattern = ".*\\.json";
        List<Path> tuplesFiles;
        try {
            tuplesFiles = listFilesMatchingPattern(tuplesPath, tuplesPattern);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        if (tuplesFiles.isEmpty()) {
            System.out.println("No tuples files found matching pattern " + tuplesPattern);
            System.exit(1);
        }

        // Map terms and labels
        String oboPath;
        if (args.length > 1) {
            oboPath = args[1];
        } else {
            oboPath = oboDir.toString();
        }
        String oboPattern = "ro.owl";
        List<Path> oboFiles;
        try {
            oboFiles = listFilesMatchingPattern(oboPath, oboPattern);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        Map<String, OntologyElementMap> ontologyElementMaps = null;
        if (oboFiles.isEmpty()) {
            System.out.println("No OBO files found matching pattern " + oboPattern);
            System.exit(2);
        } else {
            ontologyElementMaps = parseOntologyElements(oboFiles);
        }

        // Create the database and graph
        String cellKnDbName;
        if (args.length > 2) {
            cellKnDbName = args[2];
        } else {
            cellKnDbName = "Cell-KN-Ontologies";
        }
        // NEVER DO THIS: arangoDbUtilities.deleteDatabase(databaseName);
        ArangoDatabase cellKnDb = arangoDbUtilities.createOrGetDatabase(cellKnDbName);
        String cellKnGraphName;
        if (args.length > 3) {
            cellKnGraphName = args[3];
        } else {
            cellKnGraphName = "KN-Ontologies-v2.0";
        }
        // NEVER DO THIS: arangoDbUtilities.deleteGraph(db, graphName);
        ArangoGraph cellKnGraph = arangoDbUtilities.createOrGetGraph(cellKnDb, cellKnGraphName);

        // Collect vertex keys for each vertex collection to prevent constructing
        // duplicate vertices in the vertex collection
        Map<String, Set<String>> vertexKeys = new HashMap<>();

        // Collect edge keys in each edge collection to prevent constructing duplicate
        // edges in the edge collection
        Map<String, Set<String>> edgeKeys = new HashMap<>();

        // Collect all vertices and edges before inserting them into the graph for improved performace
        Map<String, ArangoVertexCollection> vertexCollections = new HashMap<>();
        Map<String, Map<String, BaseDocument>> vertexDocuments = new HashMap<>();
        Map<String, ArangoEdgeCollection> edgeCollections = new HashMap<>();
        Map<String, Map<String, BaseEdgeDocument>> edgeDocuments = new HashMap<>();

        // Read the results tuples files
        for (Path tuplesFile : tuplesFiles) {
            System.out.println("Processing tuples file " + tuplesFile);
            ArrayList<ArrayList<Node>> tuplesArrayList;
            try {
                tuplesArrayList = readJsonFile(tuplesFile.toString());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            // Construct, and update vertices
            constructVertices(tuplesArrayList, cellKnGraph, vertexKeys, vertexCollections, vertexDocuments);
            updateVertices(tuplesArrayList, ontologyElementMaps, vertexDocuments);

            // Construct, and update edges
            constructEdges(tuplesArrayList, ontologyElementMaps, cellKnGraph, edgeKeys, edgeCollections, edgeDocuments);
            updateEdges(tuplesArrayList, ontologyElementMaps, edgeDocuments);
        }
        // Insert vertices, and edges
        insertVertices(vertexCollections, vertexDocuments);
        insertEdges(edgeCollections, edgeDocuments);

        // Load the schema tuples file
        loadSchemaTuples(cellKnGraph, ontologyElementMaps);
    }
}
