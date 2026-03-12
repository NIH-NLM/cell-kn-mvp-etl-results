package gov.nih.nlm;

import com.arangodb.ArangoDatabase;
import com.arangodb.ArangoEdgeCollection;
import com.arangodb.ArangoGraph;
import com.arangodb.ArangoVertexCollection;
import com.arangodb.entity.BaseDocument;
import com.arangodb.entity.BaseEdgeDocument;
import com.arangodb.model.AqlQueryOptions;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static gov.nih.nlm.AqlQuerySetBuilder.AqlQuerySet;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInFive;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInFour;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInFourWithHierarchy;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInOne;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInThree;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInThreeWithHierarchy;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInTwo;
import static gov.nih.nlm.AqlQuerySetBuilder.getQuerySetInTwoWithHierarchy;
import static gov.nih.nlm.OntologyGraphBuilder.getDocumentCollectionName;

/**
 * Queries the fully populated ontology ArangoDB to identify all paths that connect cell set vertices inward to UBERON
 * then NCBITaxon vertices, and outward to gene then disease and drug vertices. The longest path outbound from each
 * UBERON, NCBITaxon, and MONDO nodes are included. Collects unique vertex and edge documents, then inserts them in the
 * phenotype ArangoDB.
 */
public class PhenotypeGraphBuilder {

    // Construct ArangoDB utilities
    private static final ArangoDbUtilities arangoDbUtilities = new ArangoDbUtilities();

    /**
     * Query a fully populated ontology ArangoDB to identify all paths that connect cell set vertices inward to UBERON
     * then NCBITaxon vertices, and outward to gene then disease and drug vertices.
     *
     * @param databaseName Name of database containing fully populated graph
     * @param graphName    Name of fully populated graph
     * @return All identified paths
     */
    private static List<Map> getPaths(String databaseName, String graphName) {

        List<AqlQuerySet> aqlQuerySets = new ArrayList<>();

        aqlQuerySets.add(getQuerySetInOne(graphName, "BGS"));

        aqlQuerySets.add(getQuerySetInTwo(graphName, "BMC", "BGS"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "CL", "CSD"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "CL", "GS"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "CL", "PR"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "CSD", "PUB"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "UBERON", "CHEBI"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "UBERON", "CSD"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "UBERON", "GS"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "UBERON", "NCBITaxon"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "UBERON", "PATO"));
        aqlQuerySets.add(getQuerySetInTwo(graphName, "UBERON", "PR"));

        aqlQuerySets.add(getQuerySetInTwoWithHierarchy(graphName,
                "CL",
                "NCBITaxon",
                "NCBITaxon-NCBITaxon",
                "SUB_CLASS_OF"));
        aqlQuerySets.add(getQuerySetInTwoWithHierarchy(graphName, "CL", "PATO", "PATO-PATO", "SUB_CLASS_OF"));
        aqlQuerySets.add(getQuerySetInTwoWithHierarchy(graphName, "CL", "UBERON", "UBERON-UBERON", "PART_OF"));
        aqlQuerySets.add(getQuerySetInTwoWithHierarchy(graphName, "UBERON", "GO", "GO-GO", "SUB_CLASS_OF"));

        aqlQuerySets.add(getQuerySetInThree(graphName, "CL", "GO", "NCBITaxon"));
        aqlQuerySets.add(getQuerySetInThree(graphName, "CL", "GS", "BMC"));
        aqlQuerySets.add(getQuerySetInThree(graphName, "CL", "GS", "UBERON"));

        aqlQuerySets.add(getQuerySetInThreeWithHierarchy(graphName,
                "CL",
                "GS",
                "MONDO",
                "MONDO-MONDO",
                "SUB_CLASS_OF"));

        aqlQuerySets.add(getQuerySetInFour(graphName, "CL", "GS", "PR", "CHEMBL"));
        aqlQuerySets.add(getQuerySetInFour(graphName, "CL", "GS", "MONDO", "NCBITaxon"));

        aqlQuerySets.add(getQuerySetInFourWithHierarchy(graphName, "CL", "GS", "MONDO", "HP", "HP-HP", "SUB_CLASS_OF"));

        aqlQuerySets.add(getQuerySetInFive(graphName, "CL", "GS", "RS", "CHEMBL", "MONDO"));
        aqlQuerySets.add(getQuerySetInFive(graphName, "CL", "GS", "RS", "CHEMBL", "PR"));

        ArangoDatabase db = arangoDbUtilities.createOrGetDatabase(databaseName);
        AqlQueryOptions queryOpts = new AqlQueryOptions();
        List<Map> paths = new ArrayList<>();
        for (AqlQuerySet aqlQuerySet : aqlQuerySets) {
            System.out.println(aqlQuerySet.queryStr().lines().collect(Collectors.joining()).replaceAll("\\s+", " "));
            long startTime = System.nanoTime();
            List<Map> queryPaths = db.query(aqlQuerySet.queryStr(),
                    Map.class,
                    aqlQuerySet.bindVars(),
                    queryOpts).asListRemaining();
            paths.addAll(queryPaths);
            long stopTime = System.nanoTime();
            System.out.println("Collected " + queryPaths.size() + " paths in " + (stopTime - startTime) / 1e9 + " s");
        }
        return paths;
    }

    /**
     * Collect unique vertex documents from all identified paths.
     *
     * @param paths All identified paths
     * @return Unique vertex documents
     */
    private static List<BaseDocument> getVertexDocuments(List<Map> paths) {
        System.out.println("Collecting unique vertex documents from " + paths.size() + " identified paths");
        long startTime = System.nanoTime();
        List<BaseDocument> vertexDocuments = new ArrayList<>();
        for (Map path : paths) {
            ArrayList<LinkedHashMap> vertices = (ArrayList<LinkedHashMap>) path.get("vertices");
            for (LinkedHashMap vertex : vertices) {
                BaseDocument vertexDoc = new BaseDocument(vertex);
                if (!vertexDocuments.contains(vertexDoc)) {
                    vertexDocuments.add(vertexDoc);
                }
            }
        }
        long stopTime = System.nanoTime();
        System.out.println("Collected " + vertexDocuments.size() + " unique vertex documents from " + paths.size() + " identified paths in " + (stopTime - startTime) / 1e9 + " s");
        return vertexDocuments;
    }

    /**
     * Collect unique edge documents from all identified paths.
     *
     * @param paths All identified paths
     * @return Unique edge documents
     */
    private static List<BaseEdgeDocument> getEdgeDocuments(List<Map> paths) {
        System.out.println("Collecting unique edge documents from " + paths.size() + " identified paths");
        long startTime = System.nanoTime();
        List<BaseEdgeDocument> edgeDocuments = new ArrayList<>();
        for (Map path : paths) {
            ArrayList<LinkedHashMap> edges = (ArrayList<LinkedHashMap>) path.get("edges");
            for (LinkedHashMap edge : edges) {
                BaseEdgeDocument edgeDoc = new BaseEdgeDocument(edge);
                if (!edgeDocuments.contains(edgeDoc)) {
                    edgeDocuments.add(edgeDoc);
                }
            }
        }
        long stopTime = System.nanoTime();
        System.out.println("Collected " + edgeDocuments.size() + " unique edge documents from " + paths.size() + " identified paths in " + (stopTime - startTime) / 1e9 + " s");
        return edgeDocuments;
    }

    /**
     * Insert unique vertex documents.
     *
     * @param phenotypeVertexDocuments Unique vertex documents
     * @param phenotypeGraph           Graph in phenotype database
     * @param ontologyGraph            Graph in ontology database
     */
    private static void insertVertexDocuments(List<BaseDocument> phenotypeVertexDocuments,
                                              ArangoGraph phenotypeGraph,
                                              ArangoGraph ontologyGraph) {
        System.out.println("Inserting " + phenotypeVertexDocuments.size() + " vertex documents");
        long startTime = System.nanoTime();
        Map<String, ArangoVertexCollection> phenotypeVertexCollections = new HashMap<>();
        for (BaseDocument phenotypeVertexDocument : phenotypeVertexDocuments) {
            String id = getDocumentCollectionName(phenotypeVertexDocument.getId());
            String key = phenotypeVertexDocument.getKey();
            if (!phenotypeVertexCollections.containsKey(id)) {
                phenotypeVertexCollections.put(id, arangoDbUtilities.createOrGetVertexCollection(phenotypeGraph, id));
            }
            BaseDocument ontologyVertexDocument = ontologyGraph.vertexCollection(id).getVertex(key, BaseDocument.class);
            if (phenotypeVertexCollections.get(id).getVertex(key, BaseDocument.class) == null) {
                if (ontologyVertexDocument != null) {
                    phenotypeVertexCollections.get(id).insertVertex(ontologyVertexDocument);
                } else {
                    phenotypeVertexCollections.get(id).insertVertex(phenotypeVertexDocument);
                }
            } else {
                if (ontologyVertexDocument != null) {
                    phenotypeVertexCollections.get(id).replaceVertex(key, ontologyVertexDocument);
                } else {
                    phenotypeVertexCollections.get(id).replaceVertex(key, phenotypeVertexDocument);
                }
            }
        }
        long stopTime = System.nanoTime();
        System.out.println("Inserted " + phenotypeVertexDocuments.size() + " vertex documents in " + (stopTime - startTime) / 1e9 + " s");
    }

    /**
     * Insert unique edge documents.
     *
     * @param phenotypeEdgeDocuments Unique edge documents
     * @param phenotypeGraph         Graph in phenotype database
     */
    private static void insertEdgeDocuments(List<BaseEdgeDocument> phenotypeEdgeDocuments, ArangoGraph phenotypeGraph) {
        System.out.println("Inserting " + phenotypeEdgeDocuments.size() + " edge documents");
        long startTime = System.nanoTime();
        Map<String, ArangoEdgeCollection> edgeCollections = new HashMap<>();
        for (BaseEdgeDocument edgeDocument : phenotypeEdgeDocuments) {
            String idPair = getDocumentCollectionName(edgeDocument.getId());
            String key = edgeDocument.getKey();
            String idFrom = getDocumentCollectionName(edgeDocument.getFrom());
            String idTo = getDocumentCollectionName(edgeDocument.getTo());
            if (!edgeCollections.containsKey(idPair)) {
                edgeCollections.put(idPair, arangoDbUtilities.createOrGetEdgeCollection(phenotypeGraph, idFrom, idTo));
            }
            if (edgeCollections.get(idPair).getEdge(key, BaseEdgeDocument.class) == null) {
                edgeCollections.get(idPair).insertEdge(edgeDocument);
            } else {
                edgeCollections.get(idPair).replaceEdge(key, edgeDocument);
            }
        }
        long stopTime = System.nanoTime();
        System.out.println("Inserted " + phenotypeEdgeDocuments.size() + " edge documents in " + (stopTime - startTime) / 1e9 + " s");
    }

    /**
     * Query the fully populated ontology ArangoDB to identify all paths that connect cell set vertices inward to UBERON
     * then NCBITaxon vertices, and outward to gene then disease and drug vertices. Collect unique vertex and edge
     * documents, then insert them in the phenotype ArangoDB.
     */
    public static void main(String[] args) {

        // Get all phenotype database subgraph paths in the ontology database and fully
        // populated graph
        String ontologyDatabaseName = "Cell-KN-Ontologies";
        String ontologyGraphName = "KN-Ontologies-v2.0";
        ArangoDatabase ontologyDb = arangoDbUtilities.createOrGetDatabase(ontologyDatabaseName);
        ArangoGraph ontologyGraph = arangoDbUtilities.createOrGetGraph(ontologyDb, ontologyGraphName);
        List<Map> paths = getPaths(ontologyDatabaseName, ontologyGraphName);

        // Initialize the phenotype database and subgraph
        String phenotypeDatabaseName = "Cell-KN-Phenotypes";
        String phenotypeGraphName = "KN-Phenotypes-v2.0";
        arangoDbUtilities.deleteDatabase(phenotypeDatabaseName);
        ArangoDatabase phenotypeDb = arangoDbUtilities.createOrGetDatabase(phenotypeDatabaseName);
        arangoDbUtilities.deleteGraph(phenotypeDb, phenotypeGraphName);
        ArangoGraph phenotypeGraph = arangoDbUtilities.createOrGetGraph(phenotypeDb, phenotypeGraphName);

        // Get vertex documents in the ontology graph and insert them in the phenotype graph
        List<BaseDocument> phenotypeVertexDocuments = getVertexDocuments(paths);
        insertVertexDocuments(phenotypeVertexDocuments, phenotypeGraph, ontologyGraph);

        // Get edge documents in the ontology graph and insert them in the phenotype graph
        List<BaseEdgeDocument> phenotypeEdgeDocuments = getEdgeDocuments(paths);
        insertEdgeDocuments(phenotypeEdgeDocuments, phenotypeGraph);

        // Disconnect from a local ArangoDB server instance
        arangoDbUtilities.arangoDB.shutdown();
    }
}
