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
     * @param limit        Number of cell sets to consider. If none, no limit is applied.
     * @return All identified paths
     */
    private static List<Map> getPaths(String databaseName, String graphName, int limit) {

        Map<String, Object> bindVars;
        List<Map<String, Object>> bindVarMaps = new ArrayList<>();
        List<String> queryStrings = new ArrayList<>();
        AqlQueryOptions queryOpts = new AqlQueryOptions();

        //@formatter:off

        // Consider only CS nodes
        String queryPrefix = "FOR cs IN CS ";

        // Path CS-BGS-BMC always created using NSForest results
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 2 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('BGS', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('BMC', p.vertices[2]) "
                        + "RETURN p"
        );

        // Path CS-BMC-GS always created using NSForest results
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 2 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('BMC', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('GS', p.vertices[2]) "
                        + "RETURN p"
        );

        // Path CS-BMC-CL always created using author to CL mapping results
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 2 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('BMC', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('CL', p.vertices[2]) "
                        + "RETURN p"
        );

        // Path CL-UBERON-NCBITaxon always created using author to CL mapping results,
        // and ontology contents, path UBERON-UBERON PART_OF, and path NCBITaxon-NCBITaxon SUB_CLASS_OF
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVars.put("uberonEdgeCollection", "UBERON-UBERON");
        bindVars.put("ncbiTaxonEdgeCollection", "NCBITaxon-NCBITaxon");
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 3 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('CL', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('UBERON', p.vertices[2]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('NCBITaxon', p.vertices[3]) "
                        + "LET lP1 = FIRST("
                        + "FOR n1, e1, p1 "
                        + "IN 1..64 "
                        + "OUTBOUND p.vertices[2] "
                        + "@uberonEdgeCollection "
                        + "PRUNE e1 != null AND e1.Label NOT IN [\"PART_OF\", [\"PART_OF\"]] "
                        + "FILTER p1.edges[*].Label ALL IN [\"PART_OF\", [\"PART_OF\"]] "
                        + "SORT LENGTH(p1.edges) DESC "
                        + "LIMIT 1 "
                        + "RETURN p1"
                        + ") "
                        + "LET lP2 = FIRST("
                        + "FOR n2, e2, p2 "
                        + "IN 1..64 "
                        + "OUTBOUND p.vertices[3] "
                        + "@ncbiTaxonEdgeCollection "
                        + "PRUNE e2 != null AND e2.Label NOT IN [\"SUB_CLASS_OF\", [\"SUB_CLASS_OF\"]] "
                        + "FILTER p2.edges[*].Label ALL IN [\"SUB_CLASS_OF\", [\"SUB_CLASS_OF\"]] "
                        + "SORT LENGTH(p2.edges) DESC "
                        + "LIMIT 1 "
                        + "RETURN p2"
                        + ") "
                        + "RETURN {"
                        + "vertices: FLATTEN(["
                        + "p ? p.vertices : [], "
                        + "lP1 ? lP1.vertices : [], "
                        + "lP2 ? lP2.vertices : []"
                        + "]), "
                        + "edges: FLATTEN(["
                        + "p ? p.edges : [], "
                        + "lP1 ? lP1.edges : [], "
                        + "lP2 ? lP2.edges : []"
                        + "])"
                        + "}"
        );

        // Path CL-CSD-PUB always created using author to CL mapping results
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 3 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('CL', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('CSD', p.vertices[2]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('PUB', p.vertices[3]) "
                        + "RETURN p"
        );

        // Path CL-GS always created using NSForest and author to CL mapping results
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 2 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('CL', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('GS', p.vertices[2]) "
                        + "RETURN p"
        );

        // Path CL-GS always created using NSForest and author to CL mapping results,
        // however CL-GS-MONDO may not always exist, and path MONDO-MONDO SUB_CLASS_OF
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVars.put("mondoEdgeCollection", "MONDO-MONDO");
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 3 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('CL', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('GS', p.vertices[2]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('MONDO', p.vertices[3]) "
                        + "LET lP1 = FIRST("
                        + "FOR n1, e1, p1 "
                        + "IN 1..64 "
                        + "OUTBOUND p.vertices[3] "
                        + "@mondoEdgeCollection "
                        + "PRUNE e1 != null AND e1.Label NOT IN [\"SUB_CLASS_OF\", [\"SUB_CLASS_OF\"]] "
                        + "FILTER p1.edges[*].Label ALL IN [\"SUB_CLASS_OF\", [\"SUB_CLASS_OF\"]] "
                        + "SORT LENGTH(p1.edges) DESC "
                        + "LIMIT 1 "
                        + "RETURN p1"
                        + ") "
                        + "RETURN {"
                        + "vertices: FLATTEN(["
                        + "p ? p.vertices : [], "
                        + "lP1 ? lP1.vertices : []"
                        + "]), "
                        + "edges: FLATTEN(["
                        + "p ? p.edges : [], "
                        + "lP1 ? lP1.edges : []"
                        + "])"
                        + "}"
        );

        // Path CL-GS-PR always created using NSForest, author to CL mapping, and external API results,
        // however CL-GS-PR-CHEMBL may not always exist
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 4 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('CL', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('GS', p.vertices[2]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('PR', p.vertices[3]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('CHEMBL', p.vertices[4]) "
                        + "RETURN p"
        );

        // Path CL-GS-PR always created using NSForest, author to CL mapping, and external API results,
        // however CL-GS-PR-CHEMBL-MONDO may not always exist, and path MONDO-MONDO SUB_CLASS_OF
        bindVars = new HashMap<>();
        if (limit > 0) {
            bindVars.put("limit", limit);
            queryPrefix += "LIMIT @limit ";
        }
        bindVars.put("graphName", graphName);
        bindVars.put("mondoEdgeCollection", "MONDO-MONDO");
                bindVarMaps.add(bindVars);
        queryStrings.add(
                queryPrefix
                        + "FOR v, e, p IN 5 ANY cs GRAPH @graphName "
                        + "FILTER "
                        + "IS_SAME_COLLECTION('CL', p.vertices[1]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('GS', p.vertices[2]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('PR', p.vertices[3]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('CHEMBL', p.vertices[4]) "
                        + "AND "
                        + "IS_SAME_COLLECTION('MONDO', p.vertices[5]) "
                        + "LET lP1 = FIRST("
                        + "FOR n1, e1, p1 "
                        + "IN 1..64 "
                        + "OUTBOUND p.vertices[5] "
                        + "@mondoEdgeCollection "
                        + "PRUNE e1 != null AND e1.Label NOT IN [\"SUB_CLASS_OF\", [\"SUB_CLASS_OF\"]] "
                        + "FILTER p1.edges[*].Label ALL IN [\"SUB_CLASS_OF\", [\"SUB_CLASS_OF\"]] "
                        + "SORT LENGTH(p1.edges) DESC "
                        + "LIMIT 1 "
                        + "RETURN p1"
                        + ") "
                        + "RETURN {"
                        + "vertices: FLATTEN(["
                        + "p ? p.vertices : [], "
                        + "lP1 ? lP1.vertices : []"
                        + "]), "
                        + "edges: FLATTEN(["
                        + "p ? p.edges : [], "
                        + "lP1 ? lP1.edges : []"
                        + "])"
                        + "} "
        );
        //@formatter:on

        // Execute query
        ArangoDatabase db = arangoDbUtilities.createOrGetDatabase(databaseName);
        System.out.println("Quering a fully populated ontology ArangoDB to identify paths");
        List<Map> paths = new ArrayList<>();
        for (int queryIdx = 0; queryIdx < queryStrings.size(); queryIdx++) {
            System.out.println("Query: " + queryStrings.get(queryIdx));
            List<Map> queryPaths = db.query(queryStrings.get(queryIdx),
                    Map.class,
                    bindVarMaps.get(queryIdx),
                    queryOpts).asListRemaining();
            paths.addAll(queryPaths);
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
     * @param phenotypeGraph            Graph in phenotype database
     * @param ontologyGraph             Graph in ontology database
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
        int limit = 0;
        List<Map> paths = getPaths(ontologyDatabaseName, ontologyGraphName, limit);

        // Initialize the phenotype database and subgraph
        String phenotypeDatabaseName = "Cell-KN-Phenotypes";
        String phenotypeGraphName = "KN-Phenotypes-v2.0";
        ArangoDatabase phenotypeDb = arangoDbUtilities.createOrGetDatabase(phenotypeDatabaseName);
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
