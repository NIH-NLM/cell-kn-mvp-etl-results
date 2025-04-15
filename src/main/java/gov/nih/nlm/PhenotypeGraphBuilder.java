package gov.nih.nlm;

import java.util.*;

import com.arangodb.ArangoDatabase;
import com.arangodb.ArangoEdgeCollection;
import com.arangodb.ArangoGraph;
import com.arangodb.ArangoVertexCollection;
import com.arangodb.entity.BaseDocument;
import com.arangodb.entity.BaseEdgeDocument;
import com.arangodb.model.AqlQueryOptions;

/**
 * Queries the fully populated ontology ArangoDB to identify all paths that
 * connect cell set vertices inward to UBERON then NCBITaxon vertices, and
 * outward to gene then disease and drug vertices. Collects unique vertex and
 * edge documents, then inserts them in the phenotype ArangoDB.
 */
public class PhenotypeGraphBuilder {

	// Construct ArangoDB utilities
	private static final ArangoDbUtilities arangoDbUtilities = new ArangoDbUtilities();

	/**
	 * Query a fully populated ontology ArangoDB to identify all paths that connect
	 * cell set vertices inward to UBERON then NCBITaxon vertices, and outward to
	 * gene then disease and drug vertices.
	 *
	 * @param databaseName Name of database containing fully populated graph
	 * @param graphName    Name of fully populated graph
	 * @param limit        Number of cell sets to consider. If none, no limit is
	 *                     applied.
	 * @return All identified paths
	 */
	private static List<Map> getPaths(String databaseName, String graphName, int limit) {

		// Capture bind parameters
		Map<String, Object> bindVars = new HashMap<>();
		bindVars.put("graphName", graphName);

		// Capture query options
		AqlQueryOptions queryOpts = new AqlQueryOptions();

		// Construct AQL query string
		//@formatter:off
		String queryStr = "FOR cs IN CS";
		if (limit > 0) {
			bindVars.put("limit", limit);
			queryStr += " LIMIT @limit";
		}
		// Iterate through CS collection
		queryStr += " FOR v, e, p IN 3 ANY cs GRAPH @graphName";
		// Filter for paths with 4 vertices"
		queryStr += " FILTER LENGTH(p.vertices) == 4";
		// Find only paths with correct route - CL - UBERON - NCBITaxon or CL - GS - MONDO/CHEMBL or CL - CSD - PUB
		queryStr += " AND IS_SAME_COLLECTION('CL', p.vertices[1])"
				+ " AND (IS_SAME_COLLECTION('UBERON', p.vertices[2]) OR IS_SAME_COLLECTION('GS', p.vertices[2]) OR IS_SAME_COLLECTION('CSD', p.vertices[2]))"
				+ " AND (IS_SAME_COLLECTION('NCBITaxon', p.vertices[3]) OR IS_SAME_COLLECTION('MONDO', p.vertices[3]) OR IS_SAME_COLLECTION('CHEMBL', p.vertices[3]) OR IS_SAME_COLLECTION('PUB', p.vertices[3]))"
				+ " RETURN p";
		//@formatter:on

		// Execute query
		ArangoDatabase db = arangoDbUtilities.createOrGetDatabase(databaseName);
		System.out.println("Quering a fully populated ontology ArangoDB to identify paths");
		List<Map> paths = db.query(queryStr, Map.class, bindVars, queryOpts).asListRemaining();
		return paths;
	}

	/**
	 * Collect unique vertex documents from all identified paths.
	 *
	 * @param paths All identified paths
	 * @return Unique vertex documents
	 */
	private static List<BaseDocument> getVertexDocuments(List<Map> paths) {
		System.out.println("Collecting unique vertex documents from all identified paths");
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
		return vertexDocuments;
	}

	/**
	 * Collect unique edge documents from all identified paths.
	 *
	 * @param paths All identified paths
	 * @return Unique edge documents
	 */
	private static List<BaseEdgeDocument> getEdgeDocuments(List<Map> paths) {
		System.out.println("Collecting unique edge documents from all identified paths");
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
		return edgeDocuments;
	}

	/**
	 * Insert unique vertex documents.
	 *
	 * @param vertexDocuments Unique vertex documents
	 * @param graph           Graph in phenotype database
	 */
	private static void insertVertexDocuments(List<BaseDocument> vertexDocuments, ArangoGraph graph) {
		System.out.println("Inserting " + vertexDocuments.size() + " vertex documents");
		long startTime = System.nanoTime();
		Map<String, ArangoVertexCollection> vertexCollections = new HashMap<>();
		for (BaseDocument vertexDocument : vertexDocuments) {
			String id = vertexDocument.getId().substring(0, vertexDocument.getId().indexOf("/"));
			if (!vertexCollections.containsKey(id)) {
				vertexCollections.put(id, arangoDbUtilities.createOrGetVertexCollection(graph, id));
			}
			vertexCollections.get(id).insertVertex(vertexDocument);
		}
		long stopTime = System.nanoTime();
		System.out.println("Inserted " + vertexDocuments.size() + " vertex documents in " + (stopTime - startTime) / 1e9 + " s");
	}

	/**
	 * Insert unique edge documents.
	 *
	 * @param edgeDocuments Unique edge documents
	 * @param graph         Graph in phenotype database
	 */
	private static void insertEdgeDocuments(List<BaseEdgeDocument> edgeDocuments, ArangoGraph graph) {
		System.out.println("Inserting " + edgeDocuments.size() + " edge documents");
		long startTime = System.nanoTime();
		Map<String, ArangoEdgeCollection> edgeCollections = new HashMap<>();
		for (BaseEdgeDocument edgeDocument : edgeDocuments) {
			String idPair = edgeDocument.getId().substring(0, edgeDocument.getId().indexOf("/"));
			String idFrom = idPair.substring(0, idPair.indexOf("-"));
			String idTo = idPair.substring(idPair.indexOf("-") + 1);
			if (!edgeCollections.containsKey(idPair)) {
				edgeCollections.put(idPair, arangoDbUtilities.createOrGetEdgeCollection(graph, idFrom, idTo));
			}
			edgeCollections.get(idPair).insertEdge(edgeDocument);
		}
		long stopTime = System.nanoTime();
		System.out.println("Inserted " + edgeDocuments.size() + " edge documents in " + (stopTime - startTime) / 1e9 + " s");
	}

	/**
	 * Query the fully populated ontology ArangoDB to identify all paths that
	 * connect cell set vertices inward to UBERON then NCBITaxon vertices, and
	 * outward to gene then disease and drug vertices. Collect unique vertex and
	 * edge documents, then insert them in the phenotype ArangoDB.
	 */
	public static void main(String[] args) {

		// Get all phenotype database subgraph paths in the ontology database and fully
		// populated graph
		String ontologyDatabaseName = "Cell-KN-Ontologies";
		String ontologyGraphName = "KN-Ontologies-v2.0";
		int limit = 0;
		List<Map> paths = getPaths(ontologyDatabaseName, ontologyGraphName, limit);

		// Initialize the phenotype database and subgraph
		String phenotypeDatabaseName = "Cell-KN-Phenotypes";
		String phenotypeGraphName = "KN-Phenotypes-v2.0";
		arangoDbUtilities.deleteDatabase(phenotypeDatabaseName);
		ArangoDatabase db = arangoDbUtilities.createOrGetDatabase(phenotypeDatabaseName);
		arangoDbUtilities.deleteGraph(db, phenotypeGraphName);
		ArangoGraph graph = arangoDbUtilities.createOrGetGraph(db, phenotypeGraphName);

		// Get unique vertex and edge documents in the ontology database and fully
		// populated graph
		List<BaseDocument> vertexDocuments = getVertexDocuments(paths);
		List<BaseEdgeDocument> edgeDocuments = getEdgeDocuments(paths);

		// Insert unique vertex and edge documents in the phenotype database and
		// subgraph
		insertVertexDocuments(vertexDocuments, graph);
		insertEdgeDocuments(edgeDocuments, graph);
	}
}
