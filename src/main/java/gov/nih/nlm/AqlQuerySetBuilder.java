package gov.nih.nlm;

import com.arangodb.ArangoDatabase;
import com.arangodb.model.AqlQueryOptions;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Defines sets of AQL query strings and bind variables for obtaining paths between specified nodes, including paths
 * outbound from the last node in ontologies.
 */
public class AqlQuerySetBuilder {

    // Construct ArangoDB utilities
    private static final ArangoDbUtilities arangoDbUtilities = new ArangoDbUtilities();

    /**
     * Get AQL query set to identify a path with one edge.
     *
     * @param graph Graph name
     * @param node  Node name
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInOne(String graph, String node) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("node", node);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 1 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@node, p.vertices[1])
                RETURN p
                """.formatted(graph, node);
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with two edges.
     *
     * @param graph   Graph name
     * @param nodeOne Node one name
     * @param nodeTwo Node two name
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInTwo(String graph, String nodeOne, String nodeTwo) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 2 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                RETURN p
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with two edges, including paths outbound from the last node in ontologies.
     *
     * @param graph          Graph name
     * @param nodeOne        Node one name
     * @param nodeTwo        Node two name
     * @param edgeCollection Edge collection name, usually "nodeTwo-nodeTwo"
     * @param edgeLabel      Edge label, usuall "SUB_CLASS_OF", or "PART_OF"
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInTwoWithHierarchy(String graph,
                                                            String nodeOne,
                                                            String nodeTwo,
                                                            String edgeCollection,
                                                            String edgeLabel) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        bindVars.put("edgeCollection", edgeCollection);
        bindVars.put("edgeLabel", edgeLabel);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 2 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                    LET l = FIRST(
                      FOR n1, e1, p1 IN 1..64 OUTBOUND p.vertices[2] @edgeCollection
                      PRUNE e1 != null AND e1.Label NOT IN [@edgeLabel]
                      FILTER p1.edges[*].Label ALL IN [@edgeLabel]
                      SORT LENGTH(p1.edges) DESC
                      LIMIT 1
                      RETURN p1
                    )
                RETURN {
                  vertices: FLATTEN(
                    [
                      p.vertices,
                      l ? l.vertices : []
                    ]
                  ),
                  edges: FLATTEN(
                    [
                      p.edges,
                      l ? l.edges : []
                    ]
                  )
                }
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with three edges.
     *
     * @param graph     Graph name
     * @param nodeOne   Node one name
     * @param nodeTwo   Node two name
     * @param nodeThree Node three name
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInThree(String graph, String nodeOne, String nodeTwo, String nodeThree) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        bindVars.put("nodeThree", nodeThree);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 3 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                    AND
                    IS_SAME_COLLECTION(@nodeThree, p.vertices[3])
                RETURN p
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with three edges, including paths outbound from the last node in
     * ontologies.
     *
     * @param graph          Graph name
     * @param nodeOne        Node one name
     * @param nodeTwo        Node two name
     * @param nodeThree      Node three name
     * @param edgeCollection Edge collection name, usually "nodeTwo-nodeTwo"
     * @param edgeLabel      Edge label, usuall "SUB_CLASS_OF", or "PART_OF"
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInThreeWithHierarchy(String graph,
                                                              String nodeOne,
                                                              String nodeTwo,
                                                              String nodeThree,
                                                              String edgeCollection,
                                                              String edgeLabel) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        bindVars.put("nodeThree", nodeThree);
        bindVars.put("edgeCollection", edgeCollection);
        bindVars.put("edgeLabel", edgeLabel);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 3 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                    AND
                    IS_SAME_COLLECTION(@nodeThree, p.vertices[3])
                    LET l = FIRST(
                      FOR n1, e1, p1 IN 1..64 OUTBOUND p.vertices[3] @edgeCollection
                        PRUNE e1 != null AND e1.Label NOT IN [@edgeLabel]
                        FILTER p1.edges[*].Label ALL IN [@edgeLabel]
                        SORT LENGTH(p1.edges) DESC
                        LIMIT 1
                        RETURN p1
                      )
                RETURN {
                  vertices: FLATTEN(
                    [
                      p.vertices,
                      l ? l.vertices : []
                    ]
                  ),
                  edges: FLATTEN(
                    [
                      p.edges,
                      l ? l.edges : []
                    ]
                  )
                }
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with four edges.
     *
     * @param graph     Graph name
     * @param nodeOne   Node one name
     * @param nodeTwo   Node two name
     * @param nodeThree Node three name
     * @param nodeFour  Node four name
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInFour(String graph,
                                                String nodeOne,
                                                String nodeTwo,
                                                String nodeThree,
                                                String nodeFour) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        bindVars.put("nodeThree", nodeThree);
        bindVars.put("nodeFour", nodeFour);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 4 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                    AND
                    IS_SAME_COLLECTION(@nodeThree, p.vertices[3])
                    AND
                    IS_SAME_COLLECTION(@nodeFour, p.vertices[4])
                RETURN p
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with four edges, including paths outbound from the last node in ontologies.
     *
     * @param graph          Graph name
     * @param nodeOne        Node one name
     * @param nodeTwo        Node two name
     * @param nodeThree      Node three name
     * @param nodeFour       Node four name
     * @param edgeCollection Edge collection name, usually "nodeTwo-nodeTwo"
     * @param edgeLabel      Edge label, usuall "SUB_CLASS_OF", or "PART_OF"
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInFourWithHeirarchy(String graph,
                                                             String nodeOne,
                                                             String nodeTwo,
                                                             String nodeThree,
                                                             String nodeFour,
                                                             String edgeCollection,
                                                             String edgeLabel) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        bindVars.put("nodeThree", nodeThree);
        bindVars.put("nodeFour", nodeFour);
        bindVars.put("edgeCollection", edgeCollection);
        bindVars.put("edgeLabel", edgeLabel);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 4 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                    AND
                    IS_SAME_COLLECTION(@nodeThree, p.vertices[3])
                    AND
                    IS_SAME_COLLECTION(@nodeFour, p.vertices[4])
                    LET l = FIRST(
                      FOR n1, e1, p1 IN 1..64 OUTBOUND p.vertices[4] @edgeCollection
                        PRUNE e1 != null AND e1.Label NOT IN [@edgeLabel]
                        FILTER p1.edges[*].Label ALL IN [@edgeLabel]
                        SORT LENGTH(p1.edges) DESC
                        LIMIT 1
                        RETURN p1
                    )
                RETURN {
                  vertices: FLATTEN(
                    [
                      p.vertices,
                      l ? l.vertices : []
                    ]
                  ),
                  edges: FLATTEN(
                    [
                      p.edges,
                      l ? l.edges : []
                    ]
                  )
                }
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Get AQL query set to identify a path with five edges.
     *
     * @param graph     Graph name
     * @param nodeOne   Node one name
     * @param nodeTwo   Node two name
     * @param nodeThree Node three name
     * @param nodeFour  Node four name
     * @param nodeFive  Node five name
     * @return Query set
     */
    public static AqlQuerySet getQuerySetInFive(String graph,
                                                String nodeOne,
                                                String nodeTwo,
                                                String nodeThree,
                                                String nodeFour,
                                                String nodeFive) {
        Map<String, Object> bindVars = new HashMap<>();
        bindVars.put("graph", graph);
        bindVars.put("nodeOne", nodeOne);
        bindVars.put("nodeTwo", nodeTwo);
        bindVars.put("nodeThree", nodeThree);
        bindVars.put("nodeFour", nodeFour);
        bindVars.put("nodeFive", nodeFive);
        String queryStr = """
                FOR cs IN CS
                  FOR n, e, p IN 5 ANY cs GRAPH @graph
                    FILTER
                    IS_SAME_COLLECTION(@nodeOne, p.vertices[1])
                    AND
                    IS_SAME_COLLECTION(@nodeTwo, p.vertices[2])
                    AND
                    IS_SAME_COLLECTION(@nodeThree, p.vertices[3])
                    AND
                    IS_SAME_COLLECTION(@nodeFour, p.vertices[4])
                    AND
                    IS_SAME_COLLECTION(@nodeFive, p.vertices[5])
                RETURN p
                """;
        return new AqlQuerySet(bindVars, queryStr);
    }

    /**
     * Run a few example queries.
     *
     * @param args None expected
     */
    public static void main(String[] args) {

        String database = "Cell-KN-Ontologies";
        String graph = "KN-Ontologies-v2.0";
        ArangoDatabase db = arangoDbUtilities.createOrGetDatabase(database);

        List<AqlQuerySet> aqlQuerySets = new ArrayList<>();

        // CS - BGS
        String node = "BGS";
        aqlQuerySets.add(getQuerySetInOne(graph, node));

        // CS - CL - GO
        String nodeOne = "CL";
        String nodeTwo = "GO";
        aqlQuerySets.add(getQuerySetInTwo(graph, nodeOne, nodeTwo));

        // CS - CL - GO
        nodeOne = "CL";
        nodeTwo = "GO";
        String edgeCollection = "GO-GO";
        String edgeLabel = "SUB_CLASS_OF";
        aqlQuerySets.add(getQuerySetInTwoWithHierarchy(graph, nodeOne, nodeTwo, edgeCollection, edgeLabel));

        // CS - CL - GS - PR
        nodeOne = "CL";
        nodeTwo = "GS";
        String nodeThree = "PR";
        aqlQuerySets.add(getQuerySetInThree(graph, nodeOne, nodeTwo, nodeThree));

        // CS - CL - GS - MONDO
        nodeOne = "CL";
        nodeTwo = "GS";
        nodeThree = "MONDO";
        edgeCollection = "MONDO-MONDO";
        edgeLabel = "SUB_CLASS_OF";
        aqlQuerySets.add(getQuerySetInThreeWithHierarchy(graph,
                nodeOne,
                nodeTwo,
                nodeThree,
                edgeCollection,
                edgeLabel));

        // CS - CL - GS - MONDO - CHEMBL
        nodeOne = "CL";
        nodeTwo = "GS";
        nodeThree = "MONDO";
        String nodeFour = "CHEMBL";
        aqlQuerySets.add(getQuerySetInFour(graph, nodeOne, nodeTwo, nodeThree, nodeFour));

        // CS - CL - GS - MONDO - HP
        nodeOne = "CL";
        nodeTwo = "GS";
        nodeThree = "MONDO";
        nodeFour = "HP";
        edgeCollection = "HP-HP";
        edgeLabel = "SUB_CLASS_OF";
        aqlQuerySets.add(getQuerySetInFourWithHeirarchy(graph,
                nodeOne,
                nodeTwo,
                nodeThree,
                nodeFour,
                edgeCollection,
                edgeLabel));

        // CS - CL - GS - MONDO - CHEMBL - PR
        nodeOne = "CL";
        nodeTwo = "GS";
        nodeThree = "MONDO";
        nodeFour = "CHEMBL";
        String nodeFive = "PR";
        aqlQuerySets.add(getQuerySetInFive(graph, nodeOne, nodeTwo, nodeThree, nodeFour, nodeFive));

        AqlQueryOptions queryOpts = new AqlQueryOptions();
        for (AqlQuerySet aqlQuerySet : aqlQuerySets) {
            System.out.println(aqlQuerySet.queryStr().lines().collect(Collectors.joining()).replaceAll("\\s+", " "));
            long startTime = System.nanoTime();
            List<Map> queryPaths = db.query(aqlQuerySet.queryStr(),
                    Map.class,
                    aqlQuerySet.bindVars(),
                    queryOpts).asListRemaining();
            long stopTime = System.nanoTime();
            System.out.println("Collected " + queryPaths.size() + " paths in " + (stopTime - startTime) / 1e9 + " s");
        }

        arangoDbUtilities.arangoDB.shutdown();
    }

    public record AqlQuerySet(Map<String, Object> bindVars, String queryStr) {
    }
}
