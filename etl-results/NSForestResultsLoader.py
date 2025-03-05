import argparse
import ast
from pathlib import Path

from rdflib.term import BNode, Literal, URIRef

import ArangoDbUtilities as adb
from LoaderUtilities import load_results, hyphenate, PURLBASE, RDFSBASE
from OntologyParserLoader import load_tuples_into_adb_graph, parse_obo, VALID_VERTICES

OBO_DIRPATH = Path("../data/obo")
NSFOREST_DIRPATH = Path("../data/results")


def create_tuples_from_nsforest(results):
    tuples = []
    for _, row in results.iterrows():
        uuid = row["uuid"]
        cluster_name = hyphenate(row["clusterName"])
        cluster_size = row["clusterSize"]
        binary_genes = ast.literal_eval(row["binary_genes"])
        nsforest_markers = ast.literal_eval(row["NSForest_markers"])
        cs_term = f"CS_{cluster_name}-{uuid}"
        bmc_term = f"BMC_{uuid}-NSF"
        bgc_term = f"BMC_{uuid}-BG"

        # TODO: Remove Biomarker_combination_Class
        # Biomarker_combination_Ind, INSTANCE_OF, Biomarker_combination_Class
        # ---, rdf:type, SO:0001260
        # tuples.append(
        #     (
        #         URIRef(f"{PURLBASE}/{bmc_term}"),
        #         URIRef(f"{RDFSBASE}/rdf#type"),
        #         URIRef(f"{PURLBASE}/SO_0001260"),
        #     )
        # )

        # Gene_Class, PART_OF, Biomarker_combination_Ind
        # SO:0000704, BFO:0000050, SO:0001260
        for gene in nsforest_markers:
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/GS_{gene}"),
                    URIRef(f"{PURLBASE}/BFO_0000050"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                )
            )

        # Cell_set_Ind, HAS_CHARACTERIZING_MARKER_SET, Biomarker_combination_Ind
        # ---, RO:0015004, SO:0001260
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/RO_0015004"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
            )
        )

        # Biomarker_combination_Ind, SUBCLUSTER_OF, Binary_gene_combination_Ind
        # SO:0001260, RO:0015003, SO:0001260
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{PURLBASE}/RO_0015003"),
                URIRef(f"{PURLBASE}/{bgc_term}"),
            )
        )

        # Node annotations

        # Cell_set_Ind, rdfs:label, clusterName
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Label"),
                Literal(cluster_name),
            )
        )

        # Cell_set_Ind, STATO:0000047 (count), clusterSize
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Total_cell_count"),
                Literal(str(cluster_size)),
            )
        )

        # Cell_set_Ind, -, binary_genes
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Binary_genes"),
                Literal(" ".join(binary_genes)),
            )
        )

        # Cell_set_Ind, RO:0015004 (has characterizing marker set), NSForest_markers
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Markers"),
                Literal(" ".join(nsforest_markers)),
            )
        )

        # Edge annotations for BMC terms
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{PURLBASE}/#Source_algorithm"),  # [IAO_0000064]
                Literal("NSForest-v4.0_dev"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#F_beta_confidence_score"),  # [STAT:0000663]
                Literal(str(row["f_score"])),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#PPV"),  # [STAT:0000416]
                Literal(str(row["PPV"])),
            )
        )
        # TODO: Restore when available in data
        # tuples.append(
        #     (
        #         URIRef(f"{PURLBASE}/{cs_term}"),
        #         URIRef(f"{PURLBASE}/{bmc_term}"),
        #         URIRef(f"{RDFSBASE}#Recall"),  # [STAT:0000233]
        #         Literal(str(row["recall"])),
        #     )
        # )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#TN"),  # [STAT:0000597]
                Literal(str(row["TN"])),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#TP"),  # [STAT:0000595]
                Literal(str(row["TP"])),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#FN"),  # [STAT:0000598]
                Literal(str(row["FN"])),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#FP"),  # [STAT:0000596]
                Literal(str(row["FP"])),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#Marker_count"),  # [STAT:0000047]
                Literal(str(row["marker_count"])),
            )
        )

        # Edge annotations for BGC terms
        # TODO: Restore when available in data
        # tuples.append(
        #     (
        #         URIRef(f"{PURLBASE}/{cs_term}"),
        #         URIRef(f"{PURLBASE}/{bgc_term}"),
        #         URIRef(f"{RDFSBASE}#On_target"),  # [STAT:0000047]
        #         Literal(str(row["onTarget"])),
        #     )
        # )

    return tuples


def main(parameters=None):

    parser = argparse.ArgumentParser(description="Load NSForest results")
    group = parser.add_argument_group(
        "Cell Ontology (CL)", "Version of the CL assumed loaded"
    )
    exclusive_group = group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument(
        "--test", action="store_true", help="assume the test ontology loaded"
    )
    exclusive_group.add_argument(
        "--full", action="store_true", help="assume the full ontology loaded"
    )
    parser.add_argument(
        "--label",
        default="",
        help="label to add to database_name",
    )

    if parameters is None:
        args = parser.parse_args()

    else:
        args = parser.parse_args(parameters)

    if args.test:
        db_name = "Cell-KN-v1.5"
        graph_name = "CL-Test"

    if args.full:
        db_name = "Cell-KN-v1.5"
        graph_name = "CL-Full"

    if args.label:
        db_name += f"-{args.label}"

    ro_filename = "ro.owl"
    log_filename = f"{graph_name}.log"

    print("Parse the relationship ontology")
    ro, _, _ = parse_obo(OBO_DIRPATH, ro_filename)

    print("Getting ArangoDB database and graph, and loading tuples")
    db = adb.create_or_get_database(db_name)
    adb_graph = adb.create_or_get_graph(db, graph_name)
    vertex_collections = {}
    edge_collections = {}

    nsforest_path = (
        NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-guo-2023-2025-02-22.csv"
    ).resolve()

    nsforest_results = load_results(nsforest_path).sort_values(
        "clusterName", ignore_index=True
    )

    nsforest_tuples = create_tuples_from_nsforest(nsforest_results)

    with open("NSForestResultsLoader.out", "w") as f:
        for tuple in nsforest_tuples:
            f.write(str(tuple) + "\n")

    VALID_VERTICES.add("BMC")
    VALID_VERTICES.add("CS")
    VALID_VERTICES.add("GS")

    load_tuples_into_adb_graph(
        nsforest_tuples,
        adb_graph,
        vertex_collections,
        edge_collections,
        ro=ro,
        do_update=True,
    )


if __name__ == "__main__":
    main()
