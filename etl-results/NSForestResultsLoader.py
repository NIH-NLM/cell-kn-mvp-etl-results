import ast
from pathlib import Path

from rdflib.term import BNode, Literal, URIRef

from LoaderUtilities import load_results, hyphenate, PURLBASE, RDFSBASE


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
                URIRef(f"{RDFSBASE}/rdf/#Label"),
                Literal(cluster_name),
            )
        )

        # Cell_set_Ind, STATO:0000047 (count), clusterSize
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/Total_cell_count"),
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

        # TODO: Complete
        # Edge annotations

    return tuples


if __name__ == "__main__":

    nsforest_path = Path(
        "../data/cell-kn-mvp-nsforest-results-guo-2023-2025-02-22.csv"
    ).resolve()
    nsforest_results = load_results(nsforest_path).sort_values(
        "clusterName", ignore_index=True
    )

    nsforest_tuples = create_tuples_from_nsforest(nsforest_results)
