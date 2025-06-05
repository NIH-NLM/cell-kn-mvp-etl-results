import ast
from glob import glob
import json
from pathlib import Path

import pandas as pd
from rdflib.term import Literal, URIRef

from LoaderUtilities import load_results, hyphenate, PURLBASE, RDFSBASE

NSFOREST_DIRPATH = Path("../../../data/results")
TUPLES_DIRPATH = Path("../../../data/tuples")


def create_tuples_from_nsforest(results):
    """Creates tuples from NSForest results consistent with schema
    v0.7.

    Parameters
    ----------
    results : pd.DataFrame
        DataFrame containing NSForest results

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    """
    tuples = []

    # Nodes for each cell set, marker and binary genes, and cell type
    for _, row in results.iterrows():
        uuid = row["uuid"]
        cluster_name = hyphenate(row["clusterName"])
        cluster_size = row["clusterSize"]
        binary_genes = ast.literal_eval(row["binary_genes"])
        nsforest_markers = ast.literal_eval(row["NSForest_markers"])
        cs_term = f"CS_{cluster_name}-{uuid}"
        bmc_term = f"BMC_{uuid}-NSF"
        bgc_term = f"BMC_{uuid}-BG"

        # Biomarker_combination_Ind, INSTANCE_OF, Biomarker_combination_Class
        # ---, rdf:type, SO:0001260
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}/rdf#type"),
                URIRef(f"{PURLBASE}/SO_0001260"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{PURLBASE}/SO_0001260"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("NSForest"),
            )
        )

        # Gene_Class, PART_OF, Biomarker_combination_Ind
        # SO:0000704, BFO:0000050, SO:0001260
        for gene in nsforest_markers:
            # TODO: Use gs_term?
            gs_term = f"GS_{gene}"
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/BFO_0000050"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("NSForest"),
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
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("NSForest"),
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
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{PURLBASE}/{bgc_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("NSForest"),
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

        # Biomarker_combination_Ind, RO:0015004 (has characterizing marker set), NSForest_markers
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{bmc_term}"),
                URIRef(f"{RDFSBASE}#Markers"),
                Literal(" ".join(nsforest_markers)),
            )
        )

        # Edge annotations for BMC terms
        tuples.extend(
            [
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{PURLBASE}/#Source_algorithm"),  # [IAO_0000064]
                    Literal("NSForest-v4.0_dev"),
                ),
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#F_beta_confidence_score"),  # [STAT:0000663]
                    Literal(str(row["f_score"])),
                ),
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#PPV"),  # [STAT:0000416]
                    Literal(str(row["PPV"])),
                ),
            ]
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
        tuples.extend(
            [
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#TN"),  # [STAT:0000597]
                    Literal(str(row["TN"])),
                ),
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#TP"),  # [STAT:0000595]
                    Literal(str(row["TP"])),
                ),
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#FN"),  # [STAT:0000598]
                    Literal(str(row["FN"])),
                ),
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#FP"),  # [STAT:0000596]
                    Literal(str(row["FP"])),
                ),
                (
                    URIRef(f"{PURLBASE}/{cs_term}"),
                    URIRef(f"{PURLBASE}/{bmc_term}"),
                    URIRef(f"{RDFSBASE}#Marker_count"),  # [STAT:0000047]
                    Literal(str(row["marker_count"])),
                ),
            ]
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


def main(summarize=False):
    """Load NSForest results from processing datasets corresponding to
    the Guo et al. 2023, Li et al. 2023, and Sikkema, et al. 2023
    publications, create tuples consistent with schema v0.7, and write
    the result to a JSON file. If summarizing, retain the first row
    only, and include results in output.

    Parameters
    ----------
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    None
    """
    nsforest_paths = [
        Path(p).resolve()
        for p in glob(str(NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-*.csv"))
    ]
    for nsforest_path in nsforest_paths:

        # Load NSForest results
        nsforest_results = load_results(nsforest_path).sort_values(
            "clusterName", ignore_index=True
        )
        if summarize:
            nsforest_results = nsforest_results.head(1)

        print(f"Creating tuples from {nsforest_path}")
        nsforest_tuples = create_tuples_from_nsforest(nsforest_results)
        if summarize:
            output_dirpath = TUPLES_DIRPATH / "summaries"
        else:
            output_dirpath = TUPLES_DIRPATH
        with open(
            output_dirpath / nsforest_path.name.replace(".csv", ".json"), "w"
        ) as f:
            data = {}
            if summarize:
                data["results"] = nsforest_results.to_dict()
            data["tuples"] = nsforest_tuples
            json.dump(data, f, indent=4)

        if summarize:
            break


if __name__ == "__main__":
    main(summarize=True)
    main()
