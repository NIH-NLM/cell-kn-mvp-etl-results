import ast
from pathlib import Path
import random
import string

import pandas as pd
from rdflib.term import BNode, Literal, URIRef


ALPHABET = string.ascii_lowercase + string.digits
PURLBASE = "http://purl.obolibrary.org/obo"
RDFSBASE = "http://www.w3.org/1999/02/22-rdf-syntax-ns"


def get_uuid():
    """Get an eight character random string.

    Parameters
    ----------
    None

    Returns
    -------
    An eight character random string.
    """
    return "".join(random.choices(ALPHABET, k=12))


def load_results(results_path):
    """Load NSForest results CSV file and append a UUID.

    Parameters
    ----------
    results_Path : Path
        Path of CSV file containing NSForest results

    Returns
    -------
    results : pd.DataFrame
        DataFrame containing NSForest results
    """
    results = pd.read_csv(results_path)
    if not "uuid" in results.columns:
        print(f"Add UUID column to results CSV file {results_path.name}")
        results["uuid"] = [get_uuid() for idx in results.index]
        results.to_csv(results_path)
    return results


def hyphenate(name):
    for c in [" ", "_", ","]:
        name = name.replace(c, "-").replace("--", "-")
    return name


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


def create_tuples_from_author_to_cl(results):
    tuples = []

    # Nodes for these results
    csd_term = f"CSD_{results['dataset_id'][0]}"
    pub_term = f"PUB_{results['DOI'][0]}"

    # Cell_set_dataset_Ind, SOURCE, Publication_Ind
    # IAO:0000100, dc:source, IAO:0000311
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}/dc#source"),
            URIRef(f"{PURLBASE}/{pub_term}"),
        )
    )

    # Cell_set_dataset_Ind, INSTANCE_OF, Cell_set_dataset_Class
    # IAO:0000100, rdf:type, IAO:0000100
    tuples.append(
        (
            URIRef(f"{PURLBASE}/CSD_{results['dataset_id'][0]}"),
            URIRef(f"{RDFSBASE}/rdf#type"),
            URIRef(f"{PURLBASE}/DS_{results['dataset_source'][0]}"),
        )
    )

    # Node annotations
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}#Dataset_version_id"),
            Literal(results["dataset_version_id"][0]),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}#Collection_id"),
            Literal(results["collection_id"][0]),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{csd_term}"),
            URIRef(f"{RDFSBASE}#Collection_version_id"),
            results["collection_version_id"][0],
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#PMID"),
            results["PMID"][0],
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#PMCID"),
            results["PMCID"][0],
        )
    )

    # Nodes for each cell type or cell set
    uuid_0 = results["uuid"][0]
    for _, row in results.iterrows():
        cl_term = row["cell_ontology_id"]
        uberon_term = row["uberon_entity_id"]
        author_cell_set = hyphenate(row["author_cell_set"])
        cs_term = f"CS_{author_cell_set}-{uuid_0}"

        # Cell_type_Class, PART_OF, Anatomical_structure_Class
        # CL:0000000, BFO:0000050, UBERON:0001062
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{PURLBASE}/BFO_0000050"),
                URIRef(f"{PURLBASE}/{uberon_term}"),
            )
        )

        # Cell_type_Class, HAS_EXEMPLAR_DATA, Cell_set_dataset_Ind
        # CL:0000000, RO:0015001, IAO:0000100
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{PURLBASE}/RO_0015001"),
                URIRef(f"{PURLBASE}/{csd_term}"),
            )
        )

        # Cell_set_Ind, DERIVES_FROM, Anatomical_structure_Ind
        # -, RO:0001000, UBERON:0001062
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/RO_0001000"),
                URIRef(f"{PURLBASE}/{uberon_term}-{uuid_0}"),
            )
        )

        # Cell_set_Ind, SOURCE, Cell_set_dataset_Ind
        # -, dc:source, IAO:0000100
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}/dc#source"),
                URIRef(f"{PURLBASE}/{csd_term}"),
            )
        )

        # Cell_set_Ind, COMPOSED_PRIMARILY_OF, Cell_type_Class
        # -, RO:0002473, CL:0000000
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/RO_0002473"),
                URIRef(f"{PURLBASE}/{cl_term}"),
            )
        )

        # Node annotations
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Author_cell_term"),
                Literal(row["author_cell_term"]),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{RDFSBASE}#Author_category"),
                Literal(row["author_category"]),
            )
        )

        # Edge annotations
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{RDFSBASE}#Match"),
                Literal(row["match"]),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{PURLBASE}/{cl_term}"),
                URIRef(f"{RDFSBASE}#Mapping_method"),
                Literal(row["mapping_method"]),
            )
        )

        # Nodes for each cell type and marker gene
        marker_genes = ast.literal_eval(row["NSForest_markers"])
        for gene in marker_genes:

            # Gene_Class, IS_MARKER_FOR, Cell_type_Class
            # SO:0000704, RO:0002607, CL:0000000
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/GS_{gene}"),
                    URIRef(f"{PURLBASE}/RO_0002607"),
                    URIRef(f"{PURLBASE}/{cl_term}"),
                )
            )

        # Nodes for each cell type, and marker and binary gene
        binary_genes = ast.literal_eval(row["binary_genes"])
        for gene in marker_genes + binary_genes:

            # Cell_type_Class, EXPRESSES, Gene_Class
            # CL:0000000, RO:0002292, SO:0000704
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{PURLBASE}/RO_0002292"),
                    URIRef(f"{PURLBASE}/GS_{gene}"),
                )
            )

            # Gene_Class, PART_OF, Cell_type_Class
            # SO:0000704, BFO:0000050, CL:0000000
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/GS_{gene}"),
                    URIRef(f"{PURLBASE}/BFO_0000050"),
                    URIRef(f"{PURLBASE}/{cl_term}"),
                )
            )

            # Gene_Class, EXPRESSED_IN, Anatomical_structure_Class
            # SO:0000704, RO:0002206, UBERON:0001062
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/GS_{gene}"),
                    URIRef(f"{PURLBASE}/RO_0002206"),
                    URIRef(f"{PURLBASE}/{uberon_term}"),
                )
            )

    return tuples


if __name__ == "__main__":

    nsforest_path = Path(
        "../data/cell-kn-mvp-nsforest-results-guo-2023-2025-02-22.csv"
    ).resolve()
    nsforest_results = load_results(nsforest_path).sort_values(
        "clusterName", ignore_index=True
    )
    nsforest_tuples = create_tuples_from_nsforest(nsforest_results)

    author_to_cl_path = Path(
        str(nsforest_path)
        .replace("nsforest-results", "map-author-to-cl")
        .replace("2023-2025-02-22", "2023-data-v0.4")
    )
    author_to_cl_results = load_results(author_to_cl_path).sort_values(
        "author_cell_set", ignore_index=True
    )
    author_to_cl_results = author_to_cl_results.merge(
        nsforest_results[["clusterName", "NSForest_markers", "binary_genes"]].copy(),
        left_on="author_cell_set",
        right_on="clusterName",
    )
    author_to_cl_tuples = create_tuples_from_author_to_cl(author_to_cl_results)
