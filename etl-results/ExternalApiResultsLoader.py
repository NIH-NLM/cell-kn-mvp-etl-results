from pathlib import Path

from rdflib.term import Literal, URIRef

from ExternalApiResultsFetcher import (
    RESOURCES,
    get_gene_id_to_names_map,
    get_opentargets_results,
    get_uniprot_results,
    map_gene_id_to_names,
)
from LoaderUtilities import hyphenate, PURLBASE, RDFSBASE


def create_tuples_from_opentargets(opentargets_path):

    tuples = []

    nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)

    ensp2accn = uniprot_results["ensp2accn"]

    gid2nms = get_gene_id_to_names_map()
    for gene_id in opentargets_results["gene_ids"]:
        gene_symbol = map_gene_id_to_names(gene_id, gid2nms)
        gene_term = gene_id.replace("ENSG", "ENSG_")

        # == Gene relations

        for disease in opentargets_results[gene_id]["diseases"]:

            # == Disease relations

            # Gene, IS_GENETIC_BASIS_FOR_CONDITION, Disease
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#GENETIC_BASIS_FOR"),
                    URIRef(f"{PURLBASE}/{disease['id']}"),
                )
            )

            # == Disease annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{disease['id']}"),
                        URIRef(f"{RDFSBASE}#name"),
                        Literal(disease["name"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{disease['id']}"),
                        URIRef(f"{RDFSBASE}#description"),
                        Literal(disease["description"]),
                    ),
                ]
            )

            # == Gene to Disease edge annotation

            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{PURLBASE}/{disease['id']}"),
                    URIRef(f"{RDFSBASE}#score"),
                    Literal(disease["score"]),
                )
            )

        for drug in opentargets_results[gene_id]["drugs"]:
            drug_term = drug["id"].replace("CHEMBL", "CHEMBL_")

            # == Drug_product relations

            # Gene, MOLECULARLY_INTERACTS_WITH, Drug_product
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#MOLECULARLY_INTERACTS_WITH"),
                    URIRef(f"{PURLBASE}/{drug_term}"),
                )
            )

            # Drug_product, IS_SUBSTANCE_THAT_TREATS, Disease
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{drug_term}"),
                    URIRef(f"{RDFSBASE}#IS_SUBSTANCE_THAT_TREATS"),
                    URIRef(f"{PURLBASE}/{drug['disease_id']}"),
                )
            )

            for drug_trial_id in drug["trial_ids"]:
                drug_trial_term = drug_trial_id.replace("NCT", "NCT_")

                # == Clinical_trial relations

                # Drug_product, EVALUATED_IN, Clinical_trial
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{drug_term}"),
                        URIRef(f"{RDFSBASE}#EVALUATED_IN"),
                        URIRef(f"{PURLBASE}/{drug_trial_term}"),
                    )
                )

                # == Clinical_trial annotations

                tuples.extend(
                    [
                        (
                            URIRef(f"{PURLBASE}/{drug_trial_term}"),
                            URIRef(f"{RDFSBASE}#phase"),
                            Literal(drug["trial_phase"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_trial_term}"),
                            URIRef(f"{RDFSBASE}#status"),
                            Literal(drug["trial_status"]),
                        ),
                    ]
                )
                # == Drug_product annotations

                tuples.extend(
                    [
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#name"),
                            Literal(drug["name"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#type"),
                            Literal(drug["type"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#mechanism of action"),
                            Literal(drug["action_mechanism"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#description"),
                            Literal(drug["description"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#synonyms"),
                            Literal(drug["synonyms"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#trade names"),
                            Literal(drug["trade_names"]),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#approved"),
                            Literal(drug["approved"]),
                        ),
                    ]
                )

        for interaction in opentargets_results[gene_id]["interactions"]:

            # == Interaction relations

            # Gene, GENETICALLY_INTERACTS_WITH, Gene
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#GENETICALLY_INTERACTS_WITH"),
                    URIRef(f"{PURLBASE}/{interaction['gene_b_id']}"),
                )
            )

            # == Gene to Interaction edge annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{interaction['gene_b_id']}"),
                        URIRef(f"{RDFSBASE}#evidence score"),
                        Literal(interaction["evidence_score"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{interaction['gene_b_id']}"),
                        URIRef(f"{RDFSBASE}#evidence count"),
                        Literal(interaction["evidence_count"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{interaction['gene_b_id']}"),
                        URIRef(f"{RDFSBASE}#source db"),
                        Literal(interaction["source_db"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{interaction['gene_b_id']}"),
                        URIRef(f"{RDFSBASE}#protein a"),
                        Literal(interaction["protein_a_id"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{interaction['gene_b_id']}"),
                        URIRef(f"{RDFSBASE}#protein b"),
                        Literal(interaction["protein_b_id"]),
                    ),
                ]
            )

        for pharmacogenetic in opentargets_results[gene_id]["pharmacogenetics"]:
            rs_term = pharmacogenetic["rs_id"].replace("rs", "rs_")

            # == Pharmacogenetic relations

            # Gene, HAS_QUALITY, Mutation
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#HAS_QUALITY"),
                    URIRef(f"{PURLBASE}/{rs_term}"),
                )
            )

            # Mutation, INVOLVED_IN, Variant_consequence
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{rs_term}"),
                    URIRef(f"{RDFSBASE}#INVOLVED_IN"),
                    URIRef(f"{PURLBASE}/{pharmacogenetic['variant_consequence_id']}"),
                )
            )

            for pharmacogenetic_drug in pharmacogenetic["drugs"]:
                pharmacogenetic_drug_term = pharmacogenetic_drug["id"].replace(
                    "CHEMBL", "CHEMBL_"
                )

                # Mutation, HAS_PHARMACOLOGICAL_EFFECT, Drug_product
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#HAS_PHARMACOLOGICAL_EFFECT"),
                        URIRef(f"{PURLBASE}/{pharmacogenetic_drug_term}"),
                    )
                )

            # == Pharmacogenetic annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#genotype_id"),
                        Literal(pharmacogenetic["genotype_id"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#genotype"),
                        Literal(pharmacogenetic["genotype"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#phenotype"),
                        Literal(pharmacogenetic["phenotype"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#genotype_annotation"),
                        Literal(pharmacogenetic["genotype_annotation"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#response_category"),
                        Literal(pharmacogenetic["response_category"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#evidence_level"),
                        Literal(pharmacogenetic["evidence_level"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#source"),
                        Literal(pharmacogenetic["source"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#literature"),
                        Literal(pharmacogenetic["literature"]),
                    ),
                ]
            )

            # == Variant_consequence annotations

            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{pharmacogenetic['variant_consequence_id']}"),
                    URIRef(f"{RDFSBASE}#variant_consequence_label"),
                    Literal(pharmacogenetic["variant_consequence_label"]),
                )
            )

        # == Tractability relations

        # None

        # == Gene annotations

        tuples.extend(
            [
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#symbol"),
                    Literal(gene_symbol),
                ),
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#label"),
                    Literal(opentargets_results[gene_id]["tractability"]["label"]),
                ),
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#modality"),
                    Literal(opentargets_results[gene_id]["tractability"]["modality"]),
                ),
            ]
        )

        for expression in opentargets_results[gene_id]["expression"]:
            if expression["tissue_id"][0:7] != "UBERON_":
                continue

            # == Expression relations

            # Gene, EXPRESSED_IN, Anatomical_structure
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#EXPRESSED_IN"),
                    URIRef(f"{PURLBASE}/expression['tissue_id']"),
                )
            )

            # == Gene to Anatomical_structure edge annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/expression['tissue_id']"),
                        URIRef(f"{RDFSBASE}#rna_zscore"),
                        Literal(expression["rna_zscore"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/expression['tissue_id']"),
                        URIRef(f"{RDFSBASE}#rna_value"),
                        Literal(expression["rna_value"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/expression['tissue_id']"),
                        URIRef(f"{RDFSBASE}#rna_unit"),
                        Literal(expression["rna_unit"]),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/expression['tissue_id']"),
                        URIRef(f"{RDFSBASE}#rna_level"),
                        Literal(expression["rna_level"]),
                    ),
                ]
            )


def create_tuples_from_uniprot(opentargets_path):

    tuples = []

    uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)

    ensp2accn = uniprot_results["ensp2accn"]

    for protein_id in uniprot_results["protein_ids"]:
        pass
