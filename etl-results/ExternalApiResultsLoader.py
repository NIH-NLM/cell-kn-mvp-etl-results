from pathlib import Path

from ExternalApiResultsFetcher import get_opentargets_results, RESOURCES
from LoaderUtilities import hyphenate, PURLBASE, RDFSBASE


def create_tuples_from_opentargets(opentargets_path):

    tuples = []

    nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    for gene_id in opentargets_results["gene_ids"]:
        gene_term = gene_id.replace("ENSG", "ENSG_")

        for disease in opentargets_results[gene_id]["diseases"]:

            # == Disease relations

            # Gene, IS_GENETIC_BASIS_FOR_CONDITION, Disease
            tuples.append(
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#GENETIC_BASIS_FOR",
                    f"{PURLBASE}/{disease['id']}",
                )
            )

            # == Disease annotations

            tuples.extend(
                [
                    (
                        f"{PURLBASE}/{disease['id']}",
                        f"{RDFSBASE}#name",
                        disease["name"],
                    ),
                    (
                        f"{PURLBASE}/{disease['id']}",
                        f"{RDFSBASE}#description",
                        disease["description"],
                    ),
                ]
            )

            # == Gene to Disease edge annotation

            tuples.append(
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{PURLBASE}/{disease['id']}",
                    f"{RDFSBASE}#score",
                    disease["score"],
                )
            )

        for drug in opentargets_results[gene_id]["drugs"]:
            drug_term = drug["id"].replace("CHEMBL", "CHEMBL_")

            # == Drug_product relations

            # Gene, MOLECULARLY_INTERACTS_WITH, Drug_product
            tuples.append(
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#MOLECULARLY_INTERACTS_WITH",
                    f"{PURLBASE}/{drug_term}",
                )
            )

            # Drug_product, IS_SUBSTANCE_THAT_TREATS, Disease
            tuples.append(
                (
                    f"{PURLBASE}/{drug_term}",
                    f"{RDFSBASE}#IS_SUBSTANCE_THAT_TREATS",
                    f"{PURLBASE}/{drug['disease_id']}",
                )
            )

            for drug_trial_id in drug["trial_ids"]:
                drug_trial_term = drug_trial_id.replace("NCT", "NCT_")

                # == Clinical_trial relations

                # Drug_product, EVALUATED_IN, Clinical_trial
                tuples.append(
                    (
                        f"{PURLBASE}/{drug_term}",
                        f"{RDFSBASE}#EVALUATED_IN",
                        f"{PURLBASE}/{drug_trial_term}",
                    )
                )

                # == Clinical_trial annotations

                tuples.extend(
                    [
                        (
                            f"{PURLBASE}/{drug_trial_term}",
                            f"{RDFSBASE}#phase",
                            drug["trial_phase"],
                        ),
                        (
                            f"{PURLBASE}/{drug_trial_term}",
                            f"{RDFSBASE}#status",
                            drug["trial_status"],
                        ),
                    ]
                )
                # == Drug_product annotations

                tuples.extend(
                    [
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#name",
                            drug["name"],
                        ),
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#type",
                            drug["type"],
                        ),
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#mechanism of action",
                            drug["action_mechanism"],
                        ),
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#description",
                            drug["description"],
                        ),
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#synonyms",
                            drug["synonyms"],
                        ),
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#trade names",
                            drug["trade_names"],
                        ),
                        (
                            f"{PURLBASE}/{drug_term}",
                            f"{RDFSBASE}#approved",
                            drug["approved"],
                        ),
                    ]
                )

        for interaction in opentargets_results[gene_id]["interactions"]:

            # == Interaction relations

            # Gene, GENETICALLY_INTERACTS_WITH, Gene
            tuples.append(
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#GENETICALLY_INTERACTS_WITH",
                    f"{PURLBASE}/{interaction['gene_b_id']}",
                )
            )

            # == Gene to Interaction edge annotations

            tuples.extend(
                [
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/{interaction['gene_b_id']}",
                        f"{RDFSBASE}#evidence score",
                        interaction["evidence_score"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/{interaction['gene_b_id']}",
                        f"{RDFSBASE}#evidence count",
                        interaction["evidence_count"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/{interaction['gene_b_id']}",
                        f"{RDFSBASE}#source db",
                        interaction["source_db"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/{interaction['gene_b_id']}",
                        f"{RDFSBASE}#protein a",
                        interaction["protein_a_id"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/{interaction['gene_b_id']}",
                        f"{RDFSBASE}#protein b",
                        interaction["protein_b_id"],
                    ),
                ]
            )

        for pharmacogenetic in opentargets_results[gene_id]["pharmacogenetics"]:
            rs_term = pharmacogenetic["rs_id"].replace("rs", "rs_")

            # == Pharmacogenetic relations

            # Gene, HAS_QUALITY, Mutation
            tuples.append(
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#HAS_QUALITY",
                    f"{PURLBASE}/{rs_term}",
                )
            )

            # Mutation, INVOLVED_IN, Variant_consequence
            tuples.append(
                (
                    f"{PURLBASE}/{rs_term}",
                    f"{RDFSBASE}#INVOLVED_IN",
                    f"{PURLBASE}/{pharmacogenetic['variant_consequence_id']}",
                )
            )

            for pharmacogenetic_drug in pharmacogenetic["drugs"]:
                pharmacogenetic_drug_term = pharmacogenetic_drug["id"].replace(
                    "CHEMBL", "CHEMBL_"
                )

                # Mutation, HAS_PHARMACOLOGICAL_EFFECT, Drug_product
                tuples.append(
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#HAS_PHARMACOLOGICAL_EFFECT",
                        f"{PURLBASE}/{pharmacogenetic_drug_term}",
                    )
                )

            # == Pharmacogenetic annotations

            tuples.extend(
                [
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#genotype_id",
                        pharmacogenetic["genotype_id"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#genotype",
                        pharmacogenetic["genotype"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#phenotype",
                        pharmacogenetic["phenotype"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#genotype_annotation",
                        pharmacogenetic["genotype_annotation"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#response_category",
                        pharmacogenetic["response_category"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#evidence_level",
                        pharmacogenetic["evidence_level"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#source",
                        pharmacogenetic["source"],
                    ),
                    (
                        f"{PURLBASE}/{rs_term}",
                        f"{RDFSBASE}#literature",
                        pharmacogenetic["literature"],
                    ),
                ]
            )

            # == Variant_consequence annotations

            tuples.append(
                (
                    f"{PURLBASE}/{pharmacogenetic['variant_consequence_id']}",
                    f"{RDFSBASE}#variant_consequence_label",
                    pharmacogenetic["variant_consequence_label"],
                )
            )

        # == Tractability relations

        # None

        # == Gene annotations

        tuples.extend(
            [
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#label",
                    opentargets_results[gene_id]["tractability"]["label"],
                ),
                (
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#modality",
                    opentargets_results[gene_id]["tractability"]["modality"],
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
                    f"{PURLBASE}/{gene_term}",
                    f"{RDFSBASE}#EXPRESSED_IN",
                    f"{PURLBASE}/expression['tissue_id']",
                )
            )

            # == Gene to Anatomical_structure edge annotations

            tuples.extend(
                [
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/expression['tissue_id']",
                        f"{RDFSBASE}#rna_zscore",
                        expression["rna_zscore"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/expression['tissue_id']",
                        f"{RDFSBASE}#rna_value",
                        expression["rna_value"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/expression['tissue_id']",
                        f"{RDFSBASE}#rna_unit",
                        expression["rna_unit"],
                    ),
                    (
                        f"{PURLBASE}/{gene_term}",
                        f"{PURLBASE}/expression['tissue_id']",
                        f"{RDFSBASE}#rna_level",
                        expression["rna_level"],
                    ),
                ]
            )
