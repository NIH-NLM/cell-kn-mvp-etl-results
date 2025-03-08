from pathlib import Path

from ExternalApiResultsFetcher import get_opentargets_results, RESOURCES
from LoaderUtilities import hyphenate, PURLBASE, RDFSBASE


def create_tuples_from_opentargets(opentargets_path):
    
    nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    for gene_id in opentargets_results["gene_ids"]:
        gene_term = gene_id.replace("ENSG", "ENSG_")

        diseases = opentargets_results[gene_id]["diseases"]
        for disease in diseases:
            disease_term = disease["id"]
            disease_name = disease["name"]
            disease_description = disease["description"]
            disease_score = disease["score"]
            
            # == Disease relations

            # Gene, IS_GENETIC_BASIS_FOR_CONDITION, Disease
            # ensg_id, GENETIC_BASIS_FOR, disease_term
            f"{PURLBASE}/{gene_term}"
            f"{RDFSBASE}#GENETIC_BASIS_FOR"
            f"{PURLBASE}/{disease_term}"

            # == Disease annotations

            f"{PURLBASE}/{disease_term}"
            f"{RDFSBASE}#name"
            disease_name

            f"{PURLBASE}/{disease_term}"
            f"{RDFSBASE}#description"
            disease_description

        drugs = opentargets_results[gene_id]["drugs"]
        for drug in drugs:
            # TODO: Other terms to transform?
            drug_term = drug['id'].replace("CHEMBL", "CHEMBL_")
            drug_name = drug['name']
            drug_type = drug['type']
            drug_action_mechanism = drug['action_mechanism']
            drug_description = drug['description']
            drug_synonyms = drug['synonyms']
            drug_trade_names = drug['trade_names']
            drug_disease_term = drug['disease_id']
            # drug_disease_name = drug['disease_name']
            drug_trial_phase = drug['trial_phase']
            drug_trial_status = drug['trial_status']
            drug_trial_ids = drug['trial_ids']
            drug_approved = drug['approved']

            # == Drug_product relations

            # Gene, MOLECULARLY_INTERACTS_WITH, Drug_product
            # ensg_id, MOLECULARLY_INTERACTS_WITH, chembl_id (id)
            f"{PURLBASE}/{gene_term}"
            f"{RDFSBASE}#MOLECULARLY_INTERACTS_WITH"
            f"{PURLBASE}/{drug_term}"
            
            # Drug_product, IS_SUBSTANCE_THAT_TREATS, Disease
            # chembl_id, IS_SUBSTANCE_THAT_TREATS, disease_id
            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#IS_SUBSTANCE_THAT_TREATS"
            f"{PURLBASE}/{drug_disease_term}"
            
            for drug_trial_id in drug_trial_ids:
                drug_trial_term = drug_trial_id.replace("NCT", "NCT_")

                # == Clinical_trial relations

                # Drug_product, EVALUATED_IN (in-house relation), Clinical_trial (OPMI:0004507)
                # chembl_id, EVALUATED_IN, trial_ids
                f"{PURLBASE}/{drug_term}"
                f"{RDFSBASE}#EVALUATED_IN"
                f"{PURLBASE}/{drug_trial_term}"

                # == Clinical_trial annotations

                f"{PURLBASE}/{drug_trial_term}"
                f"{RDFSBASE}#phase"
                drug_trial_phase

                f"{PURLBASE}/{drug_trial_term}"
                f"{RDFSBASE}#status"
                drug_trial_status

            # == Drug_product annotations

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#name"
            drug_name

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#type"
            drug_type

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#mechanism of action"
            drug_action_mechanism

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#description"
            drug_description

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#synonyms"
            drug_synonyms

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#trade names"
            drug_trade_names

            f"{PURLBASE}/{drug_term}"
            f"{RDFSBASE}#approved"
            drug_approved

        interactions = opentargets_results[gene_id]["interactions"]
        for interaction in interactions:
            interaction_gene_term = interaction['gene_b_id']
            interaction_evidence_score = interaction['evidence_score']
            interaction_evidence_count = interaction['evidence_count']
            interaction_source_db = interaction["source_db"]
            interaction_protein_a_id = interaction['protein_a_id']
            interaction_protein_b_id = interaction['protein_b_id']

            # == Interaction relations

            # Gene, GENETICALLY_INTERACTS_WITH (RO:0002435), Gene
            # ensg_id, GENETICALLY_INTERACTS_WITH (RO:0002435), ensg_id

            f"{PURLBASE}/{gene_term}"
            f"{RDFSBASE}#GENETICALLY_INTERACTS_WITH"
            f"{PURLBASE}/{interaction_gene_term}"

            # == Interaction relations edge annotations

            f"{PURLBASE}/{gene_term}"
            f"{PURLBASE}/{interaction_gene_term}"
            f"{RDFSBASE}#evidence score"
            interaction_evidence_score

            f"{PURLBASE}/{gene_term}"
            f"{PURLBASE}/{interaction_gene_term}"
            f"{RDFSBASE}#evidence count"
            interaction_evidence_count

            f"{PURLBASE}/{gene_term}"
            f"{PURLBASE}/{interaction_gene_term}"
            f"{RDFSBASE}#source db"
            interaction_source_db

            f"{PURLBASE}/{gene_term}"
            f"{PURLBASE}/{interaction_gene_term}"
            f"{RDFSBASE}#protein a"
            interaction_protein_a_id

            f"{PURLBASE}/{gene_term}"
            f"{PURLBASE}/{interaction_gene_term}"
            f"{RDFSBASE}#protein b"
            interaction_protein_b_id

        pharmacogenetics = opentargets_results[gene_id]["pharmacogenetics"]
        for pharmacogenetic in pharmacogenetics:
            pharmacogenetic_rs_term = pharmacogenetic["rs_id"].replace("rs", "rs_")
            pharmacogenetic_genotype_id = pharmacogenetic["genotype_id"]
            pharmacogenetic_genotype = pharmacogenetic["genotype"]
            pharmacogenetic_variant_consequence_term = pharmacogenetics["variant_consequence_id"]
            pharmacogenetic_variant_consequence_label = pharmacogenetics["variant_consequence_label"]
            pharmacogenetic_drugs = pharmacogenetics["drugs"]
            pharmacogenetic_phenotype = pharmacogenetic["phenotype"]
            pharmacogenetic_genotype_annotation = pharmacogenetic["genotype_annotation"]
            pharmacogenetic_response_category = pharmacogentic["response_category"]
            pharmacogenetic_evidence_level = pharmacogentic["evidence_level"]
            pharmacogenetic_source = pharmacogentic["source"]
            pharmacogenetic_literature = pharmacogentic["literature"]

            # == Pharmacogenetic relations

            # Gene, HAS_QUALITY, Mutation (SO:0001060)
            # ensg_id, has_rs_id, rs_id
            f"{PURLBASE}/{gene_term}"
            f"{RDFSBASE}#HAS_QUALITY"
            f"{PURLBASE}/{pharmacogenetic_rs_term}"
            
            # Mutation, INVOLVED_IN, Variant_consequence (SO:0001059)
            # rs_id, has_variant_consequence_id, variant_consequence_id
            f"{PURLBASE}/{pharmacogenetic_rs_term}"
            f"{RDFSBASE}#INVOLVED_IN"
            f"{PURLBASE}/{pharmacogenetic_variant_consequence_term}"

            for pharmacogenetic_drug in pharmacogenetic_drugs:
                pharmacogenetic_drug_term = pharmacogenetic_drug["id"].replace("CHEMBL", "CHEMBL_")

                # Mutation, HAS_PHARMACOLOGICAL_EFFECT, Drug_ID
                # rs_id, has_drugs, drugs_list
                f"{PURLBASE}/{pharmacogenetic_rs_term}"
                f"{RDFSBASE}#HAS_PHARMACOLOGICAL_EFFECT"
                f"{PURLBASE}/{pharmacogenetic_drug_term}"


            # == Pharmacogenetic annotations

    # genotype_id
    # genotype
    # phenotype_description
    # genotype_annotation
    # response_category
    # evidence_level
    # source_id
    # literature_reference_id

    # == Variant_consequence annotations

    # variant_consequence_label



    # == Tractability relations

    # None

    # Gene annotations

    # label
    # modality



    # == Expression relations

    # Gene, EXPRESSED_IN, Anatomical_structure
    # ensg_id, is_expressed, tissue_id

    # == Gene annotationss

    # rna_expression_in_anatomical_systems
    # rna_expression_in_organs

    # == Anatomical_structure annotations

    # tissue_name



    # == Edge annotations

    # ensg_id, disease_id, disease_score

    # ensg_id, tissue_id, tissue_id_rna_expression_zscore
    # ensg_id, tissue_id, tissue_id_rna_expression_value
    # ensg_id, tissue_id, tissue_id_rna_expression_unit
    # ensg_id, tissue_id, tissue_id_rna_expression_unit

