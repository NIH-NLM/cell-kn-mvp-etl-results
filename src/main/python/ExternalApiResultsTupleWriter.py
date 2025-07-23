from glob import glob
import json
from pathlib import Path
import re

from rdflib.term import Literal, URIRef

from ExternalApiResultsFetcher import (
    RESOURCES,
    NSFOREST_DIRPATH,
    HUBMAP_DIRPATH,
    get_gene_results,
    get_opentargets_results,
    get_uniprot_results,
)
from LoaderUtilities import (
    PURLBASE,
    RDFSBASE,
    get_chembl_to_pubchem_map,
    get_efo_to_mondo_map,
    get_gene_ensembl_id_to_names_map,
    load_results,
    map_chembl_to_pubchem,
    map_efo_to_mondo,
    map_gene_ensembl_id_to_names,
    map_protein_ensembl_id_to_accession,
)


NSFOREST_DIRPATH = Path("../../../data/results")
TUPLES_DIRPATH = Path("../../../data/tuples")


def get_mondo_term(disease_id, efo2mondo):
    """Return MONDO term, mapping from EFO term when necessary

    Parameters
    ----------
    disease_id : str
        Disease id which equals a MONDO or EFO term
    efo2mondo : pd.DataFrame
        DataFrame indexed by EFO containing MONDO term

    Return
    ------
    mondo_term : str
        Disease MONDO term
    """
    mondo_term = None

    if "MONDO" in disease_id:
        mondo_term = disease_id

    elif "EFO" in disease_id:
        mondo_term = map_efo_to_mondo(disease_id, efo2mondo)

    return mondo_term


def get_protein_term(protein_id, ensp2accn):
    """Map protein id to term by mapping Ensembl ids to UniProt
    accessions, if needed, and following the term naming convention
    for parsing.

    Parameters
    ----------
    protein_id : str
        Protein id provided by gget opentargets command
    ensp2accn : dict
        Mapping of Ensembl id to UniProt accession

    Returns
    -------
    protein_term : str
    """
    protein_term = None

    if "ENSP" in protein_id:
        accession = map_protein_ensembl_id_to_accession(protein_id, ensp2accn)

    else:
        accession = protein_id

    if accession is not None:
        protein_term = f"PR_{accession}"

    return protein_term


def create_tuples_from_opentargets(opentargets_path, summarize=False):
    """Creates tuples from the result of using the gget opentargets
    command to obtain resources for each unique gene id mapped from
    each gene symbol from corresponding NSForest results. If
    summnarizing, retain the first gene id only.

    Parameters
    ----------
    opentargets_path : Path
        Path to opentargets results
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    results : dict
        Dictionary containg opentargets results keyed by gene id, then
        by resource
    """
    tuples = []

    # Load the opentargets results
    nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))
    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    # Load the Gene results to obtain the UniProt names corresponding
    # to gene symbols
    _gene_path, gene_results = get_gene_results(nsforest_path)

    # Load mappings
    gid2nms = get_gene_ensembl_id_to_names_map()
    efo2mondo = get_efo_to_mondo_map()
    chembl2pubchem = get_chembl_to_pubchem_map()

    # Assign gene ids to consider
    if summarize:

        # Find a gene id with all resources, and a valid disease and interaction
        for gene_ensembl_id in opentargets_results["gene_ensembl_ids"]:

            # Find a gene id for which all resources are not empty
            is_empty = False
            for resource in RESOURCES:
                if len(opentargets_results[gene_ensembl_id][resource]) == 0:
                    is_empty = True
                    break
            if is_empty:
                continue

            # Find a valid disease
            found_disease = False
            for disease in opentargets_results[gene_ensembl_id]["diseases"]:
                if "MONDO" in disease["disease"]["id"] and disease["score"] > 0.5:
                    found_disease = True
                    break

            # Find a valid interaction
            found_interaction = False
            for interaction in opentargets_results[gene_ensembl_id]["interactions"]:
                if (
                    interaction["evidences"]
                    and interaction["evidences"][0]["evidenceScore"]
                    and interaction["evidences"][0]["evidenceScore"] > 0.5
                ):
                    found_interaction = True
                    break

            if found_disease and found_interaction:
                break

        # Consider selected gene id
        gene_ensembl_ids = [gene_ensembl_id]
        results = {}
        results[gene_ensembl_id] = opentargets_results[gene_ensembl_id]
        results[gene_ensembl_id]["symbol"] = map_gene_ensembl_id_to_names(
            gene_ensembl_id, gid2nms
        )[0]

        # Retain only the first result for each resource
        for resource in RESOURCES:
            if resource == "diseases":
                results[gene_ensembl_id][resource] = [disease]

            elif resource == "interactions":
                results[gene_ensembl_id][resource] = [interaction]

            else:
                results[gene_ensembl_id][resource] = [
                    results[gene_ensembl_id][resource][0]
                ]

    else:

        # Consider all gene ids
        gene_ensembl_ids = opentargets_results["gene_ensembl_ids"]
        results = opentargets_results

    for gene_ensembl_id in gene_ensembl_ids:

        # Map gene id to gene name
        gene_symbol = map_gene_ensembl_id_to_names(gene_ensembl_id, gid2nms)[0]

        # Follow term naming convention for parsing
        gs_term = f"GS_{gene_symbol}"  # gene_ensembl_id.replace("ENSG", "GS_")

        # == Gene relations

        for disease in results[gene_ensembl_id]["diseases"]:
            mondo_term = get_mondo_term(disease["disease"]["id"], efo2mondo)
            if mondo_term is None:
                continue
            if disease["score"] < 0.5:
                # Skip diseases with low evidence scores
                continue

            # == Disease relations

            # Gene, IS_GENETIC_BASIS_FOR_CONDITION, Disease
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{RDFSBASE}#GENETIC_BASIS_FOR"),
                    URIRef(f"{PURLBASE}/{mondo_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{mondo_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Open Targets"),
                )
            )

            # == Disease annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{mondo_term}"),
                        URIRef(f"{RDFSBASE}#Name"),
                        Literal(str(disease["disease"]["name"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{mondo_term}"),
                        URIRef(f"{RDFSBASE}#Description"),
                        Literal(str(disease["disease"]["description"])),
                    ),
                ]
            )

            # == Gene to Disease edge annotation

            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{mondo_term}"),
                    URIRef(f"{RDFSBASE}#Score"),
                    Literal(str(disease["score"])),
                )
            )

        for drug in results[gene_ensembl_id]["drugs"]:
            mondo_term = get_mondo_term(drug["diseaseId"], efo2mondo)
            if (
                mondo_term is None
                or drug["drug"]["maximumClinicalTrialPhase"] < 3
                or not drug["drug"]["isApproved"]
            ):
                continue
            # TODO: Test disease score

            # Follow term naming convention for parsing
            chembl_term = drug["drug"]["id"].replace("CHEMBL", "CHEMBL_")

            # == Drug_product relations

            # Drug_product, IS_SUBSTANCE_THAT_TREATS, Disease
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{chembl_term}"),
                    URIRef(f"{RDFSBASE}#IS_SUBSTANCE_THAT_TREATS"),
                    URIRef(f"{PURLBASE}/{mondo_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{chembl_term}"),
                    URIRef(f"{PURLBASE}/{mondo_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Open Targets"),
                )
            )

            # Drug_product, MOLECULARLY_INTERACTS_WITH, Protein
            if (
                "UniProt_name" in gene_results[gene_symbol]
                and gene_results[gene_symbol]["UniProt_name"]
            ):
                # Map gene name to protein uniprot name
                pr_term = f"PR_{gene_results[gene_symbol]['UniProt_name']}"
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#MOLECULARLY_INTERACTS_WITH"),
                        URIRef(f"{PURLBASE}/{pr_term}"),
                    )
                )
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{PURLBASE}/{pr_term}"),
                        URIRef(f"{RDFSBASE}#Source"),
                        Literal("Open Targets and UniProt"),
                    )
                )

            for indication in drug["drug"]["indications"]["rows"]:
                mondo_term = get_mondo_term(indication["disease"]["id"], efo2mondo)
                if mondo_term is None:
                    continue
                # TODO: Test disease score

                # == Indications annotations

                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Indications"),
                        Literal(mondo_term),
                    ),
                )

            for drug_trial_id in drug["ctIds"]:

                # Follow term naming convention for parsing
                nct_term = drug_trial_id.replace("NCT", "NCT_")

                # == Clinical_trial relations

                # Drug_product, EVALUATED_IN, Clinical_trial
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#EVALUATED_IN"),
                        URIRef(f"{PURLBASE}/{nct_term}"),
                    )
                )
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{PURLBASE}/{nct_term}"),
                        URIRef(f"{RDFSBASE}#Source"),
                        Literal("Open Targets"),
                    )
                )

                # == Clinical_trial annotations

                # TODO: Find another source for clinical trial data
                # tuples.extend(
                #     [
                #         (
                #             URIRef(f"{PURLBASE}/{nct_term}"),
                #             URIRef(f"{RDFSBASE}#Phase"),
                #             Literal(str(drug["trial_phase"])),
                #         ),
                #         (
                #             URIRef(f"{PURLBASE}/{nct_term}"),
                #             URIRef(f"{RDFSBASE}#Status"),
                #             Literal(str(drug["trial_status"])),
                #         ),
                #     ]
                # )

            # == Drug_product annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Name"),
                        Literal(str(drug["drug"]["name"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Target"),
                        Literal(gs_term.replace("GS_", "")),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Type"),
                        Literal(str(drug["drugType"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Mechanism_of_action"),
                        Literal(str(drug["mechanismOfAction"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Description"),
                        Literal(str(drug["drug"]["description"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Synonyms"),
                        Literal(str(drug["drug"]["synonyms"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Trade_names"),
                        Literal(str(drug["drug"]["tradeNames"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Approved"),
                        Literal(str(drug["drug"]["isApproved"])),
                    ),
                ]
            )

            pubchem_id = map_chembl_to_pubchem(
                chembl_term.replace("_", ""), chembl2pubchem
            )
            if pubchem_id:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Link_to_PubChem_record"),
                        Literal(f"pubchem.ncbi.nlm.nih.gov/compound/{pubchem_id}"),
                    )
                )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{chembl_term}"),
                    URIRef(f"{RDFSBASE}#Link_to_UniProt_ID"),
                    Literal(
                        remove_protocols(
                            gene_results[gene_symbol]["Link_to_UniProt_ID"]
                        )
                    ),
                )
            )

        # Note: Removed in NLM Cell KN MVP UX/UI specification
        # for interaction in results[gene_ensembl_id]["interactions"]:
        #     if interaction["gene_b_id"] is None:
        #         # Skip interactions missing the second protein
        #         continue
        #     if (
        #         interaction["evidence_score"] is None
        #         or interaction["evidence_score"] < 0.5
        #     ):
        #         # Skip interactions with low evidence scores
        #         continue

        #     # Map id to name
        #     gene_b_id = interaction["gene_b_id"]
        #     gene_b_symbol = map_gene_ensembl_id_to_names(gene_b_id, gid2nms)
        #     if len(gene_b_symbol) == 0:
        #         # Skip interactions with no second gene symbol
        #         continue
        #     gene_b_symbol = gene_b_symbol[0]

        #     # Follow term naming convention for parsing
        #     gs_b_term = (
        #         f"GS_{gene_b_symbol}"  # interaction["gene_b_id"].replace("ENSG", "GS_")
        #     )

        #     # == Interaction relations

        #     # Gene, GENETICALLY_INTERACTS_WITH, Gene
        #     tuples.append(
        #         (
        #             URIRef(f"{PURLBASE}/{gs_term}"),
        #             URIRef(f"{RDFSBASE}#GENETICALLY_INTERACTS_WITH"),
        #             URIRef(f"{PURLBASE}/{gs_b_term}"),
        #         )
        #     )
        #     tuples.append(
        #         (
        #             URIRef(f"{PURLBASE}/{gs_term}"),
        #             URIRef(f"{PURLBASE}/{gs_b_term}"),
        #             URIRef(f"{RDFSBASE}#Source"),
        #             Literal("Open Targets"),
        #         )
        #     )

        #     # Get protein terms, handling Ensembl ids and the term
        #     # naming convention for parsing
        #     protein_a_id = interaction["protein_a_id"]
        #     pr_a_term = get_protein_term(protein_a_id, ensp2accn)
        #     protein_b_id = interaction["protein_b_id"]
        #     pr_b_term = get_protein_term(protein_b_id, ensp2accn)

        #     # Gene, PRODUCES, Protein
        #     if pr_a_term is not None:
        #         tuples.append(
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{RDFSBASE}#PRODUCES"),
        #                 URIRef(f"{PURLBASE}/{pr_a_term}"),
        #             )
        #         )
        #         tuples.append(
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{PURLBASE}/{pr_a_term}"),
        #                 URIRef(f"{RDFSBASE}#Source"),
        #                 Literal("Open Targets"),
        #             )
        #         )

        #         for drug in results[gene_ensembl_id]["drugs"]:
        #             mondo_term = get_mondo_term(drug["disease_id"], efo2mondo)
        #             if (
        #                 mondo_term is None
        #                 or drug["trial_phase"] < 3
        #                 or not drug["approved"]
        #             ):
        #                 continue
        #             # TODO: Test disease score

        #             # Follow term naming convention for parsing
        #             chembl_term = drug["id"].replace("CHEMBL", "CHEMBL_")

        #             # == Drug_product relations

        #             # Drug_product, MOLECULARLY_INTERACTS_WITH, Protein
        #             tuples.append(
        #                 (
        #                     URIRef(f"{PURLBASE}/{chembl_term}"),
        #                     URIRef(f"{RDFSBASE}#MOLECULARLY_INTERACTS_WITH"),
        #                     URIRef(f"{PURLBASE}/{pr_a_term}"),
        #                 )
        #             )
        #             tuples.append(
        #                 (
        #                     URIRef(f"{PURLBASE}/{chembl_term}"),
        #                     URIRef(f"{PURLBASE}/{pr_a_term}"),
        #                     URIRef(f"{RDFSBASE}#Source"),
        #                     Literal("Open Targets"),
        #                 )
        #             )

        #     # == Gene to Interaction edge annotations

        #     tuples.extend(
        #         [
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{PURLBASE}/{gs_b_term}"),
        #                 URIRef(f"{RDFSBASE}#Evidence_score"),
        #                 Literal(str(interaction["evidence_score"])),
        #             ),
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{PURLBASE}/{gs_b_term}"),
        #                 URIRef(f"{RDFSBASE}#Evidence_count"),
        #                 Literal(str(interaction["evidence_count"])),
        #             ),
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{PURLBASE}/{gs_b_term}"),
        #                 URIRef(f"{RDFSBASE}#Source_DB"),
        #                 Literal(str(interaction["source_db"])),
        #             ),
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{PURLBASE}/{gs_b_term}"),
        #                 URIRef(f"{RDFSBASE}#Protein_a"),
        #                 Literal(str(protein_a_id)),
        #             ),
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{PURLBASE}/{gs_b_term}"),
        #                 URIRef(f"{RDFSBASE}#Protein_b"),
        #                 Literal(str(protein_b_id)),
        #             ),
        #         ]
        #     )

        for pharmacogenetic in results[gene_ensembl_id]["pharmacogenetics"]:
            if pharmacogenetic["variantRsId"] is None:
                # Skip pharmacogenetics missing an id
                continue

            # Follow term naming convention for parsing
            rs_term = pharmacogenetic["variantRsId"].replace("rs", "RS_")
            so_term = pharmacogenetic["variantFunctionalConsequenceId"]

            # == Pharmacogenetic relations

            # Gene, HAS_QUALITY, Mutation
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{RDFSBASE}#HAS_QUALITY"),
                    URIRef(f"{PURLBASE}/{rs_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{rs_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Open Targets"),
                )
            )

            # Mutation, INVOLVED_IN, Variant_consequence
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{rs_term}"),
                    URIRef(f"{RDFSBASE}#INVOLVED_IN"),
                    URIRef(f"{PURLBASE}/{so_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{rs_term}"),
                    URIRef(f"{PURLBASE}/{so_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Open Targets"),
                )
            )

            for pharmacogenetic_drug in pharmacogenetic["drugs"]:
                # TODO: Check drug trial phase when available
                if pharmacogenetic_drug["drugId"] is None:
                    # Skip pharmacogenetic drugs missing an id
                    continue

                # Follow term naming convention for parsing
                pharmacogenetic_chembl_term = pharmacogenetic_drug["drugId"].replace(
                    "CHEMBL", "CHEMBL_"
                )

                # Mutation, HAS_PHARMACOLOGICAL_EFFECT, Drug_product
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#HAS_PHARMACOLOGICAL_EFFECT"),
                        URIRef(f"{PURLBASE}/{pharmacogenetic_chembl_term}"),
                    )
                )
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{PURLBASE}/{pharmacogenetic_chembl_term}"),
                        URIRef(f"{RDFSBASE}#Source"),
                        Literal("Open Targets"),
                    )
                )

            # == Pharmacogenetic annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Genotype_ID"),
                        Literal(str(pharmacogenetic["genotypeId"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Genotype"),
                        Literal(str(pharmacogenetic["genotype"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Phenotype"),
                        Literal(str(pharmacogenetic["phenotypeText"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Genotype_annotation"),
                        Literal(str(pharmacogenetic["genotypeAnnotationText"])),
                    ),
                    # (
                    #     URIRef(f"{PURLBASE}/{rs_term}"),
                    #     URIRef(f"{RDFSBASE}#Response_category"),
                    #     Literal(str(pharmacogenetic["response_category"])),
                    # ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Evidence_level"),
                        Literal(str(pharmacogenetic["evidenceLevel"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Source"),
                        Literal(str(pharmacogenetic["datasourceId"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Literature"),
                        Literal(str(pharmacogenetic["literature"])),
                    ),
                ]
            )

            # == Variant_consequence annotations

            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{so_term}"),
                    URIRef(f"{RDFSBASE}#Variant_consequence_label"),
                    Literal(str(pharmacogenetic["variantFunctionalConsequence"])),
                )
            )

        # == Tractability relations

        # None

        # == Gene annotations

        # Note: Now created from NCBI Gene results
        # tuples.append(
        #     (
        #         URIRef(f"{PURLBASE}/{gs_term}"),
        #         URIRef(f"{RDFSBASE}#Symbol"),
        #         Literal(str(gene_symbol)),
        #     )
        # )

        # Note: Removed in NLM Cell KN MVP UX/UI specification
        # for tractability in results[gene_ensembl_id]["tractability"]:
        #     if tractability == {}:
        #         # Skip empty tractabilities
        #         continue
        #     tuples.extend(
        #         [
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{RDFSBASE}#Label"),
        #                 Literal(str(tractability["label"])),
        #             ),
        #             (
        #                 URIRef(f"{PURLBASE}/{gs_term}"),
        #                 URIRef(f"{RDFSBASE}#Modality"),
        #                 Literal(str(tractability["modality"])),
        #             ),
        #         ]
        #     )

        for expression in results[gene_ensembl_id]["expression"]:
            if expression["tissue"]["id"][0:7] != "UBERON_":
                # Skip expressions for tissue not in UBERON
                continue
            exp_term = expression["tissue"]["id"]

            # == Expression relations

            # Gene, EXPRESSED_IN, Anatomical_structure
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{RDFSBASE}#EXPRESSED_IN"),
                    URIRef(f"{PURLBASE}/{exp_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{exp_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("Open Targets"),
                )
            )

            # == Gene to Anatomical_structure edge annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{gs_term}"),
                        URIRef(f"{PURLBASE}/{exp_term}"),
                        URIRef(f"{RDFSBASE}#RNA_zscore"),
                        Literal(str(expression["rna"]["zscore"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gs_term}"),
                        URIRef(f"{PURLBASE}/{exp_term}"),
                        URIRef(f"{RDFSBASE}#RNA_value"),
                        Literal(str(expression["rna"]["value"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gs_term}"),
                        URIRef(f"{PURLBASE}/{exp_term}"),
                        URIRef(f"{RDFSBASE}#RNA_unit"),
                        Literal(str(expression["rna"]["unit"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gs_term}"),
                        URIRef(f"{PURLBASE}/{exp_term}"),
                        URIRef(f"{RDFSBASE}#RNA_level"),
                        Literal(str(expression["rna"]["level"])),
                    ),
                ]
            )

    return tuples, results


def create_tuples_from_gene(gene_path, summarize=False):
    """Creates tuples from the result of using the E-Utilities to
    fetch Gene data for each gene symbol in the NSForest results
    loaded from the specified path. If summnarizing, retain the first
    gene symbol only.

    Parameters
    ----------
    gene_path : Path
        Path to gene results
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    results : dict
        Dictionary containg Gene results keyed by gene symbol
    """
    tuples = []

    # Load the Gene results
    nsforest_path = Path(str(gene_path).replace("-gene.json", ".csv"))
    _gene_path, gene_results = get_gene_results(nsforest_path)

    # Assign gene symbols to consider
    if summarize:

        # Find a gene symbol for which results are not empty
        for gene_symbol in gene_results["gene_symbols"]:
            if not gene_results[gene_symbol]:
                break

        # Consider selected gene symbol
        gene_symbols = [gene_symbol]
        results = {}
        results[gene_symbol] = gene_results[gene_symbol]

    else:

        # Consider all gene symbols
        gene_symbols = gene_results["gene_symbols"]
        results = gene_results

    keys = [
        "Official_symbol",
        "Official_full_name",
        "Gene_type",
        "Link_to_UniProt_ID",
        "Organism",
        "RefSeq_gene_ID",
        "Also_known_as",
        "Summary",
        "UniProt_name",
        "mRNA_(NM)_and_protein_(NP)_sequences",
    ]
    for gene_symbol in gene_symbols:
        if not results[gene_symbol]:
            # Skip empty gene symbol
            continue
        gs_term = f"GS_{gene_symbol}"

        # == Gene relations

        # Gene, PRODUCES, Protein
        if (
            "UniProt Name" in gene_results[gene_symbol]
            and gene_results[gene_symbol]["UniProt_name"]
        ):
            # Map gene name to protein uniprot name
            pr_term = f"PR_{gene_results[gene_symbol]['UniProt Name']}"
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{RDFSBASE}#PRODUCES"),
                    URIRef(f"{PURLBASE}/{pr_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gs_term}"),
                    URIRef(f"{PURLBASE}/{pr_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("UniProt"),
                )
            )

        # == Gene annotations

        for key in keys:
            if key in gene_results[gene_symbol] and gene_results[gene_symbol][key]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{gs_term}"),
                        URIRef(f"{RDFSBASE}#{key.replace(' ', '_')}"),
                        Literal(remove_protocols(gene_results[gene_symbol][key])),
                    )
                )

    return tuples, results


def create_tuples_from_uniprot(uniprot_path, summarize=False):
    """Creates tuples from the result of using a UniProt API endpoint
    for each protein accession in the gene results corresponding to
    the specified path. If summarizing, retain the first protein
    accession only.

    Parameters
    ----------
    uniprot_path : Path
        Path to uniprot results
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    results : dict
        Dictionary containg UniProt results keyed by protein accession
    """
    tuples = []

    # Load the UniProt results
    gene_path = Path(str(uniprot_path).replace("uniprot", "gene"))
    _uniprot_path, uniprot_results = get_uniprot_results(gene_path)

    # Assign protein accessions to consider
    if summarize:

        # Find a protein accession for which results are not empty
        for protein_accession in uniprot_results["protein_accessions"]:
            if not uniprot_results[protein_accession]:
                break

        # Consider selected protein accession
        protein_accessions = [protein_accession]
        results = {}
        results[protein_accession] = uniprot_results[protein_accession]

    else:

        # Consider all protein ids
        protein_accessions = uniprot_results["protein_accessions"]
        results = uniprot_results

    keys = [
        "Protein_name",
        "UniProt_ID",
        "Gene_name",
        "Number_of_amino_acids",
        "Function",
        "Annotation_score",
        "Organism",
    ]
    for protein_accession in protein_accessions:

        # == Protein annotations

        pr_term = f"PR_{protein_accession}"
        for key in keys:
            if key in uniprot_results[protein_accession]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{pr_term}"),
                        URIRef(f"{RDFSBASE}#{key.replace(' ', '_')}"),
                        Literal(uniprot_results[protein_accession][key]),
                    )
                )

    return tuples, results


def create_tuples_from_hubmap(hubmap_path, cl_terms):
    """Creates tuples from HuBMAP data tables.

    Parameters
    ----------
    hubmap_path : Path
        Path to HuBMAP data table JSON file
    cl_terms : set(str)
        Set of all CL terms identified in all author to CL results

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    hubmap_data : dict
        Dictionary containg HuBMAP data table
    """
    tuples = []

    # Load JSON file
    with open(hubmap_path, "r") as fp:
        hubmap_data = json.load(fp)

    key_set = set(["id", "ccf_part_of"])
    anatomical_structures = hubmap_data["data"]["anatomical_structures"]
    for anatomical_structure in anatomical_structures:
        if not key_set.issubset(set(anatomical_structure.keys())):
            # Skip anatomical structure without the needed keys
            continue

        # Get the subject UBERON term
        s_uberon_term = anatomical_structure["id"].replace(":", "_")
        if "UBERON" not in s_uberon_term:
            continue

        # Get each object UBERON term
        for o_uberon_term in anatomical_structure["ccf_part_of"]:
            if "UBERON" not in o_uberon_term:
                continue
            o_uberon_term = o_uberon_term.replace(":", "_")

            # == Anatomical structure relations

            # Anatomical_structure, PART_OF, Anatomical_structure
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{s_uberon_term}"),
                    URIRef(f"{RDFSBASE}#PART_OF"),
                    URIRef(f"{PURLBASE}/{o_uberon_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{s_uberon_term}"),
                    URIRef(f"{PURLBASE}/{o_uberon_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("HuBMAP"),
                )
            )

    # Consider each cell type which has a CL term related to an UBERON
    # term
    key_set = set(["id", "ccf_located_in"])
    cell_types = hubmap_data["data"]["cell_types"]
    for cell_type in cell_types:
        if not key_set.issubset(set(cell_type.keys())):
            # Skip cell types without the needed keys
            continue

        # Get the CL term
        cl_term = cell_type["id"].replace(":", "_")
        if "CL" not in cl_term or "PCL" in cl_term:
            continue

        # Skip CL terms that do not exist in any author to CL mapping
        if cl_term not in cl_terms:
            continue

        # Get each UBERON term
        for uberon_term in cell_type["ccf_located_in"]:
            if "UBERON" not in uberon_term:
                continue
            uberon_term = uberon_term.replace(":", "_")

            # == Cell type relations

            # Cell_type, PART_OF, Anatomical_structure
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{RDFSBASE}#PART_OF"),
                    URIRef(f"{PURLBASE}/{uberon_term}"),
                )
            )
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{cl_term}"),
                    URIRef(f"{PURLBASE}/{uberon_term}"),
                    URIRef(f"{RDFSBASE}#Source"),
                    Literal("HuBMAP"),
                )
            )

    return tuples, hubmap_data


def remove_protocols(value):
    """Remove hypertext protocols.

    Parameters
    ----------
    value : any
       Any value, howerver, only strings are processed

    Returns
    -------
    value : any
       The value, if type str, with hypertext protocols removed
    """
    if isinstance(value, str):
        value = value.replace("http://", "")
        value = value.replace("https://", "")
    return value


def get_cl_terms(author_to_cl_results):
    """Create a set of clean CL terms from the given author to CL results.

    Parameters
    ----------
    author_cl_results : pd.DataFrame
        DataFrame containing author to CL results

    Returns
    -------
    set(str)
        Set of clean CL terms
    """
    return set(
        author_to_cl_results.loc[
            author_to_cl_results["cell_ontology_id"].str.contains("CL"),
            "cell_ontology_id",
        ]
        .str.replace("http://purl.obolibrary.org/obo/", "")
        .str.replace("https://purl.obolibrary.org/obo/", "")
    )


def main(summarize=False):
    """Load results from

    - using the gget opentargets command to obtain the diseases,
      drugs, interactions, pharmacogenetics, tractability, expression,
      and depmap resources for each unique gene id mapped from each
      gene symbol in the NSForest results from processing datasets
      corresponding to the Guo et al. 2023, Li et al. 2023, and
      Sikkema, et al. 2023 publications,

    - using an EBI API endpoint to obtain drug ontology data for each
      unique drug name in the opentargets results,

    - using an RxNav API endpoint for each unique drug name in the
      opentargets results to obtain the mapping from drug name to
      RXCUI, suggested spellings, prescribable drugs information, and
      drug properties, and

    - using a UniProt API endpoint for each protein id in the
      opentargets results to obtain other protein ids, descriptions,
      and comments

    Also, load data tables from HuBMAP.

    Then create tuples consistent with schema v0.7, and write the
    result to JSON files. If summarizing, retain the first gene id
    opentargets results, and protein id uniprot results only, and
    include results in output.

    Note that tuples created from HuBMAP data tables are not
    summarized.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Load results from external API fetch and create tuples
    cl_terms = None
    nsforest_paths = [
        Path(p).resolve()
        for p in glob(str(NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-*.csv"))
    ]
    for nsforest_path in nsforest_paths:

        # Collect all CL terms identified in all author to CL results
        author = re.search("results-([a-zA-Z]*)", nsforest_path.name).group(1)
        author_to_cl_path = Path(
            glob(
                str(NSFOREST_DIRPATH / f"cell-kn-mvp-map-author-to-cl-{author}-*.csv")
            )[-1]
        ).resolve()
        author_to_cl_results = load_results(author_to_cl_path)
        if cl_terms is None:
            cl_terms = get_cl_terms(author_to_cl_results)

        else:
            cl_terms = cl_terms.union(get_cl_terms(author_to_cl_results))

        opentargets_path = Path(str(nsforest_path).replace(".csv", "-opentargets.json"))

        print(f"Creating tuples from {opentargets_path}")
        opentargets_tuples, opentargets_results = create_tuples_from_opentargets(
            opentargets_path, summarize=summarize
        )
        tuples_to_load = opentargets_tuples.copy()

        gene_path = Path(str(nsforest_path).replace(".csv", "-gene.json"))

        print(f"Creating tuples from {gene_path}")

        gene_tuples, gene_results = create_tuples_from_gene(
            gene_path, summarize=summarize
        )
        tuples_to_load.extend(gene_tuples)

        uniprot_path = Path(str(nsforest_path).replace(".csv", "-uniprot.json"))

        print(f"Creating tuples from {uniprot_path}")

        uniprot_tuples, uniprot_results = create_tuples_from_uniprot(
            uniprot_path, summarize=summarize
        )
        tuples_to_load.extend(uniprot_tuples)

        if summarize:
            output_dirpath = TUPLES_DIRPATH / "summaries"

        else:
            output_dirpath = TUPLES_DIRPATH

        with open(
            output_dirpath
            / nsforest_path.name.replace("nsforest-results", "external-api").replace(
                ".csv", ".json"
            ),
            "w",
        ) as f:
            data = {}
            if summarize:
                data["opentargets"] = opentargets_results
                data["uniprot"] = uniprot_results
                data["gene"] = gene_results
            data["tuples"] = tuples_to_load
            json.dump(data, f, indent=4)

        if summarize:
            break

    # Load data from HuBMAP and create tuples
    if not summarize:
        hubmap_paths = [Path(p).resolve() for p in glob(str(HUBMAP_DIRPATH / "*.json"))]
        for hubmap_path in hubmap_paths:
            print(f"Creating tuples from {hubmap_path}")
            hubmap_tuples, _hubmap_data = create_tuples_from_hubmap(
                hubmap_path, cl_terms
            )
            with open(TUPLES_DIRPATH / f"hubmap-{hubmap_path.name}", "w") as f:
                data = {}
                data["tuples"] = hubmap_tuples
                json.dump(data, f, indent=4)


if __name__ == "__main__":
    main(summarize=True)
    main()
