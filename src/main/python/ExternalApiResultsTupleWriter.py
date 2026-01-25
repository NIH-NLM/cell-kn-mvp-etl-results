from glob import glob
import json
from pathlib import Path

from rdflib.term import Literal, URIRef

from ExternalApiResultsFetcher import (
    CELLXGENE_PATH,
    OPENTARGETS_RESOURCES,
    OPENTARGETS_PATH,
    GENE_PATH,
    UNIPROT_PATH,
    HUBMAP_DIRPATH,
)
from LoaderUtilities import (
    DEPRECATED_TERMS,
    PURLBASE,
    RDFSBASE,
    collect_results_sources_data,
    get_chembl_to_pubchem_map,
    get_efo_to_mondo_map,
    get_gene_ensembl_id_to_names_map,
    get_gene_entrez_id_to_names_map,
    get_gene_name_to_entrez_ids_map,
    map_chembl_to_pubchem,
    map_efo_to_mondo,
    map_gene_ensembl_id_to_names,
    map_gene_entrez_id_to_names,
    map_gene_name_to_entrez_ids,
    map_protein_ensembl_id_to_accession,
)

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

    if mondo_term in DEPRECATED_TERMS:
        print(f"Warning: MONDO term {mondo_term} deprecated")
        return None

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


def create_tuples_from_cellxgene(summarize=False):
    """Creates tuples from the result of using the CELLxGENE curation
    API to fetch metadata for dataset version ids. If summnarizing,
    retain the first dataset only.

    Parameters
    ----------
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    results : dict
        Dictionary containg cellxgene results keyed by dataset version
        id
    """
    tuples = []

    # Load the cellxgene results
    print(f"Loading cellxgene results from {CELLXGENE_PATH}")
    with open(CELLXGENE_PATH, "r") as fp:
        cellxgene_results = json.load(fp)

    # Assign datasets to consider
    if summarize:

        # Consider the first dataset only
        results = {}
        for dataset_version_id, dataset_metadata in cellxgene_results.items():
            results[dataset_version_id] = dataset_metadata

    else:

        # Consider all datasets
        results = cellxgene_results

    # Create tuples for each dataset
    for dataset_version_id, dataset_metadata in results.items():
        csd_term = f"CSD_{dataset_version_id}"
        pub_term = f"PUB_{dataset_metadata['Link_to_publication']}"

        # Cell_set_dataset_Ind, SOURCE, Publication_Ind
        # IAO:0000100, dc:source, IAO:0000311
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{csd_term}"),
                URIRef(f"{RDFSBASE}/dc#Source"),
                URIRef(f"{PURLBASE}/{pub_term}"),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{csd_term}"),
                URIRef(f"{PURLBASE}/{pub_term}"),
                URIRef(f"{RDFSBASE}#Source"),
                Literal("Manual Mapping"),
            )
        )

        # CSD node annotations
        keys = [
            "Link_to_publication",
            "Link_to_CELLxGENE_collection",
            "Link_to_CELLxGENE_dataset",
            "Dataset_name",
            "Number_of_cells",
            "Organism",
            "Tissue",
            "Disease_status",
            "Collection_ID",
            "Collection_version_ID",
            "Dataset_ID",
            "Dataset_version_ID",
            "Zenodo/Nextflow_workflow/Notebook",
        ]
        for key in keys:
            value = cellxgene_results[dataset_version_id][key]
            if isinstance(value, str):
                value = value.replace("https://", "")
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{csd_term}"),
                    URIRef(f"{RDFSBASE}#{key.replace(' ', '_')}"),
                    Literal(value),
                )
            )

    return tuples, results


def create_tuples_from_opentargets(summarize=False):
    """Creates tuples from the result of using the Open Targets
    Platform GraphQL API to obtain resources for gene Ensembl ids. If
    summnarizing, retain the first gene Ensembl id only.

    Parameters
    ----------
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
    print(f"Loading opentargets results from {OPENTARGETS_PATH}")
    with open(OPENTARGETS_PATH, "r") as fp:
        opentargets_results = json.load(fp)

    # Load the Gene results to obtain the UniProt names corresponding
    # to gene names
    print(f"Loading gene results from {GENE_PATH}")
    with open(GENE_PATH, "r") as fp:
        gene_results = json.load(fp)

    # Load mappings
    gene_ensembl_id_to_names = get_gene_ensembl_id_to_names_map()
    gene_name_to_entrez_ids = get_gene_name_to_entrez_ids_map()
    efo2mondo = get_efo_to_mondo_map()
    chembl2pubchem = get_chembl_to_pubchem_map()

    # Assign gene ids to consider
    if summarize:

        # Find a gene id with all resources, and a valid disease and interaction
        gene_ensembl_id = None
        for gene_ensembl_id in opentargets_results["gene_ensembl_ids"]:

            # Find a gene id for which all resources are not empty
            is_empty = False
            for resource in OPENTARGETS_RESOURCES:
                if len(opentargets_results[gene_ensembl_id][resource]) < 3:
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
        results[gene_ensembl_id]["name"] = map_gene_ensembl_id_to_names(
            gene_ensembl_id, gene_ensembl_id_to_names
        )[0]

        # Retain the first three results for each resource
        for resource in OPENTARGETS_RESOURCES:
            results[gene_ensembl_id][resource] = results[gene_ensembl_id][resource][0:3]

    else:

        # Consider all gene ids
        gene_ensembl_ids = opentargets_results["gene_ensembl_ids"]
        results = opentargets_results

    for gene_ensembl_id in gene_ensembl_ids:

        # Map gene Ensembl id to gene name and Entrez id
        gene_name = map_gene_ensembl_id_to_names(
            gene_ensembl_id, gene_ensembl_id_to_names
        )
        if gene_name == []:
            continue
        else:
            gene_name = gene_name[0]
        gene_entrez_id = map_gene_name_to_entrez_ids(gene_name, gene_name_to_entrez_ids)
        if gene_entrez_id == []:
            continue
        else:
            gene_entrez_id = gene_entrez_id[0]

        # Follow term naming convention for parsing
        gs_term = f"GS_{gene_name}"  # gene_ensembl_id.replace("ENSG", "GS_")

        # == Gene relations

        for disease in results[gene_ensembl_id]["diseases"]:
            mondo_term = get_mondo_term(disease["disease"]["id"], efo2mondo)
            if mondo_term is None:
                continue
            if disease["score"] < 0.5:
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
                or drug["drug"]["hasBeenWithdrawn"]
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
                "UniProt_name" in gene_results[gene_entrez_id]
                and gene_results[gene_entrez_id]["UniProt_name"]
            ):
                # Map gene name to protein uniprot name
                pr_term = f"PR_{gene_results[gene_entrez_id]['UniProt_name']}"
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

            if drug["drug"]["indications"]:
                for indication in drug["drug"]["indications"]["rows"]:
                    mondo_term = get_mondo_term(indication["disease"]["id"], efo2mondo)
                    if mondo_term is None or indication["maxPhaseForIndication"] < 4:
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

                # None

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
                    (
                        URIRef(f"{PURLBASE}/{chembl_term}"),
                        URIRef(f"{RDFSBASE}#Withdrawn"),
                        Literal(str(drug["drug"]["hasBeenWithdrawn"])),
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
                            gene_results[gene_entrez_id]["Link_to_UniProt_ID"]
                        )
                    ),
                )
            )

        for pharmacogenetic in results[gene_ensembl_id]["pharmacogenetics"]:
            if pharmacogenetic["variantRsId"] is None:
                continue

            # Follow term naming convention for parsing
            rs_term = pharmacogenetic["variantRsId"].replace("rs", "RS_")
            so_term = pharmacogenetic["variantFunctionalConsequenceId"]
            if so_term in DEPRECATED_TERMS:
                print(f"Warning: SO term {so_term} deprecated")

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

        for expression in results[gene_ensembl_id]["expression"]:
            if expression["tissue"]["id"][0:7] != "UBERON_":
                continue
            exp_term = expression["tissue"]["id"]
            if exp_term in DEPRECATED_TERMS:
                print(f"Warning: Expression term {exp_term} deprecated")

            # == Expression relations

            # Gene, EXPRESSED_IN, Anatomical_structure
            # NOTE: Removed to resolve issue 105
            # tuples.append(
            #     (
            #         URIRef(f"{PURLBASE}/{gs_term}"),
            #         URIRef(f"{RDFSBASE}#EXPRESSED_IN"),
            #         URIRef(f"{PURLBASE}/{exp_term}"),
            #     )
            # )
            # tuples.append(
            #     (
            #         URIRef(f"{PURLBASE}/{gs_term}"),
            #         URIRef(f"{PURLBASE}/{exp_term}"),
            #         URIRef(f"{RDFSBASE}#Source"),
            #         Literal("Open Targets"),
            #     )
            # )

            # == Gene to Anatomical_structure edge annotations

            # NOTE: Removed to resolve issue 105
            # tuples.extend(
            #     [
            #         (
            #             URIRef(f"{PURLBASE}/{gs_term}"),
            #             URIRef(f"{PURLBASE}/{exp_term}"),
            #             URIRef(f"{RDFSBASE}#RNA_zscore"),
            #             Literal(str(expression["rna"]["zscore"])),
            #         ),
            #         (
            #             URIRef(f"{PURLBASE}/{gs_term}"),
            #             URIRef(f"{PURLBASE}/{exp_term}"),
            #             URIRef(f"{RDFSBASE}#RNA_value"),
            #             Literal(str(expression["rna"]["value"])),
            #         ),
            #         (
            #             URIRef(f"{PURLBASE}/{gs_term}"),
            #             URIRef(f"{PURLBASE}/{exp_term}"),
            #             URIRef(f"{RDFSBASE}#RNA_unit"),
            #             Literal(str(expression["rna"]["unit"])),
            #         ),
            #         (
            #             URIRef(f"{PURLBASE}/{gs_term}"),
            #             URIRef(f"{PURLBASE}/{exp_term}"),
            #             URIRef(f"{RDFSBASE}#RNA_level"),
            #             Literal(str(expression["rna"]["level"])),
            #         ),
            #     ]
            # )

    return tuples, results


def create_tuples_from_gene(summarize=False):
    """Creates tuples from the result of using the E-Utilities to
    fetch Gene data for gene names. If summnarizing, retain the first
    gene name only.

    Parameters
    ----------
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    tuples : list(tuple(str))
        List of tuples (triples or quadruples) created
    results : dict
        Dictionary containg Gene results keyed by gene name
    """
    tuples = []

    # Load the Gene results
    print(f"Loading gene results from {GENE_PATH}")
    with open(GENE_PATH, "r") as fp:
        gene_results = json.load(fp)

    # Load mappings
    gene_entrez_id_to_names = get_gene_entrez_id_to_names_map()

    # Assign gene names to consider
    if summarize:

        # Find a gene name for which results are not empty
        gene_entrez_id = None
        for gene_entrez_id in gene_results["gene_entrez_ids"]:
            if len(gene_results[gene_entrez_id]) > 0:
                break

        # Consider selected gene name
        gene_entrez_ids = [gene_entrez_id]
        results = {}
        results[gene_entrez_id] = gene_results[gene_entrez_id]

    else:

        # Consider all gene names
        gene_entrez_ids = gene_results["gene_entrez_ids"]
        results = gene_results

    keys = [
        "Gene_ID",
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
    for gene_entrez_id in gene_entrez_ids:
        if not results[gene_entrez_id]:
            continue

        # Map gene Entrez id to gene name
        gene_name = map_gene_entrez_id_to_names(
            gene_entrez_id, gene_entrez_id_to_names
        )[0]
        gs_term = f"GS_{gene_name}"

        # == Gene relations

        # Gene, PRODUCES, Protein
        if (
            "UniProt_name" in gene_results[gene_entrez_id]
            and gene_results[gene_entrez_id]["UniProt_name"]
        ):
            # Map gene name to protein uniprot name
            pr_term = f"PR_{gene_results[gene_entrez_id]['UniProt_name']}"
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
            if (
                key in gene_results[gene_entrez_id]
                and gene_results[gene_entrez_id][key]
            ):
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{gs_term}"),
                        URIRef(f"{RDFSBASE}#{key.replace(' ', '_')}"),
                        Literal(remove_protocols(gene_results[gene_entrez_id][key])),
                    )
                )

    return tuples, results


def create_tuples_from_uniprot(summarize=False):
    """Creates tuples from the result of using a UniProt API endpoint
    for protein accessions. If summarizing, retain the first protein
    accession only.

    Parameters
    ----------
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
    print(f"Loading uniprot results from {UNIPROT_PATH}")
    with open(UNIPROT_PATH, "r") as fp:
        uniprot_results = json.load(fp)

    # Assign protein accessions to consider
    if summarize:

        # Find a protein accession for which results are not empty
        for protein_accession in uniprot_results["protein_accessions"]:
            if len(uniprot_results[protein_accession]) > 0:
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


def create_tuples_from_hubmap(hubmap_path, cl_terms, summarize=False):
    """Creates tuples from HuBMAP data tables.

    Parameters
    ----------
    hubmap_path : Path
        Path to HuBMAP data table JSON file
    cl_terms : set(str)
        Set of all CL terms identified in all author to CL results
    summarize : bool
        Flag to summarize results, or not

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
    found_anatomical_structure_relation = False
    for anatomical_structure in anatomical_structures:
        if not key_set.issubset(set(anatomical_structure.keys())):
            continue

        # Get the subject UBERON term
        s_uberon_term = anatomical_structure["id"].replace(":", "_")
        if "UBERON" not in s_uberon_term:
            continue
        if s_uberon_term in DEPRECATED_TERMS:
            print(f"Warning: UBERON term {s_uberon_term} deprecated")

        # Get each object UBERON term
        for o_uberon_term in anatomical_structure["ccf_part_of"]:
            if "UBERON" not in o_uberon_term:
                continue
            o_uberon_term = o_uberon_term.replace(":", "_")
            if o_uberon_term in DEPRECATED_TERMS:
                print(f"Warning: UBERON term {o_uberon_term} deprecated")

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
            if summarize:
                found_anatomical_structure_relation = True
                break

        if summarize and found_anatomical_structure_relation:
            break

    # Consider each cell type which has a CL term related to an UBERON
    # term
    key_set = set(["id", "ccf_located_in"])
    cell_types = hubmap_data["data"]["cell_types"]
    found_cell_type_anatomical_structure_relation = False
    for cell_type in cell_types:
        if not key_set.issubset(set(cell_type.keys())):
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
            if uberon_term in DEPRECATED_TERMS:
                print(f"Warning: UBERON term {uberon_term} deprecated")

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
            if summarize:
                found_cell_type_anatomical_structure_relation = True
                break

        if summarize and found_cell_type_anatomical_structure_relation:
            break

    if summarize:
        hubmap_data = {}
        hubmap_data["cell_type"] = cell_type
        hubmap_data["anatomical_structure"] = anatomical_structure

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


def main(summarize=False):
    """Load results from:

    - Using the Open Targets Platform GraphQL API to obtain the
      diseases, drugs, interactions, pharmacogenetics, tractability,
      expression, and depmap resources for each unique gene Ensembl id

    - Using the E-Utilities to fetch Gene data for each unique gene
      Entrez id

    - Using a UniProt API endpoint for each protein accession in the
      gene results

    Also, load data tables from HuBMAP.

    Then create tuples consistent with schema v0.7, and write the
    result to JSON files. If summarizing, retain the first gene id
    opentargets results, and protein id uniprot results only, and
    include results in output.

    Note that tuples created from HuBMAP data tables are not
    summarized.

    Parameters
    ----------
    summarize : bool
        Flag to summarize results, or not

    Returns
    -------
    None
    """
    # Collect paths to all NSForest results, and author cell set to CL
    # term mappings identified in the results sources. Collect the
    # dataset_version_ids used for creating the NSForest results
    # paths. Collect the unique gene names, Ensembl identifiers, and
    # Entrez identifiers corresponding to all NSForet results.
    (
        _nsforest_paths,
        _silhouette_paths,
        _author_to_cl_paths,
        _dataset_version_ids,
        cl_terms,
        _gene_names,
        _gene_ensembl_ids,
        _gene_entrez_ids,
    ) = collect_results_sources_data()

    print(f"Creating tuples from {CELLXGENE_PATH}")
    cellxgene_tuples, cellxgene_results = create_tuples_from_cellxgene(
        summarize=summarize
    )
    tuples_to_load = cellxgene_tuples.copy()

    print(f"Creating tuples from {OPENTARGETS_PATH}")
    opentargets_tuples, opentargets_results = create_tuples_from_opentargets(
        summarize=summarize
    )
    tuples_to_load.extend(opentargets_tuples)

    print(f"Creating tuples from {GENE_PATH}")
    gene_tuples, gene_results = create_tuples_from_gene(summarize=summarize)
    tuples_to_load.extend(gene_tuples)

    print(f"Creating tuples from {UNIPROT_PATH}")
    uniprot_tuples, uniprot_results = create_tuples_from_uniprot(summarize=summarize)
    tuples_to_load.extend(uniprot_tuples)

    if summarize:
        output_dirpath = TUPLES_DIRPATH / "summaries"

    else:
        output_dirpath = TUPLES_DIRPATH

    with open(output_dirpath / "cell-kn-mvp-external-api-results.json", "w") as f:
        data = {}
        if summarize:
            data["cellxgene"] = cellxgene_results
            data["opentargets"] = opentargets_results
            data["uniprot"] = uniprot_results
            data["gene"] = gene_results
        data["tuples"] = tuples_to_load
        json.dump(data, f, indent=4)

    # Load data from HuBMAP and create tuples
    hubmap_paths = [Path(p).resolve() for p in glob(str(HUBMAP_DIRPATH / "*.json"))]
    for hubmap_path in hubmap_paths:
        print(f"Creating tuples from {hubmap_path}")
        hubmap_tuples, hubmap_data = create_tuples_from_hubmap(
            hubmap_path, cl_terms, summarize=summarize
        )
        with open(output_dirpath / f"hubmap-{hubmap_path.name}", "w") as f:
            data = {}
            if summarize:
                data["hubmap"] = hubmap_data
            data["tuples"] = hubmap_tuples
            json.dump(data, f, indent=4)
        if summarize:
            break


if __name__ == "__main__":
    main(summarize=True)
    main()
