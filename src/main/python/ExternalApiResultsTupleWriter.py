from glob import glob
import json
from pathlib import Path

from rdflib.term import Literal, URIRef

from ExternalApiResultsFetcher import (
    RESOURCES,
    get_opentargets_results,
    get_uniprot_results,
)
from LoaderUtilities import (
    PURLBASE,
    RDFSBASE,
    get_gene_id_to_names_map,
    map_protein_id_to_accession,
    map_gene_id_to_names,
)

NSFOREST_DIRPATH = Path("../../../data/results")
TUPLES_DIRPATH = Path("../../../data/tuples")


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
        accession = map_protein_id_to_accession(protein_id, ensp2accn)

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

    # Load the UniProt results to obtain the Ensembl id to UniProt
    # accession mapping
    _uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)
    ensp2accn = uniprot_results["ensp2accn"]

    # Assign gene ids to consider
    gid2nms = get_gene_id_to_names_map()
    if summarize:

        # Find a gene id with all resources, and a valid disease and interaction
        for gene_id in opentargets_results["gene_ids"]:

            # Find a gene id for which all resources are not empty
            is_empty = False
            for resource in RESOURCES:
                if len(opentargets_results[gene_id][resource]) == 0:
                    is_empty = True
                    break
            if not is_empty:

                # Find a valid disease
                found_disease = False
                for disease in opentargets_results[gene_id]["diseases"]:
                    if "MONDO" in disease["id"] and disease["score"] > 0.5:
                        found_disease = True
                        break

                # Find a valid interaction
                found_interaction = False
                for interaction in opentargets_results[gene_id]["interactions"]:
                    if (
                        interaction["evidence_score"] is not None
                        and interaction["evidence_score"] > 0.5
                    ):
                        found_interaction = True
                        break

                if found_disease and found_interaction:
                    break

        # Consider selected gene id
        gene_ids = [gene_id]
        results = {}
        results[gene_id] = opentargets_results[gene_id]
        results[gene_id]["symbol"] = map_gene_id_to_names(gene_id, gid2nms)[0]

        # Retain only the first result for each resource
        for resource in RESOURCES:
            if resource == "diseases":
                results[gene_id][resource] = [disease]
            elif resource == "interactions":
                results[gene_id][resource] = [interaction]
            else:
                results[gene_id][resource] = [results[gene_id][resource][0]]

    else:

        # Consider all gene ids
        gene_ids = opentargets_results["gene_ids"]
        results = opentargets_results

    for gene_id in gene_ids:

        # Map id to name
        gene_symbol = map_gene_id_to_names(gene_id, gid2nms)[0]

        # Follow term naming convention for parsing
        gene_term = f"GS_{gene_symbol}"  # gene_id.replace("ENSG", "GS_")

        # == Gene relations

        for disease in results[gene_id]["diseases"]:
            if "EFO" in disease["id"]:
                # Skip EFO terms
                continue

            if disease["score"] < 0.5:
                # Skip diseases with low evidence scores
                continue

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
                        URIRef(f"{RDFSBASE}#Name"),
                        Literal(str(disease["name"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{disease['id']}"),
                        URIRef(f"{RDFSBASE}#Description"),
                        Literal(str(disease["description"])),
                    ),
                ]
            )

            # == Gene to Disease edge annotation

            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{PURLBASE}/{disease['id']}"),
                    URIRef(f"{RDFSBASE}#Score"),
                    Literal(str(disease["score"])),
                )
            )

        for drug in results[gene_id]["drugs"]:
            if "EFO" in drug["disease_id"]:
                # Skip EFO terms
                continue

            # Follow term naming convention for parsing
            drug_term = drug["id"].replace("CHEMBL", "CHEMBL_")

            # == Drug_product relations

            # Drug_product, IS_SUBSTANCE_THAT_TREATS, Disease
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{drug_term}"),
                    URIRef(f"{RDFSBASE}#IS_SUBSTANCE_THAT_TREATS"),
                    URIRef(f"{PURLBASE}/{drug['disease_id']}"),
                )
            )

            for drug_trial_id in drug["trial_ids"]:

                # Follow term naming convention for parsing
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
                            URIRef(f"{RDFSBASE}#Phase"),
                            Literal(str(drug["trial_phase"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_trial_term}"),
                            URIRef(f"{RDFSBASE}#Status"),
                            Literal(str(drug["trial_status"])),
                        ),
                    ]
                )

                # == Drug_product annotations

                tuples.extend(
                    [
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Name"),
                            Literal(str(drug["name"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Type"),
                            Literal(str(drug["type"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Mechanism_of_action"),
                            Literal(str(drug["action_mechanism"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Description"),
                            Literal(str(drug["description"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Synonyms"),
                            Literal(str(drug["synonyms"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Trade_names"),
                            Literal(str(drug["trade_names"])),
                        ),
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#Approved"),
                            Literal(str(drug["approved"])),
                        ),
                    ]
                )

        for interaction in results[gene_id]["interactions"]:
            if interaction["gene_b_id"] is None:
                # Skip interactions missing the second protein
                continue
            if (
                interaction["evidence_score"] is None
                or interaction["evidence_score"] < 0.5
            ):
                # Skip interactions with low evidence scores
                continue

            # Map id to name
            gene_b_id = interaction["gene_b_id"]
            gene_b_symbol = map_gene_id_to_names(gene_b_id, gid2nms)[0]

            # Follow term naming convention for parsing
            gene_b_term = (
                f"GS_{gene_b_symbol}"  # interaction["gene_b_id"].replace("ENSG", "GS_")
            )

            # == Interaction relations

            # Gene, GENETICALLY_INTERACTS_WITH, Gene
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#GENETICALLY_INTERACTS_WITH"),
                    URIRef(f"{PURLBASE}/{gene_b_term}"),
                )
            )

            # Get protein terms, handling Ensembl ids and the term
            # naming convention for parsing
            protein_a_id = interaction["protein_a_id"]
            protein_a_term = get_protein_term(protein_a_id, ensp2accn)
            protein_b_id = interaction["protein_b_id"]
            protein_b_term = get_protein_term(protein_b_id, ensp2accn)

            # Gene, PRODUCES, Protein
            if protein_a_term is not None:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{RDFSBASE}#PRODUCES"),
                        URIRef(f"{PURLBASE}/{protein_a_term}"),
                    )
                )

                for drug in results[gene_id]["drugs"]:
                    if "EFO" in drug["disease_id"]:
                        # Skip EFO terms
                        continue

                    # == Drug_product relations

                    # Drug_product, MOLECULARLY_INTERACTS_WITH, Protein
                    tuples.append(
                        (
                            URIRef(f"{PURLBASE}/{drug_term}"),
                            URIRef(f"{RDFSBASE}#MOLECULARLY_INTERACTS_WITH"),
                            URIRef(f"{PURLBASE}/{protein_a_term}"),
                        )
                    )

            # == Gene to Interaction edge annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{gene_b_term}"),
                        URIRef(f"{RDFSBASE}#Evidence_score"),
                        Literal(str(interaction["evidence_score"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{gene_b_term}"),
                        URIRef(f"{RDFSBASE}#Evidence_count"),
                        Literal(str(interaction["evidence_count"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{gene_b_term}"),
                        URIRef(f"{RDFSBASE}#Source_db"),
                        Literal(str(interaction["source_db"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{gene_b_term}"),
                        URIRef(f"{RDFSBASE}#Protein_a"),
                        Literal(str(protein_a_id)),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{gene_b_term}"),
                        URIRef(f"{RDFSBASE}#Protein_b"),
                        Literal(str(protein_b_id)),
                    ),
                ]
            )

        for pharmacogenetic in results[gene_id]["pharmacogenetics"]:
            if pharmacogenetic["rs_id"] is None:
                # Skip pharmacogenetics missing an id
                continue

            # Follow term naming convention for parsing
            rs_term = pharmacogenetic["rs_id"].replace("rs", "RS_")
            variant_consequence_term = pharmacogenetic[
                "variant_consequence_id"
            ].replace("SO:", "SO_")

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
                    URIRef(f"{PURLBASE}/{variant_consequence_term}"),
                )
            )

            for pharmacogenetic_drug in pharmacogenetic["drugs"]:
                if pharmacogenetic_drug["id"] is None:
                    # Skip pharmacogenetic drugs missing an id
                    continue

                # Follow term naming convention for parsing
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
                        URIRef(f"{RDFSBASE}#Genotype_id"),
                        Literal(str(pharmacogenetic["genotype_id"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Genotype"),
                        Literal(str(pharmacogenetic["genotype"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Phenotype"),
                        Literal(str(pharmacogenetic["phenotype"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Genotype_annotation"),
                        Literal(str(pharmacogenetic["genotype_annotation"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Response_category"),
                        Literal(str(pharmacogenetic["response_category"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Evidence_level"),
                        Literal(str(pharmacogenetic["evidence_level"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{rs_term}"),
                        URIRef(f"{RDFSBASE}#Source"),
                        Literal(str(pharmacogenetic["source"])),
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
                    URIRef(f"{PURLBASE}/{variant_consequence_term}"),
                    URIRef(f"{RDFSBASE}#Variant_consequence_label"),
                    Literal(str(pharmacogenetic["variant_consequence_label"])),
                )
            )

        # == Tractability relations

        # None

        # == Gene annotations

        tuples.append(
            (
                URIRef(f"{PURLBASE}/{gene_term}"),
                URIRef(f"{RDFSBASE}#Symbol"),
                Literal(str(gene_symbol)),
            )
        )

        for tractability in results[gene_id]["tractability"]:
            if tractability == {}:
                # Skip empty tractabilities
                continue
            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{RDFSBASE}#Label"),
                        Literal(str(tractability["label"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{RDFSBASE}#Modality"),
                        Literal(str(tractability["modality"])),
                    ),
                ]
            )

        for expression in results[gene_id]["expression"]:
            if expression["tissue_id"][0:7] != "UBERON_":
                # Skip expressions for tissue not in UBERON
                continue

            # == Expression relations

            # Gene, EXPRESSED_IN, Anatomical_structure
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#EXPRESSED_IN"),
                    URIRef(f"{PURLBASE}/{expression['tissue_id']}"),
                )
            )

            # == Gene to Anatomical_structure edge annotations

            tuples.extend(
                [
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{expression['tissue_id']}"),
                        URIRef(f"{RDFSBASE}#RNA_zscore"),
                        Literal(str(expression["rna_zscore"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{expression['tissue_id']}"),
                        URIRef(f"{RDFSBASE}#RNA_value"),
                        Literal(str(expression["rna_value"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{expression['tissue_id']}"),
                        URIRef(f"{RDFSBASE}#RNA_unit"),
                        Literal(str(expression["rna_unit"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{expression['tissue_id']}"),
                        URIRef(f"{RDFSBASE}#RNA_level"),
                        Literal(str(expression["rna_level"])),
                    ),
                ]
            )

    return tuples, results


def create_tuples_from_uniprot(opentargets_path, summarize=False):
    """Creates tuples from the result of using a UniProt API endpoint
    for each protein id in the opentargets results corresponding to
    the specified path to obtain other protein ids, descriptions, and
    comments. If summnarizing, retain the first protein id only.

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
        Dictionary containg UniProt results keyed by protein id
    """
    tuples = []

    # Load the UniProt results, and assign the Ensembl id to UniProt
    # accession mapping
    _uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)
    ensp2accn = uniprot_results["ensp2accn"]

    # Assign protein ids to consider
    if summarize:

        # Find a protein id for which results are not empty
        for protein_id in uniprot_results["protein_ids"]:
            if uniprot_results[protein_id] != [] and "ENSP" in protein_id:
                break

        # Consider selected protein id
        protein_ids = [protein_id]
        results = {}
        results[protein_id] = uniprot_results[protein_id]
        results[protein_id]["accession"] = get_protein_term(
            protein_id, ensp2accn
        ).replace("PR_", "")
    else:

        # Consider all protein ids
        protein_ids = uniprot_results["protein_ids"]
        results = uniprot_results

    for protein_id in protein_ids:
        if results[protein_id] == []:
            # Skip empty proteins
            continue

        # == Protein annotations

        # Get protein term, handling Ensembl ids and the term naming
        # convention for parsing
        protein_term = get_protein_term(protein_id, ensp2accn)
        if protein_term is None:
            # Skip unmappable protein ids
            continue
        if "proteinDescription" in results[protein_id]:
            if "recommendedName" in results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Recommended_name"),
                        Literal(
                            str(
                                results[protein_id]["proteinDescription"][
                                    "recommendedName"
                                ]["fullName"]["value"]
                            )
                        ),
                    )
                )
            if "alternativeNames" in results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Alternative_name"),
                        Literal(
                            str(
                                results[protein_id]["proteinDescription"][
                                    "alternativeNames"
                                ][0]["fullName"]["value"]
                            )
                        ),
                    )
                )
            if "submissionNames" in results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Submission_name"),
                        Literal(
                            str(
                                results[protein_id]["proteinDescription"][
                                    "submissionNames"
                                ][0]["fullName"]["value"]
                            )
                        ),
                    )
                )
            if "cdAntigenNames" in results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#CD_antigen_name"),
                        Literal(
                            str(
                                results[protein_id]["proteinDescription"][
                                    "cdAntigenNames"
                                ][0]["value"]
                            )
                        ),
                    )
                )
        if "comments" in results[protein_id]:
            if results[protein_id]["comments"] != []:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Comment"),
                        Literal(
                            str(results[protein_id]["comments"][0]["texts"][0]["value"])
                        ),
                    ),
                )

    return tuples, results


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

    Then create tuples consistent with schema v0.7, and write the
    result to a JSON file. If summarizing, retain the first gene id
    opentargets results, and protein id uniprot results only, and
    include results in output.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    nsforest_paths = [
        Path(p).resolve()
        for p in glob(str(NSFOREST_DIRPATH / "cell-kn-mvp-nsforest-results-*.csv"))
    ]
    for nsforest_path in nsforest_paths:
        opentargets_path = Path(str(nsforest_path).replace(".csv", "-opentargets.json"))

        print(f"Creating tuples from {opentargets_path}")
        opentargets_tuples, opentargets_results = create_tuples_from_opentargets(
            opentargets_path, summarize=summarize
        )
        uniprot_tuples, uniprot_results = create_tuples_from_uniprot(
            opentargets_path, summarize=summarize
        )
        tuples_to_load = opentargets_tuples.copy()
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
            data["tuples"] = tuples_to_load
            json.dump(data, f, indent=4)

        if summarize:
            break


if __name__ == "__main__":
    main(summarize=True)
    main()
