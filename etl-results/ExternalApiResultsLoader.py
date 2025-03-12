import argparse
import json
from pathlib import Path

from rdflib.term import Literal, URIRef

import ArangoDbUtilities as adb
from ExternalApiResultsFetcher import (
    RESOURCES,
    get_opentargets_results,
    get_uniprot_results,
)
from LoaderUtilities import (
    PURLBASE,
    RDFSBASE,
    get_gene_id_to_names_map,
    hyphenate,
    map_protein_id_to_accession,
    map_gene_id_to_names,
)
from OntologyParserLoader import load_tuples_into_adb_graph, parse_obo, VALID_VERTICES

OBO_DIRPATH = Path("../data/obo")
NSFOREST_DIRPATH = Path("../data/results")
TUPLES_DIRPATH = Path("../data/tuples")


def get_protein_term(protein_id, ensp2accn):

    protein_term = None

    if "ENSP" in protein_id:
        accession = map_protein_id_to_accession(protein_id, ensp2accn)

    else:
        accession = protein_id

    if accession is not None:
        protein_term = f"PR_{accession}"

    return protein_term


def create_tuples_from_opentargets(opentargets_path):

    tuples = []

    nsforest_path = Path(str(opentargets_path).replace("-opentargets.json", ".csv"))

    opentargets_path, opentargets_results = get_opentargets_results(nsforest_path)

    uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)

    ensp2accn = uniprot_results["ensp2accn"]

    gid2nms = get_gene_id_to_names_map()
    for gene_id in opentargets_results["gene_ids"]:
        gene_symbol = map_gene_id_to_names(gene_id, gid2nms)
        gene_term = gene_id.replace("ENSG", "GS_")

        # == Gene relations

        for disease in opentargets_results[gene_id]["diseases"]:
            if "EFO" in disease["id"]:
                continue
            elif disease["score"] < 0.5:
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

        for interaction in opentargets_results[gene_id]["interactions"]:
            if interaction["gene_b_id"] is None:
                continue
            if (
                interaction["evidence_score"] is None
                or interaction["evidence_score"] < 0.5
            ):
                continue
            gene_b_term = interaction["gene_b_id"].replace("ENSG", "GS_")

            # Gene, GENETICALLY_INTERACTS_WITH, Gene
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#GENETICALLY_INTERACTS_WITH"),
                    URIRef(f"{PURLBASE}/{gene_b_term}"),
                )
            )

            # Gene, PRODUCES, Protein
            protein_a_term = get_protein_term(interaction["protein_a_id"], ensp2accn)
            if protein_a_term is not None:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{RDFSBASE}#PRODUCES"),
                        URIRef(f"{PURLBASE}/{protein_a_term}"),
                    )
                )
            protein_b_term = get_protein_term(interaction["protein_b_id"], ensp2accn)
            if protein_b_term is not None:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{RDFSBASE}#PRODUCES"),
                        URIRef(f"{PURLBASE}/{protein_b_term}"),
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
                        Literal(str(interaction["protein_a_id"])),
                    ),
                    (
                        URIRef(f"{PURLBASE}/{gene_term}"),
                        URIRef(f"{PURLBASE}/{gene_b_term}"),
                        URIRef(f"{RDFSBASE}#Protein_b"),
                        Literal(str(interaction["protein_b_id"])),
                    ),
                ]
            )

        for pharmacogenetic in opentargets_results[gene_id]["pharmacogenetics"]:
            if pharmacogenetic["rs_id"] is None:
                continue
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
                    continue
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
            [
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{RDFSBASE}#Symbol"),
                    Literal(str(gene_symbol)),
                )
            ]
        )

        for tractability in opentargets_results[gene_id]["tractability"]:
            if tractability == {}:
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

        for expression in opentargets_results[gene_id]["expression"]:
            if expression["tissue_id"][0:7] != "UBERON_":
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

    return tuples


def create_tuples_from_uniprot(opentargets_path):

    tuples = []

    uniprot_path, uniprot_results = get_uniprot_results(opentargets_path)

    ensp2accn = uniprot_results["ensp2accn"]

    for protein_id in uniprot_results["protein_ids"]:

        if uniprot_results[protein_id] == []:
            continue

        protein_term = get_protein_term(protein_id, ensp2accn)
        if protein_term is None:
            continue

        if "proteinDescription" in uniprot_results[protein_id]:
            if "recommendedName" in uniprot_results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Recommended_name"),
                        Literal(
                            str(
                                uniprot_results[protein_id]["proteinDescription"][
                                    "recommendedName"
                                ]["fullName"]["value"]
                            )
                        ),
                    )
                )

            if "alternativeNames" in uniprot_results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Alternative_name"),
                        Literal(
                            str(
                                uniprot_results[protein_id]["proteinDescription"][
                                    "alternativeNames"
                                ][0]["fullName"]["value"]
                            )
                        ),
                    )
                )

            if "submissionNames" in uniprot_results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Submission_name"),
                        Literal(
                            str(
                                uniprot_results[protein_id]["proteinDescription"][
                                    "submissionNames"
                                ][0]["fullName"]["value"]
                            )
                        ),
                    )
                )

            if "cdAntigenNames" in uniprot_results[protein_id]["proteinDescription"]:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#CD_antigen_name"),
                        Literal(
                            str(
                                uniprot_results[protein_id]["proteinDescription"][
                                    "cdAntigenNames"
                                ][0]["value"]
                            )
                        ),
                    )
                )

        if "comments" in uniprot_results[protein_id]:
            if uniprot_results[protein_id]["comments"] != []:
                tuples.append(
                    (
                        URIRef(f"{PURLBASE}/{protein_term}"),
                        URIRef(f"{RDFSBASE}#Comment"),
                        Literal(
                            str(
                                uniprot_results[protein_id]["comments"][0]["texts"][0][
                                    "value"
                                ]
                            )
                        ),
                    ),
                )

    return tuples


def main(parameters=None):

    parser = argparse.ArgumentParser(description="Load external API results")
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

    # print("Getting ArangoDB database and graph, and loading tuples")
    # db = adb.create_or_get_database(db_name)
    # adb_graph = adb.create_or_get_graph(db, graph_name)
    # vertex_collections = {}
    # edge_collections = {}

    for author in ["guo", "li", "sikkema"]:

        opentargets_path = (
            NSFOREST_DIRPATH
            / f"cell-kn-mvp-nsforest-results-{author}-2023-2025-02-22-opentargets.json"
        ).resolve()

        opentargets_tuples = create_tuples_from_opentargets(opentargets_path)
        uniprot_tuples = create_tuples_from_uniprot(opentargets_path)

        tuples_to_load = opentargets_tuples.copy()
        tuples_to_load.extend(uniprot_tuples)

        with open(TUPLES_DIRPATH / f"ExternalApiResultsLoader-{author}.json", "w") as f:
            results = {}
            results["tuples"] = tuples_to_load
            json.dump(results, f, indent=4)

        # VALID_VERTICES.add("BMC")
        # VALID_VERTICES.add("CHEMBL")
        # VALID_VERTICES.add("CS")
        # VALID_VERTICES.add("CSD")
        # VALID_VERTICES.add("DOID")
        # VALID_VERTICES.add("DS")
        # VALID_VERTICES.add("GS")
        # VALID_VERTICES.add("HP")
        # VALID_VERTICES.add("PUB")
        # VALID_VERTICES.add("RS")
        # VALID_VERTICES.add("SO")

        # load_tuples_into_adb_graph(
        #     tuples_to_load,
        #     adb_graph,
        #     vertex_collections,
        #     edge_collections,
        #     ro=ro,
        #     do_update=True,
        # )


if __name__ == "__main__":
    main()
