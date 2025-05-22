import ast
from glob import glob
import json
from pathlib import Path
import re

from rdflib.term import Literal, URIRef

from E_Utilities import get_data_for_pmid
from LoaderUtilities import load_results, hyphenate, PURLBASE, RDFSBASE

NSFOREST_DIRPATH = Path("../../../data/results")
TUPLES_DIRPATH = Path("../../../data/tuples")


def create_tuples_from_author_to_cl(results):
    """Creates tuples from manual author cell set to CL term mapping
    consistent with schema v0.7.

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

    # Nodes for these results
    csd_term = f"CSD_{results['dataset_id'][0]}"
    pub_term = f"PUB_{results['DOI'][0].replace('/', '-')}"

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
            Literal(results["collection_version_id"][0]),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#PMID"),
            Literal(str(results["PMID"][0])),
        )
    )
    tuples.append(
        (
            URIRef(f"{PURLBASE}/{pub_term}"),
            URIRef(f"{RDFSBASE}#PMCID"),
            Literal(str(results["PMCID"][0])),
        )
    )
    data = get_data_for_pmid(results["PMID"][0])
    for key in data.keys():
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{pub_term}"),
                URIRef(f"{RDFSBASE}#{key.capitalize()}"),
                Literal(data[key]),
            )
        )

    # Nodes for each cell type or cell set
    uuid_0 = results["uuid"][0]
    for _, row in results.iterrows():
        uuid = row["uuid"]
        cl_term = row["cell_ontology_id"]
        uberon_term = row["uberon_entity_id"]
        author_cell_set = hyphenate(row["author_cell_set"])
        cs_term = f"CS_{author_cell_set}-{uuid}"

        # Cell_type_Class, PART_OF, Anatomical_structure_Class
        # CL:0000000, BFO:0000050, UBERON:0001062
        tuples.append(
            (
                URIRef(f"{cl_term}"),
                URIRef(f"{PURLBASE}/BFO_0000050"),
                URIRef(f"{uberon_term}"),
            )
        )

        # Cell_type_Class, HAS_EXEMPLAR_DATA, Cell_set_dataset_Ind
        # CL:0000000, RO:0015001, IAO:0000100
        tuples.append(
            (
                URIRef(f"{cl_term}"),
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
                URIRef(f"{uberon_term}-{uuid_0}"),
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
                URIRef(f"{cl_term}"),
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
                URIRef(f"{cl_term}"),
                URIRef(f"{RDFSBASE}#Match"),
                Literal(row["match"]),
            )
        )
        tuples.append(
            (
                URIRef(f"{PURLBASE}/{cs_term}"),
                URIRef(f"{cl_term}"),
                URIRef(f"{RDFSBASE}#Mapping_method"),
                Literal(row["mapping_method"]),
            )
        )

        # Nodes for each cell type and marker gene
        marker_genes = ast.literal_eval(row["NSForest_markers"])
        for gene in marker_genes:
            gene_term = f"GS_{gene}"

            # Gene_Class, IS_MARKER_FOR, Cell_type_Class
            # SO:0000704, RO:0002607, CL:0000000
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{PURLBASE}/RO_0002607"),
                    URIRef(f"{cl_term}"),
                )
            )

        # Nodes for each cell type, and marker and binary gene
        binary_genes = ast.literal_eval(row["binary_genes"])
        for gene in marker_genes + binary_genes:
            gene_term = f"GS_{gene}"

            # Cell_type_Class, EXPRESSES, Gene_Class
            # CL:0000000, RO:0002292, SO:0000704
            tuples.append(
                (
                    URIRef(f"{cl_term}"),
                    URIRef(f"{PURLBASE}/RO_0002292"),
                    URIRef(f"{PURLBASE}/{gene_term}"),
                )
            )

            # Gene_Class, PART_OF, Cell_type_Class
            # SO:0000704, BFO:0000050, CL:0000000
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{PURLBASE}/BFO_0000050"),
                    URIRef(f"{cl_term}"),
                )
            )

            # Gene_Class, EXPRESSED_IN, Anatomical_structure_Class
            # SO:0000704, RO:0002206, UBERON:0001062
            tuples.append(
                (
                    URIRef(f"{PURLBASE}/{gene_term}"),
                    URIRef(f"{PURLBASE}/RO_0002206"),
                    URIRef(f"{uberon_term}"),
                )
            )

    return tuples


def main(summarize=False):
    """Load manual author cell set to CL term mapping for the NSForest
    results from processing datasets corresponding to the Guo et
    al. 2023, Li et al. 2023, and Sikkema, et al. 2023 publications,
    create tuples consistent with schema v0.7, and write the result to
    a JSON file. If summarizing, retain the first row only, and
    include results in output.

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

        # Map NSForest results filename to manual author cell set to
        # CL term mapping filename, then load mapping results,
        # dropping "uuid" column in order to merge "uuid" column from
        # NSForest results
        author = re.search("results-([a-zA-Z]*)", nsforest_path.name).group(1)
        author_to_cl_path = Path(
            glob(
                str(NSFOREST_DIRPATH / f"cell-kn-mvp-map-author-to-cl-{author}-*.csv")
            )[-1]
        ).resolve()
        author_to_cl_results = (
            load_results(author_to_cl_path)
            .sort_values("author_cell_set", ignore_index=True)
            .drop(columns=["uuid"])
        )

        # Merge NSForest results with manual author cell set to CL
        # term mapping since author cell sets may not align exactly
        author_to_cl_results = author_to_cl_results.merge(
            nsforest_results[
                ["clusterName", "NSForest_markers", "binary_genes", "uuid"]
            ].copy(),
            left_on="author_cell_set",
            right_on="clusterName",
        )

        print(f"Creating tuples from {author_to_cl_path}")
        author_to_cl_tuples = create_tuples_from_author_to_cl(author_to_cl_results)
        if summarize:
            output_dirpath = TUPLES_DIRPATH / "summaries"
        else:
            output_dirpath = TUPLES_DIRPATH
        with open(
            output_dirpath / author_to_cl_path.name.replace(".csv", ".json"),
            "w",
        ) as f:
            data = {}
            if summarize:
                data["results"] = author_to_cl_results.to_dict()
            data["tuples"] = author_to_cl_tuples
            json.dump(data, f, indent=4)

        if summarize:
            break


if __name__ == "__main__":
    main(summarize=True)
    main()
