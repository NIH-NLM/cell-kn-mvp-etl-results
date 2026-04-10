from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path("../../src/main/python").resolve()))

from ExternalApiResultsFetcher import get_opentargets_results
from LoaderUtilities import get_efo_to_mondo_map, map_efo_to_mondo

efo2mondo = get_efo_to_mondo_map()

all_genes_human = pd.read_csv("all_genes_human.csv")["x"].to_list()

with open("fetch_all_genes_human.csv", "w") as fp:
    fp.write(
        ",".join(
            [
                "ensembl_id",
                "approvedSymbol",
                "drugId",
                "phase",
                "diseaseId",
                "\n",
            ]
        )
    )

    n_gn = 1000
    for i_gn in range(0, len(all_genes_human), n_gn):
        opentargets_results = get_opentargets_results(
            all_genes_human[i_gn : i_gn + n_gn],
            opentargets_path=Path(__file__).resolve().parent
            / f"opentargets-{i_gn:05d}.json",
        )
        for ensembl_id in opentargets_results["gene_ensembl_ids"]:
            for drug in opentargets_results[ensembl_id]["drugs"]:
                diseaseId = drug["diseaseId"]
                if "EFO" in diseaseId:
                    diseaseId = map_efo_to_mondo(diseaseId, efo2mondo)
                if not diseaseId or "MONDO" not in diseaseId:
                    continue
                fp.write(
                    ",".join(
                        [
                            ensembl_id,
                            drug["approvedSymbol"],
                            drug["drugId"],
                            str(drug["phase"]),
                            diseaseId,
                            "\n",
                        ]
                    )
                )
