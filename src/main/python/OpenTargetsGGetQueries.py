gget_queries = {
    "target": {
        "purpose": "Obtain target attributes",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query target($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    dbXrefs {
        id
        source
    }
    proteinIds {
      id
      source
    }
    transcriptIds
    approvedSymbol
    approvedName
    associatedDiseases {
      count
      rows {
        score
        disease {
          id
          description
          dbXRefs
          name
        }
      }
    }
    knownDrugs {
      count
      rows {
        phase
        status
        drugId
        drugType
        diseaseId
        approvedSymbol
        ctIds
        approvedName
        mechanismOfAction
        drug {
          id
          description
          maximumClinicalTrialPhase
          isApproved
          synonyms
          tradeNames
          name
          indications {
            count
            rows {
              disease {
                id
                name
                description
              }
              maxPhaseForIndication
              references {
                source
                ids
              }
            }
          }
        }
      }
    }
    interactions {
      count
      rows {
        count
        score
        sourceDatabase
        targetA {
          proteinIds {
            id
          }
          id
          approvedSymbol
        }
        intABiologicalRole
        speciesA {
          taxonId
        }

        targetB {
          proteinIds {
            id
          }
          id
          approvedSymbol
        }
        intBBiologicalRole
        speciesB {
          taxonId
        }
        evidences {
          pubmedId
          evidenceScore
        }
      }
    }
    pharmacogenomics {
      variantRsId
      genotypeId
      genotype
      variantFunctionalConsequenceId
      variantFunctionalConsequence {
        label
      }
      drugs {
        drugId
        drugFromSource
        drug {
          name
        }
      }
      phenotypeText
      genotypeAnnotationText
      pgxCategory
      isDirectTarget
      evidenceLevel
      datasourceId
      literature
      haplotypeFromSourceId
      targetFromSourceId
      studyId
      datatypeId
      phenotypeFromSourceId
      variantId
      haplotypeId
    }
    tractability {
      label
      modality
      value
    }
    expressions {
      tissue {
        id
        label
        anatomicalSystems
        organs
      }
      rna {
        zscore
        value
        unit
        level
      }
    }
    depMapEssentiality {
      screens {
        depmapId
        expression
        geneEffect
        cellLineName
        diseaseCellLineId
        diseaseFromSource
        mutation
      }
      tissueId
      tissueName
    }
  }
}
""",
    },
    "diseases": {
        "purpose": "Duplicate and extend gget opentagets -r diseases command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query diseases($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    associatedDiseases {
      count
      rows {
        score
        disease {
          id
          description
          dbXRefs
          name
        }
      }
    }
  }
}
""",
    },
    "drugs": {
        "purpose": "Duplicate and extend gget opentagets -r drugs command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query diseases($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    knownDrugs {
      count
      rows {
        phase
        status
        drugId
        drugType
        diseaseId
        approvedSymbol
        ctIds
        approvedName
        mechanismOfAction
        drug {
          id
          description
          maximumClinicalTrialPhase
          isApproved
          synonyms
          tradeNames
          name
          indications {
            count
            rows {
              maxPhaseForIndication
              references {
                source
                ids
              }
            }
          }
        }
      }
    }
  }
}
""",
    },
    "interactions": {
        "purpose": "Duplicate and extend gget opentagets -r interactions command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query interactions($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    interactions {
      count
      rows {
        count
        score
        sourceDatabase
        targetA {
          proteinIds {
            id
          }
          id
          approvedSymbol
        }
        intABiologicalRole
        speciesA {
          taxonId
        }

        targetB {
          proteinIds {
            id
          }
          id
          approvedSymbol
        }
        intBBiologicalRole
        speciesB {
          taxonId
        }
        evidences {
          pubmedId
          evidenceScore
        }
      }
    }
  }
}
""",
    },
    "pharmacogenetics": {
        "purpose": "Duplicate and extend gget opentagets -r pharmacogenetics command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query pharmacogenetics($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    pharmacogenomics {
      variantRsId
      genotypeId
      genotype
      variantFunctionalConsequenceId
      variantFunctionalConsequence {
        label
      }
      drugs {
        drugId
        drugFromSource
        drug {
          name
        }
      }
      phenotypeText
      genotypeAnnotationText
      pgxCategory
      isDirectTarget
      evidenceLevel
      datasourceId
      literature
      haplotypeFromSourceId
      targetFromSourceId
      studyId
      datatypeId
      phenotypeFromSourceId
      variantId
      haplotypeId
    }
  }
}
""",
    },
    "tractability": {
        "purpose": "Duplicate and extend gget opentagets -r tractability command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query tractability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    tractability {
      label
      modality
      value
    }
  }
}
""",
    },
    "expression": {
        "purpose": "Duplicate and extend gget opentagets -r expression command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query expression($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    expressions {
      tissue {
        id
        label
        anatomicalSystems
        organs
      }
      rna {
        zscore
        value
        unit
        level
      }
    }
  }
}
""",
    },
    "depmap": {
        "purpose": "Duplicate and extend gget opentagets -r depmap command",
        "variables": {
            "ensemblId": "ENSG00000169252",
        },
        "query_string": """
query depmap($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    depMapEssentiality {
      screens {
        depmapId
        expression
        geneEffect
        cellLineName
        diseaseCellLineId
        diseaseFromSource
        mutation
      }
      tissueId
      tissueName
    }
  }
}
""",
    },
}
