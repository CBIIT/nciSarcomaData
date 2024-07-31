library(rcellminer) # Load dependency functions 
library(nciSarcomaData) # Loads data

# LOAD DATA ----
## Data
sarcomaAct <- exprs(getAct(nciSarcomaData::drugData))
sarcomaExp <- getAllFeatureData(nciSarcomaData::molData)[["exp"]] # "mir" also exists for miRNA data

## Annotations
sarcomaSampleAnnot <- getSampleData(nciSarcomaData::molData)
sarcomaDrugAnnot <- getFeatureAnnot(nciSarcomaData::drugData)[["drug"]]

## Check data loaded correctly
dim(sarcomaAct)
# [1] 291  71
dim(sarcomaExp)
# [1] 18432    71

## Check metadata loaded correctly
dim(sarcomaSampleAnnot)
# [1] 71  6
dim(sarcomaDrugAnnot)
# [1] 291   7

write.csv(sarcomaDrugAnnot, "./inst/extdata/drug/export_annotations.csv")

