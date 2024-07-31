library(rcellminer)
library(stringr)
tmpEnv <- new.env()
##
load("data/molData.RData", envir = tmpEnv)
#load("data/drugData.RData", envir = tmpEnv)

nciSarcomaMiame <- tmpEnv$molData@sampleData

## nciSclcESetList=list()
nciSarcomaESetList=tmpEnv$molData@eSetList

#Sarcoma0Exp <- rcellminer::getAllFeatureData(tmpEnv$molData)[["exp"]]
SarcomaExp <- exprs(tmpEnv$molData@eSetList$exp)

gene_exp_info = tmpEnv$molData@eSetList$exp@featureData@data
dim(SarcomaExp); dim(gene_exp_info) #  17804    77  # 17804    13
cells <- colnames(SarcomaExp)
stopifnot(identical(rownames(SarcomaExp),rownames(gene_exp_info)))

## first add an ID column to gene info

gene_exp_info = cbind(ID=rownames(gene_exp_info),gene_exp_info)

##  updated information on cell lines
## miame : sample information
minfo=read.csv("inst/extdata/drug/export_cell_line_curated.csv",stringsAsFactors = F,row.names = 1)
dim(minfo) # 71 - 5

minfo=minfo[cells,]
stopifnot(identical(colnames(SarcomaExp),rownames(minfo)))
#--------------------------------------------------------------------------------------------------
# Make MIAME object
#--------------------------------------------------------------------------------------------------
minfo$OncoTree2[which(minfo$OncoTree2=="")]=NA
minfo$OncoTree3[which(minfo$OncoTree3=="")]=NA
minfo$OncoTree4[which(minfo$OncoTree4=="")]=NA

sarcomaMiame <- new("MIAME", name="NCI/DTP Sarcoma project", lab="NCI/DTP",
                 samples=list(Name = rownames(minfo),
                              TissueType = minfo$panel,
                              OncoTree1  = minfo$OncoTree1,
                              OncoTree2  = minfo$OncoTree2,
                              OncoTree3  = minfo$OncoTree3,
                              OncoTree4  = minfo$OncoTree4))
#--------------------------------------------------------------------------------------------------
# Make and save MolData object
#--------------------------------------------------------------------------------------------------

expData <- ExpressionSet(SarcomaExp)
stopifnot(identical(rownames(SarcomaExp), rownames(gene_exp_info)))
featureData(expData) <- new("AnnotatedDataFrame", data=gene_exp_info)
# ok

## start


nciSarcomaESetList[["exp"]] <- expData


molData <- new("MolData", eSetList = nciSarcomaESetList, sampleData = sarcomaMiame)

save(molData, file = "data/molData.RData")

## update sample info in  drug data

load("data/drugData.RData", envir = tmpEnv)

actData = tmpEnv$drugData@act
repeatActData = tmpEnv$drugData@repeatAct
stopifnot(identical(colnames(actData),cells))

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = sarcomaMiame)

save(drugData, file = "data/drugData.RData")

##




