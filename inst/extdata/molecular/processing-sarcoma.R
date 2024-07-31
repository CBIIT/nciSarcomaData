
setwd("./inst/extdata/molecular")
vtitle=as.character(read.csv("AromaTrExAnnotate.csv",header=F, stringsAsFactors = F)[1,])

m=read.csv("AromaTrExAnnotate.csv",stringsAsFactors = F)
dim(m) # 18690 - 83
colnames(m)=vtitle

length(m[,1]) # 18690
length(unique(m[,1])) # 18690

length(m$GeneSymbol) # 18690
length(unique(m$GeneSymbol)) # 17805

length(m$EntrezID) # 18690
length(unique(m$EntrezID)) # 17805

nogeneid=which(m$EntrezID==-1); length(nogeneid) # 258
nogenen=which(is.na(m$GeneSymbol)); length(nogenen) # 258
stopifnot(identical(nogeneid,nogenen))

m2=m[-nogeneid,]; dim(m2) # 18432 83
length(unique(m2$EntrezID)) # 17804 (removed -1)

sort(table(m2$EntrezID),decreasing = T)

# 441056    728096     24150    122183    653219      4097      6133      6171      9085    339044    414060
# 7         6         5         5         5         4         4         4         4         4         4

i=which(m2$EntrezID==441056) ## DUX4L4 PSEUDO GENE CHR 4 AND 10
m2[i,]

## proceed like rnaseq
ID=m2$GeneSymbol; length(unique(ID)) # 17804
vdup=unique(ID[duplicated(ID)]); length(vdup) #565

# fix issue for non unique gene
i=which(ID %in% vdup); length(i); # 1193
length(unique(ID[i])) # 565

ID[i]=paste0(ID[i],":",m2$TranscriptClusterID[i])
length(unique(ID)) # 18432 OK
rownames(m2)=ID

expdata=m2[,2:71]
rownames(expdata)=ID
dim(expdata) # 18432 - 70

metadata=m2[,c(1,72:83)]
rownames(metadata)=ID
dim(metadata) # 18432 - 13
##
# 65 cell lines with activites, need to check also micro rna cell lines and do jointure!
##
mirann=read.csv("AnnotatedmiRNA.csv",stringsAsFactors = F,sep = "\t")
length(mirann$GeneName); length(unique(mirann$GeneName)) # 800 -800

mirdata=read.delim("miRNAOnlyDESeqAdjustedCountLog2Plus1.tsv",sep="\t",stringsAsFactors = F)
dim(mirdata)

# spread data horizontally
library(tidyr)
library(dplyr)
mymirdata= mirdata %>% spread(CellLineName,Expression)
dim(mymirdata) # 800 - 59
rownames(mymirdata)=mymirdata[,1]
mymirdata=mymirdata[,-1]

rownames(mirann)=mirann$GeneName
mirann=mirann[rownames(mymirdata),]

stopifnot(identical(rownames(mymirdata), mirann$GeneName))
mirann=mirann[,c(4,1:3,5:9)]
colnames(mirann)[1]="ID"

cell.exp=colnames(expdata); cell.mir=colnames(mymirdata)
length(cell.exp); length(cell.mir) # 70, 58

length(intersect(cell.exp,cell.mir)) # 58 OK
setdiff(cell.exp,cell.mir)
# [1] "CHA-59"         "CHLA-25"        "CHLA-258"       "CHLA-32"        "ES-2"           "ES-6"
# [7] "KHOS-240S"      "KHOS-312H"      "Rh28 PX11/LPAM" "Hs 706.T"       "Hs 729"         "Hs 132.T"
# to add to mirRNA data + A549/ATCC >> 71 cell lines
##
cell.drugs=as.character(read.delim("../drug/sarcoma_drug_activity.txt",stringsAsFactors = F,header = F)[1,-1])
# 65 cell lines

length(intersect(cell.exp,cell.drugs)) # 64 minus 1

setdiff(cell.drugs,cell.exp)
# A549/ATCC to add in expression data >> total 71 cell lines

setdiff(cell.exp,cell.drugs)
# "hMSC"    "NHAC-kn" "NHDF"    "NHOst"   "SkMC"    "JJ012" to add to drug data >> 71
## ---- update expdata and mirna data

cell.all = sort(union(cell.exp,cell.drugs)); length(cell.all)
expdata[,"A549/ATCC"]=NA
dim(expdata)
expdata=expdata[,cell.all]

#
vnew=c("CHA-59","CHLA-25","CHLA-258","CHLA-32","ES-2","ES-6","KHOS-240S","KHOS-312H","Rh28 PX11/LPAM",
       "Hs 706.T","Hs 729","Hs 132.T","A549/ATCC")
mymirdata[,vnew]=NA
dim(mymirdata)
mymirdata=mymirdata[,cell.all]

stopifnot(identical(colnames(expdata),colnames(mymirdata)))

## miame : sample information
minfo=read.csv("../drug/export_cell_line.csv",stringsAsFactors = F,row.names = 1)
dim(minfo) # 71 - 4

minfo=minfo[cell.all,]
stopifnot(identical(colnames(expdata),rownames(minfo)))
#--------------------------------------------------------------------------------------------------
# Make MIAME object
#--------------------------------------------------------------------------------------------------
minfo$OncoTree2[which(minfo$OncoTree2=="")]=NA
minfo$OncoTree3[which(minfo$OncoTree3=="")]=NA

sarcomaMiame <- new("MIAME", name="NCI/DTP Sarcoma project", lab="NCI/DTP",
                 samples=list(Name = rownames(minfo),
                              TissueType = minfo$panel,
                              OncoTree1  = minfo$OncoTree1,
                              OncoTree2  = minfo$OncoTree2,
                              OncoTree3  = minfo$OncoTree3,
                              OncoTree4  = replicate(71,NA)))
#--------------------------------------------------------------------------------------------------
# Make and save MolData object
#--------------------------------------------------------------------------------------------------
library(rcellminer)
expData <- ExpressionSet(as.matrix(expdata))
stopifnot(identical(rownames(expdata), rownames(metadata)))
featureData(expData) <- new("AnnotatedDataFrame", data=metadata)
# ok
mirData <- ExpressionSet(as.matrix(mymirdata))
stopifnot(identical(rownames(mymirdata), rownames(mirann)))
featureData(mirData) <- new("AnnotatedDataFrame", data=mirann)

## start
sarcomaESetList <- list()

sarcomaESetList[["exp"]] <- expData
sarcomaESetList[["mir"]] <- mirData

molData <- new("MolData", eSetList = sarcomaESetList, sampleData = sarcomaMiame)

save(molData, file = "../../../data/molData.RData")

## add drug data
## update rcellminerUtilsCDB package




