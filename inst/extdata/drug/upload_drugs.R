library(dplyr)
library(tidyr)
library(rcellminer)

setwd("./inst/extdata/drug")

cells=read.csv("export_cell_line.csv")
drugs=read.csv("export_compound.csv")
drugAnnot=read.csv("export_annotations.csv", row.names = 1)
experims=read.csv("export_experiment.csv",stringsAsFactors = F)
activity=read.csv("export_experiment_curve.csv",stringsAsFactors = F)

dim(cells);dim(drugs);dim(experims);dim(activity)
# [1] 71 5
# [1] 445   3
# [1] 33503     4
# [1] 33503     7

act_detail=merge(activity,experims,by.x=1,by.y=1); dim(act_detail)
# write.csv(act_detail,"Sarcoma_drug_experiment_activity.csv",row.names = F)
act_light=act_detail[,c(1,6,7,8,10)]
act_light2=act_light[order(act_light$nsc,act_light$cell_line),]

avg_act <-  act_light2 %>%
            group_by(nsc,cell_line,flag_ic50) %>%
            summarize(mean = mean(log_ic50, na.rm = TRUE))

dim(avg_act) #28502 - 4

## checking
stat_act <- avg_act %>%
  group_by(nsc,cell_line) %>%
  summarize(count = n())
# 1 group have 3 others have 2 averges
length(which(stat_act$count>1))
#[1] 226

# mark duplicates, keep the one with flag "=" if found

ok=replicate(nrow(avg_act),"y")
curr.nsc=avg_act$nsc[1] ; curr.cell=avg_act$cell_line[1]

for (k in 2:nrow(avg_act)){
    if ( avg_act$nsc[k]==curr.nsc & avg_act$cell_line[k]==curr.cell) {
      if (avg_act$flag_ic50[k] %in% c("<",">")) {
        ok[k]="no"
      } else
        { ok[k-1]="no" }
    } else {
      curr.nsc = avg_act$nsc[k]
      curr.cell = avg_act$cell_line[k]
    }
}

i=which(ok=="no"); length(i) # 227 # OK

results=cbind(avg_act,OK=ok)
write.csv(results,"avg_activities2.csv",row.names = F)
# now remove the i index
avg_act_filtered= avg_act[-i,c(1,2,4)]
dim(avg_act_filtered) # 28275 3

# can use tidyr spread to create repAct matrix , advantage speed
testresult= avg_act_filtered %>% spread(cell_line,mean)
dim(testresult) # 445 - 66

avg_act_filtered_flag= avg_act[-i,c(1,2,3)]
dim(avg_act_filtered_flag) # 28275 3
testresultflag= avg_act_filtered_flag %>% spread(cell_line,flag_ic50)
dim(testresultflag) # 445 - 66

stopifnot(identical(colnames(testresult),colnames(testresultflag)))
stopifnot(identical(testresult[,1],testresultflag[,1]))

# new : remove the control A549 from the stats
nb.eq2=apply(data.frame(testresultflag[,-4]),1,function(x) length(which(as.character(x)=="=")))
nb.more2=apply(data.frame(testresultflag[,-4]),1,function(x) length(which(as.character(x)==">")))
nb.less2=apply(data.frame(testresultflag[,-4]),1,function(x) length(which(as.character(x)=="<")))

sd2 = apply(data.frame(testresult[,c(2:3,5:66)]),1,function(x) sd(x,TRUE))
rg2 = apply(data.frame(testresult[,c(2:3,5:66)]),1,function(x) max(x,na.rm=T)-min(x,na.rm=T))
stats_results2=cbind(testresult[,1],nb.eq2=nb.eq2,nb.more2=nb.more2,nb.less2=nb.less2,sd2=sd2,range=rg2)
write.table(stats_results2,"sarcoma_drug_filtering_no_A549.txt",sep="\t",row.names = F)

length(which(nb.eq2==0))
# 68
IQR(nb.eq2) # 57
quantile(nb.eq2)
# 0%  25%  50%  75% 100%
# 0    3   39   60   64

# FILTERING: use both nb.eq2 > 10and then sd2 > 0.2
index=which(nb.eq2>10 & sd2 > 0.2); length(index) # 291
final=as.data.frame(testresult[index,])
rownames(final)=as.character(final$nsc)
final=final[,-1]
dim(final) # 291 - 65
write.table(final,"sarcoma_drug_activity.txt",sep="\t",row.names = F)

# NEED TO APPLY NEGATIVE TO LOG_IC50 !!!
final2 = final * -1

## stop here
vnew =c("hMSC","NHAC-kn","NHDF","NHOst","SkMC","JJ012") # to add to drug data >> 71
## retrieve cell line names from molecular data
cell.all=nciSarcomaData::molData@sampleData@samples$Name

final2[,vnew]=NA
final2=final2[,cell.all]

#-----------------------------------------------------
# metadata + no repeat activity
rownames(drugs)=as.character(drugs$nsc)
fdrugs=drugs[rownames(final2),]; dim(fdrugs)

# Final check that fdrugs (final data) matches processing order of final2
stopifnot(identical(rownames(final2),rownames(fdrugs)))

#--------------------------------------------------------------------------------------------------
# GET DRUG ANNOTATION TABLE.
#--------------------------------------------------------------------------------------------------
# merge metadata from nci60 and sarcoma compound metadata
metadrugs=merge(fdrugs,drugAnnot,by.x=c("nsc", "cas", "drug_name"),by.y=c("nsc", "cas", "name"), all.x = T)
dim(metadrugs) # 290 -15
setdiff(rownames(fdrugs),rownames(drugAnnot))

kk=c(2,41,43,117,118,287)

metadrugs$drug_name=as.character(metadrugs$drug_name)

annot291=metadrugs
annot291$pubchem_id[which(annot291$pubchem_id==-1)]=NA

rownames(annot291)=as.character(annot291$nsc)
annot291=annot291[rownames(final2),]
setdiff(rownames(final2),rownames(annot291))
setdiff(rownames(annot291),rownames(final2))
stopifnot(identical(rownames(final2), rownames(annot291)))
colnames(annot291)=c("NSC","CAS","NAME","FDA_STATUS","MOA","PUBCHEM_ID","SMILES")

stopifnot(identical(colnames(final2), colnames(exprs(nciSarcomaData::molData@eSetList$exp))))

sarcomaMiame <- nciSarcomaData::molData@sampleData

actData <- ExpressionSet(as.matrix(final2))
featureData(actData) <- new("AnnotatedDataFrame", data=annot291)

repeatActData <- actData

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = sarcomaMiame)

save(drugData, file = "../../../data/drugData.RData")
