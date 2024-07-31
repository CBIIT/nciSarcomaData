setwd("./inst/extdata/drug")

cells=read.csv("export_cell_line.csv")
drugs=read.csv("export_compound.csv")
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
#
library(dplyr)

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
#
i=which(ok=="no"); length(i) # 227 # OK

results=cbind(avg_act,OK=ok)
write.csv(results,"avg_activities2.csv",row.names = F)
# now remove the i index
avg_act_filtered= avg_act[-i,c(1,2,4)]
dim(avg_act_filtered) # 28275 3

# can use tidyr spread to create repAct matrix , advantage speed
library(tidyr)
library(dplyr)
testresult= avg_act_filtered %>% spread(cell_line,mean)
dim(testresult) # 445 - 66

avg_act_filtered_flag= avg_act[-i,c(1,2,3)]
dim(avg_act_filtered_flag) # 28275 3
testresultflag= avg_act_filtered_flag %>% spread(cell_line,flag_ic50)
dim(testresultflag) # 445 - 66

stopifnot(identical(colnames(testresult),colnames(testresultflag)))
stopifnot(identical(testresult[,1],testresultflag[,1]))

# stats for filtering nsc
# nb.eq=apply(data.frame(testresultflag),1,function(x) length(which(as.character(x)=="=")))
# nb.more=apply(data.frame(testresultflag),1,function(x) length(which(as.character(x)==">")))
# nb.less=apply(data.frame(testresultflag),1,function(x) length(which(as.character(x)=="<")))
#
# sd = apply(data.frame(testresult[,2:66]),1,function(x) sd(x,TRUE))
#
# stats_results=cbind(testresult[,1],nb.eq=nb.eq,nb.more=nb.more,nb.less=nb.less,sd=sd)
#
# write.table(testresult,"sarcoma_drug_activity.txt",sep="\t",row.names = F)
# write.table(testresultflag,"sarcoma_drug_flag.txt",sep="\t",row.names = F)
# write.table(stats_results,"sarcoma_drug_filtering.txt",sep="\t",row.names = F)

# new : remove the control A549 from the stats
nb.eq2=apply(data.frame(testresultflag[,-4]),1,function(x) length(which(as.character(x)=="=")))
nb.more2=apply(data.frame(testresultflag[,-4]),1,function(x) length(which(as.character(x)==">")))
nb.less2=apply(data.frame(testresultflag[,-4]),1,function(x) length(which(as.character(x)=="<")))

sd2 = apply(data.frame(testresult[,c(2:3,5:66)]),1,function(x) sd(x,TRUE))
rg2 = apply(data.frame(testresult[,c(2:3,5:66)]),1,function(x) max(x,na.rm=T)-min(x,na.rm=T))
stats_results2=cbind(testresult[,1],nb.eq2=nb.eq2,nb.more2=nb.more2,nb.less2=nb.less2,sd2=sd2,range=rg2)
write.table(stats_results2,"sarcoma_drug_filtering_no_A549.txt",sep="\t",row.names = F)

# rg63 = apply(data.frame(testresult[,c(2:3,5:24,26:66)]),1,function(x) max(x,na.rm=T)-min(x,na.rm=T))
# length(which(rg63<0.5)) # 83 should be 100 ??

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

# check if all drugs are public in nci60 # now action to do now since all are published in Sarcoma ---
nci60drug=rownames(exprs(rcellminerData::drugData@act))
nci60drugpriv=rownames(exprs(rcellminerDataInt::drugData@act))

length(intersect(nci60drug,rownames(final))) # 128 public
length(intersect(nci60drugpriv,rownames(final))) #  263 both
chkgenes=setdiff(intersect(nci60drugpriv,rownames(final)),intersect(nci60drug,rownames(final)))
write.csv(chkgenes,"~/chkgenes.csv")
length(chkgenes) # 135 ALMOST ALL ARE PUBLIC EXCEPT ONE 761386
# plus 28 new NSC ??!!!
setdiff(rownames(final),nci60drugpriv)
# [1] "95100"  "727859" "751249" "751549" "754350" "754355" "754358" "754364" "754367" "755766" "755774" "755927"
# [13] "756644" "756658" "756661" "756874" "756875" "757333" "758490" "759674" "760842" "761432" "762381" "762598"
# [25] "763932" "764237" "765262" "766887
## end check ------------------------------------------------------------------------------------------


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

stopifnot(identical(rownames(final2),rownames(fdrugs)))
# merge metadata from nci60 and scarcoma compound metada
# #
# metanci60=rcellminerDataInt::drugData@act@featureData@data
# dim(metanci60)


#--------------------------------------------------------------------------------------------------
# GET DRUG ANNOTATION TABLE.
#--------------------------------------------------------------------------------------------------
if (!require(RMySQL))
  stop("Package RMySQL must be installed")

pswd="EpmbvC16"
con = dbConnect(dbDriver("MySQL"),
                host = "discovery.nci.nih.gov",
                user = "GBGdtb",
                port = 1521,
                password = pswd)

rs <- dbSendQuery(con, statement = paste("use", "nci60public2_2"))


query2 <- "Select  nsc, cas, name, testing_status, mechanism, pubchem_id, smiles, confidential_flag, failure_reason, prefix  from drug"
rs <- dbSendQuery(con, statement = query2)

drugAnnot <- fetch(rs, n = -1)
dbDisconnect(con)
n=dim(drugAnnot)[1] # 141046
pnsc=c()
for (i in 1:n){
  if (drugAnnot$prefix[i]=="S") pnsc=c(pnsc,drugAnnot$nsc[i]) else
    pnsc=c(pnsc,paste0(drugAnnot$prefix[i],drugAnnot$nsc[i]))
}

drugAnnot$pnsc = pnsc


# comments NSCs are not unique ...!!!! ADDed field "prefix"
length(drugAnnot$nsc)
# 141046
length(unique(drugAnnot$nsc))
# 136287
length(unique(drugAnnot$pnsc)) # 141046
rownames(drugAnnot)=pnsc


# merge metadata from nci60 and scarcoma compound metada

metadrugs=merge(fdrugs,drugAnnot,by.x=0,by.y=0)
dim(metadrugs) # 290 -15
setdiff(rownames(fdrugs),rownames(drugAnnot))
# "751549"

metadrugs=merge(fdrugs,drugAnnot,by.x=0,by.y=0, all.x=T)
dim(metadrugs) # 291 -15
colnames(metadrugs)
# [1] "Row.names"         "nsc.x"             "cas.x"             "drug_name"         "nsc.y"
# [6] "cas.y"             "name"              "testing_status"    "mechanism"         "pubchem_id"
# [11] "smiles"            "confidential_flag" "failure_reason"    "prefix"            "pnsc"

kk=c(2,41,43,117,118,287)

metadrugs$drug_name=as.character(metadrugs$drug_name)
metadrugs$drug_name[kk]=metadrugs$name[kk]

annot291=metadrugs[,c(2,3,4,8,9,10,11)]
annot291$pubchem_id[which(annot291$pubchem_id==-1)]=NA

rownames(annot291)=as.character(annot291$nsc.x)
annot291=annot291[rownames(final2),]
stopifnot(identical(rownames(final2), rownames(annot291)))
colnames(annot291)=c("NSC","CAS","NAME","FDA_STATUS","MOA","PUBCHEM_ID","SMILES")
##
library(rcellminer)


stopifnot(identical(colnames(final2), colnames(exprs(nciSarcomaData::molData@eSetList$exp))))

sarcomaMiame <- nciSarcomaData::molData@sampleData

actData <- ExpressionSet(as.matrix(final2))
featureData(actData) <- new("AnnotatedDataFrame", data=annot291)

repeatActData <- actData

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = sarcomaMiame)

save(drugData, file = "../../../data/drugData.RData")

# update rcellminercdbUtilsCDB with drug synonyms
