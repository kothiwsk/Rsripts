#!/usr/bin/Rscript



## Collect arguments
args<-commandArgs(TRUE)

## Default settings when no arguments passed
if(length(args) < 2){
 args<-c("--help")
}

if("--help" %in% args){
 cat("
     This script calculates area under the curve from *.data files output by bcl neural network in the result folder

     Arguments:
     --args1=filename   - string, file containing bioactivity data
     --args2=filename   - string, file containing protein uniprot id
     --args3=filename   - string, file containing 
     --help             - print this text
    
     Example:
     ./calculated_area_under_curve.R --arg1=kinase_list.txt --arg2=result_label5.obj --arg3=selectivity_data.txt
     \n\n
     ")
 q(save="no")
}

## Parse arguments(expect the form --arg=value)
parseArgs<-function(x) strsplit(sub("^--","",x),"=")
argsDF<-as.data.frame(do.call("rbind",parseArgs(args)))
argsL<-as.list(as.character(argsDF$V2))
names(argsL)<-argsDF$V1

bioactivity<-argsL$arg1
proteins<-argsL$arg2

bioactivity_data<-read.table(bioactivity,header=T,sep="\t",fill=T)
protein_ids<-read.table(proteins,header=FALSE)
names(bioactivity_data)
head(bioactivity_data,nrows=10)
molecules<-unique(bioactivity_data$CMPD_CHEMBLID)
proteins<-unique(protein_ids$V1)

protein_string<-paste(proteins,collapse=' ')
frame<-data.frame("molecule",protein_string)
write.table(frame,"molecular_activity.txt",sep=" ",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
## find activity of molecule of interest for proteins listed in  proteins 
find_activity <- function(molecule, bioactivity_data,proteins){
 subset_data<-subset(bioactivity_data,CMPD_CHEMBLID == molecule)
 
# get activity values
 subset_data_2=subset_data[match(proteins,subset_data$PROTEIN_ACCESSION),]$PCHEMBL_VALUE
 subset_data_2[!is.na(subset_data_2)] <- 1
 subset_data_2[is.na(subset_data_2)] <- 0
 activity<-paste(subset_data_2,collapse=' ')

 data<-data.frame(molecule,activity,stringsAsFactors=FALSE)
 write.table(data,"molecular_activity.txt",sep=" ",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
}

for(i in 1:(length(molecules))) find_activity(molecules[i],bioactivity_data,proteins)
