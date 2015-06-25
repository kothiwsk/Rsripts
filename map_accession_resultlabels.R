#!/usr/bin/Rscript

## Collect arguments
args<-commandArgs(TRUE)

## Default settings when no arguments passed
if(length(args) < 3){
 args<-c("--help")
}

if("--help" %in% args){
 cat("
     This script calculates area under the curve from *.data files output by bcl neural network in the result folder

     Arguments:
     --args1=filename   - string, filename containing all cases
     --args2=filename   - string, filename containing indices to indicate which cases are used
     --args3=filename   - string, filename containing cases with accession numbers
     --args4=filename   - string, filename containing predicted activities
     --args5=filename   - string, filename containing blind dataset activities
     --args6=filename   - string, filename containing cutoffs for different cases
     --help             - print this text
    
     Example:
     Rscript ~/scripts/Rscripts/map_accession_resultlabels.R --args1=kinase_list.txt --args2=result_label5.obj --args3=peterson_proteins_accession.txt --args4=predictions.txt --args5=molecular_activity.txt --args6=log_merge.txt
	(files in /home/kothiwsk/ddr_project/kinase_panel/Zarrinkar_NatBiotech/peterson_dataset/blind_dataset_validation/result5_indices)
     \n\n
     ")
 q(save="no")
}

## Parse arguments(expect the form --arg=value)
parseArgs<-function(x) strsplit(sub("^--","",x),"=")
argsDF<-as.data.frame(do.call("rbind",parseArgs(args)))
argsL<-as.list(as.character(argsDF$V2))

names(argsL)<-argsDF$V1

kinases_file<-argsL$args1
obj_file<-argsL$args2
accession_file<-argsL$args3
predicted_file<-argsL$args4
blind_file<-argsL$args5
cutoff_file<-argsL$args6


kinases<-read.table(kinases_file,header=FALSE,sep="\t")
objects<-read.table(obj_file,header=FALSE)
accession<-read.table(accession_file,header=FALSE,sep="\t",stringsAsFactors=FALSE)
predictions<-read.csv(predicted_file,header=FALSE)
experimental<-read.table(blind_file,header=TRUE,sep=" ")
cutoffs.lines<-scan(cutoff_file,skip=143,nlines=1,sep='\n',what="character")
#head(kinases,n=20)
#head(objects,n=20)
#head(accession,n=20)
#head(predictions)
#head(experimental)

## get cutoffs for all the cases
string1<-strsplit(cutoffs.lines,"Result",fixed=TRUE)
cutoff_data<-data.frame(index=numeric(),cutoff=numeric(),stringsAsFactors=FALSE)
for(i in 2:(length(string1[[1]]))){
 value<-strsplit(string1[[1]][i],"cutoff ",fixed=TRUE)[[1]][2]
 test<-gsub("\t","",value,fixed=TRUE)
 cutoff_data[i-1,]<-list(i-1,test)
}
## get subset of molecules common between experimental and prediction files
experimental_subset<-experimental[match(gsub(" *$","",predictions[,c("V1")]),experimental$molecule),]
experimental_subset<-experimental_subset[complete.cases(experimental_subset),]
experimental_subset<-experimental_subset[!duplicated(experimental_subset),]
prediction_subset<-predictions[match(experimental_subset$molecule,gsub(" *$","",predictions$V1)),]
library(pROC)
## get mappings
objects_correctedindex<-objects[,1]+1
kinases_used<-kinases$V1[objects_correctedindex]
subset_accession<-accession[match(kinases_used,accession$V1),]$V2
mol_names<-experimental_subset[,c("molecule")]
for(i in 1:length(subset_accession)){
 protein_name<-subset_accession[i]
 if(!is.na(protein_name)){
  common_name<-accession[which(accession$V2==protein_name),]$V1
  pred<-prediction_subset[,i+1]
  exp<-experimental_subset[,c(protein_name)]
  cur_kinase<-data.frame(mol_names,pred,exp)
  cur_kinase<-cur_kinase[with(cur_kinase,order(-pred)),]
  test<-transform(cur_kinase,x = as.numeric(pred) > as.numeric(cutoff_data$cutoff[i]))
  test$x<-test$x*1


  rocval<-auc(controls=test$x[test$exp==0],cases=test$x[test$exp==1])
  cat(common_name,rocval[1],"\n")
  #print(common_name)
  #print(test)
 }
}
