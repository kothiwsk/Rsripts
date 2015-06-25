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
     --args1=filename   - string, peterson kinase list
     --args2=filename   - string, selectivity data for peterson
     --args3=filename   - string, common kinases between peterson and zarrinkar  datasets
     --args4=filename   - string, predicted output file
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
selectivity_file<-argsL$args2
overlap_file<-argsL$args3
predicted_file<-argsL$args4
#blind_file<-argsL$args5
#cutoff_file<-argsL$args6


kinases<-read.table(kinases_file,header=FALSE,sep="\t")
selectivity<-read.table(selectivity_file,header=FALSE)
overlap<-read.table(overlap_file,header=FALSE,sep="\t",stringsAsFactors=FALSE)
predictions<-readLines(predicted_file)
#experimental<-read.table(blind_file,header=TRUE,sep=" ")
#cutoffs.lines<-scan(cutoff_file,skip=143,nlines=1,sep='\n',what="character")
#head(kinases,n=20)
#head(objects,n=20)
#head(accession,n=20)
#head(predictions)
#head(experimental)

zarrinkar_output<-function(suffix,data,name){
 files=""
 for(item in data){
  suffix_sdf<-paste(item,".sdf",sep="")
  files<-paste(files,suffix_sdf,sep=" ")
 }

 command_line<-paste("cat ",files," > ",name,"/",suffix,".sdf",sep="")
 system(command_line)
}

for(name in kinases$V1){
 kinase_idx<-which(kinases$V1==name)
 peterson_actives<-selectivity[which(selectivity[[kinase_idx+1]]==1),]$V1
 files=""
 if( !name %in% overlap$V1) next
 if( length(peterson_actives) < 1) next
 for(item in peterson_actives){
  file_name<-strsplit(item,"b/",fixed=TRUE)[[1]][2]
  file_name<-gsub(".csv",".sdf",file_name,fixed=TRUE)
  files<-paste(files,file_name,sep=" ")
 }
 dir_name<-chartr("/ ","__",name)

 system(paste("mkdir ",dir_name,sep=""))
 command_line<-paste("cat ",files," > ",dir_name,"/","peterson_actives.sdf",sep="")
 system(command_line)

 line_number<-(which(grepl(name,predictions,fixed=TRUE)==TRUE))
 if(length(line_number)>1){
  for(j in line_number){
   check<-strsplit(predictions[j],"\"",fixed=TRUE)[[1]][2]
   if(check==name){
    line_number<-j
    break
   }
  }
 }
 cur_kinase_predictions<-scan(predicted_file,skip=line_number,nlines=55,sep="\n",what="character")
 predicted_dataframe<-data.frame(matrix(ncol=5,nrow=0))
 for(i in 2:length(cur_kinase_predictions)){
  value<-unlist(strsplit(cur_kinase_predictions[i]," +"))
  predicted_dataframe[i-1,]<-value
 }
 colnames(predicted_dataframe)<-c("idx","mol_name","prob","exp","pred")
 zarrinkar_actives<-predicted_dataframe[which(predicted_dataframe$exp==1),]$mol_name
 zarrinkar_tp<-predicted_dataframe[which(predicted_dataframe$exp==1 & predicted_dataframe$pred == 1),]$mol_name
 zarrinkar_fp<-predicted_dataframe[which(predicted_dataframe$exp==0 & predicted_dataframe$pred == 1),]$mol_name
 zarrinkar_fn<-predicted_dataframe[which(predicted_dataframe$exp==1 & predicted_dataframe$pred == 0),]$mol_name

 if(length(zarrinkar_actives) > 0) zarrinkar_output(deparse(substitute(zarrinkar_actives)),zarrinkar_actives,dir_name)
 if(length(zarrinkar_tp) > 0) zarrinkar_output(deparse(substitute(zarrinkar_tp)),zarrinkar_tp,dir_name)
 if(length(zarrinkar_fp) > 0) zarrinkar_output(deparse(substitute(zarrinkar_fp)),zarrinkar_fp,dir_name)
 if(length(zarrinkar_fn) > 0) zarrinkar_output(deparse(substitute(zarrinkar_fn)),zarrinkar_fn,dir_name)
}
