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
     --args1=filename   - string, filename containing all output node identifier on seperate lines
     --args2=filename   - string, filename containing indices of output nodes for which AUC needs to be calculated, with indices on new lines
     --args3=filename   - string, filename containing activity data of training datapoints on new lines
     --help             - print this text
    
     Example:
     ./calculated_area_under_curve.R --arg1=kinase_list.txt --arg2=er=result_label5.obj --arg3=selectivity_data.txt
     \n\n
     ")
 q(save="no")
}

## Parse arguments(expect the form --arg=value)
parseArgs<-function(x) strsplit(sub("^--","",x),"=")
argsDF<-as.data.frame(do.call("rbind",parseArgs(args)))
argsL<-as.list(as.character(argsDF$V2))

names(argsL)<-argsDF$V1
library(MESS)

kinases<-argsL$args1
obj<-argsL$args2
binding_data<-argsL$args3
kinase_list<-read.table(kinases,header=FALSE,sep="\t")
objects<-read.table(obj,header=FALSE)
data = list.files(pattern="*.data$")
bdata<-read.table(binding_data,header=FALSE,sep=" ")

calculate_AUC <- function(name, data,kinase,obj,binding_data){
 test<-strsplit(name,"_")
 test<-strsplit(test[[1]][4],".",fixed=TRUE)
 test<-strsplit(test[[1]][1],"l",fixed=TRUE)[[1]][2]
 #index in filename starts from 0 but objects when read by R, begin at line 1.
 kiname<-kinase$V1[obj$V1[strtoi(test)+1]+1]
 # first column in binding data is CAS number so 0th kinase starts at 2
 nth_kinase<-binding_data[obj$V1[strtoi(test)+1]+2]
  
 value = format(round(auc(data$FPR,data$TPR),2))
 data<-data.frame(kiname,value,sum(nth_kinase),stringsAsFactors=FALSE)
 write.table(data,"area_under_curve.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
}

for(i in 1:(length(data))) calculate_AUC(data[i],read.table(data[i],header=TRUE),kinase_list,objects,bdata)






#args<-commandArgs(TRUE)
#print (args)
#avge_min<- args[1]
#output_file<- args[2]
## load the necessary libraries
#library(nlme)
#
## set the output file
#png(output_file, 400, 300)
##get(getOption("device"))()
#
##first plot
#data<-read.table(avge_min,header=TRUE)
#
#par(mar=c(4,4,1.5,0)+.1)
#plot(x=data$bin,y=data$bcl,ylim=c(0,3),lty=1,xlab='number of rotatable bonds',ylab='average minimum rmsd (in Angstrom)',lwd=0.75,type='l',xaxp=c(0,length(data$bcl),length(data$bcl)),cex.lab=1.2, font.axis=1.5, las=2, xaxt='n')
#par(new=T)
#plot(x=data$bin,y=data$moe,ylim=c(0,3),lty=2,xaxt='n',axes=F,ylab='',lwd=1,type='l',xaxt='n',yaxt='n',xlab='')
#par(new=T)
#plot(x=data$bin,y=data$frog,ylim=c(0,3),lty=3,xaxt='n',axes=F,ylab='',lwd=1,type='l',xaxt='n',yaxt='n',xlab='')
#par(new=T)
#plot(x=data$bin,y=data$confab,ylim=c(0,3),lty=4,xaxt='n',axes=F,ylab='',lwd=1,type='l',xaxt='n',yaxt='n',xlab='')
#par(new=T)
#plot(x=data$bin,y=data$confimport,ylim=c(0,3),lty=5,xaxt='n',axes=F,ylab='',lwd=1,type='l',xaxt='n',yaxt='n',xlab='')
#axis(1,at=1:4,lab=c("0-3","4-6","6-10",">10"))
#legend("topleft",c("bcl","confab","confimport","frog","moe"),lty=c(1,4,5,3,2),lwd=c(0.75,0.75,0.75,0.75,0.75),title=("method"))
#graphics.off()
## close the output file
#
## unload the libraries
#detach("package:nlme")
