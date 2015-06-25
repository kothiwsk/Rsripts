
getaverage<-function(file){
 curfile<-read.table(filename,header=TRUE)
 curfile[1]<-NULL 
 cmean<-colMeans(curfile)
 mean(cmean)
}

directories<-dir()[file.info(dir())$isdir]
mean_data<-data.frame(matrix(ncol=6,nrow=length(directories)))
for(i in 1:length(directories)){
 d<-directories[i]
 files=list.files(path=paste("./",d,sep=""),pattern="*.txt")
 kinase_col<-list(d)
 for(f in files){
  filename<-paste("./",d,"/",f,sep="")
  kinase_col<-c(kinase_col,getaverage(filename))
 }
 mean_data[i,]<-kinase_col
}
print(mean_data)
