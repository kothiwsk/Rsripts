files=list.files(pattern="CHEMBL*")

for(f in files){
 x<-readLines(f)
 chemblid<-strsplit(f,".",fixed=TRUE)[[1]][1]
 chemblid<-paste("M  END\n> <CHEMBLID>\n",chemblid,sep="")
 y<-gsub("M  END",paste(chemblid,"\n",sep=""),x)
 
 cat(y,file=paste(f,".added",sep=""),sep="\n")
}
