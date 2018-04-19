args<-commandArgs()


myfile<-args[5]
mymethod<-args[6]
myoutput<-args[7]
min_miRNA<-as.numeric(args[8])
min_gene<-as.numeric(args[9])

mydata<-read.delim(header=F,stringsAsFactors=F,myfile)
mydata[mydata$V2==0] <- NA
mydata$V2 <- replace(mydata$V2, is.na(mydata$V2), min(na.omit(mydata$V2)))

mydata[mydata$V4==0] <- NA
mydata$V4 <- replace(mydata$V4, is.na(mydata$V4), min(na.omit(mydata$V4)))



result<-data.frame()
for(linea in 1:nrow(mydata)){
  
  miRNA<-mydata[linea,1]
  gene<-mydata[linea,3]
  RPKM_m<-as.numeric(unlist(strsplit(mydata[linea,2], ",")))
  RPKM_m[RPKM_m==0]<-min_miRNA
  
  RPKM_g<-as.numeric(unlist(strsplit(mydata[linea,4], ",")))
  RPKM_g[RPKM_g==0]<-min_gene
  
  mycor<-cor.test(log(RPKM_m),log(RPKM_g),method = mymethod )

  result[linea,1]<-miRNA
  result[linea,2]<-gene
  result[linea,3]<-mycor$estimate[[1]]
  result[linea,4]<-mycor$p.value
}

colnames(result)<-c("miRNA","GeneName","R","P-value")


result<-result[order(result$`P-value`),]
write.table(result,row.names = F,file=myoutput,sep="\t")

#result$FDR<-p.adjust(result$`P-value`,n=nrow(result))
#result<-result[order(result$FDR),]
#write.table(result[result$`P-value`<=0.05,],row.names = F,file=myoutput,sep="\t")
