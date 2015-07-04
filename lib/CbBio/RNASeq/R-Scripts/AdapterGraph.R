AdapterGraph<-function(dir,file, name){
  #Setting the directory
  workingDir<-dir
  setwd(workingDir)
  #Reading the stats file from Reaper
  data<-read.table(file, header=TRUE, sep="\t")
  #Obtaining the size of the reads and the correspondant number of reads
  readsizes<-data[,1]
  readvalues<-data[,2]
  #Saving the plot in a pdf file in the provided directory
  pdf(file=paste("Reads_histogram_of_", name,".pdf",sep=""), paper="a4")
  #Setting graph parameters
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.8, cex.main=1, las=1, mar=c(5,5,4,1), mgp=c(3,1,0), oma=c(2, 2, 1, 1))
  #Plotting the bar graph
  barplot(readvalues, main=paste("Reads size distribution of ", name, " after adapter removal", sep=""), xlab="Read size", ylab="", names.arg=readsizes, ylim=c(0,(max(readvalues)+(max(readvalues)/5))), col="yellowgreen", beside=TRUE, space=0)
  mtext("Number of Reads", side=2, line=4, las=0, cex=0.8)
  #Closing the file
  dev.off()
}