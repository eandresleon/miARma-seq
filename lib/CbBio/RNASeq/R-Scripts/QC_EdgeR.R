#########################################################################  
#	QC_EdgeR R function			                          	 									#
#																		                                    #
#	Created at Computational Biology and Bioinformatics Group (CbBio)	    #
#	Institute of Biomedicine of Seville. IBIS (Spain)					            #
#	Copyright (c) 2016 IBIS. All rights reserved.						              #
#	mail : miarmaseq-devel@cbbio.es 								                      #
#########################################################################

#Printing the function information and  help message
write("
      
      #########################################################################  
      QC_EdgeR R function			                          	 									
      
      Created at Computational Biology and Bioinformatics Group (CbBio)	    
      Institute of Biomedicine of Seville. IBIS (Spain)					            
      Copyright (c) 2016 IBIS. All rights reserved.						              
      mail : miarmaseq-devel@cbbio.es 								                      
      #########################################################################
      
      
      Get help typing QC_EdgeR_help()\n",stderr())

#Defining help function
QC_EdgeR_help<-function(){
  write("
        QC_EdgeR is a R function which takes a tab file with the number of the reads from htseq-count analysis and a target file 
        with the experimental conditions of the samples to analyze the quality of the samples. For this purpose QC_EdgeR genrates a pdf 
        report with differents plots: boxplots and density plots with the distribution of the reads in each sample with 
        and without normalization, MDS and PCA plots and heatmaps of the samples with the more expressed elements(genes, miRNAs...)
        and between the samples. 
        
        QC_EdgeR takes 10  arguments:
        
        Mandatory arguments:
        [projectdir] Path of the directory where will be saved the QC report.
        [dir] the path of the directory which contains the files. This directory will be configured as working directory for R
        [file] name of the tab file which contains the number of reads from the htseq analysis
        [targetfile] Path of the tabulated file which contains the experimental condiction of each sample. 
        First column must contain the names to be used to the plots, and the next columns the condition of each factor. 
        QC_EdgeR works with 1 or 2 factors.
        [label] Name of the experiment to appear in the title of the plots and in the name of the pdf file results
        [filter] This value refers to filter processing in the reads (Should be \"yes\" or \"no\").
        
        Optional arguments:
        [cpmvalue] Cutoff for the counts per million value to be used in filter processing (1 cpm by default).
        [repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)", stderr())
}

QC_EdgeR<-function(projectdir,dir,file,targetfile,label,filter, cpmvalue=1, repthreshold=2,normethod="TMM", bcvvalue=0.4){
  
  #Checking the mandatory parameters
  if(missing(projectdir) | missing(dir) | missing(file) | missing(targetfile) | missing(label) | missing(filter)){
    QC_EdgeR_help()
    stop("QC_EdgeR Error: Missing mandatory argument")
  }
  
  #########################################################################
  #1- IMPORTING AND FORMATING THE DATA
  #########################################################################
  
  #installing need packages
  list.of.packages <- c("ggplot2", "Rcpp","edgeR","NOISeq","gplots","lattice","genefilter","RColorBrewer")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)){
                
  	if(R.version$minor>5){
  		if (!requireNamespace("BiocManager"))
  		install.packages("BiocManager")
  		BiocManager::install(new.packages,
  		update = F, 
  		ask = F)
  	}
  	else{
  		source("http://bioconductor.org/biocLite.R")
  		biocLite(new.packages,
  	                   suppressUpdates=T,
  	                   suppressAutoUpdate=T,
  					   ask=F)
  	}
  }
  
  
  #Loading the needed packagge
  require(edgeR)
  require(ggplot2)
  require(gplots)
  require(genefilter)
  
  check_parameter<-function(data){
    cont<-0
    for(a in 1:length(data)){
      if(data[[a]]>0){
        cont=cont+1
      }
    }
    return(cont)
  }
  
  #Printing the date and information of the proccess
  time<-Sys.time()
  message<-paste(time, " Starting quality control analysis of ", label,"\n", sep="")
  cat(message, sep="")
  
  #Importing the data
  workingDir<-dir
  setwd(workingDir)
  filepath<-file.path(workingDir,file)
  #added check.names=F suggested by Guillaume Noell
  data <- read.table(filepath, header=TRUE, sep="\t",check.names=F)
  
  data<-data[,sort(colnames(data))]
  
  #########################################################################
  #2-EDGER ANALYSIS-BOXPLOT, DENSITY PLOT AND CLUSTERING
  #########################################################################
  
  #Importing targets
  targets<-readTargets(targetfile)
  targets<-targets[order(targets$Filename),]
  
  options(warn=-1)
  #Checking the dimensions of the targets
  #For 1 dimension target file must contain 2 columns: names and condition and for 2 dimensions must contain 3: names, condition 1 and condition 2
  if(length(colnames(targets))==3){
    group<-factor(targets[,3])
  }else if(length(colnames(targets))==4){
    group<-factor(paste(targets[,3],targets[,2],sep="_"))
  }else{
    stop("Target file has an incorrect format. Please see documentation to get information of the correct format of this file")
  }
  #Merging the conditions 
  cbind(targets,group=group)
  
  #Obtaining the number of replicates
  df<- cbind(targets,group=group)
  levelsnames<-levels(group)
  repnumber<-NA
  for(a in 1:length(levels(group))){
    repnumber[a] <- length(rownames(subset(df, group ==levelsnames[a])))
  }  
  #Generating the comparison matrix
  design<-model.matrix(~0+group)
  colnames(design)<-levels(group)
  
  #DGEList element creation
  dge<-DGEList(counts=data,group=group)
  
  
  #Filtering the counts
  if(filter=="yes"){
    keep<-rowSums(cpm(dge)>cpmvalue) >= repthreshold
    dge<-dge[keep, ]
    #Recalculating the library size
    dge$samples$lib.size <-colSums(dge$counts)
  }
  
  #Opening pdf file
  filename<-paste(label, "QC_Report.pdf", sep="_")
  filepath<-file.path(projectdir,filename)
  #Saving the path of the result file to return it.
  resultspaths<-NA
  resultspaths<- filepath
  pdf(file=filepath, paper="a4")
  #Setting the plot features
  par(mfrow=c(1,2), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.8, cex.main=0.8)
  #Setting the colors plot for each replicate
  mycolors<-colorRampPalette(c("steelblue","indianred1","yellowgreen"))(length(levels(group)))
  boxcol<-NA
  counta<-1
  for (a in 1:length(levels(group))){
    for(b in 1:repnumber[a]){
      boxcol[counta]<-mycolors[a]
      counta<-counta+1
    }
  }
  #Obtaining the names of the samples
  samplenames<-as.character(targets[,2])
  #Boxplot of the samples
  boxplot(log2(data), main=paste("Boxplot of ",label," samples",sep=""),las=2, names=samplenames, ylab="log2(counts)", xlab="")
  #Normalization of the samples
  dgenorm<-calcNormFactors(dge,method=normethod)
  #Boxplot of the normalized samples
  boxplot(log2(dgenorm$counts), col=boxcol, main=paste("Boxplot of ", label , " normalized samples",sep=""),las=2, names=samplenames, ylab="log2(counts)", xlab="")
  
  
  #Density plots of the number of reads 
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.8)
  plot(density(log10(data[,1])),col=boxcol[1], main=paste("Density plot of number of reads of ", label, sep=""), lwd=1.5, xlab="log10(counts)")
  lineformat<-NA
  countb<-1
  for (a in 1:length(levels(group))){
    for(b in 1:repnumber[a]){
      lineformat[countb]<-b
      countb<-countb+1
    }
  }
  number<-0
  for (a in 2:length(colnames(data))){
    number<-check_parameter(data[,a])
    if(number>1){
      lines(density(log10(data[,a])),col=boxcol[a], lwd=1.5, lty=lineformat[a])
    }
  }
  legend("topright", samplenames, col=boxcol, lty=lineformat, lwd=4, cex=0.8)
  
  #Density plot of normalized counts data
  plot(density(log10(dgenorm$counts[,1])),col=boxcol[1], main=paste("Density plot of normalized and filtered number of reads of ", label, sep=""), lwd=1.5, xlab="log10(counts)")
  for (a in 2:length(colnames(data))){
    number<-check_parameter(dgenorm$counts[,a])
    if(number>1){
      lines(density(log10(dgenorm$counts[,a])),col=boxcol[a], lwd=1.5, lty=lineformat[a])
    }
  }
  legend("topright", samplenames, col=boxcol, lty=lineformat, lwd=4, cex=0.8)
  
  if(length(targets$Name)>2){
    #Clustering analysis plot
    par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
    pr.hc.c<- hclust(na.omit(dist(t(data))))
    plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of ", label, sep=""), labels=samplenames,col=boxcol)
    #Normalized clustering analysis plot 
    pr.hc.c<- hclust(na.omit(dist(t(dgenorm$counts))))
    plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of Normalized samples of ", label, sep=""), labels=samplenames,col=boxcol)
    
    #########################################################################
    #3. MDS PLOT (EDGER) AND PCA ANALYSIS
    #########################################################################
    
    #Exploring the data by MDS plot
    labels<-colnames(targets)
    par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.8, cex.main=0.8)
    if(length(colnames(targets))==2){
      plotMDS.default(dgenorm, labels=samplenames, col=boxcol, main=paste("MDS plot of ",label, " samples",sep=""), cex=0.75)
    }else if(length(colnames(targets))==3){
      plotMDS.default(dgenorm, labels=samplenames, col=boxcol, main=paste("MDS plot of ",label, " samples",sep=""), xlab=paste("logFC ",labels[2],sep=""), ylab=paste("logFC ",labels[3], sep=""), cex=0.75)
    }
  }
  
  #3.1 PCA plot
  data.PC = prcomp(t(as.matrix(dgenorm)))
  plot.default(data.PC$x,col=boxcol,main=paste("PCA plot of ",label, " samples",sep=""))
  text(data.PC$x[,1],data.PC$x[,2],labels=samplenames, cex=0.7, pos=1)
  
  #########################################################################
  #4-HEATMAP PLOT
  #########################################################################
  
  rsd <- rowSds(as.matrix(dgenorm))
  sel <- order(rsd, decreasing=TRUE)[1:250]
  #hmcol= colorRampPalette(brewer.pal(9, "YlGnBu"))(100) 
  #heatmap(na.omit(as.matrix(dgenorm[sel,])),margins=c(10,8),main="Heatmap 250 most DE entities",cexRow=0.5,cexCol=0.5,labCol=samplenames, col=hmcol)
  heatmap(na.omit(as.matrix(dgenorm[sel,])),margins=c(10,8),main="Heatmap 250 most DE entities",cexRow=0.01,cexCol=0.5,labCol=samplenames)  
  dev.off()
  #Printing the date and information of the proccess
  time<-Sys.time()
  message<-paste(time, " The file ", filename," has been generated in ",projectdir,"\n", sep="")
  cat(message, sep="")
  return(resultspaths)
}
