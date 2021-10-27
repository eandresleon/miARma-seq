#########################################################################  
#  DE_EdgeR R function  		                           									#
#																		                                    #
#	Created at Computational Biology and Bioinformatics Group (CbBio)	    #
#	Institute of Biomedicine of Seville. IBIS (Spain)					            #
#	Copyright (c) 2016 IBIS. All rights reserved.						              #
#	mail : miarmaseq-devel@cbbio.es 								                      #
#########################################################################

#Printing the function information and  help message
write("

#########################################################################  
DE_EdgeR R function    	                           									
																		                                    
Created at Computational Biology and Bioinformatics Group (CbBio)	    
Institute of Biomedicine of Seville. IBIS (Spain)					            
Copyright (c) 2016 IBIS. All rights reserved.						              
mail : miarmaseq-devel@cbbio.es 								                      
#########################################################################

Get help typing DE_EdgeR_help()\n",stderr())

#Defining help function
DE_EdgeR_help<-function(){
  write("
DE_EdgeR function calls a R function called DE_EdgeR which takes a tab file with the number of the reads from htseq-count analysis, 
a target file with the experimental conditions of the samples and the contrast file with the contrast to analyze the differential 
expression between the defined samples. This function outputs a pdf file with the descriptive plots of the analysis and a tab file 
(xls extension to easy exploration) for condition evaluated. This function accepts one or two factors experiments.
   
DE_EdgeR takes 12  arguments:

Mandatory arguments:
  [projectdir] Path of the directory where will be saved the QC report.
  [dir] the path of the directory which contains the files. This directory will be configured as working directory for R
  [file] Path of the tab file which contains the number of reads from the htseq analysis
  [targetfile] Path of the tabulated file which contains the experimental condiction of each sample. 
  First column must contain the names to be used to the plots, and the next columns the condition of each factor. 
  QC_EdgeR works with 1 or 2 factors.
  [label] Name of the experiment to appear in the title of the plots and in the name of the pdf file results.
  [contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. 
  The syntax of the contrast will be: name_of_contrast=contrast to evaluate.Any type of contrast can be done but 
  condition name must be one of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
  [filter] This value refers to filter processing in the reads (Should be \"yes\" or \"no\")
  
Optional arguments:
  [cpmvalue] Cutoff for the counts per million value to be used in filter processing (1 cpm by default).
  [repthreshold] Number of replicates that have to contains at least a defined number of reads per million to perform the filter (2 replicates by default)
  [normethod] Normalization method. It can be one of \"TMM\" (default), \"RLE\", \"upperquartile\" or \"none\" (no normalization).
  [replicates] Value to indicate if replicates samples are present in the analysis. It can be \"yes\" (by default) or \"no\".
  [bcvvalue] Value for the common BCV (square- root-dispersion) in experiments without replicates. Standard values from well-controlled 
  experiments are 0.4 for human data (by default), 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
  [rpkm] Value to indicate if a file with RPKM values should be created. It can be \"yes\" (by default) or \"no\".
  [cpm] Value to indicate if a file with CPM values should be created. It can be \"yes\" (by default) or \"no\".

  Example:
    
    DE_EdgeR(projectdir=\"/Project\",dir=\"/Project/htseq-results\",file=\"htseqresults.tab\",targetfile=\"/Project/target.txt\",label=\"Hypoxia\",contrastfile=\"/Project/contrast.txt\", filter=\"yes\", repthreshold=2, cpmvalue=5, normethod=\"TMM\", replicates=\"yes\", bcvvalue=0.4)", stderr())
  
}


DE_EdgeR<-function(projectdir,dir,file,targetfile,label,contrastfile, filter, cpmvalue=1, repthreshold=2, normethod="TMM", replicates="yes", bcvvalue=0.4,rpkm=T,cpm=F,file_size=F){
  
  #Checking the mandatory parameters
  if(missing(projectdir) | missing(dir) | missing(file) | missing(targetfile) | missing(label) | missing(filter) | missing(contrastfile)){
    DE_EdgeR_help()
    stop("DE_EdgeR Error: Missing mandatory argument")
  }
  
  
  #########################################################################
  #1- IMPORTING AND FORMATING THE DATA
  #########################################################################
  
  #installing need packages
	list.of.packages <- c("ggplot2", "Rcpp","edgeR","NOISeq","gplots","lattice","genefilter","RColorBrewer","ggrepel")
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
  library("genefilter")
  library("ggplot2")
  library("ggrepel")
  #Printing the date and information of the proccess
  time<-Sys.time()
  message<-paste(time, " Starting Differential expression analysis with EdgeR of ", label,"\n", sep="")
  cat(message, sep="")
  
  #Importing the data
  workingDir<-projectdir
  setwd(workingDir)
  
  #added check.names=F suggested by Guillaume Noell
  data <- read.table(file, header=TRUE, sep="\t",check.names=F)
  data<-data[,sort(colnames(data))]
  

  #Importing targets
  targets<-readTargets(targetfile,row.names="Filename")
  targets<-targets[order(targets$Filename),]

  #########################################################################
  #2- FILTERING AND NORMALIZATION
  #########################################################################
   
  #Checking the dimensions of the targets
  #For 1 dimension target file must contain 2 columns: names and condition and for 2 dimensions must contain 3: names, condition 1 and condition 2
  if(length(colnames(targets))==3){
    group<-factor(targets[,3])
  }else if(length(colnames(targets))==4){
    group<-factor(paste(targets[,3],targets[,4],sep="_"))
  }else{
    stop("Target file has an incorrect format. Please see documentation to get information of the correct format of this file")
  }
  #Merging the conditions
  cbind(targets,group=group)
  
  #DGEList element creation
  dge<-NULL
  if(rpkm == "yes"){
    gene.length<-read.table(file_size,header=T)
    idx<-match(rownames(data),gene.length$Gene)
    results_counts<-gene.length[idx,]
    results_counts[is.na(results_counts$Length),"Length"]<-0
    
    dge<-DGEList(counts=data,genes=results_counts,group = group)
    
  }else{
   dge<-DGEList(counts=data,group=group)
  }
  
  #Filtering the counts
  if(filter=="yes"){
    keep<-rowSums(cpm(dge)>cpmvalue) >= repthreshold
    dge<-dge[keep, ]
    #Recalculating the library size
    dge$samples$lib.size <-colSums(dge$counts)
  }
  
  #Normalization of the samples
  dgenorm<-calcNormFactors(dge, method=normethod)
  if(cpm=="yes"){
    write(paste("[",Sys.time(),"]"," Calculating normalized CPMs",sep=""), stderr())
    mycpm<-cpm(dgenorm,normalized.lib.sizes=T,log=F)
    results_cpm<-file.path(projectdir,paste(label,"_EdgeR_normalized_reads.tsv", sep=""))
    write.table(mycpm,file=results_cpm,sep="\t",col.names = NA,row.names = T)
  }
  if(rpkm=="yes"){
    write(paste("[",Sys.time(),"]"," Calculating RPKMs",sep=""), stderr())
    y_rpkm<-rpkm(dgenorm,dge$genes$Length)
    write(paste("[",Sys.time(),"]"," RPKMs done!",sep=""), stderr())
    
    results_rpkm<-file.path(projectdir,paste(label,"_EdgeR_RPKM.tsv", sep=""))
    write.table(y_rpkm,results_rpkm,col.names=NA,row.names=T,sep="\t")
  }
  #########################################################################
  #3. DISPERSION ESTIMATION AND DIFFERENTIAL EXPRESSION (DE) ANALYSIS
  ######################################################################### 
  
  #Making contrasts
  #Reading contrast file
  contrastconds<-read.delim(file=contrastfile,col.names=NA)
  
  #Generating plot file
  plotpaths<-NA
  filename<-paste(label, "EdgeR_plots.pdf", sep="_")
  filepath<-file.path(projectdir,filename)
  plotpaths[1]<- filepath
  pdf(file=filepath, paper="a4")
  filepaths<-NA
  #Analysis with replicates
  if(tolower(replicates)=="yes"){
    #Dispersion estimation and DE analysis differs according to to the number of factors. 
    filepaths<-NA
    if(length(colnames(targets))==3){
      #For 1 factor qCML method is recommended to estimate the dispersion and exact test to perform DE analysis
      dgenorm <- estimateCommonDisp(dgenorm)
      dgenorm <- estimateTagwiseDisp(dgenorm)
      
      #Plotting the tagwise dispersions against log2-CPM
      par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.8, cex.main=1.2)
      plotBCV(dgenorm, col.tagwise="forestgreen",pch=19,col.common="indianred1", cex=0.3, main="Biological Variation Plot")
      # Mean-Variance Plot: Four things are shown in the plot: the raw variances of the counts (grey dots), the variances using the tagwise dispersions (light blue dots), 
      #the variances using the common dispersion (solid blue line), and the variance = mean a.k.a. poisson variance (solid black line).  
      meanVarPlot <- suppressWarnings(plotMeanVar( dgenorm ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ))
      
      #DE analysis for each condition
      for (a in 1:length(rownames(contrastconds))){
        #Obtaining the name of the contrast and the conditions to make the contrast
        contrastsname<-strsplit(as.character(contrastconds[a,1]), "=")[[1]]
        contrastval<-strsplit(as.character(contrastsname[2]), "-")[[1]]
        #Performing the DE analysis for each pair
        et <- exactTest(dgenorm,pair=c(contrastval[2],contrastval[1]))
        #Extracting the statistical data order by p-value
        top<-topTags(et, n=nrow(et))
        if(rpkm=="yes"){
          rownames(top$table)<-top$table$Gene
          top$table<-top$table[,grep("Gene",colnames(top$table),invert = T)]
        }
        #Writting the data in a tab file with xls extension
        resultsfile<-file.path(projectdir,paste(label,"_EdgeR_results_",contrastsname[1],".xls", sep=""))
        write.table(top, file=resultsfile, sep = "\t", col.names = NA , qmethod = "double")
        #Plotting the tagwise log-fold-changes against log-cpm
        de <- decideTestsDGE(et, p=0.05, adjust="BH")
        detags <- rownames(dgenorm)[as.logical(de)]
        plotSmear(et, de.tags=detags, col="forestgreen", main=paste("Expression plot of ", contrastsname[1], sep=""), cex=0.3)
        abline(h = c(-2, 2), col = "indianred")
        #Saving the path of the results file to return it to the function
        filepaths[a]<-resultsfile
        
		    #volcano plot
        data<-top$table
        data$threshold = data$FDR < 0.05
        
        data[data$FDR > 0.05,"threshold"]<-"NoDE"
        data[data$FDR <= 0.05,"threshold"]<-"DE"
        
        de_data<-data[data$FDR<=0.05,]
        
        selected_FC<-rbind(head(de_data[order(de_data$logFC,decreasing = T),], 10),head(de_data[order(de_data$logFC,decreasing = F),],10))
        g = ggplot(data=data, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
          geom_point(alpha=0.4, size=1.75) +
          xlab("log2 fold change") + ylab("-log10 FDR") +
          geom_text_repel(data=selected_FC, aes(label=rownames(selected_FC)),colour="black",size=3) #adding text for the top 20 genes
        
        plot(g)

        #bar plots
        
        plot_barplot<-function(mydata,rpkm){
          genes<-rownames(mydata)
          FC<-mydata$logFC
          for(gene in genes){
            FC<-round(mydata[rownames(mydata)==gene,]$logFC,2)
            boxplot(split(as.numeric(rpkm[gene,]),contrastval),col=c("salmon","lightblue"),main=paste("Relative ", gene, " expression (logFC:",FC,")",sep=""), ylab="log(RPKM)")
          }
        }
        
        y_rpkm[y_rpkm==0] <- NA
        y_rpkm <- replace(y_rpkm, is.na(y_rpkm), min(na.omit(y_rpkm)))
        
        plot_barplot(unique(selected_FC),log(y_rpkm))
        
        #Printing the date and information of the proccess
        time<-Sys.time()
        message<-paste(time, " The file ", resultsfile ," has been generated in ", projectdir," for ",contrastsname[1]," comparison\n", sep="")
        cat(message, sep="")
      }
    }else if(length(colnames(targets))==4){
      #Generating the comparison matrix
      design<-model.matrix(~0+group)
      colnames(design)<-levels(group)
      
      #For 2 factors GLM method is used to estimate dispersion.
      dgenorm<-estimateGLMCommonDisp(dgenorm,design)
      dgenorm<-estimateGLMTrendedDisp(dgenorm,design)
      dgenorm<-estimateGLMTagwiseDisp(dgenorm,design)
      
      #Plotting the tagwise dispersions against log2-CPM
      par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.8, cex.main=1.2)
      plotBCV(dgenorm, col.tagwise="forestgreen",pch=19,col.common="indianred1", cex=0.3, main="Biological Variation Plot")
      # Mean-Variance Plot: Four things are shown in the plot: the raw variances of the counts (grey dots), the variances using the tagwise dispersions (light blue dots), 
      #the variances using the common dispersion (solid blue line), and the variance = mean a.k.a. poisson variance (solid black line).  
      meanVarPlot <- suppressWarnings(plotMeanVar( dgenorm ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ))
      
      # Fitting to the linear model
      fit<-glmFit(dgenorm,design)
      
      #Obtaining the conditions
      contrastvector<-as.vector(contrastconds$NA.)
      
      #Contrast DEG calculation for each experimental condition
      filepaths<-NA
      for (a in 1:(length(contrastvector))){
        #Obtaining the name of the contrast
        contrastsname<-strsplit(as.character(contrastconds[a,1]), "=")[[1]]
        #Generating the contrast
        my.contrasts<-makeContrasts(contrast=contrastvector[a],levels=design)
        lrt<-glmLRT(fit, contrast=my.contrasts)
        #Estadistical data extraction
        top<-topTags(lrt, n=nrow(lrt))
        if(rpkm=="yes"){
          rownames(top$table)<-top$table$Gene
          top$table<-top$table[,grep("Gene",colnames(top$table),invert = T)]
        }
        #Saving the data in a xls file
        resultsfile<-file.path(projectdir,paste(label,"_EdgeR_results_",contrastsname[1],".xls", sep=""))
        write.table(top, file=resultsfile, sep = "\t", col.names = NA, qmethod = "double")
        filepaths[a]<- resultsfile
        #Plotting the tagwise log-fold-changes against log-cpm
        de <- decideTestsDGE(lrt, p=0.05, adjust="BH")
        detags <- rownames(dgenorm)[as.logical(de)]
        plotSmear(lrt, de.tags=detags, col="forestgreen", main=paste("Expression plot of ", contrastsname[1], sep=""), cex=0.3)
        abline(h = c(-2, 2), col = "indianred")
        
        #volcano plot
		
        de_data<-top$table
        de_data$threshold = de_data$FDR < 0.05
        
        de_data[de_data$FDR > 0.05,"threshold"]<-"NoDE"
        de_data[de_data$FDR <= 0.05,"threshold"]<-"DE"
        
        selected_FC<-rbind(head(de_data[order(de_data$logFC,decreasing = T) & de_data$FDR<=0.05,], 5),head(de_data[order(de_data$logFC,decreasing = F) & de_data$FDR<=0.05,],5))
        g = ggplot(data=de_data, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
          geom_point(alpha=0.4, size=1.75) +
          xlab("log2 fold change") + ylab("-log10 FDR") +
          geom_text_repel(data=selected_FC, aes(label=rownames(selected_FC)),colour="black",size=3) #adding text for the top 20 genes
        
        plot(g)
        
        #Printing the date and information of the proccess
        time<-Sys.time()
        message<-paste(time, " The file ", resultsfile ," has been generated in ", projectdir," for ",contrastsname[1]," comparison\n", sep="")
        cat(message, sep="")
      }
    }else{
      stop("Target file has an incorrect format. Please see documentation to get information of the correct format of this file")
    }
    dev.off()
  }else if(tolower(replicates)=="no"){
    #If there is no replicates:
      bcv <- bcvvalue
      #DE analysis for each condition
      for (a in 1:length(rownames(contrastconds))){
        #Obtaining the name of the contrast and the conditions to make the contrast
        contrastsname<-strsplit(as.character(contrastconds[a,1]), "=")[[1]]
        contrastval<-strsplit(as.character(contrastsname[2]), "-")[[1]]
        #Performing the DE analysis for each pair
        et <- exactTest(dgenorm,pair=c(contrastval[2],contrastval[1]), dispersion=bcv^2)
        #Extracting the statistical data order by p-value
        top<-topTags(et, n=nrow(et))
        #Writting the data in a tab file with xls extension
        resultsfile<-file.path(projectdir,paste(label,"_EdgeR_results_",contrastsname[1],".xls", sep=""))
        write.table(top, file=resultsfile, sep = "\t", col.names = NA , row.names = TRUE, qmethod = "double")
        #Plotting the tagwise log-fold-changes against log-cpm
        de <- decideTestsDGE(et, p=0.05, adjust="BH")
        detags <- rownames(dgenorm)[as.logical(de)]
        plotSmear(et, de.tags=detags, col="forestgreen", main=paste("Expression plot of ", contrastsname[1], sep=""), cex=0.3)
        abline(h = c(-2, 2), col = "indianred")
        #Saving the path of the results file to return it to the function
        filepaths[a]<-resultsfile
        
        #volcano plot
        de_data<-top$table
        de_data$threshold = de_data$FDR < 0.05
        
        de_data[de_data$FDR > 0.05,"threshold"]<-"NoDE"
        de_data[de_data$FDR <= 0.05,"threshold"]<-"DE"
        
        selected_FC<-rbind(head(de_data[order(de_data$logFC,decreasing = T) & de_data$FDR<=0.05,], 5),head(de_data[order(de_data$logFC,decreasing = F) & de_data$FDR<=0.05,],5))
        g = ggplot(data=de_data, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
          geom_point(alpha=0.4, size=1.75) +
          xlab("log2 fold change") + ylab("-log10 FDR") +
          geom_text_repel(data=selected_FC, aes(label=rownames(selected_FC)),colour="black",size=3) #adding text for the top 20 genes
        
        plot(g)
        #Printing the date and information of the proccess
        time<-Sys.time()
        message<-paste(time, " The file ", resultsfile ," has been generated in ", projectdir," for ",contrastsname[1]," comparison\n", sep="")
        cat(message, sep="")
      }
  }else{
    stop("Replicates has an incorrect value. This argument only accepts yes or no")
  }
  resultspathsdef<-c(plotpaths,filepaths)
  return(resultspathsdef)
}
