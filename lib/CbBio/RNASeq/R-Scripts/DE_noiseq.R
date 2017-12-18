#########################################################################  
#  DE_noiseq R function			                           									#
#																		                                    #
#	Created at Computational Biology and Bioinformatics Group (CbBio)	    #
#	Institute of Biomedicine of Seville. IBIS (Spain)					            #
#	Copyright (c) 2016 IBIS. All rights reserved.						              #
#	mail : miarmaseq-devel@cbbio.es 								                      #
#########################################################################


#Printing the function information and  help message
write("

######################################################################### 
  DE_noiseq R function  		                           									
																		                                    
	Created at Computational Biology and Bioinformatics Group (CbBio)	    
	Institute of Biomedicine of Seville. IBIS (Spain)					            
	Copyright (c) 2016 IBIS. All rights reserved.						              
	mail : miarmaseq-devel@cbbio.es 								                      
#########################################################################
      
      
Get help typing DE_noiseq_help()\n",stderr())

#Defining help function
DE_noiseq_help<-function(){
  write("
DE_noiseq function takes a tab file with the number of the reads from htseq-count analysis, a target file with the experimental conditions 
of the samples and the contrast file with the contrast to analyze the differential expression between the defined samples. 
This function outputs a pdf file with the descriptive plots of the analysis and two tab files (xls extension to easy exploration) 
for condition evaluated. The Noiseq_results tab file contains the whole results data and the Noiseq_DE tab file contains only the data of the 
differentially expressed elements.
        
DE_noiseq takes 22 arguments:
        
Mandatory arguments:
  [projectdir] Path of the directory where the output files will be saved.
  [dir] Path of the directory which contains the file with the reads. This directory will be configured as working directory for R
  [file] Name of the tab file which contains the number of reads from the htseq analysis
  [targetfile] Path of the target file. This is a tabulated file which contains the experimental condiction of each sample. 
  The names of the samples defined in this file will be used to the plots.
  [label] Name of the experiment to appear in the title of the plots and in the name of the results files
  [filter] This value refers to filter processing in the reads (Should be \"yes\" or \"no\")
  [contrastfile] Path of the contrast file.This file has one column with the contrasts user want to evaluate. 
  The syntax of the contrast will be: name_of_contrast=condition1-condition2.Any type of contrast can be done but 
  condition name must be one of the conditions present in targets file. Also contrast must differ of 0 (ie: cond=WT-WT is not allowed)
        
Optional arguments
  [filtermethod] Method that will be used to filter proccess. Method must be one of 1,2 or 3 for CPM method, Wilcoxon method and Proportion test, respectively (by default CPM method)
  [cpmvalue] Cutoff for the counts per million value to be used in filter processing with methods 1 and 3 (1 cpm by default).
  [cutoffvalue] Cutoff for the coefficient of variation per condition to be used in filter processing with method 1 (in percentage, 100 by default).
  [lengthfile] Path of the file with the information about the length of each element. This file must contain the name of the element and the legnth of this element under the name.
  [gcfile] Path of the file with the information about the gc content of each element. This file must contain the name of the element and the gc content (%) of this element under the name.
  [biotypefile] Path of the file with the information about the biotype of each element. This file must contain the name of the element and the biotype of this element under the name.
  [chromsfile] Path of the file with the genomic location information of each element. This file must contain 4 columns corresponding to: the name of the element, the number of chromosome, and the coordinates of the start and end of the element.
  [normethod] Normalization method. It can be one of \"rpkm\" (default), \"uqua\" (upper quartile), \"tmm\" (trimmed mean of M) or \"n\" (no normalization).
  [replicatevalue] type of replicates to be used. It can be Technical, biological or none. By default, technical replicates option is chosen.
  [kvalue] Counts equal to 0 are replaced by k. By default, k = 0.5.
  [lcvalue] Length correction is done by dividing expression by length^lc. By default, lc = 0.
  [pnrvalue] Percentage of the total reads used to simulated each sample when no replicates are available. By default, pnr = 0.2.
  [nssvalue] Number of samples to simulate for each condition (nss>= 2). By default, nss = 5. 
  [vvalue] Variability in the simulated sample total reads. By default, v = 0.02. Sample total reads is computed as a random value from a uniform distribution in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]
  [qvalue] Probability of differential expression. By default=0.8
        
DE_noiseq(projectdir=\"/Project\",dir=\"/Project/htseq-results\",file=\"htseqresults.tab\",targetsfile=\"/Project/target.txt\",label=\"Hypoxia\",filter=\"yes\", contrastfile=\"/Project/contrast.txt\", filtermethod=3, cpmvalue=5, normethod=\"rpkm\", qvalue=0.9)", stderr())
}

DE_noiseq<-function(projectdir,dir,file,targetsfile,label,filter,contrastfile,lenghtfile=NULL, gcfile=NULL, biotypefile=NULL, chromsfile=NULL, filtermethod=1,cpmvalue=1,cutoffvalue=100, normethod="rpkm",replicatevalue="technical",kvalue=0.5,lcvalue=0,pnrvalue=0.2,nssvalue=5,vvalue=0.02,qvalue=0.8,rpkm=F,cpm=F,file_size=F){
  
  #Checking the mandatory parameters
  if(missing(projectdir) | missing(dir) |  missing(file) | missing(targetsfile) | missing(label) | missing(filter) | missing(contrastfile)){
    DE_noiseq_help()
    stop("DE_noiseq Error: Missing mandatory argument")
  }
  
  #Loading the required packages
  require("NOISeq")
  
  #Printing the date and information of the proccess
  time<-Sys.time()
  message<-paste(time, "Starting Differential expression analysis with Noiseq of ", label,"\n", sep="")
  cat(message, sep="")
  
  #Importing the data
  workingDir<-dir
  setwd(workingDir)
  #added check.names=F suggested by Guillaume Noell
  data<-read.table(file=file, sep="\t", header=T,check.names=F)
  #Checking any additional data provided to use in data exploration
  if(!missing(lenghtfile)){
    mylength<-read.table(file=lenghtfile, sep="\t", header=T)
  }else{
    mylength<-NULL
  }
  if(!missing(gcfile)){
    mygc<-read.table(file=gcfile, sep="\t", header=T)
  }else{
    mygc<-NULL
  }
  if(!missing(biotypefile)){
    mybiotype<-read.table(file=biotypefile, sep="\t", header=T)
  }else{
    mybiotype<-NULL
  }
  if(!missing(chromsfile)){
    mychroms<-read.table(file=chromsfile, sep="\t", header=T)
  }else{
    mychroms<-NULL
  }
  
  if(cpm=="yes"){
    write(paste("[",Sys.time(),"]"," Calculating normalized CPMs",sep=""), stderr())
    mycpm<-tmm(data,long=1000,lc=0,k=0)
    results_cpm<-file.path(projectdir,paste(label,"_NOISeq_normalized_reads.tsv", sep=""))
    write.table(mycpm,file=results_cpm,sep="\t",col.names = NA,row.names = T)
  }
  if(rpkm=="yes"){
    gene.length<-read.table(file_size,header=T)
    idx<-match(rownames(data),gene.length$Gene)
    results_counts<-gene.length[idx,]
    results_counts[is.na(results_counts$Length),"Length"]<-0
    
    write(paste("[",Sys.time(),"]"," Calculating RPKMs",sep=""), stderr())
    mycpm<-rpkm(data,long=results_counts$Length,lc=1,k=0)
    
    results_cpm<-file.path(projectdir,paste(label,"_NOISeq_RPKM.tsv", sep=""))
    write.table(mycpm,file=results_cpm,sep="\t",col.names = NA,row.names = T)
  }
  #Generating factors for experimental conditions
  target<-read.table(file=targetsfile, sep="\t", header=T,row.names=1)
  #Generating factors using targets. Only one factor is allowed
  myfactors <- data.frame(Factor1=target[,2])
  samplenames<-target[,1]
  
  #Filtering the counts
  if(filter=="yes"){
    #Printing the date and information of the proccess
    time<-Sys.time() 
    message<-paste(time, " Filtering the reads of ", label," with Noiseq filter method ", filtermethod, "\n", sep="")
    cat(message, sep="")
    #Checking the method selected to filter the data
    if (filtermethod=="3"){
      myfilt = filtered.data(data, factor = myfactors$Factor1, norm = FALSE, method = 3, cpm = cpmvalue)
    }else if(filtermethod=="2"){
      myfilt = filtered.data(data, factor = myfactors$Factor1, norm = FALSE, method = 2)
    }else if(filtermethod=="1"){
      myfilt = filtered.data(data, factor = myfactors$Factor1, norm = FALSE, method = 1,  cpm = cpmvalue, cv.cutoff=cutoffvalue)
    }else{
      #If the value of filtermethod is != 1,2 or 3, CPM method will be used to filter
      time<-Sys.time() 
      message<-paste(time, " Filtering the reads of ", label," with CPM filter method because the filter method=", filtermethod, " is not an allowed value\n", sep="")
      cat(message, sep="")
      myfilt= filtered.data(data, factor = myfactors$Factor1, norm = FALSE, method =1,  cpm = cpmvalue, cv.cutoff=cutoffvalue)
    }
    #Merging the data and making the ExpressionSet element
    mydata <- readData(data = myfilt, factors = myfactors, length = mylength, biotype = mybiotype, chromosome = mychroms, gc = mygc)
  }else{
    #Merging the data and making the ExpressionSet element
    mydata <- readData(data = data, factors = myfactors, length = mylength, biotype = mybiotype, chromosome = mychroms, gc = mygc)
  }
  #Loading the contrastfile
  contrast<-read.table(file=contrastfile, sep="\t", header=T)
  
  #Performing the DE analysis for each contrast
  filepaths<-NA
  counta<-1
  for (a in 1:length(rownames(contrast))){
    #Obtaining the name of the contrast and the conditions to make the contrast
    contrastsname<-strsplit(as.character(contrast[a,1]), "=")[[1]]
    contrastval<-strsplit(as.character(contrastsname[2]), "-")[[1]]
    #Generating the QCReport for these conditions
    QCname<-file.path(projectdir,paste(label,"_QCReport_NOISeq_",contrastsname[1],".pdf", sep=""))
    filepaths[counta]<-QCname
    #QCreport(mydata, file=QCname, samples=c(contrastval[1],contrastval[2]), factor="Factor1")
    #Computation of differential expression between two experimental conditions from read count data
    mynoiseq <- noiseq(mydata, norm = normethod, factor = "Factor1", conditions=c(contrastval[1],contrastval[2]),  replicates = replicatevalue, k=kvalue, pnr = pnrvalue, nss = nssvalue, v = vvalue, lc = lcvalue)
    mynoiseqname<-file.path(projectdir,paste(label,"_NOISeq_results_",contrastsname[1],".xls", sep=""))
    filepaths[counta+1]<-mynoiseqname
    write.table(mynoiseq@results, file=mynoiseqname, sep = "\t", col.names = NA , row.names = TRUE, qmethod = "double")
    #Obtaining the Differentially expressed elements with q probability
    mynoiseq.deg = degenes(mynoiseq, q = qvalue, M = NULL)
    mynoiseq.degname<-file.path(projectdir,paste(label,"_NOISeq_DE_of_",contrastsname[1],"_with_q_", qvalue,".xls", sep=""))
    filepaths[counta+2]<-mynoiseq.degname
    #write.table(mynoiseq.deg, file=mynoiseq.degname, sep = "\t", col.names = NA , row.names = TRUE, qmethod = "double")
    #Generating the DE plots for the condition analyzed
    plotsname<-file.path(projectdir,paste(label, "_NOISeq_DE_plots_",contrastsname[1],".pdf", sep=""))
    filepaths[counta+3]<-plotsname  
    pdf(file=plotsname, paper="a4")
    par(mfrow=c(1,1), col.main="midnightblue", col.lab="midnightblue", col.axis="midnightblue", bg="white", fg="midnightblue", font=2, cex.axis=0.8, cex.main=1.2)
    DE.plot(mynoiseq, q = qvalue, graphic = "expr", log.scale = TRUE, col="forestgreen", main=paste("Expression plot of ",contrastsname[1], sep=""))
    DE.plot(mynoiseq, q = qvalue, graphic = "MD", col="forestgreen", main=paste("MD plot of ",contrastsname[1], sep=""))
    if(!missing(chromsfile)){
      DE.plot(mynoiseq, q = qvalue, graphic = "chrom", chromosomes = NULL, log.scale = TRUE, join = FALSE, main=paste("Manhattan plot of ",contrastsname[1], sep=""))
      if(!missing(biotypefile)){
        DE.plot(mynoiseq, q = qvalue, graphic = "distr", chromosomes = NULL, log.scale = TRUE, join = FALSE, main=paste("Distribution of DE features of ",contrastsname[1], sep=""))
      }
    }
    #Printing the date and information of the proccess
    time<-Sys.time()
    message<-paste(time, " The following files have been generated in ", projectdir," for ",contrastsname[1]," comparison:\n" ,QCname ,"\n ",mynoiseqname, "\n", mynoiseq.degname,"\n", plotsname, "\n", sep="")
    cat(message, sep="")
    dev.off()
    counta<-counta+4
  }
  return(filepaths)
}
