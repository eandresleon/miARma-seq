#########################################################################  
#  F_Analysis R function  		                           									#
#																		                                    #
#	Created at Computational Biology and Bioinformatics Group (CbBio)	    #
#	Institute of Biomedicine of Seville. IBIS (Spain)					            #
#	Copyright (c) 2016 IBIS. All rights reserved.						              #
#	mail : miarmaseq-devel@cbbio.es 								                      #
#########################################################################

#Printing the function information and  help message
write("

#########################################################################  
F_Analysis R function    	                           									
																		                                    
Created at Computational Biology and Bioinformatics Group (CbBio)	    
Institute of Biomedicine of Seville. IBIS (Spain)					            
Copyright (c) 2016 IBIS. All rights reserved.						              
mail : miarmaseq-devel@cbbio.es 								                      
#########################################################################

Get help typing F_Analysis_help()\n",stderr())

#Defining help function
F_Analysis_help<-function(){
  write("
F_Analysis function calls a R function called F_Analysis which takes a tab file with the number of the reads from htseq-count analysis, 
a target file with the experimental conditions of the samples and the contrast file with the contrast to analyze the differential 
expression between the defined samples. This function outputs a pdf file with the descriptive plots of the analysis and a tab file 
(xls extension to easy exploration) for condition evaluated. This function accepts one or two factors experiments.
   
F_Analysis takes 4 arguments:

Mandatory arguments:
  [projectdir] Path of the directory where will be saved the QC report.
  [up] Path of the tab file which contains the number of reads from the htseq analysis
  [down] Path of the tabulated file which contains the experimental condiction of each sample. 
  [universe] Path of the tabulated file which contains the experimental condiction of each sample. 
  [organism] Path of the tabulated file which contains the experimental condiction of each sample. 
  [method] Path of the tabulated file which contains the experimental condiction of each sample. 
  [seq_id] Type of entitie: gene_id or transcript_id
  [label] label

  Example:
    
    F_Analysis(projectdir=\"./Project\",up=\"up\",down=\"down\",universe=\"universe\",organism=\"human\",method=\"edger\",seq_id=\"transcript_id\",label=\"TSA\")", stderr())
  
}

F_Analysis<-function(projectdir,up,down,universe,organism,method,seq_id,mydataset="hsapiens_gene_ensembl",label){
  
  #Checking the mandatory parameters
  if(missing(projectdir) | missing(up) | missing(down) | missing(universe) | missing(organism) | missing(method) | missing(seq_id) | missing(label)){
    F_Analysis_help()
    stop("F_Analysis Error: Missing mandatory argument")
  }
  
  #########################################################################
  #1- IMPORTING AND FORMATING THE DATA
  #########################################################################
  
  #installing need packages
	list.of.packages <- c("goseq","biomaRt","Hmisc","geneLenDataBase","GO.db","org.Hs.eg.db")
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
	library(goseq)
	require(biomaRt)
	require(Hmisc)
	
	#Printing the date and information of the proccess
	time<-Sys.time()
	message<-paste(time, " Starting Functional analysis with goseq ", organism,"\n", sep="")
	cat(message, sep="")
	
	#Importing the data
	workingDir<-projectdir
	setwd(workingDir)

	up_file<-read.table(up)
	up_file<-as.vector(up_file$V1)
	up<-unique(up_file)
	
	down_file<-read.table(down)
	down_file<-as.vector(down_file$V1)
	down<-unique(down_file)
	
	universe_file<-read.table(universe)
	universe_file<-as.vector(universe_file$V1)
	universe<-unique(universe_file)
	
	mapping_table<-NA
	genes_up<-NA
	genes_down<-NA
	
	
	if((tolower(organism) == "human") || (tolower(organism) == "man") || (tolower(organism) == "homo sapiens")){
		mydataset="hsapiens_gene_ensembl"
		mybuild="hg19"
	}
	if((tolower(organism) == "mouse") || (tolower(organism) == "mice") || (tolower(organism) == "mus musculus")){
		mydataset="mmusculus_gene_ensembl"
		mybuild="mm9"
	}
	if((tolower(organism) == "rat") || (tolower(organism) == "rattus") || (tolower(organism) == "rattus norvegicus")){
		mydataset="rnorvegicus_gene_ensembl"
		mybuild="rn5"
	}
	
	#ensembl = useMart("ensembl",dataset=mydataset)
	ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset=mydataset, host="grch37.ensembl.org")
	
	if(tolower(seq_id)=="transcript_id" ){
		mapping_table<-getBM(
		  attributes=c('ensembl_gene_id','ensembl_transcript_id'), 
		  filters = 'ensembl_transcript_id', values=universe, 
		  mart=ensembl,
		  uniqueRows=T
		)
		#up
		idx<-match(up,mapping_table$ensembl_transcript_id)
		up_entrez<-mapping_table[idx,1]
		genes_up<-as.integer(unique(mapping_table$ensembl_gene_id) %nin% up_entrez)
		names(genes_up)<-unique(mapping_table$ensembl_gene_id)
		#down
		idx<-match(down,mapping_table$ensembl_transcript_id)
		down_entrez<-mapping_table[idx,1]
		genes_down<-as.integer(unique(mapping_table$ensembl_gene_id)  %nin% down_entrez)
		names(genes_down)<-unique(mapping_table$ensembl_gene_id)
  	    #de
  	    all_de<-c(up,down)
  	    idx<-match(all_de,mapping_table$ensembl_transcript_id)
  	    de_entrez<-mapping_table[idx,1]
  	    genes_de<-as.integer(unique(mapping_table$ensembl_gene_id)  %nin% de_entrez)
  	    names(genes_de)<-unique(mapping_table$ensembl_gene_id)
		
	} else if(tolower(seq_id)=="gene_name" ){
	  mapping_table<-getBM(
	    attributes=c('ensembl_gene_id','external_gene_name'), 
	    filters = 'external_gene_name', values=universe, 
	    mart=ensembl,
	    uniqueRows=T
	  )
	  #up
	  idx<-match(up,mapping_table$external_gene_name)
	  up_entrez<-mapping_table[idx,1]
	  genes_up<-as.integer(unique(mapping_table$ensembl_gene_id) %nin% up_entrez)
	  names(genes_up)<-unique(mapping_table$ensembl_gene_id)
	  
	  #down
	  idx<-match(down,mapping_table$external_gene_name)
	  down_entrez<-mapping_table[idx,1]
	  genes_down<-as.integer(unique(mapping_table$ensembl_gene_id)  %nin% down_entrez)
	  names(genes_down)<-unique(mapping_table$ensembl_gene_id)
	  
	  #de
	  all_de<-c(up,down)
	  idx<-match(all_de,mapping_table$external_gene_name)
	  de_entrez<-mapping_table[idx,1]
	  genes_de<-as.integer(unique(mapping_table$ensembl_gene_id)  %nin% de_entrez)
	  names(genes_de)<-unique(mapping_table$ensembl_gene_id)
	  table(genes_de)
	  
	}else{
		genes_up<-as.integer(unique(universe) %nin% up)
		names(genes_up)<-unique(universe)
		genes_down<-as.integer(unique(universe) %nin% down)
		names(genes_down)<-unique(universe)
	}
	
	pwf_up=try(nullp(genes_up,mybuild,"ensGene",plot.fit=F),silent=T)
	go_up<-try(goseq(pwf_up,mybuild,'ensGene',method="Wallenius",test.cats=c("GO:CC", "GO:BP", "GO:MF","KEGG")),silent=T)
	file_up<-paste(method,"_", label, "_FAnalysys_Up.xls",sep="")
	
	pwf_dw=try(nullp(genes_down,mybuild,"ensGene",plot.fit=F),silent=T)
	go_dw<-try(goseq(pwf_dw,mybuild,'ensGene',method="Wallenius",test.cats=c("GO:CC", "GO:BP", "GO:MF","KEGG")),silent=T)
	file_down<-paste(method,"_", label, "_FAnalysys_Down.xls",sep="")
	
	pwf_de=try(nullp(genes_de,mybuild,"ensGene",plot.fit=F),silent=T)
	go_de<-try(goseq(pwf_de,mybuild,'ensGene',method="Wallenius",test.cats=c("GO:CC", "GO:BP", "GO:MF","KEGG")),silent=T)
	file_de<-paste(method,"_", label, "_FAnalysys_DE.xls",sep="")
	
	
	pathways<-mapPathwayToName()
	if( (!is.null(nrow(go_up)))&(!is.null(nrow(go_dw)))&(!is.null(nrow(go_de)))){
		#up
		go_up[is.na(go_up$ontology),"ontology"] <- "KEGG"
		kegg_up<-go_up[go_up$ontology=="KEGG",]
		go_up<-go_up[-grep("KEGG",go_up$ontology),]
		kegg_up$term<-pathways[match(kegg_up$category,pathways$path),2]
		write.table(file=file_up,rbind(go_up,kegg_up),sep="\t",col.names=T,row.names=F)
		# down
		go_dw[is.na(go_dw$ontology),"ontology"] <- "KEGG"
		kegg_dw<-go_dw[go_dw$ontology=="KEGG",]
		go_dw<-go_dw[-grep("KEGG",go_dw$ontology),]
		kegg_dw$term<-pathways[match(kegg_dw$category,pathways$path),2]
		write.table(file=file_down,rbind(go_dw,kegg_dw),sep="\t",col.names=T,row.names=F)
		# de
		go_de[is.na(go_de$ontology),"ontology"] <- "KEGG"
		kegg_de<-go_de[go_de$ontology=="KEGG",]
		go_de<-go_de[-grep("KEGG",go_de$ontology),]
		kegg_de$term<-pathways[match(kegg_de$category,pathways$path),2]
		write.table(file=file_de,rbind(go_de,kegg_de),sep="\t",col.names=T,row.names=F)
		
		resultsfiles<-NA
		resultsfiles<-(c(file_up,file_down,file_de))
		return(resultsfiles);
	}	else{
		return(NA)
	}
}

mapPathwayToName <- function() {
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, sep="")
  
  pathway_id_name <- data.frame()
  cont<-0
  for (line in readLines(pathway_list_REST_url)) {
    cont<-cont+1
    tmp <- strsplit(line, "\t")[[1]]
    pathway_id <-gsub("\\(?[a-z|:,.]+","",tmp[1],perl=T)
    pathway_name <- tmp[2]
    pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
    pathway_id_name[cont, 1] = pathway_id
    pathway_id_name[cont, 2] = pathway_name
    
  }
  names(pathway_id_name) <- c("path","pathway_name")
  pathway_id_name
}
