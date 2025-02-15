#ml R/4.3.2
# cd /data/rodriguesrr/Fontana_Tosti_Mediterranean/PBMC_RNAseq/deseq2
# Rscript /data/rodriguesrr/scripts/R/GSEA_analyses_Feb2024.R Timepoint_Diet_Week_00_MD_vs_Week_00_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_08_MD_vs_Week_08_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_16_MD_vs_Week_16_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_08_C_vs_Week_00_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_16_C_vs_Week_08_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_08_MD_vs_Week_00_MD-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_16_MD_vs_Week_08_MD-lfcUnshrunk-deseq2-results-geneNames.csv

args = commandArgs(trailingOnly=TRUE)

# take the deseq2-results.csv file and loop over the different analysis
# alternatively, take the lfcUnshrunk-deseq2-results-geneNames.csv

# 	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	GeneID	GeneName
# ENSMUSG00000007777.9	109.4207621	0.066114938	0.236740472	0.279271799	0.780036247	0.99987983	ENSMUSG00000007777.9	0610009B22Rik

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("msigdb", "clusterProfiler", "DOSE"))

library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(dplyr)
library(stringr)
library(enrichplot)
#library(DOSE)


# lfcUnshrunk-deseq2-results-geneNames.csv
get_file = function(infile){
  indf = read.csv(infile, header=T)
  #print(head(indf))
  
  id_names = str_split_fixed(indf[, "gene_id"], "_", 2)
  colnames(id_names) = c("EnsemblID", "GeneName")
  #print(head(id_names))
  
  indf = cbind(indf, id_names)
  rownames(indf) = indf[, "gene_id"]
  print(head(indf))
  
  return(indf)
}





do_everything = function(infile){
        outsig = get_file(infile)
        
        #### Ordered gene list with log2FoldChange for GSEA analyses
        genes = outsig[, grep("log2FoldChange", colnames(outsig), value=T)]
        names(genes) <- outsig$GeneName
        genes<-na.omit(genes)
        genes = sort(genes , decreasing = TRUE)
        #print(head(genes))
        #print(tail(genes))

        
        #### for ORA add some pval or padj cutoff
        #sigGenes <- outsig[ which(outsig[, grep("pvalue", colnames(outsig), value=T)] < 0.05), "EnsemblID"] %>% unique(.)
        #genes = sigGenes
        
        ##msigbr to download all genes and their GO IDs #Ref:http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
        msigDB <-  msigdbr::msigdbr(species = "Homo sapiens")

        
        for(d in c("REACTOME", "GOBP", "KEGG", "HALLMARK", "CP")){
            print(paste("#---------", d, "---------#"))
            gseaResCsv = paste(infile, d, "GSEA.csv", sep="_")
            gseaResPdf = paste(infile, d, "GSEA-plot.pdf", sep="_")
            em <- geneEnrichments(genes, backgroundGenes = NULL, qval = 1, database = d, method = "GSEA", msigDB = msigDB, ofilename = gseaResPdf)
            print(dim(em))
            #print(head(em))
            #write.csv(em, gseaResCsv)
        }


}


# https://forum.posit.co/t/column-names-as-variables-in-dplyr-select-v-filter/140188  # colName="gs_subcat"  get({{colName}})
get_DB = function(msigDB, subcat , pattern){
  db <- msigDB %>% dplyr::filter(gs_subcat == subcat) %>% 
    mutate(gs_name = str_remove(gs_name, pattern)) %>% 
    dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% 
    mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
    mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% 
    dplyr::select(!c(gs_exact_source))
    
  #print(db)
  
  return(db)
}


#-----#
# works fine. do not merge with above function for gs_subcat
#----#
get_hallmark_DB = function(msigDB){
  db <- msigDB %>% dplyr::filter(gs_cat=="H") %>% 
    mutate(gs_name = str_remove(gs_name, "HALLMARK_")) %>% 
    dplyr::select(gs_name, gs_id, gene_symbol) %>% 
    mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
    mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
    dplyr::select(!c(gs_id))
  
  #print(db)
    
  return(db)
}



generate_single_and_multiple_GSEA_plots = function(geneList, pathway, db, ofilename){
    
    
    pathres = as.data.frame(pathway@result)
    idx =  which(pathres[, "qvalue"] <= 0.1, arr.ind = T)
    
    if(length(idx) == 0) {
      return(0)
    } else{
      print(idx)
      print(head(pathres$ID[idx]))
      titl = NULL
      
      if(length(idx) == 1){
          titl = pathway$Description[1]
      } else{
          titl = "Pathways qval <= 0.1"
      }
      p = gseaplot2(pathway, geneSetID = idx, title = titl)  #, pvalue_table = TRUE, base_size = 7)
      
      # works for one plot at a time
      # (1)
        # convert the db to list where the name is the term and the entries are the genes
        # db_for_fgsea = split(db$gene_symbol, db$gs_name)
        # p = plotEnrichment(pathway=db_for_fgsea[["R-HSA-72662:Activation of the mrna upon binding of the cap binding complex and eifs and subsequent binding to 43s"]],
        #                    stats = geneList)
      # (2)
        # p = gseaplot2(pathway, geneSetID = 1, title = pathway$Description[1])
      

      pdf(ofilename)
      print(p)
      dev.off()
      
      return(1)
    }

}




run_enrichment = function(sigGenes, backgroundGenes = NULL, qval = 1, msigDB = msigDB,
                            db, method, targetedGO = NULL, ofilename){
    pathway = NULL
    backgroundGenes <- unique(backgroundGenes)
    pathwayData = NULL

    if(method =="ORA"){ 
      sigGenes <- unique(sigGenes)
      pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=db, universe = backgroundGenes) 
    }

    if(method =="GSEA"){ 
      pathway <- clusterProfiler::GSEA(sigGenes, TERM2GENE=db, pvalueCutoff = qval, pAdjustMethod = "BH")
    }

    print(pathway)
    

    if(is.null(pathway)==TRUE){
      pathwayData <- as.data.frame(matrix(nrow=0, ncol = 3))
    } else{
      generate_single_and_multiple_GSEA_plots(sigGenes, pathway, db, ofilename)  # rankedlist, gseaResult, database
      pathwayData <- pathway@result %>% as.data.frame(.) %>% dplyr::filter(qvalue <= qval)
    }

    return(pathwayData)

}





#### Enrichment function
geneEnrichments <- function(sigGenes, backgroundGenes = NULL, qval = 1, msigDB = msigDB,
                            database=c("KEGG", "REACTOME", "GOBP", "HALLMARK", "CP"), method = c("ORA","GSEA"), targetedGO = NULL, ofilename = "outfile.pdf") {

  db_to_use = NULL
  if (database == "CP") {
    db_to_use = get_DB(msigDB, "CP", "____")
  } else if (database == "KEGG") {
    db_to_use = get_DB(msigDB, "CP:KEGG", "KEGG_")
  } else if (database == "REACTOME"){
    db_to_use = get_DB(msigDB, "CP:REACTOME", "REACTOME_")
  } else if (database == "GOBP"){ # C5 (ontology gene sets, 16008 gene sets) ; GO:BP (GO biological process, 7647 gene sets)
    db_to_use <- get_DB(msigDB, "GO:BP", "GOBP_")
  } else if (database == "HALLMARK"){
    db_to_use <- get_hallmark_DB(msigDB)
  }
  
  pathwayData = run_enrichment(sigGenes, backgroundGenes = NULL, qval, msigDB, db_to_use, method, targetedGO, ofilename)

  return(pathwayData)
}


# infile = "../deseq2/Timepoint_Diet_Week_00_MD_vs_Week_00_C-lfcUnshrunk-deseq2-results-geneNames.csv"
for(infile in args){
  do_everything(infile)
}
