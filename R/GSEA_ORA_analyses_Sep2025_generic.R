
# How to run clusterprofiler if:
    #----------------#
    # (A) deseq2 code written by RRR to analyze the RSEM based counts from: 
    #----------------#
    # NRR (RNA-seek (v1.7.1)) RNAseq pipeline
        # /data/rodriguesrr/Fontana_Tosti_Mediterranean/exfoliome_data/README.txt
        # /data/rodriguesrr/Fontana_Tosti_Mediterranean/exfoliome_data/RNAseek_hg38_out/raw_fastq_from_CCRSF/DEG_ALL/RSEM.isoforms.expected_count.all_samples.txt
    # or coming from the CCRSF seq facility.
        # /data/rodriguesrr/Fontana_Tosti_Mediterranean/PBMC_RNAseq/deseq2/run-deseq2.R 

        # take the deseq2-results.csv file and loop over the different analysis
        # alternatively, take the lfcUnshrunk-deseq2-results-geneNames.csv
        # 	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	GeneID	GeneName
        # ENSMUSG00000007777.9	109.4207621	0.066114938	0.236740472	0.279271799	0.780036247	0.99987983	ENSMUSG00000007777.9	0610009B22Rik


    # example usage of enrichment analysis using clusterprofiler:
        #ml R/4.3.2
        # cd /data/rodriguesrr/Fontana_Tosti_Mediterranean/PBMC_RNAseq/deseq2
        # Rscript /data/rodriguesrr/scripts/R/GSEA_analyses_Feb2024.R Timepoint_Diet_Week_00_MD_vs_Week_00_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_08_MD_vs_Week_08_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_16_MD_vs_Week_16_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_08_C_vs_Week_00_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_16_C_vs_Week_08_C-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_08_MD_vs_Week_00_MD-lfcUnshrunk-deseq2-results-geneNames.csv Timepoint_Diet_Week_16_MD_vs_Week_08_MD-lfcUnshrunk-deseq2-results-geneNames.csv
    #-------#-------#



    #----------------#
    # (B) counts from analysis using Olink stats package: 
    #----------------#
    # /data/MicrobiomeCore/Maria-Clavijo-Salomon/analysis/run_script.R 

    # example usage of enrichment analysis using clusterprofiler:
        # ml R/4.4.3
        # cd /data/MicrobiomeCore/Maria-Clavijo-Salomon/analysis/parametric/fgsea
        # Rscript /data/rodriguesrr/scripts/R/GSEA_analyses_Mar2025_Olink.R ../CpG/all_CpG_NR_TT_W3_W1_PatID_sig.pval_34_sig.adj.pval_19_.csv ../CpG/all_CpG_R_TT_W3_W1_PatID_sig.pval_22_sig.adj.pval_7_.csv ../CpG/all_CpG_W1_TT_R_NR__sig.pval_7_sig.adj.pval_2_.csv ../CpG/all_CpG_W3_TT_R_NR__sig.pval_6_sig.adj.pval_0_.csv ../FMT/all_FMT_C2_TT_R_NR__sig.pval_4_sig.adj.pval_0_.csv ../FMT/all_FMT_C3_TT_R_NR__sig.pval_8_sig.adj.pval_0_.csv ../FMT/all_FMT_NR_TT_C2_Pre_PatID_sig.pval_1_sig.adj.pval_0_.csv ../FMT/all_FMT_NR_TT_C3_Pre_PatID_sig.pval_29_sig.adj.pval_1_.csv ../FMT/all_FMT_Pre_TT_R_NR__sig.pval_29_sig.adj.pval_0_.csv ../FMT/all_FMT_R_TT_C2_Pre_PatID_sig.pval_37_sig.adj.pval_0_.csv ../FMT/all_FMT_R_TT_C3_Pre_PatID_sig.pval_1_sig.adj.pval_0_.csv ../vsHD/all_C2_NR_vs_HD_TT_NR_Ctrl__sig.pval_27_sig.adj.pval_2_.csv ../vsHD/all_C2_R_vs_HD_TT_R_Ctrl__sig.pval_32_sig.adj.pval_16_.csv ../vsHD/all_C3_NR_vs_HD_TT_NR_Ctrl__sig.pval_29_sig.adj.pval_13_.csv ../vsHD/all_C3_R_vs_HD_TT_R_Ctrl__sig.pval_26_sig.adj.pval_15_.csv ../vsHD/all_Pre_NR_vs_HD_TT_NR_Ctrl__sig.pval_18_sig.adj.pval_0_.csv ../vsHD/all_Pre_R_vs_HD_TT_R_Ctrl__sig.pval_40_sig.adj.pval_32_.csv
    #-------#-------#



    #----------------#
    # (C) the deseq2 results from nfcore pipeline (differentialabundance 1.5.0) using counts from STAR_SALMON (e.g. rnaseq 3.21.0)
    #----------------#
    # to be implemented #
    # /data/MicrobiomeCore/CCRSF_DME_Globus/NovaSeq_X_Plus/2025_07_11_metatranscript_198_SE_samps_NAS_CS039538_Flowcell_22V7C2LT4/nf-core-analysis/filter_input_files_as_per_contrast_then_nf_deg_v1.5.0_gsea_gprofiler.R
    # /data/MicrobiomeCore/CCRSF_DME_Globus/NovaSeq_X_Plus/2025_07_11_metatranscript_198_SE_samps_NAS_CS039538_Flowcell_22V7C2LT4/nf-core-analysis/manually_run_each_analysis_separately/FMT_Pre_vs_Resp1yr_Adv_PD1_mel_Responder/gsea/tables/differential/FMT_Pre_vs_Resp1yr_Adv_PD1_mel_Responder.deseq2.results.tsv
        # gene_id baseMean        log2FoldChange  lfcSE   pvalue  padj
        # ENSG00000000003 28.73054        0.09902214      0.59188281      0.8013004       0.9008476

    # /data/MicrobiomeCore/CCRSF_DME_Globus/NovaSeq_X_Plus/2025_07_11_metatranscript_198_SE_samps_NAS_CS039538_Flowcell_22V7C2LT4/nf-core-analysis/manually_run_each_analysis_separately/results/enrichment/RvsNR_noNRchange.csv
        #gene_id,gene_name,RvsNR_noNRchange
        #ENSG00000123349,PFDN5,Up

    # example usage of enrichment analysis using clusterprofiler:
        # ml R/4.5.0
        # Rscript /data/rodriguesrr/scripts/R/GSEA_ORA_analyses_Sep2025_generic.R ../RvsNR_noNRchange.csv
    #-------#-------#





# res = msigdbr_collections()

# print(res, n =25)

# ---- deprecated as of July 2025 ---- #
# # A tibble: 23 × 3
##    gs_cat gs_subcat         num_genesets.
#    <chr>  <chr>                    <int>
#  1 C1     ""                         299
#  2 C2     "CGP"                     3384
#  3 C2     "CP"                        29  # msigdbr::msigdbr(species = "Homo sapiens", category="C2", subcategory="CP") is 4197
#  4 C2     "CP:BIOCARTA"              292
#  5 C2     "CP:KEGG"                  186
#  6 C2     "CP:PID"                   196
#  7 C2     "CP:REACTOME"             1615
#  8 C2     "CP:WIKIPATHWAYS"          664
#  9 C3     "MIR:MIRDB"               2377
# 10 C3     "MIR:MIR_Legacy"           221
# 11 C3     "TFT:GTRD"                 518
# 12 C3     "TFT:TFT_Legacy"           610
# 13 C4     "CGN"                      427
# 14 C4     "CM"                       431
# 15 C5     "GO:BP"                   7658
# 16 C5     "GO:CC"                   1006
# 17 C5     "GO:MF"                   1738
# 18 C5     "HPO"                     5071
# 19 C6     ""                         189
# 20 C7     "IMMUNESIGDB"             4872
# 21 C7     "VAX"                      347
# 22 C8     ""                         700
# 23 H      ""                          50


# ---- Sep 2025 uses the following ---- #
# # A tibble: 25 × 4
#    gs_collection gs_subcollection  gs_collection_name               num_genesets
#    <chr>         <chr>             <chr>                                   <int>
#  1 C1            ""                "Positional"                              302
#  2 C2            "CGP"             "Chemical and Genetic Perturbat…         3538
#  3 C2            "CP"              "Canonical Pathways"                       19  # msigdbr::msigdbr(species = "Homo sapiens", collection="C2", subcollection="CP") is 579
#  4 C2            "CP:BIOCARTA"     "BioCarta Pathways"                       292
#  5 C2            "CP:KEGG_LEGACY"  "KEGG Legacy Pathways"                    186
#  6 C2            "CP:KEGG_MEDICUS" "KEGG Medicus Pathways"                   658
#  7 C2            "CP:PID"          "PID Pathways"                            196
#  8 C2            "CP:REACTOME"     "Reactome Pathways"                      1787
#  9 C2            "CP:WIKIPATHWAYS" "WikiPathways"                            885
# 10 C3            "MIR:MIRDB"       "miRDB"                                  2377
# 11 C3            "MIR:MIR_LEGACY"  "MIR_Legacy"                              221
# 12 C3            "TFT:GTRD"        "GTRD"                                    505
# 13 C3            "TFT:TFT_LEGACY"  "TFT_Legacy"                              610
# 14 C4            "3CA"             "Curated Cancer Cell Atlas gene…          148
# 15 C4            "CGN"             "Cancer Gene Neighborhoods"               427
# 16 C4            "CM"              "Cancer Modules"                          431
# 17 C5            "GO:BP"           "GO Biological Process"                  7583
# 18 C5            "GO:CC"           "GO Cellular Component"                  1042
# 19 C5            "GO:MF"           "GO Molecular Function"                  1855
# 20 C5            "HPO"             "Human Phenotype Ontology"               5748
# 21 C6            ""                "Oncogenic Signature"                     189
# 22 C7            "IMMUNESIGDB"     "ImmuneSigDB"                            4872
# 23 C7            "VAX"             "HIPC Vaccine Response"                   347
# 24 C8            ""                "Cell Type Signature"                     866
# 25 H             ""                "Hallmark"                                 50


args = commandArgs(trailingOnly=TRUE)
#print(args)






#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("msigdb", "clusterProfiler", "DOSE"))

library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(msigdbr) # packageVersion("msigdbr") # ‘25.1.1’ updated in July 23 2025
library(clusterProfiler)
library(fgsea)
library(dplyr)
library(stringr)
library(enrichplot)
#library(DOSE)


# output file after performing stats on the olink data
get_file_Olink_results = function(infile){
  indf = read.csv(infile, header=T)
  print(head(indf))
  
  
  #id_names = str_split_fixed(indf[, "gene_id"], "_", 2)
  #colnames(id_names) = c("EnsemblID", "GeneName")
  #print(head(id_names))
  
  #indf = cbind(indf, id_names)
  indf = indf[, c("UniProt", "estimate")]
  colnames(indf) = c("UniProt", "log2FoldChange")
  #rownames(indf) = indf[, "UniProt"]
  print(head(indf))

  outsig = map_UniProt_to_GeneName(infile, indf, "/data/MicrobiomeCore/Maria-Clavijo-Salomon/DAVID_id_conversion/Uniprot_to_Official_GeneSymbol.txt")

  return(outsig)
}


get_file_deg_list = function(infile){
  indf = read.csv(infile, header=T)
  print(head(indf))

  colnames(indf) = c("Ensembl", "Symbol", "FCdir")
  return(indf[, c("Symbol", "FCdir")])
}


map_UniProt_to_GeneName = function(infile, indf, mapper_infile){
    mapdf = read.delim(mapper_infile, header=T, check.names=F, sep="\t")

    mergeddf = merge(indf, mapdf, by=1)
    print(head(mergeddf))
    write.csv(mergeddf, paste0(basename(infile), "_w_GeneSymbol.csv"))

    return(mergeddf)
}



do_everything = function(infile, type = "Olink", method = "GSEA"){
        outsig = NULL

        if(type == "Olink"){
          outsig = get_file_Olink_results(infile)
          # continues below to the GSEA code
        } else if(type == "Ensembl_Symbol_FCdir"){
          outsig = get_file_deg_list(infile)
          # exists out to go ORA enrichment
        }

        ##msigbr to download all genes and their GO IDs #Ref:http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
        msigDB <-  msigdbr::msigdbr(species = "Homo sapiens")


        if(method == "GSEA"){
              #### Ordered gene list with log2FoldChange for GSEA analyses
              genes = outsig[, "log2FoldChange"]
              names(genes) <- outsig$To #outsig$GeneName
              genes<-na.omit(genes)
              genes = sort(genes , decreasing = TRUE)
              #print(head(genes))
              #print(tail(genes))

              
              #### for ORA add some pval or padj cutoff
              #sigGenes <- outsig[ which(outsig[, grep("pvalue", colnames(outsig), value=T)] < 0.05), "EnsemblID"] %>% unique(.)
              #genes = sigGenes
              
              ##### this is not tested in Sep 2025 when the ORA function was written and related code was updated
              for(d in c("REACTOME", "GOBP", "KEGG", "HALLMARK", "CP", "IMMUNESIGDB", "CellTypeSig")){
                  print(paste("#---------", d, "---------#"))
                  gseaResCsv = paste(basename(infile), d, "GSEA.csv", sep="_")
                  gseaResPdf = paste(basename(infile), d, "GSEA-plot.pdf", sep="_")
                  em <- geneEnrichments(genes, backgroundGenes = NULL, qval = 1, database = d, method = "GSEA", msigDB = msigDB, ofilename = gseaResPdf)
                  print(dim(em))
                  print(head(em))
                  write.csv(em, gseaResCsv)
              }
          }




          if(method == "ORA"){
              for(d in c("Cell_marker_2.0_Human", "REACTOME", "GOBP", "KEGG_LEGACY", "KEGG_MEDICUS", "HALLMARK", "WIKIPATHWAYS", "BIOCARTA", "IMMUNESIGDB", "CellTypeSig")){
                print(paste("#---------", d, "---------#"))
                  
                for (direction in c("up", "down", "all")){

                    # get the specific list of genes.
                    gList = outsig[, "Symbol"]
                    if(direction != "all"){
                        gList = outsig[ which( tolower(outsig[, "FCdir"]) == direction ),  "Symbol"]
                    }
                     
                    enrichResCsv = paste(basename(infile), d, direction, "ORA.csv", sep="_")
                    
                    em <- geneEnrichments(gList, backgroundGenes = NULL, qval = 1, database = d, method = "ORA", msigDB = msigDB, ofilename = NULL)
                    print(dim(em))
                    print(head(em))
                    write.csv(em, enrichResCsv)
                }
              }

          }

}



# https://forum.posit.co/t/column-names-as-variables-in-dplyr-select-v-filter/140188  # colName="gs_subcat"  get({{colName}})
get_DB = function(msigDB, subcat , pattern){
  db <- msigDB %>% dplyr::filter(gs_subcollection == subcat) %>% 
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
  db <- msigDB %>% dplyr::filter(gs_collection=="H") %>% 
    mutate(gs_name = str_remove(gs_name, "HALLMARK_")) %>% 
    dplyr::select(gs_name, gs_id, gene_symbol) %>% 
    mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
    mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
    dplyr::select(!c(gs_id))
  
  #print(db)
    
  return(db)
}



get_celltype_DB = function(msigDB){
  db <- msigDB %>% dplyr::filter(gs_collection=="C8") %>% 
    #mutate(gs_name = str_remove(gs_name, "HALLMARK_")) %>% 
    dplyr::select(gs_name, gs_id, gene_symbol) %>% 
    mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
    mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
    dplyr::select(!c(gs_id))
  
  #print(db)
    
  return(db)
}



generate_single_and_multiple_GSEA_plots = function(geneList, pathway, db, ofilename){
    
    
    pathres = as.data.frame(pathway@result)
    idx =  which(pathres[, "qvalue"] <= 0.15, arr.ind = T)
    
    if(length(idx) == 0) {
      return(0)
    } else{
      print(idx)
      print(head(pathres$ID[idx]))
      titl = NULL
      
      if(length(idx) == 1){
          titl = pathway$Description[1]
      } else{
          titl = "Pathways qval <= 0.15"
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
      if(method =="GSEA"){
        generate_single_and_multiple_GSEA_plots(sigGenes, pathway, db, ofilename)  # rankedlist, gseaResult, database
      }
      pathwayData <- pathway@result %>% as.data.frame(.) %>% dplyr::filter(qvalue <= qval)
    }

    return(pathwayData)

}




#### Enrichment function
geneEnrichments <- function(sigGenes, backgroundGenes = NULL, qval = 1, msigDB = msigDB,
                            database=c("KEGG_LEGACY", "KEGG_MEDICUS", "REACTOME", "GOBP", "HALLMARK", "CP", "WIKIPATHWAYS", "BIOCARTA", "IMMUNESIGDB", "CellTypeSig"), method = c("ORA","GSEA"), targetedGO = NULL, ofilename = "outfile.pdf") {

  db_to_use = NULL
  if (database == "CP") {
    db_to_use = get_DB(msigDB, "CP", "____")
  } else if (database == "WIKIPATHWAYS") {
    db_to_use = get_DB(msigDB, "CP:WIKIPATHWAYS", "WP_")
  } else if (database == "BIOCARTA") {
    db_to_use = get_DB(msigDB, "CP:BIOCARTA", "____")
  } else if (database == "KEGG_LEGACY") {
    db_to_use = get_DB(msigDB, "CP:KEGG_LEGACY", "KEGG_")
  } else if (database == "KEGG_MEDICUS") {
    db_to_use = get_DB(msigDB, "CP:KEGG_MEDICUS", "KEGG_")
  } else if (database == "REACTOME"){
    db_to_use = get_DB(msigDB, "CP:REACTOME", "REACTOME_")
  } else if (database == "GOBP"){ # C5 (ontology gene sets, 16008 gene sets) ; GO:BP (GO biological process, 7647 gene sets)
    db_to_use <- get_DB(msigDB, "GO:BP", "GOBP_")
  } else if (database == "HALLMARK"){
    db_to_use <- get_hallmark_DB(msigDB)
  } else if (database == "IMMUNESIGDB"){ 
    db_to_use <- get_DB(msigDB, "IMMUNESIGDB", "____")
  } else if (database == "CellTypeSig"){
    db_to_use <- get_celltype_DB(msigDB)
  } else if (database == "Cell_marker_2.0_Human"){
    library(readxl)
    df = read_excel("/data/MicrobiomeCore/databases/CellMarker_2.0/Cell_marker_Human.xlsx")
    db_to_use <- df %>% dplyr::select(cell_name, marker)
  } 
    
  pathwayData = run_enrichment(sigGenes, backgroundGenes = NULL, qval, msigDB, db_to_use, method, targetedGO, ofilename)

  return(pathwayData)
}


# infile = "../deseq2/Timepoint_Diet_Week_00_MD_vs_Week_00_C-lfcUnshrunk-deseq2-results-geneNames.csv"
print(args)
for(infile in args){
  do_everything(infile, type = "Ensembl_Symbol_FCdir", method = "ORA")
}

print("Done!")