
library(JAMS)
source("/data/rodriguesrr/scripts/R/fix_qiime2_taxonomy.R")


create_JAMS_LKT_SExpObj_from_qiime2_16S = function(map_file_tsv, taxonomy_file_tsv, raw_count_table_file_tsv){

    #Load metadata
    phenotable <- fread(file = map_file_tsv, data.table = FALSE)
    rownames(phenotable) <- phenotable$SampleID
    phenotable$Sample <- phenotable$SampleID
    phenotable$SampleID <- NULL
    
    phenotable[is.na(phenotable)] <- "N_A"
    
    print(head(phenotable))


    #Load taxonomy table
    taxtable <- fix_qiime2_taxonomy(taxonomy_file_tsv)
    print(head(taxtable))

    tt <- taxtable
    print(head(tt))


    #Load counts
    ctstable <- read.delim(file = raw_count_table_file_tsv, row.names=1, header=T, check.names=F)
    print(head(ctstable))

    #Maintain same order as cts
    tt <- tt[rownames(ctstable), ]
    print(head(tt))


    #Call ASVs as non-redundant LKTs
    tt$LKT <- paste(tt$LKT, rownames(tt), sep="_ASV_")
    print(head(tt))


    rownames(tt) <- tt$LKT
    rownames(ctstable) <- tt$LKT

    
    pheno <- phenotable
    if(setequal(rownames(phenotable), colnames(ctstable))){
          # do nothing
    } else{
          print("Samples missing in either map or count file(s). Keeping samples mentioned in map file.")
          print(setdiff(rownames(phenotable), colnames(ctstable)))
          
          print(setdiff(colnames(ctstable), rownames(phenotable)))
          
          ctstable = ctstable[, rownames(phenotable)]
          print(dim(ctstable))
          
    }

    pheno <- pheno[colnames(ctstable), ]
    print(dim(pheno))
    
    cdict <- make_colour_dictionary(pheno = pheno, class_to_ignore = "N_A", colour_of_class_to_ignore = "#bcc2c2", colour_table = NULL, shuffle = FALSE)

    
    #Make SummarizedExperiment object
    expvec <- list()
    expvec$LKT <- make_SummarizedExperiment_from_tables(pheno = pheno, counttable = ctstable, featuretable = tt, analysisname = "LKT")



      #--- ignore ---#
      if(FALSE){
            #Create ASV to LKT dictionary
            asv2lkt <- data.frame(ASV = rownames(tt), LKT = tt$LKT, stringsAsFactors = FALSE)

            #Add bells and whistles to experiment object
            metadata(expvec$LKT)$asv2lkt <- asv2lkt
            ctable <- plyr::rbind.fill(cdict)
            ctable <- ctable[!duplicated(ctable$Name), ]
            rownames(ctable) <- ctable$Name
            for(anl in names(expvec)){
                metadata(expvec[[anl]])$ctable <- ctable
            }
      }



    #Export
    saveRDS(expvec$LKT, file = "RR.ExpObj.rds")

    #testing_ = read_rds("ExpObj.rds")

    #BECAUSE THERE ARE MORE ASVs THAN LKTs, YOU HAVE TO GLOM BY GENUS (SPECIES)
    #DONT FORGET TO glomby="Species" IN EVERY PLOTTING FUNCTION OR THE JAMS PLOTTING FUNCTIONS WILL FAIL

    return(expvec$LKT)
}
