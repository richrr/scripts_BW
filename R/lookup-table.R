args = commandArgs(trailingOnly=TRUE)


#cd /data/rodriguesrr/Koltsova/data/data-info/IL17-rnaseq-liv-ile-aor-metag-feces
#ml R/3.6.0
#ml R/4.3.0
#Rscript lookup-table.R /data/rodriguesrr/Koltsova/analysis/Nov2021/JAMS/JAMSbeta/batch_correct_Project_IL17/Project_IL17_Metadata_and_Relabund_PPM_light_2021-12-27.xlsx

#R

library(readxl)


infile = args[1]

shits = excel_sheets(infile)
fture_shits = grep("_featuretable", shits, value=T)
fture_shits



ftureBIGDF = NULL
for(f in fture_shits){

  if(f %in% c("FeatType_featuretable", "resfinder_featuretable")){
      next
  }

  
    print(f)
    df = data.frame(read_excel(infile, sheet=f))
    #print(head(df))
    print(dim(df))
    
    # https://localcoder.org/in-r-use-gsub-to-remove-all-punctuation-except-period#solution_3
    # https://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html
    df[, "row.names"] =  gsub("(?![._-])[[:punct:][:space:]]", "_", df[, "row.names"], perl=T) # "-What_-_.-"


    if(f %in% c("LKT_featuretable")){
        df = df[, c("row.names", "LKT")]
    } else if(f %in% c("Product_featuretable")){
        df[, "row.names"] = paste0("PROD_", df[, "row.names"])
        df = df[, c("row.names", "Accession")]
    } else {  # if(f %in% c("ECNumber_featuretable", "Pfam_featuretable", "SUPERFAMILY_featuretable", "TIGRFAM_featuretable", "Interpro_featuretable", "GO_featuretable","MetaCyc_featuretable"))
        df = df[, c("row.names", "Description")]
    } 
    
    colnames(df) = c("ID", "Value")

    #print(head(df[, 1]))

    
     if(is.null(nrow(ftureBIGDF))){
       ftureBIGDF = df
     } else{
       ftureBIGDF = rbind(ftureBIGDF, df)
     }
   }


   #print(head(ftureBIGDF))
   print(dim(ftureBIGDF))
   write.table(ftureBIGDF, "lookup-table.tsv", quote=F, row.names=F, sep="\t")
