args = commandArgs(trailingOnly=TRUE)

# read excel file
# for ppm sheets (sanitize name; prefix for each namespace)



#cd /data/rodriguesrr/Koltsova/data/data-info/IL17-rnaseq-liv-ile-aor-metag-feces
#ml R/3.6.0
#ml R/4.3.0
# Rscript /data/rodriguesrr/scripts_BW/R/prep-shotgun-data.R /data/rodriguesrr/Koltsova/analysis/Nov2021/JAMS/JAMSbeta/batch_correct_Project_IL17/Project_IL17_Metadata_and_Relabund_PPM_light_2021-12-27.xlsx

#R


library(readxl)


infile = args[1]

shits = excel_sheets(infile)
PPM_shits = grep("_PPM", shits, value=T)
PPM_shits

BIGDF = NULL
for(f in PPM_shits){

  if(f %in% c("FeatType_PPM", "resfinder_PPM")){
      next
  }


  print(f)
  df = data.frame(read_excel(infile, sheet=f))
  print(head(df))
  print(dim(df))
  
  # https://localcoder.org/in-r-use-gsub-to-remove-all-punctuation-except-period#solution_3
  # https://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html
  df[, "row.names"] =  gsub("(?![._-])[[:punct:][:space:]]", "_", df[, "row.names"], perl=T) # "-What_-_.-"
  
  if(f %in% c("Product_PPM")){
      df[, "row.names"] = paste0("PROD_", df[, "row.names"])
  }

  print(head(df[, 1]))

  
  if(is.null(nrow(BIGDF))){
    BIGDF = df
  } else{
    BIGDF = rbind(BIGDF, df)
  }
}

print(dim(BIGDF))

write.csv(BIGDF, "concatanated-PPM-files.csv", quote=F, row.names=F)
