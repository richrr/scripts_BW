args <- commandArgs(trailingOnly = TRUE)


file <- (args[1])
ofile = paste0(file, "-file-for-Multi-var-Cox_Regression.txt")

CoxReg_parsed_result_files = args[-1]
print(CoxReg_parsed_result_files)

pvalueCol =  "Score_(logrank)_test P-value"  # previously just called as "P-value"

BIGDF = NULL
for(f in CoxReg_parsed_result_files){
    df = read.delim(f, header=T, check.names=F)
    #print(head(df))
    
    sdf = NULL
    if(any(grepl("Cutpoint", colnames(df)))){  # from continuous variables with cutpoint
        sdf = df[which(df[, pvalueCol] <= 0.05), c("LKT", pvalueCol, "Cutpoint")]
    } else if(any(grepl("biomarker", colnames(df)))){  # categorical
        sdf = df[which(df[, pvalueCol] <= 0.05), c("LKT", pvalueCol)]
        sdf$Cutpoint = "category"    # 1.5        
    }

    #print(head(sdf))
    
    if(is.null(nrow(BIGDF))){
        BIGDF = sdf
    } else{
        BIGDF = rbind(BIGDF, sdf)
    }
    
}

library(dplyr)
BIGDFsorted <- arrange(BIGDF, across(pvalueCol)) # although the order doesn't matter, the most signif items will be provided first to the regression model
head(BIGDF)
head(BIGDFsorted)




# keep only biomarkers that were significant in the univariate analysis as the infile for the multivariate analysis
table <- read.table(file, sep = "\t", header = TRUE, na.strings=c("", "NA", "NaN", "N_A"))
datadf <- data.frame(table)

#head(datadf)

biomarkers = unique(intersect(BIGDFsorted[, "LKT"], colnames(datadf)))

cutpoints = BIGDFsorted[which(BIGDFsorted[, "LKT"] %in% biomarkers), "Cutpoint"]


keepdf = datadf[, c(biomarkers, "dmfs_event", "dmfs_time")]
write.table(keepdf, ofile, row.names=F, sep="\t", quote=F)




cmd= paste("Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Multi-var-Cox_Regression.R ./", ofile, paste0(biomarkers, collapse=','), paste0(cutpoints, collapse=',') )

print("Run the following command")
cmd

print("Done")
q()
