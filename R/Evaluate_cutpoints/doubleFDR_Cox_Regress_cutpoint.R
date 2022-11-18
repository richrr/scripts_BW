args = commandArgs(trailingOnly=TRUE)
library(stringr)

# ml R/4.0.5
# cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/ev-cutpoints/june-15-2021/pitts
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/doubleFDR_Cox_Regress_cutpoint.R lkt-pfs/analysis-results


#### the below can be done using BH FDR (default) or by Storey FDR
# for lkt in dirs
  # for cutp in cutps
    # get pvals
  # calc FDR and get the lowest FDR
  # use the above q vals per lkt
# calc FDR


fdrMethod = "BH"
if(length(args)>1 && tolower(args[2]) == "storey"){
  library(qvalue)
  fdrMethod = "Storey"
}


BIGDF = NULL
for (lkt in list.dirs(path = args[1], full.names = TRUE, recursive = TRUE)){
  if(grepl("LKT", lkt)){
    print(lkt)
    
    cutps = list.files(path = lkt, pattern = "-Results.csv", full.names = T)
    #print(cutps)
    if(length(cutps) < 1){
      next
    }
    cutpdf = c()
    
    for(cutp in cutps){
        df = t(read.csv(cutp, header=T, row.names=1))
        tdf = df[, c("HR", "CI", "pvalue"),drop=F]
        
        if(is.null(nrow(cutpdf))){
          cutpdf = as.data.frame(tdf)
        } else{
          cutpdf = rbind(cutpdf, tdf)
        }
    }
    
    # calculate FDR
    FDR = NULL
    pv = as.numeric(as.character(cutpdf[,"pvalue"]))
    if(fdrMethod == "Storey"){
      #print(min(pv))
      #print(max(pv))
      # return pv as it is if only one entry or if min=max pval
      if(length(pv) == 1 || min(pv) == max(pv)){
          FDR = pv
      } else{
          lambdaf = seq(min(pv),(max(pv)-0.0000001),0.01)
          # /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/ev-cutpoints/june-15-2021/houst/lkt-pfs/analysis-results/Tue_Jun_15_2021_07_25_06\ PM/LKT__s__Megasphaera_Unclassified/*csv
          # # put max value as additional entry to avoid the error  "If length of lambda greater than 1, you need at least 4 values."
          if(length(lambdaf)<4){  
            lambdaf = c(lambdaf, max(pv)) 
          }
          Storey.qobj <- qvalue(p = pv, lambda=lambdaf, pi0.meth = "bootstrap")
          FDR <- unlist(Storey.qobj[3])
      }
    } else{ # BH
      FDR = p.adjust(pv,method="fdr")
    } 
    cutpdf = cbind(cutpdf, FDR)
    
    # name of LKT
    name = str_split(lkt, "/")[[1]]
    name = name[length(name)]
    # if multiple pvalues give same qvalue, this will pick the first pvalue (even if it is higher than some other pvalue)
    #min.cutp.qval = cutpdf[which.min(cutpdf$FDR),] 
    # so pick the lowest pvalue, since it is bound to have lowest qvalue
    min.cutp.qval = cutpdf[which.min(cutpdf$pvalue),]
    rownames(min.cutp.qval) = name
    BIGDF = rbind(BIGDF, min.cutp.qval)
  }

}


# calculate FDR
FDRlkt = NULL
pvl = as.numeric(as.character(BIGDF[,"FDR"]))
if(fdrMethod == "Storey"){
  #print(min(pvl))
  #print(max(pvl))
  # not expected to happen, but return pvl as it is if only one entry or if min=max pval
  if(length(pvl) == 1 || min(pvl) == max(pvl)){
      FDR = pvl
  } else{    
    Storeyl.qobj <- qvalue(p = pvl, lambda=seq(min(pvl),(max(pvl)-0.0000001),0.01), pi0.meth = "bootstrap")
    FDRlkt <- unlist(Storeyl.qobj[3])
  }
} else{ # BH
  FDRlkt = p.adjust(pvl,method="fdr")
}

BIGDF = cbind(BIGDF, FDRlkt)

colnames(BIGDF) = c("HR", "CI", "pvalue", paste(fdrMethod,"FDR1 (qval across cutpoints per lkt)"), paste(fdrMethod,"FDR2 (qval from FDR1 across lkts )"))
print(BIGDF)

write.csv(BIGDF, paste0("double-", fdrMethod,"-FDR-output.csv"),  quote=F)

print("Done!")
q()
