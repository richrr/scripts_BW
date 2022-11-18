args = commandArgs(trailingOnly=TRUE)


library(reshape2)

phead = function(indf){
  print(head(indf))
}

long_to_wide = function(f){
  df = read.csv(f, header=T, check.names=F)
  sdf = df[, c("partner1", "partner2", "combinedCoefficient")]
  #phead(sdf)

  data_wide <- dcast(sdf, partner1 ~ partner2, value.var="combinedCoefficient")
  # replace correlation coeff "NA" with 0 so we can work it (sort/filter/etc) in excel.
  data_wide[is.na(data_wide)] <- 0
  #phead(data_wide)

  write.csv(data_wide, paste0(f,"-wide.csv"), row.names=F, quote=F)

  colnames(data_wide) = paste(strsplit(f,"/")[[1]][1], colnames(data_wide), sep="~")
  phead(data_wide)
  return(data_wide)
}


BIGDF = NULL
for(f in args){
  print(f)
  wdf = long_to_wide(f)

  if(is.null(dim(BIGDF))){
    BIGDF = wdf
  } else{
    BIGDF = merge(BIGDF, wdf, all=T, by=1)
  }
}

# replace correlation coeff "NA" with 0 so we can work it (sort/filter/etc) in excel.
BIGDF[is.na(BIGDF)] <- 0
colnames(BIGDF)[1] = 'ID'



############
# add taxonomy
############
taxonomy = read.delim("/data/MicrobiomeCore/analysis/stephanie_20190528/allruns_q2d2/exported/taxonomy.tsv", header=T, check.names=F, sep="\t")
head(taxonomy)
ASVID = paste0("ASV_", taxonomy[,1])
taxonomy = cbind(ASVID, taxonomy)
BIGDF = merge(BIGDF, taxonomy, by=1, all.x = T)
head(BIGDF)

write.table(BIGDF, "corr-summary-wide.tsv", quote=F, row.names=F, sep='\t')



############
# add comp summary
############
compdf = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/stool/microbe-meta-analysis-comp-summary-table.tsv", sep='\t', header=T, check.names=F)
#head(compdf)

BIGDF = merge(BIGDF, compdf, by=1, all.x = T)
head(BIGDF)

write.table(BIGDF, "corr-summary-wide-w-comp-summary.tsv", quote=F, row.names=F, sep='\t')


print("Done!")
