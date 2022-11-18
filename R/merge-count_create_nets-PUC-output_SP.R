args = commandArgs(trailingOnly=TRUE)

library(stringr)


library(reshape2)

phead = function(indf){
  print(head(indf))
}

ptail = function(indf){
  print(tail(indf))
}


keep_cols = function(f){
  sdf = read.delim(f, header=T, check.names=F,sep="\t")
  sdf$ID = paste0(sdf$partner1 , "<==>", sdf$partner2)
  phead(sdf)

  sdf = sdf[, !(names(sdf) %in% c("partner1","partner2", "Feature ID", "Taxon","Confidence"))]
  sdf <- sdf[,c(ncol(sdf),1:ncol(sdf)-1)]

  prefix = strsplit(f,"/")[[1]][1]
  print(prefix)
  colnames(sdf) = paste(prefix, colnames(sdf), sep="~")

  #phead(sdf)
  #ptail(sdf)

  return(sdf)
}



BIGDF = NULL
for(f in args){
  print(f)
  wdf = keep_cols(f)

  if(is.null(dim(BIGDF))){
    BIGDF = wdf
  } else{
    BIGDF = merge(BIGDF, wdf, all=T, by=1)
  }
}

# replace correlation coeff "NA" with 0 so we can work it (sort/filter/etc) in excel.
BIGDF[is.na(BIGDF)] <- 0
colnames(BIGDF)[1] = 'ID'


# make p1 and p2
pair = str_split( BIGDF$ID ,"<==>")
pairs = t(as.data.frame(pair))
colnames(pairs) = c("partner1","partner2")
rownames(pairs) = BIGDF$ID
#head(pairs)

pos = apply(BIGDF[,grep("combinedCoefficient", colnames(BIGDF), value=T)], 1, function(x) sum(x>0))
neg = apply(BIGDF[,grep("combinedCoefficient", colnames(BIGDF), value=T)], 1, function(x) sum(x<0))
tmp = cbind(pos, neg)
#head(tmp)

abs_diff = abs(tmp[,"pos",drop=F] - tmp[,"neg",drop=F])
colnames(abs_diff) = c("abs_diff")
#head(abs_diff)

abs_median = apply(BIGDF[,grep("combinedCoefficient", colnames(BIGDF), value=T)], 1, function(x) abs(median(x[x!=0])))
#head(abs_median)

tmp = cbind(pairs, tmp)
tmp = cbind(tmp, abs_diff)
tmp = cbind(tmp, abs_median)
BIGDF = cbind(tmp, BIGDF)
head(BIGDF)


#add tax and comp results

############
# add taxonomy
############
taxonomy = read.delim("/data/MicrobiomeCore/analysis/stephanie_20190528/allruns_q2d2/exported/taxonomy.tsv", header=T, check.names=F, sep="\t")
head(taxonomy)
ASVID = paste0("ASV_", taxonomy[,1])
taxonomy = cbind(ASVID, taxonomy)
BIGDF = merge(BIGDF, taxonomy, by=1, all.x = T)
head(BIGDF)

write.table(BIGDF, "BIG-corr-summary-long.tsv", quote=F, row.names=F, sep='\t')



############
# add comp summary
############
# all tissues
compdf = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/microbe-meta-analysis-comp-summary-table-consis_fc.tsv", sep='\t', header=T, check.names=F)
#head(compdf)

BIGDF1 = merge(BIGDF, compdf, by=1, all.x = T)
head(BIGDF1)

write.table(BIGDF1, "BIG-corr-summary-long-w-comp-summary-consis_fc.tsv", quote=F, row.names=F, sep='\t')

print("Done!")
