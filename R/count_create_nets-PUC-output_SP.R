args = commandArgs(trailingOnly=TRUE)

library(stringr)


library(reshape2)

phead = function(indf){
  print(head(indf))
}

ptail = function(indf){
  print(tail(indf))
}

keep_PUC = function(f){
  df = read.csv(f, header=T, check.names=F)

  sdf = df[-which(grepl("Percentage of Unexpected Correlation=", df$partner1)  |   df$partner1 == 'x') ,
                 c("partner1", "partner2", "combinedCoefficient", "PUC")]
  sdf$ID = paste0(sdf$partner1 , "<==>", sdf$partner2)
  sdf = sdf[, c("ID", "partner1", "partner2", "combinedCoefficient", "PUC")]
  sdf = sdf[which(sdf$PUC == 1), c("ID", "combinedCoefficient")]


  #prefix = strsplit(f,"/")[[1]][1]  # when only running stool
  prefix = paste(strsplit(f,"/")[[1]][1], strsplit(f,"/")[[1]][3], sep="/")
  #print(prefix)
  colnames(sdf) = paste(prefix, colnames(sdf), sep="~")

  #phead(sdf)
  #ptail(sdf)

  return(sdf)
}



BIGDF = NULL
for(f in args){
  print(f)
  wdf = keep_PUC(f)

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


# sum positive and negative correlations across groups
rownames_for_df = BIGDF$ID
rownames(BIGDF) = BIGDF$ID
BIGDF = BIGDF[,-1]

pos = apply(BIGDF, 1, function(x) sum(x>0))
neg = apply(BIGDF, 1, function(x) sum(x<0))
tmp = cbind(pos, neg)
#head(tmp)


tmp = cbind(pairs, tmp)
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

write.table(BIGDF, "corr-summary-long.tsv", quote=F, row.names=F, sep='\t')



############
# add comp summary
############
# all tissues
compdf = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/microbe-meta-analysis-comp-summary-table-consis_fc.tsv", sep='\t', header=T, check.names=F)
#head(compdf)

BIGDF1 = merge(BIGDF, compdf, by=1, all.x = T)
head(BIGDF1)

write.table(BIGDF1, "corr-summary-long-w-comp-summary-consis_fc.tsv", quote=F, row.names=F, sep='\t')



# ran only for stool
if(FALSE){
compdf = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/stool/microbe-meta-analysis-comp-summary-table.tsv", sep='\t', header=T, check.names=F)
#head(compdf)

BIGDF1 = merge(BIGDF, compdf, by=1, all.x = T)
head(BIGDF1)

write.table(BIGDF1, "corr-summary-long-w-comp-summary.tsv", quote=F, row.names=F, sep='\t')



compdf = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/stool/microbe-meta-analysis-comp-summary-table-prefisher-cut.tsv", sep='\t', header=T, check.names=F)
#head(compdf)

BIGDF2 = merge(BIGDF, compdf, by=1, all.x = T)
head(BIGDF2)

write.table(BIGDF2, "corr-summary-long-w-comp-summary-prefisher-cut.tsv", quote=F, row.names=F, sep='\t')
}

print("Done!")
