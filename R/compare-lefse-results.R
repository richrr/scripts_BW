args = commandArgs(trailingOnly=TRUE)


### Note: 
note = "Note: The class count only counts non-empty values. It does not check whether the classes belong to the same 'group' and if the lda score direction is same. Please manually check the results and edit column names before sending to client."

BIGDF = NULL
for(file in args){
  print(file)
  tmpdf = read.delim(file, header=T, check.names=F)
  tmpdf = tmpdf[,-2] # remove "log average"
  colnames(tmpdf) = paste(file, colnames(tmpdf), sep="__")
  colnames(tmpdf)[1] = "ID"
  #print(head(tmpdf))
  
  if(is.null(nrow(BIGDF))){
    BIGDF = tmpdf
  } else{
    BIGDF = merge(BIGDF, tmpdf, by=1, all=T)
  }
  
}
rownames(BIGDF) = BIGDF[,"ID"]
#head(BIGDF)
dim(BIGDF)
write.table(BIGDF, "merged.lefse.tsv", row.names=F, sep="\t")

lda.eff.size.cols = grep("LDA effect size", colnames(BIGDF), value=T)
d = BIGDF[,c(lda.eff.size.cols)]
#head(d)
d_new <- d[apply(d, 1, function(x) any(!is.na(x))), ]
#head(d_new)
sdf = BIGDF[rownames(d_new), ]

class.cols = grep("class$", colnames(sdf), value=T)
class.counts = apply(sdf[, class.cols], 1, function(x) sum((!is.na(x)) & (x != "")))
#head(class.counts)

# the treatment group ids might be different in different comparisons.
#same.class.counts = apply(sdf[, class.cols], 1, function(x) {length(unique(x))==1})
#head(same.class.counts)

sdf = cbind(sdf, class.counts)
#sdf = cbind(sdf, class.counts, same.class.counts)
head(sdf)
write.table(sdf, "lda.eff.size.w.class.counts.merged.lefse.tsv", row.names=F, sep="\t")

print("************************************************")
print(note)
print("************************************************")
