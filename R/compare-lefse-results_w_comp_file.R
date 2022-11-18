args = commandArgs(trailingOnly=TRUE)
library(stringr)


# refer /data/Goldszmid_Lab/analysis/collab_evangelos/JAMS/JAMSbeta/manual_plots/TransNet/analys/lefse/make-dirs-do-metanalysis-lefse.R to see how to build the commands to run this script


#### the last arg is which are the treatment and control groups
treat_control_groups = tail(args, n=1)
splitted = str_split(treat_control_groups, "___vs___")[[1]]
treat_groups = str_split(splitted[1], "~~AND~~")[[1]]
control_groups = str_split(splitted[2], "~~AND~~")[[1]]
#treat_groups
#control_groups


# everything except last
files_to_be_merged = head( args, -1)
#files_to_be_merged



### Note: 
note = "Note: The class count only counts non-empty values. It does not check whether the classes belong to the same 'group' and if the lda score direction is same. Please manually check the results and edit column names before sending to client."

BIGDF = NULL
for(file in files_to_be_merged){
  print(file)
  tmpdf = read.delim(file, header=T, check.names=F)
  tmpdf = tmpdf[,-2] # remove "log average"
  
  # removing the file path prefix
  splitted_path = unlist(strsplit(file, "/"))
  nfile = tail(splitted_path, n=1)
  #print(nfile)

  colnames(tmpdf) = paste(colnames(tmpdf), nfile, sep="__") # paste(nfile, colnames(tmpdf), sep="__") #
  #print(colnames(tmpdf))
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

class.cols = grep("^class__", colnames(sdf), value=T)  # grep("class$", colnames(sdf), value=T)
class.counts = apply(sdf[, class.cols], 1, function(x) sum((!is.na(x)) & (x != "")))
#head(class.counts)

# the treatment group ids might be different in different comparisons.
same.treat.class.counts = apply(sdf[, class.cols], 1, function(x) {sum(x %in% treat_groups)})
#head(same.treat.class.counts)

# the control group ids might be different in different comparisons.
same.control.class.counts = apply(sdf[, class.cols], 1, function(x) {sum(x %in% control_groups)})
#head(same.control.class.counts)

sdf = cbind(sdf, class.counts, same.treat.class.counts, same.control.class.counts)
head(sdf)
write.table(sdf, "lda.eff.size.w.class.counts.merged.lefse.tsv", row.names=F, sep="\t")

print("************************************************")
print(note)
print("************************************************")
