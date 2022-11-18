# can take several files as input
library(stringr)

args = commandArgs(trailingOnly=TRUE)

BIGDF = ''
for(infile in args)
{
  print(infile)
  df = read.delim(infile, header=F, stringsAsFactors=F, check.names=F, sep="\t")
  df = df[-which(df$V1 == "x"),,drop=F]

  cdf = str_split(df$V1,":")
  pairs = t(as.data.frame(cdf))
	colnames(pairs) = c("ID","Value")
  pairs[,"Value"] = gsub(" ", "", pairs[,"Value"])

  colnames(pairs) = c("ID", infile)
  rownames(pairs) = pairs[,"ID"]
  pairs = pairs[,-1,drop=F]
  #print(pairs)

  if(is.null(dim(BIGDF))){
    BIGDF = pairs
  } else{
    #BIGDF = merge(BIGDF, pairs, all=T, by=1)
    BIGDF = cbind(BIGDF, pairs)
  }

}

BIGDF = cbind(rownames(BIGDF),BIGDF)
head(BIGDF)
write.table(BIGDF, "net-stats-summary.tsv", quote=F, row.names=F, sep='\t')
