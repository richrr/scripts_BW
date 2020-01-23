args = commandArgs(trailingOnly=TRUE)


infile = args[1] # the tab-delimited otu table to be filtered (typically ...../exported_201910071340/otu_table.txt). The first line is skipped.
keepSampsFile = args[2]  # the tab separated map file containing which samples are to be kept. ALWAYS the first column
outPrefix = args[3]


data = read.delim(infile, header=2, skip=1, as.is = T, check.names = F, row.names=1, stringsAsFactors=F)
head(data)
dim(data)


keep = read.delim(keepSampsFile, header=T, as.is = T, check.names = F, row.names=1, stringsAsFactors=F)
head(keep)
dim(keep)


keepData = data[, rownames(keep)]
head(keepData)
dim(keepData)

write.csv(keepData, paste0(outPrefix, "_otu_table.csv"), quote=F)
