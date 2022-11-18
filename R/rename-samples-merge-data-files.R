args = commandArgs(trailingOnly=TRUE)

# map file [args1] has the new column ID [args2] which is to be used to (subset and) rename the samples in datafile [args3]
# add missing columns with NA in the datafiles so they have the same number and names of columns, then do rbind
#     the above step is similar to tranpose, merge all, transpose (basically keeps all columns from both data files)
# write output to file [args5]

map1 = read.delim(args[1], header=T, check.names=F, as.is=T, colClasses = c("character"))
head(map1)

#map2 = read.delim(args[2], header=T, check.names=F, as.is=T, colClasses = c("character"))
#head(map2)

newColName = args[2]

data1 = read.csv(args[3], check.names=F, header=T, row.names=1)
data1 = data1[, as.vector(unlist(map1[, "ID"]))]
#head(data1)

colnames(data1) = as.vector(unlist(map1[, newColName]))
ID = rownames(data1)
data1 = cbind(ID, data1)
head(data1)


data2 = read.csv(args[4], check.names=F, header=T)
head(data2)

dim(data1)
dim(data2)


# fill in non-overlapping columns with NAs
data1[setdiff(names(data2), names(data1))] <- NA
data2[setdiff(names(data1), names(data2))] <- NA


#head(data2)
dim(data1)
dim(data2)


outdf = rbind(data1, data2)
dim(outdf)

write.csv(outdf, args[5], row.names=F, quote=F)
