# clean up any pre-existing comma and punct/space chars in the first column of file
#### optional: un/comment the add prefix to row ids
# then convert tsv to csv
# output the lookup table

# can take several files as input

args = commandArgs(trailingOnly=TRUE)

for(infile in args)
{
  print(infile)
  df = read.delim(infile, header=T, check.names=F, row.names=1)

  ID =  gsub("(?![._-])[[:punct:][:space:]]", "_", rownames(df), perl=T) # "-What_-_.-"  # added on May 12 2025

  #ID = gsub(",", ".", rownames(df), fixed=TRUE) # commented on May 12 2025
   #ID = paste0("ASV_", ID)
  outdf = cbind(ID, df)
  write.csv(outdf, paste0(infile, ".csv"), quote=F, row.names=F)

  lookuptable = cbind(ID, rownames(df))
  colnames(lookuptable) = c("ID", "Original_ID")
  write.table(lookuptable, paste0(infile, ".lookup-table.tsv"), sep = "\t", quote=T, row.names=F)

}
