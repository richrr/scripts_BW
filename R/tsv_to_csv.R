# clean up any pre-existing comma in the first column of file
#### optional: un/comment the add prefix to row ids
# then convert tsv to csv

# can take several files as input

args = commandArgs(trailingOnly=TRUE)

for(infile in args)
{
  print(infile)
  df = read.delim(infile, header=T, check.names=F, row.names=1)

  ID = gsub(",", ".", rownames(df), fixed=TRUE)
   #ID = paste0("ASV_", ID)
  outdf = cbind(ID, df)

  write.csv(outdf, paste(infile, ".csv", sep=''), quote=F, row.names=F)
}
