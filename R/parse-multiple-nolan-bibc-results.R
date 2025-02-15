args = commandArgs(trailingOnly=TRUE)

# adapted from /data/rodriguesrr/scripts/R/parse-multiple-bibc-results.R

# always run it in the directory where you have the bibc results.
# cd /data/rodriguesrr/PD1_Milano/Sep-2024/Longitudinal_Multiomics_Project_LMP/analys/comp_corr/DCR_Y_merge_time/p1_cp1_cfdr1/per_analysis/net/bibc
# ml R/4.4.1
# Rscript /data/rodriguesrr/scripts/R/parse-multiple-nolan-bibc-results.R edges_degree_BiBC_gene_micro.tsv all-bibc-results.csv

patern = args[1] # e.g. "-results.txt"   # tab delimited
outfile = args[2] # e.g. "all-bibc-results.csv"

files = list.files(path = ".", pattern = patern, recursive = T)

df = read.delim(files[1], header=T, sep="\t")
#df = df[, c(1,3)]  # drop Node_degrees
colnames(df) = c("ID", "Node_degree", files[1])
#head(df)
#tail(df)

for(f in files[-1]){
    #print(f)
    tmp_df = read.delim(f, header=T, sep="\t")
    tmp_df = tmp_df[, c(1,3)]  # drop Node_degrees
    colnames(tmp_df) = c("ID", f)
    #head(tmp_df)
    df = merge(df, tmp_df, by="ID")
}

#head(df)
#tail(df)
# To remove both (NAs and empty)   # https://stackoverflow.com/questions/6437164/removing-empty-rows-of-a-data-file-in-r
df <- df[!apply(is.na(df) | df == "", 1, all),]
#head(df)
#tail(df)
write.csv(df, paste0(outfile,".degree.BiBCscore.csv"), row.names=F)

# drop Node_degrees
df = df[ , -which(names(df) %in% c("Node_degree"))]
# give NA a BiBC score of 0 
df[is.na(df)] <- 0
#head(df)
write.csv(df, outfile, row.names=F)


relativize = function(infile){
        df = read.csv(infile, header=T, check.names=F, row.names=1)
        #reldf = df/sum(df) # this is incorrect since it divides by the sum of the matrix
        reldf = prop.table(as.table(as.matrix(df)), 2)
        # this and most of the code in TransNet keeps quote=F, so if the id names have
        # comma it could be problematic.
        # this part removes ',' from the id names, so that the csv files can be used
        # later without running into issues.
        ID = gsub(",", ".", rownames(df), fixed=TRUE)
        reldf = cbind(ID, reldf)
        write.csv(reldf, paste(infile, ".rel.csv", sep=''), quote=F, row.names=F)
        return(reldf)
}


reldf = relativize(outfile)

# higher values will give lower ranks with the method (alternatively directly get rank(-x))
# the ties method of min matches excel
reldf[,-1] <- apply(reldf[,-1], 2, function (x) {rank(1/rank(as.numeric(x) , ties.method = "min") , ties.method = "min")})
write.csv(reldf, paste0(outfile,".rel.ranked.csv"), quote=F, row.names=F)


# convert the ranks to perc by dividing rank by max rank
percrankdf = reldf
percrankdf[,-1] = apply(reldf[,-1], 2, function (x) {(as.numeric(x)*100)/max(as.numeric(x))})

MinPercRank = apply(percrankdf[,-1], 1, function (x) min(as.numeric(x)))

FreqRankTop25Perc = apply(percrankdf[,-1], 1, function (x) sum(as.numeric(x) < 25))

percrankdf = cbind(percrankdf, MinPercRank, FreqRankTop25Perc)

head(percrankdf)

write.csv(percrankdf, paste0(outfile,".perc.ranked.csv"), quote=F, row.names=F)
