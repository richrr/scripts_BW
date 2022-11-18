#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# usage: /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/TransNet/analysis/corr/merged/ENR.ER.LR/p1_cp1_cfdr1/per_analysis
# Rscript /data/rodriguesrr/scripts/R/split_per_combnation_omics_type.R test.csv LKT,ENSG,METACYC,PHENO,INTERPRO,MICROBEGENE
# split 50 million corrs with 95G RAM in 6hrs 


library(stringr)
library(gtools)

df = read.csv(args[1], as.is=T, header=F)
#head(df)

pair = str_split( df$V1 ,"<==>")
pairs = t(as.data.frame(pair))
colnames(pairs) = c("partner1","partner2")
head(pairs)

uniq_omics = unlist(str_split(args[2], ","))
uniq_omics


perms = permutations(n=length(uniq_omics),r=2,v=uniq_omics,repeats.allowed=T)
#perms

sections = list()

#Name: pheno_w_pheno: pheno_w_pheno
#Name: pheno_w_metacyc: pheno_w_metacyc, metacyc_w_pheno

seen_before = c()

for(r in 1:nrow(perms)){
    #print(perms[r,])
    if(perms[r, 1] == perms[r, 2]){
        #print("Within")
        key = paste(c(perms[r, 1], perms[r, 2]),collapse='_w_')
        seen_before = c(seen_before, key)
        sections[[key]] = c(key)
    } else{
        #print("Between")
        key1 = paste(c(perms[r, 1], perms[r, 2]),collapse='_w_')
        key2 = paste(c(perms[r, 2], perms[r, 1]),collapse='_w_')
        if(key1 %in% seen_before || key2 %in% seen_before){
          #pass
        } else{
          seen_before = c(seen_before, key1, key2)
          sections[[key1]] = c(key1, key2)
        }
    }
}


get_rows_w_partial_match = function(vals){
  pattern = paste0("^", unlist(str_split(vals, "_w_"))) # to search for start of string
  print(pattern)
  sdf = pairs[ which( grepl(pattern[1], pairs[,"partner1"]) & 
                    grepl(pattern[2], pairs[,"partner2"]) ), ,drop=F ]
  return(sdf)
}


for (name in names(sections)) {
    print(name)
    vals = sections[[name]]
    print(vals)
    
    edges_in_this_set = ''
    if(length(vals) == 1){
        edges_in_this_set = get_rows_w_partial_match(vals)
    } else if(length(vals) == 2){
        edges_in_this_set = rbind(get_rows_w_partial_match(vals[1]), get_rows_w_partial_match(vals[2]))
    }
    outres = paste(as.vector(edges_in_this_set[,1]),as.vector(edges_in_this_set[,2]),sep="<==>")
    #print(outres)
    write(outres , paste0(name,".txt"))
}

print("Done")
