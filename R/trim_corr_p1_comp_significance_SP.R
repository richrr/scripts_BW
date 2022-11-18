args = commandArgs(trailingOnly=TRUE)

infile = args[1]
micro_meta = args[2]



df = read.delim(infile, check.names=F, header=T, sep="\t")
head(df)

mm = read.delim(micro_meta, check.names=F, header=T, sep="\t", as.is=T)
head(mm)

filter_non_sig_microbes = function(patt1, patt2){

    #print(patt1)
    #print(patt2)


    subdf = df[, c(patt1, patt2)]

    res = apply(subdf, 1, function(x){
            ifelse((!is.na(x[2]) && (x[2] < 0.1)), x[1], 0)
    })
    return(res)

}


outdf = df[, c(1,2)]

# for each corr and comp, keep corr only for dam w fisher < 0.1
for(entry in 1:nrow(mm)){
  organ = mm[entry, 1]
  corr = mm[entry, 2]
  comp = mm[entry, 3]

  patt1 = paste(organ, corr, "combinedCoefficient", sep="~")
  patt2 = paste0( paste(organ, comp, sep='/'), " ~ Analys1~combinedPvalue")

  res = filter_non_sig_microbes(patt1, patt2)
  outdf = cbind(outdf, res)
  colnames(outdf)[ncol(outdf)] = paste0(patt1, "~DAM")
}


head(outdf)



pos_DAM = apply(outdf[, -c(1,2)], 1, function(x) sum(x>0))
neg_DAM = apply(outdf[, -c(1,2)], 1, function(x) sum(x<0))
tmp = cbind(pos_DAM, neg_DAM)
head(tmp)

abs_diff_DAM = abs(tmp[,"pos_DAM",drop=F] - tmp[,"neg_DAM",drop=F])
colnames(abs_diff_DAM) = c("abs_diff_DAM")
#head(abs_diff_DAM)

abs_median_DAM = apply(outdf[,grep("combinedCoefficient", colnames(outdf), value=T)], 1, function(x) abs(median(x[x!=0])))
abs_median_DAM[is.na(abs_median_DAM)] = 0
#head(abs_median_DAM)

tmp = cbind(outdf, tmp)
tmp = cbind(tmp, abs_diff_DAM)
BIGDF = cbind(tmp, abs_median_DAM)
head(BIGDF)






############
# add taxonomy
############
taxonomy = read.delim("/data/MicrobiomeCore/analysis/stephanie_20190528/allruns_q2d2/exported/taxonomy.tsv", header=T, check.names=F, sep="\t")
head(taxonomy)
ASVID = paste0("ASV_", taxonomy[,1])
taxonomy = cbind(ASVID, taxonomy)
BIGDF = merge(BIGDF, taxonomy, by=1, all.x = T)
head(BIGDF)



############
# add comp summary
############
# all tissues
compdf = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/microbe-meta-analysis-comp-summary-table-consis_fc.tsv", sep='\t', header=T, check.names=F)
#head(compdf)

BIGDF1 = merge(BIGDF, compdf, by=1, all.x = T)
head(BIGDF1)



#############
# freq
##############
df2 = read.delim("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/freq-summary-w-tax.tsv", header=T, check.names=F)
rownames(df2) = as.vector(df2[,"ID"])
head(df2)

p1 = as.vector(BIGDF1[,"partner1"])

freqdf = df2[p1,]
head(freqdf)

outdff = cbind(BIGDF1, freqdf)

write.table(outdff, paste0(infile,"-DAM.txt"), quote=F, row.names=F, sep='\t')


print("Done!")
