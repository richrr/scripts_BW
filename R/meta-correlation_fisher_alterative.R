
## from nolan jun 17 2020.
#The website that I used to check my calculated pvalues is here: http://www.meta-mar.com/corr
#I have included the input that I used for the code, the output for the code, and the code itself.
#I use the "metacor" function from the "meta" package #for the new pvalues (the alternatives to Fisher).
#Here is the documentation for that package https://cran.r-project.org/web/packages/meta/meta.pdf
##

# Perform meta-correlation of multiple rho correlation values. Takes as input a csv
# with the names of studies/treatment groups as column headers, the FIRST row
# (excluding the header) contains the sample size in each group, and all following
# rows with correlations and their corresponding rho values

# EXAMPLE
#                                                                study1_rho           study2_rho         study3_rho           study4_rho
# sample_size                                                     9                    9                  7                    8
# ASV101<==>ASV1                                                  0.038635506         -0.512450039       -0.142857143         -0.730552020
# ASV105<==>ASV1                                                 -0.352734816         -0.031167254        0.535714286          0.428571429
# ASV105<==>ASV101                                                0.065552829          0.316771269       -0.392857143          0.107786364

rm(list = ls())
library(meta)

###########################################
# For testing purposes
#
# df <- data.frame(group1=as.numeric(),
#                  group2=as.numeric(),
#                  group3=as.numeric(),
#                  group4=as.numeric())
# df[1,] = c(0.781512605042017,1,0.214285714285714,0.214285714285714)
# df[2,] = c(12,10,12,12)
# row.names(df) = c("cor", "ss")
#
# newdf = data.frame(t(df))
#
# metacor(cor = newdf$cor, n =newdf$ss, studlab = row.names(newdf), sm = "ZCOR")
#
# result$pval.fixed
###########################################

data = read.csv("only_corrs_in_net_w_samp_size.csv", stringsAsFactors = F, header = T, row.names = 1, check.names = F)

data = data.frame(t(data), check.names = F)

pvals = list()


for (i in 2:length(data)){
  print(colnames(data)[i])
  print(paste(i, " out of ", length(data) - 1, " pvals have been calculated"))
  res = metacor(cor = data[,i], n = data[,1], studlab = row.names(data), sm = "ZCOR")
  print(res$pval.fixed)
  pvals[colnames(data)[i]] = res$pval.fixed
}

pvaldf = t(data.frame(pvals, check.names = F))
colnames(pvaldf) = "meta_corr_pval"

data_transp = t(data)

output = merge(data_transp, pvaldf, by = "row.names")
rownames(output) = output[,1]
output = output[,-1]

write.csv(output, "correlations_w_meta_corr_pvalue.csv")
