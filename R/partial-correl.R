# install.packages("ppcor")
# library(ppcor)
# install.packages("corpcor")
library("corpcor")


args = commandArgs(trailingOnly=TRUE)


# args:
#1: this needs samples as columns and variables as rows. The script transposes it to variables as column and samples as rows.
#2: the variables (e.g. genes, microbes) to keep. provide it as one column file without header
#3: the samples to keep. provide it as coma separated string. missing values are not allowed so forcing to compile the list of samples
    ## rather than give a map file and column to pick
#4: output-file-prefix

# output
# this does not return pvalues


indf = read.csv(args[1], header=T, check.names=F, row.names=1)
#head(indf)

data = t(indf)
#head(data)
dim(data)



degs = as.vector(unlist(read.csv(args[2], header=F, check.names=F, as.is=T)))
head(degs)



# Missing values are not allowed.
Samps = strsplit(args[3], ",")[[1]]
Samps



###########
# upper triangle to 3 column dataframe
###########
UTritoDFNames = function(X){
  ind <- which(upper.tri(X, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(X)
  out = data.frame(row = nn[[1]][ind[, 1]],
           col = nn[[2]][ind[, 2]],
           val = X[ind])
  return(out)
}




dfSamps = data[Samps, degs]
head(dfSamps)
dim(dfSamps)

# from ppcor
#resSamps = pcor(dfSamps)

# from corpcor # https://cran.r-project.org/web/packages/corpcor/corpcor.pdf
# partial correlations (fast and recommend way)
# pcr1 = pcor.shrink(X)
# other possibilities to estimate partial correlations
# pcr2 = cor2pcor( cor.shrink(X) )



resSamps =  pcor.shrink(dfSamps)
coeffSamps = UTritoDFNames(resSamps)
head(coeffSamps)
write.table(coeffSamps, args[4], quote=F, row.names=F, sep='\t')

q()
