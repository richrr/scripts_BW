######Supplemental Material to
######State-of-the-art normalizations improve NMR-based metabolomics analysis
######Stefanie M. Kohl, Matthias S. Klein, Peter J. Oefner, Rainer Spang, Wolfram Gronwald
###### 11306_2011_350_MOESM2_ESM

######In the Following, the R code for the normalizations used in the paper is given.


######The non-normalized data matrix of feature intensities is stored in a matrix called "data"
######Each row represents a feature, each column a sample



# R.3.6.1 from linuxbrew has openblas as default which throws error during normalization. ERROR; return code from pthread_create() is 22
# the source of the error appears to be openblas >= 0.3.4. You can see in your R session history that you are using openblas 0.3.5.
# use the default R (3.6) from biowulf (/usr/local/apps/R/3.6/3.6.0/bin/R)

#usage:
#ml R
#/usr/local/apps/R/3.6/3.6.0/bin/Rscript /data/rodriguesrr/scripts/R/multiple-normalization-methods_bw.R --file no-singletons-otu_table.csv.rel.million.csv --symbolColumnName ID --q --output quantile/16s_rel_mill_ --logbase 2 --unlog
#module unload R


#=================================================================================
# checks and install the R library
# loads it
#=================================================================================
checkInstallAndLoadRlibs = function(lib, pathRlibs){
	if(! lib %in% rownames(installed.packages(lib.loc = pathRlibs))){
		install.packages(lib, lib=pathRlibs, repos='http://cran.us.r-project.org')
	}
	require(lib, , lib.loc = pathRlibs, character.only = TRUE)
}

if(FALSE){
libraries <- c("argparser", "data.table")
uname = Sys.getenv(c("USER"))
pathRlibs = paste0("/data/",uname, "/R_LIBS")
print(pathRlibs)
ifelse(!dir.exists(pathRlibs), dir.create(pathRlibs), FALSE)

for(lib in libraries)
{
	checkInstallAndLoadRlibs(lib, pathRlibs)
}


# r < v3.5
#install.packages(c("devtools"), lib=pathRlibs, repos='http://cran.us.r-project.org')
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("Biobase","preprocessCore"))



tryCatch({
	library(affy)
	#library(preprocessCore)
	}, error = function (e) {
	 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib="~/R/3.6/library", repos='http://cran.us.r-project.org') # use the above libpath if needed.
	BiocManager::install(c("affy", "preprocessCore"))
})

}


library(argparser)
library(affy)
library(data.table)


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Perform all or the requested normalization.')
p <- add_argument(p, "--file", help="file which has non-normalized data of feature intensities. Each row represents a feature, each column a sample", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--output", help="output file", default="./normalizations/res_")
p <- add_argument(p, "--symbolColumnName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID") # 'UniqueID'

# explicitly mention atleast one of the following
p <- add_argument(p, "--all", help="Run all Normalizations", flag=TRUE)
# OR
p <- add_argument(p, "--pq", help="Probabilistic Quotient Normalization", flag=TRUE)
p <- add_argument(p, "--cl", help="Cyclic Loess Normalization. It takes too long to run", flag=TRUE)
p <- add_argument(p, "--c", help="Contrast Normalization. Tends to give Error in sqrt(k)", flag=TRUE)
p <- add_argument(p, "--q", help="Quantile Normalization", flag=TRUE)
p <- add_argument(p, "--lb", help="Linear Baseline Normalization", flag=TRUE)
p <- add_argument(p, "--lw", help="Li-Wong Normalization", flag=TRUE)
p <- add_argument(p, "--cs", help="Cubic Spline Normalization", flag=TRUE)
p <- add_argument(p, "--a", help="Auto Scaling", flag=TRUE)
p <- add_argument(p, "--p", help="Pareto Scaling", flag=TRUE)
p <- add_argument(p, "--v", help="Variance Stabilization Normalization. It takes long time to run", flag=TRUE)
p <- add_argument(p, "--logbase", help="calc log (1+x) using the base", default=0) # allowed: 0 (no log), 1 (e), 2, 10
p <- add_argument(p, "--unlog", help="unlog the data then minus 1", flag=TRUE) # default no unlog




argv <- parse_args(p)
#print (p)

outputFile = argv$output
symbleColumnName =  argv$symbolColumnName
logbase = argv$logbase

# set all normalization flags to TRUE
if(argv$all){
    print("Requested ALL normalizations, however, not running loess and contrast normalization!")
    argv$pq = T
    argv$q = T
    argv$lb = T
    argv$lw = T
    argv$cs = T
    argv$a = T
    argv$p = T
    argv$v = T # takes long time to run but lesser than loess

    #argv$cl = T ## takes too long
    #argv$c = T ## tends to give Error in sqrt(k) : non-numeric argument to mathematical function in maffy.normalize

}


# create if outdir doesn't exist:
odir =  strsplit(outputFile, "/")[[1]]
odir = paste(odir[1:(length(odir) -1)], collapse='/')
print(odir)
dir.create(odir, recursive = TRUE)


# read the file and create a matrix
expressionData = read.csv( argv$file, header=TRUE, check.names=FALSE, na.strings=c("","na","NA", "Na", "NaN"), row.names= symbleColumnName)
data = data.matrix(expressionData)


if(logbase != 0){
    if(logbase == 1) {
          data = log(data + 1)  # using the default base e i.e. exp(1)
    } else {
          data = log(data + 1, logbase)
    }
}


pow <- function(x=10, y=6) {
   # function to print x raised to the power y
   result <- x^y
   return(result)
}


unlog_minus_1 <- function(InData, base){
	if(base != 0){
	    if(base == 1) {
	          InData = apply(InData, 2, function(x) {(exp(x) - 1)})    # using the base e
	    } else {
	          InData = apply(InData, 2, function(x) {(pow(base, x) - 1)})
	    }
	}
	return(InData)
}



# output to file
send_to_write = function(matrix, ofile, delim=','){
    # this keeps the header of row names empty
    #write.table(matrix, ofile, sep=delim, col.names = NA, row.names = TRUE)


    df = setDT(data.frame(matrix,check.names=F), keep.rownames = TRUE)[]
    setnames(df, 1, symbleColumnName)
    write.table(df, ofile, sep=delim, row.names = F)
}




#####---Quantile Normalization-----########################################
if(argv$q){
        print("Running Quantile")
	normalize.quantile <- get("normalize.quantiles", en=asNamespace("affy"))
	quantile.data <- normalize.quantile(data)
	#quantile.data <- normalize.quantiles(data) #preprocessCore

	# if unlog
	if(argv$unlog){
		quantile.data = unlog_minus_1(quantile.data, logbase)
	}

	# add row and column names
	quantile.data = data.frame(quantile.data)
	row.names(quantile.data)=row.names(data)
	colnames(quantile.data)=colnames(data)

	send_to_write(quantile.data, paste(outputFile, ".txt", sep='quantile'))
}



if(FALSE){

#####---Probabilistic Quotient Normalization-----##########################
if(argv$pq){
        print("Running PQN")
	reference <- apply(data,1,median)
	quotient <- data/reference
	quotient.median <- apply(quotient,2,median)
	pqn.data <- t(t(data)/quotient.median)
	send_to_write(pqn.data, paste(outputFile, ".txt", sep='pqn'))
}



#####---Linear Baseline Normalization-----#################################
if(argv$lb){
        print("Running Linear Baseline")
	linear.baseline <- apply(data,1,median)
	baseline.mean <- mean(linear.baseline)
	sample.means <- apply(data,2,mean)
	linear.scaling <- baseline.mean/sample.means
	linear.baseline.data <- t(t(data)*linear.scaling)
	send_to_write(linear.baseline.data, paste(outputFile, ".txt", sep='linear.baseline'))
}



#####---Li-Wong Normalization-----#########################################
if(argv$lw){
        print("Running Li-Wong")
	#---First step: Find baseline sample
	average.intensity <- apply(data,2,mean)
	median.number <- round(ncol(data)/2 + 0.1)				#R has an add way of rounding.
												#the additional 0.1 ensures that it rounds properly
	ordering <- order(average.intensity)
	median.sample.number <- ordering[median.number]
	median.sample <- data[,median.sample.number]

	#---Apply normalization
	liwong.data=vector()
	for(i in 1:ncol(data)){
		liwong.model <- normalize.invariantset(data=data[,i],
						ref=median.sample,
						prd.td=c(0.003,0.007))		#the threshold of the rank-invariant set might need to be adjusted from case to case
		liwong.sample <- predict(liwong.model$n.curve$fit,		#chosen such that the rank-invariant set it sufficiently large
						data[,i])
		liwong.data <- cbind(liwong.data,liwong.sample$y)
	}

	# add column names
	liwong.data = data.frame(liwong.data)
	colnames(liwong.data)=colnames(data)

	send_to_write(liwong.data, paste(outputFile, ".txt", sep='liwong'))
}



#####---Cubic Spline Normalization-----####################################
if(argv$cs){
        print("Running Spline")
	spline.data <- normalize.qspline(data,
					samples=0.02,
					target=apply(data,1,mean))

	# add row and column names
	spline.data = data.frame(spline.data)
	row.names(spline.data)=row.names(data)
	colnames(spline.data)=colnames(data)

	send_to_write(spline.data, paste(outputFile, ".txt", sep='spline'))
}



#####---Auto Scaling-----##################################################
if(argv$a){
        print("Running Auto")
	centered.data <- data - apply(data,1,mean)
	scaling.auto <- apply(data,1,sd)
	auto.data <- centered.data/scaling.auto
	send_to_write(auto.data, paste(outputFile, ".txt", sep='auto'))
}



#####---Pareto Scaling-----################################################
if(argv$p){
        print("Running Pareto")
	centered.data <- data - apply(data,1,mean)
	scaling.pareto <- sqrt(apply(data,1,sd))
	pareto.data <- centered.data/scaling.pareto
	send_to_write(pareto.data, paste(outputFile, ".txt", sep='pareto'))
}



#####---Variance Stabilization Normalization (VSN)-----####################
if(argv$v){
        print("Running VSN")
	library(vsn)									#load package unless it is already loaded
	vsn.model <- vsn2(data)
	vsn.data <- predict(vsn.model,data)
	send_to_write(vsn.data, paste(outputFile, ".txt", sep='vsn'))
}



#####---Cyclic Loess Normalization-----####################################
if(argv$cl){
        print("Running Loess")
	loess.data <- normalize.loess(data,
						subset=1:nrow(data),
						epsilon=10^-2,
						maxit=2,
						log.it=FALSE,
						verbose=TRUE,
						span=0.75,
						family.loess="gaussian")
	send_to_write(loess.data, paste(outputFile, ".txt", sep='loess'))
}



#####---Contrast Normalization-----########################################
if(argv$c){
        print("Running Contrast")
	#---First adaption: Make the data matrix non-negative
	smallvalue <- function(x, threshold=1e-11){				#threshold was chosen such that it is sufficiently smaller than the data
		for(i in 1:length(x)){
			if(!x[i]>0)	x[i] <- threshold
		}
	}

	nonnegative.data=smallvalue(data)

	#---Apply normalization
	maffy.data <- maffy.normalize(nonnegative.data,
						subset=1:nrow(nonnegative.data),
						span=0.75,
						verbose=TRUE,
						family="gaussian",
						log.it=FALSE)

	#---Second adaption: Subtract 10% Quantile from each sample
	subtract <- function(x){
		t(t(x)-apply(x,2,quantile,0.1))
	}

	contrast.data <- subtract(maffy.data)
	send_to_write(contrast.data, paste(outputFile, ".txt", sep='contrast'))
}

}


print("DONE")
