# use division w raw or log transf then substraction #

#cd Stephanie_Prescott/12KLMN/analys/pheno/comp/poolExpts/
#Rscript /data/rodriguesrr/scripts/R/meta-stats_ref_raw_div_vs_log_substr_bw.R --files /data/rodriguesrr/Stephanie_Prescott/06-08-2020/data/12KLMN-pheno-data.csv /data/rodriguesrr/Stephanie_Prescott/06-08-2020/data/12KLMN-pheno-data.csv /data/rodriguesrr/Stephanie_Prescott/06-08-2020/data/12KLMN-pheno-data.csv /data/rodriguesrr/Stephanie_Prescott/06-08-2020/data/12KLMN-pheno-data.csv --mapFiles 12K-pheno-map.txt 12L-pheno-map.txt 12M-pheno-map.txt 12N-pheno-map.txt --lists /data/rodriguesrr/Stephanie_Prescott/06-08-2020/data/12KLMN-pheno-list.txt --mapFileSampleIdColName ID --AnalysToDoList 12KLMN-pheno-analys_file_comp.xlsx --comparMethod tt --dataFileSymbolColName ID --output ./comp/tt_


#written for R/3.6.1. might work on 3.6.0


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

libraries <- c("argparser", "hash", "psych", "gtools", "corpcor", "reshape", "doParallel", "readxl", "tools", "stringr")
uname = Sys.getenv(c("USER"))
pathRlibs = paste0("/data/",uname, "/R_LIBS")
print(pathRlibs)
ifelse(!dir.exists(pathRlibs), dir.create(pathRlibs), FALSE)

for(lib in libraries)
{
	checkInstallAndLoadRlibs(lib, pathRlibs)
}


p <- arg_parser('Pool experiments (dataets) by normalizing the samples to reference and do one test')
p <- add_argument(p, "--files", help="data files", nargs=Inf) # required; but written as optional format so I explicitly mention the "--files"
p <- add_argument(p, "--mapFiles", help="metadata for the data files", nargs=Inf) # required; but written as optional format so I explicitly mention it
p <- add_argument(p, "--output", help="output file", default="./pooled_datasets_")

p <- add_argument(p, "--comparMethod", help="method to compare A vs B", default="tt") # allowed: "tt" (Welch t test), "mw" (Mann-Whitney U Test), limma

##### this script assumes for >=0 that the data used for analyses is always logged (either before running the script or inside the script)
p <- add_argument(p, "--logbase", help="calc log using the base", default=2) # allowed: 0 (no log. be careful since this may mean that the data was logged outside this script), 1 (e), 2, 10
# -1 would mean keep it raw (the data was NOT log transformed outside NOR should it be inside the script.)

p <- add_argument(p, "--dataFileSymbolColName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID")

# required; but written as optional format so I explicitly mention the following in cmd line
p <- add_argument(p, "--lists", help="lists (of measurement labels e.g. taxa names, phenotypics tests) to do analysis for", nargs=1)
p <- add_argument(p, "--AnalysToDoList", help="Analysis To Do List. 5 column excel file or tab-delimited file", default="AnalysToDoList.txt")
p <- add_argument(p, "--mapFileSampleIdColName", help="Column in map file in which to search the sample IDs.", default="#SampleID") # use the columns in the mapping file to identify the samples ids to perform analysis
p <- add_argument(p, "--pairedInfoColumn", help="Column in map file in which to search the ids for pairing samples. This info will help to identify, for example, samples from same mice", default="Mouse") # Column in map file in which to search the mouse ids. This info will help to identify samples from same mice, useful when doing paired comparisons




#==================================================================================================================
# get the supplied arguments
#==================================================================================================================
argv <- parse_args(p)
#print (p)


# check data and map
if(length(argv$files)!= length(argv$mapFiles))
{
  print("number of arguments for --files and --mapFiles should match.")
  quit()
}


if(length(argv$files) ==1 || length(argv$mapFiles) == 1  )
{
  print("At least 2 files are required. Give --files file1 file2 ... in cmd or just run perform-analyses_bw.R instead")
  quit()
}



outputFile = argv$output
#pvalThreshold = argv$pvalThreshold
#foldchThreshold = argv$foldchThreshold
#combinedPvalueCutoff = argv$combPvalCutoff
#combinedFDRCutoff = argv$combFDRCutoff

# create if outdir doesn't exist:

odir =  strsplit(outputFile, "/")[[1]]
odir = paste(odir[1:(length(odir) -1)], collapse='/')
print(odir)
ifelse(!dir.exists(odir), dir.create(odir, recursive = TRUE), "Outdir exists!")

listFile = argv$lists
list = read.csv(listFile,header = FALSE,sep="\t")

genes = as.vector(list[!duplicated(list[,1]),1] )	#delete duplicated array probes
#genes = genes[order(genes)] # ascending sort

logbase = argv$logbase

#foldchVar = 'FolChMedian'
#if(argv$foldchMean){foldchVar = "FoldChange"}
#outputFile = paste(outputFile , foldchVar, '_', sep='')



phead = function(indf){
  print(head(indf))
}

ptail = function(indf){
  print(tail(indf))
}

#----------------------
# -1 (raw), 0 (logged outside), 1 (log base e inside script), 2 (log base 2 inside script), 10 (log base 10 inside script)
#----------------------
logIt = function(expressionData){
	if(logbase > 0){
    if(logbase == 1) {
          expressionData = log(expressionData + 1)  # using the default base e i.e. exp(1)
    } else {
          expressionData = log(expressionData + 1, logbase)
    }
	}
	return(expressionData)
}

attach_map_info = function(d, m){
	data = t(read.csv(d, header = T, check.names=FALSE, na.strings=c("","na","NA", "Na", "NaN"), row.names=1))
	#print("####### Before LogIt ########")
	#phead(data)
	#print(paste0("####### Log base ", logbase, " ########"))
	data = logIt(data)
	#phead(data)


	ID = rownames(data)
	data = cbind(ID, data)
	#phead(data)

	map = read.delim(m,header = TRUE,sep="\t",check.names=FALSE) #, row.names=1)
	#phead(map)

	new = merge(map, data, by=1)
  #phead(new)

	#print(levels(new$Treatment))
	return(new)
}


df = list()
#------------------------------------------------------------------
# read the files one after the other:
#------------------------------------------------------------------
for( numb in 1:length(argv$files)){
	df[[numb]] = attach_map_info(argv$files[numb], argv$mapFiles[numb])
	phead(df[[numb]])
	write.csv(df[[numb]], paste0(outputFile, numb, "w-map.csv"), quote=F, row.names=F)
}


#------------------------------------------------------------------
# read analyses file:
#------------------------------------------------------------------
if(file_ext(argv$AnalysToDoList) == 'xls' ||  file_ext(argv$AnalysToDoList) == 'xlsx'){
	print("Excel version!")
	AnalysToDoList = read_excel(argv$AnalysToDoList) # tibble
	#print(AnalysToDoList)
	AnalysToDoList = as.data.frame(AnalysToDoList)
} else if (file_ext(argv$AnalysToDoList) == 'tsv' ||  file_ext(argv$AnalysToDoList) == 'txt'){
	print("TSV version!")
	AnalysToDoList = read.delim(argv$AnalysToDoList, header = T)
	AnalysToDoList[AnalysToDoList == ""] <- NA  # change empty strings to NA
} else {
	print ("Need a 5 column excel or tab-delimited file.")
	q()
}

# remove commented Analysis
AnalysToDoList = AnalysToDoList[!grepl("#", AnalysToDoList[,1]),]
#print(AnalysToDoList)



# allows multiple cores if the multi core flag is true
#allowed_cores = 1
#if(argv$multicore){
#    allowed_cores = argv$cores
#}

#cl <<- makeCluster(allowed_cores, type="FORK")
##registerDoParallel(cl)  # this function is used to register the parallel backend with the foreach package. not needed if not using foreach


#------------------------------------------------------------------
# allows uneven row cbind
# https://stackoverflow.com/questions/7962267/cbind-a-dataframe-with-an-empty-dataframe-cbind-fill
#------------------------------------------------------------------
cbind.fill <- function(...){
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x)
        rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


#------------------------------------------------------------------
# normalize to median of reference group
#------------------------------------------------------------------
normalize_to_categ2 = function(indf, g, ColHead, categ1, categ2){
		subdf_categ1 = indf[which(indf[, ColHead] == categ1), c(ColHead, g)]

		# reference
		subdf_categ2 = indf[which(indf[, ColHead] == categ2), c(ColHead, g)]
		median2 = median(as.numeric(as.character(subdf_categ2[, g])) , na.rm=T)
		print(paste(categ2, "median", median2))

		send = 0
		if(logbase >= 0){  # log was done either outside or inside the script
				send1 = as.numeric(as.character(subdf_categ1[, g])) - median2
				send2 = as.numeric(as.character(subdf_categ2[, g])) - median2
				send = cbind.fill(send1, send2)
		} else{ # if no log, i.e. raw
				send1 = as.numeric(as.character(subdf_categ1[, g])) / median2
				send2 = as.numeric(as.character(subdf_categ2[, g])) / median2
				send = cbind.fill(send1, send2)
		}
		colnames(send) = c(categ1, categ2)
		return(send)

}


#------------------------------------------------------------------
# loop over all datasets for a single gene
#------------------------------------------------------------------
PoolExpts = function(g, df, ColHead, categ1, categ2, comparMethod, indx){

		BIGDF = NULL
		for(sdf in df){
				keep = sdf[which(sdf[,ColHead] == categ1 | sdf[,ColHead] == categ2), c(ColHead, g)]
				print(keep)

				# get median per dataset, and rbind to the next df
				tmp_df = normalize_to_categ2(keep, g, ColHead, categ1, categ2)

				if(is.null(dim(BIGDF))){
						BIGDF = tmp_df
				} else{
						BIGDF = rbind(BIGDF, tmp_df)
				}
		}
		print(BIGDF)
		write.csv(BIGDF, paste(outputFile, g, ColHead, categ1, "vs", categ2, "pooled.csv", sep="_"), quote=F, row.names=F)
		return(BIGDF)
}


#------------------------------------------------------------------
#http://stackoverflow.com/questions/9455167/converting-data-frame-factors-with-data-matrix-weird-results
#else as.numeric gives weird numbers instead of actual values
#------------------------------------------------------------------
numericizeVector = function(c){
    nc = sapply(c, function(x) as.numeric(as.character(x)))
    return(nc)
}


#------------------------------------------------------------------
# t test with default: two-sided, paired = FALSE, var.equal = FALSE
## use equal variance
#------------------------------------------------------------------
run.comp.tests = function(ldf, comparMethod){

	c1 = as.vector(ldf[,1])
	c2 = as.vector(ldf[,2])

	c1 =  numericizeVector(c1)
	#print(c1)
	c2 = numericizeVector(c2)
	#print(c2)

	outLine = ''
	p = ''

	mean_c1 = NA
	mean_c2 = NA

	median_c1 = NA
	median_c2 = NA

	# there are at least 3 samples with values (i.e. without NA)
	if(sum(!is.na(c1)) >= 3){
			mean_c1 = mean(c1, na.rm=TRUE)
			median_c1 = median(c1, na.rm=TRUE)
	}

	if(sum(!is.na(c2)) >= 3){
			mean_c2 = mean(c2, na.rm=TRUE)
			median_c2 = median(c2, na.rm=TRUE)
	}


	fold_change_mean = NA
	fold_change_median = NA

	if(sum(!is.na(c1)) < 3 || sum(!is.na(c2)) < 3){
			outLine = as.matrix(c(NA, mean_c1, mean_c2, fold_change_mean, median_c1 , median_c2 , fold_change_median, NA))
	} else{  #if(is.na(pairedd) || pairedd == "")

			if(logbase >= 0) { # some form of log transform
					fold_change_mean = mean_c1-mean_c2
					fold_change_median = median_c1-median_c2
			 } else { # not log transformed
			  	fold_change_mean = mean_c1/mean_c2
			  	fold_change_median = median_c1/median_c2
			 }

			if(comparMethod == 'tt'){
					# check if it works or throws an error (e.g. "data are essentially constant")
					obj<-try(t.test(c1,c2), silent=TRUE)
					if (is(obj, "try-error")) {
							# when error
							outLine = as.matrix(c(NA, mean_c1, mean_c2, fold_change_mean, median_c1 , median_c2 , fold_change_median, NA))
					} else {
							p = t.test(c1,c2)
							outLine = as.matrix(c(p$method, mean_c1, mean_c2, fold_change_mean, median_c1 , median_c2 , fold_change_median, p$p.value))
					}
			} else if(comparMethod == 'mw'){
					p = wilcox.test(c1,c2)
					outLine = as.matrix(c(p$method, mean_c1, mean_c2, fold_change_mean, median_c1 , median_c2 , fold_change_median, p$p.value)) # since this test does not explicitly give the means in the "estimate"
			}
	}
	return(outLine)

}


#------------------------------------------------------------------
# for each gene, compare groups after pooling datasets
#------------------------------------------------------------------
calculatePoolComparison = function(df, ColHead, categ1, categ2, comparMethod, indx){
	print(ColHead)
	print(categ1)
	print(categ2)
	print(comparMethod)
	print(indx)

	all_gene_res = ''
	for(g in genes){
			per_g_pooled_df = PoolExpts(g, df, ColHead, categ1, categ2, comparMethod, indx)
			per_g_res = t(run.comp.tests(per_g_pooled_df, comparMethod))
			#print(per_g_res)

			if(is.null(dim(all_gene_res))){
					all_gene_res = per_g_res
			} else{
					all_gene_res = rbind(all_gene_res, per_g_res)
			}
	}
	# drop the method column
	all_gene_res = all_gene_res[,-1]
	colnames(all_gene_res) = c(
								paste0(ColHead, "||", categ1, " Mean"), paste0(ColHead, "||", categ2, " Mean"),
								paste0(ColHead, "||", categ1, "_vs_", ColHead, "||", categ2, " FolChMean"),
								paste0(ColHead, "||", categ1, " Median"), paste0(ColHead, "||", categ2, " Median"),
								paste0(ColHead, "||", categ1, "_vs_", ColHead, "||", categ2, " FolChMedian"),
								paste0(ColHead, "||", categ1, "_vs_", ColHead, "||", categ2, " ", comparMethod, " p-value")
	 			)
	colnames(all_gene_res) = paste(	paste0("Analys ", indx), colnames(all_gene_res))
	ID = genes
	all_gene_res = cbind(ID, all_gene_res)
	print(head(all_gene_res))
	return(all_gene_res)
}


COLHEADGRPSEP = '_||_'
n=0
compout = c()

#==================================================================================================================
# perform analysis
#==================================================================================================================
for(indx in 1:nrow(AnalysToDoList)){
    cfull = AnalysToDoList[indx,]
		#print(cfull)
		ColHead = cfull[1]
		ColHead = toString(ColHead[1,1])
		#print(ColHead)

		c = cfull[-1]
		print(c)

    if(c[1] != "" && c[2] != "" && is.na(c[3]) && is.na(c[4])){
				print("1")
				categ1= toString(c[1,1])
				categ2= toString(c[1,2])

				#print(ColHead)
				#print (categ1)
				#print (categ2)

				# calculate comparison
        outEachExperiment = calculatePoolComparison(df, ColHead, categ1, categ2, argv$comparMethod, indx)
				#phead(outEachExperiment)

        if (length(compout)==0){
					compout <<- outEachExperiment
        }else{
					compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
        } # end else

    }

		write.csv(compout, paste(outputFile, "comp.csv", sep="_"), quote=F, row.names=F)
}

q()
