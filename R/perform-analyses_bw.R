##### to do:
## partial correlation
## parallelize method:
##   - analysis is sequential. parallelization currently only implemented for correlation (since those calcs take the longest) calc.
##   - next, we can implement for parallelizing different analysis or comparisons or other correlations

## write out the per analysis results to tmp folder in lscracth to clear buffer and reduce the size of the big (>10M data frame). test the speed w and w/o writing. allow writing option based on the flag or prod(number of correlations and corr analysis).


# usage: Rscript /data/rodriguesrr/scripts/R/perform-analyses_bw.R --help
# cd /data/rodriguesrr/TransNetTest
# Rscript /data/rodriguesrr/scripts/R/perform-analyses_bw.R ./data/abundanceTables/rnaseq-pool_RQ.limf.pheno.lipid.csv --lists ./data/lists/liver-fdr15-subset.txt ./data/lists/liver-fdr15-subset.txt --mapFile data/map/rnaseq-pool-diet.map.txt --mapFileSampleIdColName SampleID --AnalysToDoList data/analys_files/t2-diet-comp-corr.txt --comparMethod mw --correlMethod sp --dataFileSymbolColName ID --pairedInfoColumn TestPair

# Rscript /data/rodriguesrr/scripts/R/perform-analyses_bw.R ./data/abundanceTables/rnaseq-pool_RQ.limf.pheno.lipid.csv --lists ./data/lists/liver-fdr15.txt ./data/lists/liver-fdr15.txt --mapFile data/map/rnaseq-pool-diet.map.txt --mapFileSampleIdColName SampleID --AnalysToDoList data/analys_files/t2-diet-comp-corr.txt --comparMethod mw --correlMethod sp --dataFileSymbolColName ID --pairedInfoColumn TestPair

# Rscript /data/rodriguesrr/scripts/R/perform-analyses_bw.R ./data/abundanceTables/rnaseq-pool_RQ.limf.pheno.lipid.csv --lists ./data/lists/liver-fdr15.txt ./data/lists/liver-fdr15.txt --mapFile data/map/rnaseq-pool-diet.map.txt --mapFileSampleIdColName SampleID --AnalysToDoList data/analys_files/t2-diet-comp-corr.xlsx --comparMethod mw --correlMethod sp --dataFileSymbolColName ID --pairedInfoColumn TestPair



# Allowed analyses (conditions):
# add 4th column with "correlation" to indicate whatever operation is to be done, needs to be done on correlations
# 1. A vs B: 2 non-empty columns, last two columns empty [A	B empty empty]
# 2. delta A vs delta B: 3 non-empty columns, 3rd column delta, 4th column empty [A	B	delta	empty]  # needs 4 groups: two comma separated from A and two comma separated from B
# 3. A vs B paired:  3 non-empty columns, 3rd column paired, 4th column empty [A	B	paired	empty]
# 4. Correl A: 2 non-empty columns, 4th column correlation [A	empty	empty	correlation]
# 5. Correl delta A: 3 non-empty columns, 2nd column empty, 3rd column delta, 4th column correlation [A	empty	delta	correlation] # needs 2 groups: two comma separated from A
# 6. Correl A vs. Correl B: 3 non-empty columns, 4th column correlation [A	B	empty	correlation]
# 7. Correl delta A vs. Correl delta B: 4 non-empty columns, 3rd column delta, 4th column correlation [A	B	delta	correlation] # needs 4 groups: two comma separated from A and two comma separated from B
# 8. Correl A vs. Correl B paired: 4 non-empty columns, 3rd column paired, 4th column correlation [A	B	paired	correlation]
## no need to do paired in 8
# 9. Timelag Correl A and B (where A and B are different time points): 4 non-empty columns, 3rd column timelag, 4th column correlation [A	B	timelag	correlation]
## std. correlations: g1T1-g2T1 or g1T2-g2T2
## calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc.
## for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene

#written for R/3.6


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



# load the functions from a different file to keep code neater
source("/data/rodriguesrr/scripts/R/comp-correl-delta-paired-util-functs_bw.R")
#source("/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/get_base_stats_low_express_probes.R")



#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('perform-analyses_bw.R Calculates comparison and correlation analysis between lists')
p <- add_argument(p, "expressionDataFile", help="csv file containing values (measurements) to be used for analyses") # required

# truly optional
p <- add_argument(p, "--output", help="output folder(s) and file prefix", default="./perf.analys/result_")
p <- add_argument(p, "--dataFileSymbolColName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID")

# required; but written as optional format so I explicitly mention the following in cmd line
p <- add_argument(p, "--lists", help="lists (of measurement labels e.g. taxa names, phenotypics tests) to do analysis for", nargs=2)
p <- add_argument(p, "--AnalysToDoList", help="Analysis To Do List. 5 column excel file or tab-delimited file", default="AnalysToDoList.txt")
p <- add_argument(p, "--mapFile", help="tab delimited mapFile", default="mapFile.txt") # use the mapping file to identify the samples group affiliation to perform analysis


#p <- add_argument(p, "--mapColumns", nargs=2, help="Columns in map file in which to search the sample group associations", default=c("#SampleID","Experiment")) # use the columns in the mapping file to identify the samples group affiliation to perform analysis
# this changed to taking only one column, the other column will be in the analysis to do list
p <- add_argument(p, "--mapFileSampleIdColName", help="Column in map file in which to search the sample IDs.", default="#SampleID") # use the columns in the mapping file to identify the samples ids to perform analysis


p <- add_argument(p, "--pairedInfoColumn", help="Column in map file in which to search the ids for pairing samples. This info will help to identify, for example, samples from same mice", default="Mouse") # Column in map file in which to search the mouse ids. This info will help to identify samples from same mice, useful when doing paired comparisons

# optional methods
p <- add_argument(p, "--comparMethod", help="method to compare A vs B", default="limma") # allowed: "tt" (Welch t test), "mw" (Mann-Whitney U Test), limma
p <- add_argument(p, "--correlMethod", help="method to correlate A vs B", default="pearson") # allowed: pearson, "kendall", or "spearman", can be abbreviated.
p <- add_argument(p, "--logbase", help="calc log using the base", default=2) # allowed: 0 (no log. be careful since this may mean that the data was logged outside this script), 1 (e), 2, 10


p <- add_argument(p, "--background", help="the background level usually obtained from the raw files. this is used to find low expression probes", default=0) # this needs to be removed and given as a flag


p <- add_argument(p, "--stringtorelativize", help="string that indicates the entries to be relativized. Anything is allowed; however, e.g., use OTU, SEED, ENSG, ENSMUSG, etc.", default="OTU") # see ensembl prefixs here https://www.uea.ac.uk/documents/429378/432070/Ensemble%2BPresentation.pdf/b5a33c43-1c4c-4677-b588-1a7cac8db622


# flags
p <- add_argument(p, "--multicore", help="run in multicore/multiprocessors", flag=TRUE)
p <- add_argument(p, "--cores", help="use these many cores", default=8, type="numeric") # do not use more than 8. overhead of parallelization gives little speed up above 8 cores.
p <- add_argument(p, "--warnings", help="print warnings", flag=TRUE)
p <- add_argument(p, "--noCorrelationsRequested", help="no correlations calc. requested. Use this for preliminary analysis or if you have a huge (e.g. >10K input elements as rows)", flag=TRUE)
p <- add_argument(p, "--skipDiffInputsCheck", help="does not quit if diff input check fails", flag=TRUE) # avoid using this
p <- add_argument(p, "--writeCorrPerAnalysis", help="write files per analysis for correlations when requested. This will be set when the number of correlations per group is greater than 5 million", flag=TRUE)




#==================================================================================================================
# get the supplied arguments
#==================================================================================================================
argv <- parse_args(p)
#print (p)

expressionDataFile = argv$expressionDataFile

# check lists
if(length(argv$lists)!= 2)
{
  print("2 arguments for Lists are required. Give --lists list1 list2 in cmd")
  quit()
}
list1File = argv$lists[1]
list2File = argv$lists[2]

logbase = argv$logbase
backgrond = argv$background
outputFile = argv$output
symbleColumnName =  argv$dataFileSymbolColName

# create if outdir doesn't exist:
odir =  strsplit(outputFile, "/")[[1]]
odir = paste(odir[1:(length(odir) -1)], collapse='/')
print(odir)
ifelse(!dir.exists(odir), dir.create(odir, recursive = TRUE), "Outdir exists!")


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



mapColumnHeaders = unique(AnalysToDoList[,1])


# the samples in samplIdCol belong to the groups in exptCol
#samplIdCol = argv$mapColumns[1]
#exptCol = argv$mapColumns[2]

# now changed to only ask for column of sample ids in map file
samplIdCol = argv$mapFileSampleIdColName


# allows multiple cores if the multi core flag is true
allowed_cores = 1
if(argv$multicore){
    allowed_cores = argv$cores
}

cl <<- makeCluster(allowed_cores, type="FORK")
#registerDoParallel(cl)


RunFreqAnalys = 0

#==================================================================================================================
# create a group_sampleid dictionary
#==================================================================================================================
mapFile=read.delim(argv$mapFile, header = TRUE, sep="\t", check.names=FALSE)



#======================================
if(FALSE){ # list (exptcol) of list (group k) of vector (samples)
	MAPLIST = list()
	for(exptCol in mapColumnHeaders){

		tmpList = list()
		for(k in levels(mapFile[, exptCol])){
			vals = as.vector(mapFile[which(mapFile[,exptCol]==k), samplIdCol])
			tmpList[[k]] = vals
		}

		MAPLIST[[exptCol]] = tmpList
	}
	#print(MAPLIST)

	for(k in names(MAPLIST)){
	    #print(k)

			for(v in MAPLIST[k]){
					#print(v) # prints the names and elements

					for (n in names(v)){
						print(n)
						print(v[[n]])

						str_r_l = paste(k, n, length(v[[n]]), "values", paste(as.vector(v[[n]]), collapse=", "), sep="---")
						print(str_r_l)
					}
			}
	}
}
#======================================

COLHEADGRPSEP = '_||_'  #__<=-=>__
# key is exptcol + group k
# value is the vector of samples
dict = hash()
for(exptCol in mapColumnHeaders){
	# for unique items in the experiment column
	for(k in levels(mapFile[, exptCol])){
		# extract the rows whose col value in exptCol matches k
		# assign the sample names (i.e. values in samplIdCol) to key k
		vals = as.vector(mapFile[which(mapFile[,exptCol]==k), samplIdCol])
		expK = paste(exptCol, k, sep=COLHEADGRPSEP)
		dict[expK] = vals
	}
}
#print(dict)



#==================================================================================================================
# check the number and order of samples assigned to groups in dict. MAY be as per the order in which the samples are added (order in the mapping file)
#==================================================================================================================
for(k in keys(dict)){
    str_r_l = paste(k, length(dict[[k]]), "values", paste(as.vector(dict[[k]]), collapse=", "), sep="---")
    print(str_r_l)
}



#==================================================================================================================
# check if the order of samples in the two categories is correct; affects any time there is delta, paired, timelag
#==================================================================================================================
checkPairing = function(catg1, catg2, dict, mapFile, pairedInfoColumn, samplIdCol){
		s1 = as.vector(dict[[catg1]])
		s2 = as.vector(dict[[catg2]])

		tmp_mapFile = as.matrix(mapFile)
		row.names(tmp_mapFile) = mapFile[,samplIdCol]
		#print(tmp_mapFile)
		m1 = as.vector(tmp_mapFile[s1, pairedInfoColumn])
		print(m1)
		m2 = as.vector(tmp_mapFile[s2, pairedInfoColumn])
		print(m2)

		bool = all.equal(m1,m2)

		if(bool == "TRUE"){
		    print("Passed pairing check")
		}else{
		    print("XXXXXXXXXXXXXX----Failed pairing check----XXXXXXXXXXXXXX")
		    print(bool)
		    quit()
		}
		return(1)
}



for(indx in 1:nrow(AnalysToDoList)){
    cfull = AnalysToDoList[indx,]

		ColHead = cfull[1]
		ColHead = toString(ColHead[1,1])
		c = cfull[-1]

		#print (cfull)
		#print(ColHead)
		#print(c)

		if( !is.na(c[3]) && (c[3] == "paired" || c[3] == "delta" || c[3] == "timelag")){
			print(c)

      if(c[1] != "" && c[2] != "" && c[3] == "delta" && is.na(c[4]) ){
        print("2")

        categA= toString(c[1,1])
        categB= toString(c[1,2])

				cA = strsplit(categA, ",")[[1]]
        categ1 = paste( ColHead, cA[1] , sep=COLHEADGRPSEP)
        categ2 = paste( ColHead, cA[2] , sep=COLHEADGRPSEP)
				checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

				cB = strsplit(categB, ",")[[1]]
        categ3 = paste( ColHead, cB[1] , sep=COLHEADGRPSEP)
        categ4 = paste( ColHead, cB[2] , sep=COLHEADGRPSEP)
        checkPairing(categ3, categ4, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

				next
      } # end delta if


      if(c[1] != "" && c[2] != "" && c[3] == "paired" && is.na(c[4])){
        print("3")

				categ1= paste( ColHead, toString(c[1,1]) , sep=COLHEADGRPSEP)
        categ2= paste( ColHead, toString(c[1,2]) , sep=COLHEADGRPSEP)
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

				next
      } # end paired if


      if(c[1] != "" && is.na(c[2]) && c[3] == "delta" && c[4] == "correlation"){
        print("5")

        categA= toString(c[1,1])
        cA = strsplit(categA, ",")[[1]]
        categ1 = paste( ColHead, cA[1] , sep=COLHEADGRPSEP)
        categ2 = paste( ColHead, cA[2] , sep=COLHEADGRPSEP)
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

				next
      } # end delta correl if


      if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == "correlation"){
        print("7")

        categA= toString(c[1,1])
        categB= toString(c[1,2])

        cA = strsplit(categA, ",")[[1]]
        categ1 = paste( ColHead, cA[1] , sep=COLHEADGRPSEP)
        categ2 = paste( ColHead, cA[2] , sep=COLHEADGRPSEP)
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

        cB = strsplit(categB, ",")[[1]]
				categ3 = paste( ColHead, cB[1] , sep=COLHEADGRPSEP)
        categ4 = paste( ColHead, cB[2] , sep=COLHEADGRPSEP)
        checkPairing(categ3, categ4, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

				next
      } # end delta correl comparison if


      if(c[1] != "" && c[2] != "" && c[3] == "timelag" && c[4] == "correlation"){
        print("9")

        # std. correlations: g1T1-g2T1 or g1T2-g2T2
        # calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc.
        # for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene
        categ1= paste( ColHead, toString(c[1,1]) , sep=COLHEADGRPSEP)
        categ2= paste( ColHead, toString(c[1,2]) , sep=COLHEADGRPSEP)
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

				next
      } # end timelag if


    } # end outer if

}# end for



#	output filename
pairFile = "./testpair.txt"
#==================================================================================================================
# generate or load existing pairlist to avoid re-generating pairs(needs time) (for correlations and comparison of correlations)
#==================================================================================================================
list1 = read.csv(list1File,header = FALSE,sep="\t")
list2 = read.csv(list2File,header = FALSE,sep="\t")

genes1 = as.vector(list1[!duplicated(list1[,1]),1] )	#delete duplicated array probes
genes1 = genes1[order(genes1)] # ascending sort
genes2 = as.vector(list2[!duplicated(list2[,1]),1] )
genes2 = genes2[order(genes2)]

#print(genes1)
#print(genes2)


LargeNumbCorrsForCalc = FALSE
pairs = ''
if(!argv$noCorrelationsRequested){
  if (!file.exists(pairFile)){
	#if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then generate unique pairs
	if(identical(genes1, genes2)){
		pairs = t(combn(genes1,2))[,2:1]
			# use below for testing purposes and/or when calculating partial correlation # has same gene pairs (e.g. gene1 gene1)
			#pairs = combinations(length(unique(genes1)), 2, genes1, repeats.allowed=TRUE) # from https://gist.github.com/randy3k/10015496 for testing purposes
	}else{
		pairs = expand.grid(genes1,genes2)  # does genes1 x genes2 pairs ## not sure about this-> this has all the above & their opposite (e.g. gene1 gene2 and gene2 gene1) & same gene pairs (e.g. gene1 gene1)
	}
	#write.csv(pairs,pairFile,row.names=FALSE)
	#file.remove("./testpair.txt")
  }else{
	pairs = read.csv(pairFile)
  }

	if(nrow(pairs) > (5*(10^6))){
		LargeNumbCorrsForCalc = TRUE
	}
}


if(! argv$writeCorrPerAnalysis){
	argv$writeCorrPerAnalysis = FALSE
}

# if requested or if the number of correlations is > 5million; write corr per analysis.
if(argv$writeCorrPerAnalysis || LargeNumbCorrsForCalc){
	print("Writing corr results per analysis.")
	argv$writeCorrPerAnalysis = TRUE
}



#==================================================================================================================
# load expression data
#==================================================================================================================
expressionData = read.csv(expressionDataFile,header = TRUE,check.names=FALSE,na.strings=c("","na","NA", "Na", "NaN"))
expressionData = expressionData[!duplicated(expressionData[,symbleColumnName]),]    #delete duplicated measurements (genes, taxa labels, phenotypic labels, etc.)
rownames(expressionData) = expressionData[,symbleColumnName]

if(logbase != 0){
    if(logbase == 1) {
          expressionData[, !names(expressionData) %in% c(symbleColumnName)] = log(expressionData[, !names(expressionData) %in% c(symbleColumnName)] + 1)  # using the default base e i.e. exp(1)
    } else {
          expressionData[, !names(expressionData) %in% c(symbleColumnName)] = log(expressionData[, !names(expressionData) %in% c(symbleColumnName)] + 1, logbase)
    }
}

write.csv(expressionData,paste(c(outputFile,"log", logbase, "transform-output.csv"),collapse='-'),row.names=FALSE)


#==================================================================================================================
# get some basic statistics about the entire expression data
#==================================================================================================================
#get_base_stats_exp_data(expressionDataFile, backgrond, outputFile, sample_perc_threshold=10)




# which elements from list1 are not present in list2
checkDiffInputs = function(list1, list2, search, inn){
    diff = setdiff(list1, list2)
    diff = diff[!diff %in% c('DUMMY')]  # remove known tmp values
    if(length(diff) > 0){
           strr = paste("\n\nThe following ids in the", search ,"file are not present in the", inn, "file\n", sep=' ')
           cat(strr)
           print(diff)
           if(!argv$skipDiffInputsCheck){
               quit()
           }
    }
    #return(diff)
}


#==================================================================================================================
# check if all ids from the different files are consistent. This is especially to check typos or different cases ("ABX" vs "Abx")
#==================================================================================================================
# partners file and expt file to check gene name
checkDiffInputs(genes1, as.vector(expressionData[,symbleColumnName]), "partners", "expression")
checkDiffInputs(genes2, as.vector(expressionData[,symbleColumnName]), "partners", "expression")

# mapping file and expt file to check sample name
checkDiffInputs(as.vector(mapFile[, samplIdCol]), as.vector(colnames(expressionData)), "mapping", "expression")


# analysis file and mapping file to check group names
tmp_a1 = sort(as.vector(AnalysToDoList[,2]))
tmp_a2 = sort(as.vector(AnalysToDoList[,3]))
tmp_a12 = unique(sort(append(tmp_a1, tmp_a2)))

# elements contain ','
csv_tmp_a12 = grep("," , tmp_a12, value = TRUE)
# split based on ','
no_csv_tmp_a12 = unlist(strsplit(csv_tmp_a12, ','))
# remove elements with ',' and replace with the vector obtained by splitting on ','
uniq_analys_elems = append(setdiff( tmp_a12, csv_tmp_a12 ), no_csv_tmp_a12)

# elements contains ';' in case of anova like analysis
semcol_tmp_a12 = grep(";" , uniq_analys_elems, value = TRUE)
# split based on ';'
no_semcol_tmp_a12 = unlist(strsplit(semcol_tmp_a12, ';'))

# the above may contain elements with "_vs_"
vs_no_semcol_tmp_a12 = grep("_vs_" , no_semcol_tmp_a12, value = TRUE)
# split based on '_vs_'
no_vs_no_semcol_tmp_a12 = unlist(strsplit(vs_no_semcol_tmp_a12, '_vs_'))
# remove elements with '_vs_' and replace with the vector obtained by splitting on '_vs_'
uniq_analys_elems_tmp = append(setdiff( no_semcol_tmp_a12, vs_no_semcol_tmp_a12 ), no_vs_no_semcol_tmp_a12)


# remove elements with ';' and replace with the vector obtained by splitting on ';'
uniq_analys_elems = append(setdiff( uniq_analys_elems, semcol_tmp_a12 ), uniq_analys_elems_tmp)
#print(uniq_analys_elems)

CHeads = unique(sort(as.vector(AnalysToDoList[,1])))
CHeadsSampsVec = unique(as.vector(unlist(mapFile[, CHeads], use.names=FALSE)))
checkDiffInputs(unique(uniq_analys_elems[uniq_analys_elems != ""]), CHeadsSampsVec, "analysis", "mapping")



calcFC = function(v1 , v2){
    v1 = as.numeric(as.vector(v1))
    v2 = as.numeric(as.vector(v2))

    FC = c()
    for(i in 1:length(v1)){
        if(is.na(v1[i]) || is.na(v2[i])){
            FC = c(FC, NA)
        }else{
            fc = as.numeric(v1[i])/as.numeric(v2[i])
            FC = c(FC, fc)
        }
    }

    return(FC)
}



# write corr output to file per analysis.
write_to_File_CorrPerAnalysis = function(corroutPerA, idx){
	write.csv(corroutPerA,paste0(outputFile,'Analys-',idx,"-corr-output.csv"),row.names=FALSE)
}


n=0 # total number of analysis # this may not be accurate anymore since I stopped these counts for correl calc in parallel mode to avoid wasting time in writing to std out

corrout = c()
compout = c()
#folchout = c()



#==================================================================================================================
# perform analysis
#==================================================================================================================
for(indx in 1:nrow(AnalysToDoList)){
    #c = as.vector(AnalysToDoList[indx,])
		cfull = AnalysToDoList[indx,]

		ColHead = cfull[1]
		ColHead = toString(ColHead[1,1])
		c = cfull[-1]

    #print(c)
    if(c[1] != "" && c[2] != "" && is.na(c[3]) && is.na(c[4])){
				print("1")
				categ1= toString(c[1,1])
				categ2= toString(c[1,2])

        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
					genes = genes1
				}else{                              # else take union
					genes = union(genes1,genes2)
				} # end else

				#print(expressionData)
				#print(ColHead)
				#print (categ1)
				#print (categ2)

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)

				# calculate comparison
        outEachExperiment = calculateComparison(genes, expressionData, categ1, categ2, dict, argv$comparMethod, c[3], indx)
				#print(head(outEachExperiment))

        if (length(compout)==0){
					compout <<- outEachExperiment
        }else{
					compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
        } # end else

        # calculate fold change
        #FolCh =  calcFC(outEachExperiment[,5],outEachExperiment[,6]) # median to calc FC
        #FolCh = t(rbind(outEachExperiment[,1], t(FolCh)))
        #colnames(FolCh) = c(colnames(outEachExperiment)[1],paste(categ1,"vs",categ2,"FoldChange",sep=" "))

        #if (length(folchout)==0){
					#folchout <<- FolCh
        #}else{
					#folchout <<- merge(folchout,FolCh,sort=FALSE,by=1)
        #} # end else

				RunFreqAnalys = 1
				next

    } else if(c[1] != "" && c[2] != "" && c[3] == "delta" && is.na(c[4])){
				print("2")
				categA= toString(c[1,1])
				categB= toString(c[1,2])

				cA = strsplit(categA, ",")[[1]]
				categ1 = cA[1]
				categ2 = cA[2]

				cB = strsplit(categB, ",")[[1]]
				categ3 = cB[1]
				categ4 = cB[2]

				if(categ1 == "" || categ2 == "" || categ3 == "" || categ4 == ""){
				    print("Error in number of categories. Give categ1,categ2\tcateg3,categ4")
				    next
				}

				genes = ''
				# if gene list1 is the same as gene list2 then generate unique pairs
				if(identical(genes1, genes2)){  # sum(genes1 != genes2)==0
						genes = genes1
				}else{                              # else take union
	    			#print("Reached here")
	    			genes = union(genes1,genes2)
				} # end else

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)
				categ3 = paste(ColHead, categ3, sep=COLHEADGRPSEP)
				categ4 = paste(ColHead, categ4, sep=COLHEADGRPSEP)

				# calculate comparison
				outEachExperiment = calculateComparisonDelta(genes, expressionData, categ1, categ2, categ3, categ4, dict, argv$comparMethod, indx)
				#print(head(outEachExperiment))

				if (length(compout)==0){
						compout <<- outEachExperiment
				}else{
						compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
				} # end else

        # calculate fold change
        #FolCh = outEachExperiment[,c(1,7)]
        #colnames(FolCh) = c(colnames(outEachExperiment)[1],paste(categA,"vs",categB,"FoldChange",sep=" "))

        #if (length(folchout)==0){
	    		#folchout <<- FolCh
        #}else{
	    		#folchout <<- merge(folchout,FolCh,sort=FALSE,by=1)
        #} # end else

				RunFreqAnalys = 1
				next

    }else if(c[1] != "" && c[2] != "" && c[3] == "paired" && is.na(c[4])){
        print("3")
        categ1= toString(c[1,1])
        categ2= toString(c[1,2])

        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
					genes = genes1
				}else{                              # else take union
	    		genes = union(genes1,genes2)
				} # end else

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)

        # calculate comparison
        outEachExperiment = calculateComparison(genes, expressionData, categ1, categ2, dict, argv$comparMethod, c[3], indx)
				#print(head(outEachExperiment))

        if (length(compout)==0){
	    		compout <<- outEachExperiment
        }else{
	    		compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
        } # end else

        # calculate fold change
        #FolCh =  calcFC(outEachExperiment[,5],outEachExperiment[,6]) # median to calc FC
        #FolCh = t(rbind(outEachExperiment[,1], t(FolCh)))
        #colnames(FolCh) = c(colnames(outEachExperiment)[1],paste(categ1,"vs",categ2,"FoldChange",sep=" "))

        #if (length(folchout)==0){
	    		#	folchout <<- FolCh
        #}else{
	    		#	folchout <<- merge(folchout,FolCh,sort=FALSE,by=1)
        #} # end else

				RunFreqAnalys = 1
				next

    }else if(c[1] != "" && is.na(c[2]) && is.na(c[3]) && c[4] == "correlation"){
        print("4")
        categ= toString(c[1,1])
				categ = paste(ColHead, categ, sep=COLHEADGRPSEP)

        # calculate correlations

				start_time <- Sys.time()

        tmp_tmp = ''
	    	if(argv$multicore){
					tmp_tmp = calculateCorrelation_parallel(pairs, expressionData, categ, dict, argv$correlMethod, indx)
	    	} else {
					tmp_tmp = calculateCorrelation(pairs, expressionData, categ, dict, argv$correlMethod, indx)
	    	}

				end_time <- Sys.time()
				timetaken = end_time - start_time
				print("Time Taken:")
				print(timetaken)

	   		outEachExperiment = tmp_tmp
				#print(head(outEachExperiment))

				if(argv$writeCorrPerAnalysis){
					write_to_File_CorrPerAnalysis(outEachExperiment, indx)
				} else if (length(corrout)==0){
					corrout <<- outEachExperiment
				}else{
					corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
				} # end else

        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
	    		genes = genes1
				}else{                              # else take union
	    		genes = union(genes1,genes2)
				} # end else

				# calculate partial correlations only for the given list of genes (all pairs of the list)
				# instead of calc PC with the 35K genes which will be lot of (unnecessary) computation
				if(FALSE){
				outEachExperiment = calculatePartialCorrelation(pairs, genes, expressionData, categ, dict, argv$correlMethod, indx)
				if (length(corrout)==0){
				  corrout <<- outEachExperiment
				}else{
				    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
				} # end else
				# the above does not have gene1-gene2 and gene2-gene1
				}

				next

    }else if(c[1] != "" && is.na(c[2]) && c[3] == "delta" && c[4] == "correlation"){
				print("5")
				categA= toString(c[1,1])
				cA = strsplit(categA, ",")[[1]]
				categ1 = cA[1]
				categ2 = cA[2]

				if(categ1 == "" || categ2 == ""){
				    print("Error in number of categories. Give categ1,categ2")
				    next
				}

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)

        # calculate correlations
        outEachExperiment = calculateDeltaCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, indx)
				#print(head(outEachExperiment))

				if(argv$writeCorrPerAnalysis){
					write_to_File_CorrPerAnalysis(outEachExperiment, indx)
				} else if (length(corrout)==0){
						corrout <<- outEachExperiment
				}else{
						corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
				} # end else

				next

    }else if(c[1] != "" && c[2] != "" && is.na(c[3]) && c[4] == "correlation"){
				print("6")
				# compare correlations of gene pairs
				categ1= toString(c[1,1])
				categ2= toString(c[1,2])

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)

				# calculate correlations
				outEachExperiment = calculateComparisCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, c[3], indx)
				#print(head(outEachExperiment))

				if(argv$writeCorrPerAnalysis){
					write_to_File_CorrPerAnalysis(outEachExperiment, indx)
				} else if (length(corrout)==0){
						corrout <<- outEachExperiment
				}else{
						corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
				} # end else

				next

    }else if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == "correlation"){
				print("7")
				categA= toString(c[1,1])
				categB= toString(c[1,2])

				cA = strsplit(categA, ",")[[1]]
				categ1 = cA[1]
				categ2 = cA[2]

				cB = strsplit(categB, ",")[[1]]
				categ3 = cB[1]
				categ4 = cB[2]

				if(categ1 == "" || categ2 == "" || categ3 == "" || categ4 == ""){
				    print("Error in number of categories. Give categ1,categ2\tcateg3,categ4")
				    next
				}

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)
				categ3 = paste(ColHead, categ3, sep=COLHEADGRPSEP)
				categ4 = paste(ColHead, categ4, sep=COLHEADGRPSEP)

				# calculate correlations
				outEachExperiment = calculateComparisDeltaCorrelation(pairs, expressionData, categ1, categ2, categ3, categ4, dict, argv$correlMethod, indx)
				#print(head(outEachExperiment))

				if(argv$writeCorrPerAnalysis){
					write_to_File_CorrPerAnalysis(outEachExperiment, indx)
				} else if (length(corrout)==0){
					corrout <<- outEachExperiment
				}else{
					corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
				} # end else

				next

    }else if(c[1] != "" && c[2] != "" && c[3] == "paired" && c[4] == "correlation"){
        print("8")

    }else if(c[1] != "" && c[2] != "" && c[3] == "timelag" && c[4] == "correlation"){
				print("9")
				# std. correlations: g1T1-g2T1 or g1T2-g2T2
				# calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc.
				# for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene
				categ1= toString(c[1,1])
				categ2= toString(c[1,2])

				categ1 = paste(ColHead, categ1, sep=COLHEADGRPSEP)
				categ2 = paste(ColHead, categ2, sep=COLHEADGRPSEP)

				# calculate correlations
				outEachExperiment = calculateTwoCategCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, indx)
				#print(head(outEachExperiment))

				if(argv$writeCorrPerAnalysis){
					write_to_File_CorrPerAnalysis(outEachExperiment, indx)
				} else if (length(corrout)==0){
						corrout <<- outEachExperiment
				}else{
						corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
				} # end else

				next

    }else{
        print("Unknown")
    }
}


if(argv$warnings){print(warnings())}


if(! argv$writeCorrPerAnalysis){
	write.csv(corrout,paste(outputFile,"output.csv",sep='corr-'),row.names=FALSE)
} else {
	# merge the diff corr out files
	CorrOFiles = list.files(path = odir, pattern="Analys-.*-corr-output.csv", full.names=T)
	corrout = read.csv(CorrOFiles[1], header = TRUE, check.names=FALSE, na.strings=c("","na","NA", "Na", "NaN"))

	print(paste0("Deleting corr analysis file: ", CorrOFiles[1]))
	system(paste0("rm -f ", CorrOFiles[1]))

	for (cfile in CorrOFiles[-1]) {
    cdftmp = read.csv(cfile, header = TRUE, check.names=FALSE, na.strings=c("","na","NA", "Na", "NaN"))
		corrout <<- merge(corrout,cdftmp,sort=FALSE,by=1)

		print(paste0("Deleting corr analysis file: ", cfile))
		system(paste0("rm -f ", cfile))
	}

	write.csv(corrout,paste(outputFile,"output.csv",sep='corr-'),row.names=FALSE)
}

write.csv(compout,paste(outputFile,"output.csv",sep='comp-'),row.names=FALSE)
#write.csv(folchout,paste(outputFile,"output.csv",sep='folch-'),row.names=FALSE)

print("Finished performing the requested analyses.")


if(argv$multicore){
	#stop cluster
	stopCluster(cl)
}


# if there is comparison of any sort only then run the following commands #
if(RunFreqAnalys==1){

	print("Performing the frequency analyses.")

	fodir =  strsplit(argv$output, "/")[[1]]
	freqdir = paste(odir, "/freq/", fodir[length(fodir)], sep='')
	freqcmd= paste("Rscript /data/rodriguesrr/scripts/R/calc-frequency_bw.R", argv$expressionDataFile, "--lists", list1File, list1File, "--mapFile", argv$mapFile, "--mapFileSampleIdColName", samplIdCol, "--AnalysToDoList", argv$AnalysToDoList, "--comparMethod",  argv$comparMethod, "--correlMethod", argv$correlMethod, "--symbolColumnName", symbleColumnName, "--pairedInfoColumn", argv$pairedInfoColumn, "--output", freqdir, "--stringtorelativize", argv$stringtorelativize, sep=" ")

	print(freqcmd)
	system(freqcmd, wait=TRUE)

	print("Finished performing the requested frequency analyses.")

}

print("You may want to run merge-comp-correl-results-files_bw.R next.")
