#library(argparser)
#library(stringr)


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

libraries <- c("argparser", "stringr")
uname = Sys.getenv(c("USER"))
pathRlibs = paste0("/data/",uname, "/R_LIBS")
print(pathRlibs)
ifelse(!dir.exists(pathRlibs), dir.create(pathRlibs), FALSE)

for(lib in libraries)
{
	checkInstallAndLoadRlibs(lib, pathRlibs)
}



### for (i) FC and (ii) correlation, we can use any dataset (of the multiple datasets) since we
### (i) calc. correlations only on consistent genes and (ii) calc PUC only on consistent pairs.
### using median values "USES" all datasets.


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Extract pairs which pass the FDR and PUC criteria across datasets for a single correlation (group)')
p <- add_argument(p, "--file", help="file which has merged correlations from multiple expts", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--consistent", help="file which has list of consistent elements (mostly pairs)", nargs=1) # required; but written as optional format so I explicitly mention the "--consistent"
p <- add_argument(p, "--group", help="use the correlation for 'group'. Generate the network for this group") # required; but written as optional format so I explicitly mention the arg
#p <- add_argument(p, "--groups", help="use the correlations for 'groups'. Calc. PUC using both these groups (states). Default only uses one state as mentioned in --group", nargs=2) # not implemented
p <- add_argument(p, "--foldchange", help="file containing foldchange to be used in analysis") ### use the merged file with consistent genes across datasets
p <- add_argument(p, "--output", help="output file", default="./build_netw_")
p <- add_argument(p, "--indivPvalCutoff", help="individualPvalueCutoff", default=0.3, type="numeric") # 0.05
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.05, type="numeric") # 0.05
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.15, type="numeric") #
p <- add_argument(p, "--foldchthresh", help="fold change threshold", default=0, type="numeric") # default is logged data so check whether greater than 0. use 1 for unlog data
p <- add_argument(p, "--logbase", help="calc log using the base", default=2) # allowed: 0 (no log), 1 (e), 2, 10

p <- add_argument(p, "--foldchMean", help="use fold change from mean", flag=TRUE)  # default fold change is median
p <- add_argument(p, "--foldchGeoMean", help="use fold change from mean", flag=TRUE)  # default fold change is median
#p <- add_argument(p, "--transposeOutput", help="tranpose the results so the rows are analysis and columns are gene or pairs", flag=TRUE)  # easier for downstream grep of required analysis, do not use when you expect many genes or pairs since it might truncate when you open in excel or librecalc due ot limited columns

p <- add_argument(p, "--noPUC", help="do not calc. PUC, so you do not need to use fold change", flag=TRUE)

p <- add_argument(p, "--analysisfc", help="Partial header (Analysis number) to be used for median fold change calc", default="ALL") #e.g. "Analys 1 " ; avoids cutting the required columns as input for the script  # make sure there is space after the number to avoid selecting analysis 11, etc.

p <- add_argument(p, "--analysiscorr", help="Partial header (Analysis number) to be used for selecting the columns for combined pval and median correlations", default="ALL") #e.g. "Analys 1 " ; avoids cutting the required columns as input for the script   # make sure there is space after the number to avoid selecting analysis 11, etc.

p <- add_argument(p, "--numbDataFromFCfile", help="use fold change file to figure out the number of datasets", flag=TRUE)  # default is figure out from the pvalue cols of merged col file.
						# ue this argument in very very rare cases, where the number of datasets for the corr and comparisons are different.
						# e.g. the corr is calculated by pooling 2 expts of say 5 samples into 1 dataset of 10 samples, but the comp file is for 2 datasets (expts)


argv <- parse_args(p)
#print (p)


if(length(argv$file) < 1)
{
  print("At least 1 file is required. Give --file file ... in cmd")
  quit()
}

outputFile = argv$output

individualPvalueCutoff = argv$indivPvalCutoff
combinedPvalueCutoff = argv$combPvalCutoff
combinedFDRCutoff = argv$combFDRCutoff
FoldChangeFile = argv$foldchange
search_group = argv$group
search_group_santized = gsub('+', '_', search_group, fixed=TRUE)
#print(search_group_santized)
search_group_santized = gsub("\\", '', search_group_santized, fixed=TRUE)
#print(search_group_santized)

foldchthresh = argv$foldchthresh
logbase = argv$logbase

foldchVar = 'FolChMedian'
if(argv$foldchMean){foldchVar = "FoldChange"}
if(argv$foldchGeoMean){foldchVar = "FolChGeoMean"}
outputFile = paste(outputFile , foldchVar, '_', sep='')


noPUC = argv$noPUC
analysisfc = argv$analysisfc
analysiscorr = argv$analysiscorr

if(analysisfc != "ALL" || analysiscorr != "ALL"){
  outputFile = paste(outputFile , "FC", analysisfc, "Corr", analysiscorr, '_', sep='')
  outputFile = gsub(' ', '_', outputFile)
}

networkFile = paste(c(outputFile, search_group_santized, "indiv-pval", individualPvalueCutoff ,"comb-pval", combinedPvalueCutoff, "comb-fdr", combinedFDRCutoff, ".csv"), collapse='_')


#------------------------------------------------------------------------------------------------
# check whether the first "number of files" columns are the same and remove them after setting row names
#------------------------------------------------------------------------------------------------
remove_redundant_columns = function(this_df, total_numb_input_files){
         l_identical_inp = as.character(this_df[, 1])


         # if there are more than 1 datasets
         if(total_numb_input_files>1){
             # create a data frame of all datasets
             for(z in 2:total_numb_input_files){
                    l_identical_inp = cbind(l_identical_inp, as.character(this_df[,z]))
             }

             if(all(apply(l_identical_inp, 2, identical, l_identical_inp[, 1]))){
                    print("All are equal")
                    rownames(this_df) <- this_df[, 1]                           ## set column 1 as rownames
                    this_df <- this_df[, -c(1:total_numb_input_files)]          ## remove the first "total numb of input files" columns
             } else {
                    print(l_identical_inp)
                    print("All are NOT equal")
                    quit()
             }
          } else{ # the rare case where you only have one dataset
          	rownames(this_df) <- this_df[, 1]
          	this_df <- this_df[, -1, drop=F] # june 21 2020
          }

      return(this_df)
}


#------------------------------------------------------------------------------------------------
# calc PUC (%) at each FDR threshold and plot
#------------------------------------------------------------------------------------------------
calc_PUC_at_thresholds = function(PUCoutfile, df, str='all'){

    #write.csv(df, paste("tmp", str, "csv", sep='.'))
    combined_fdr = c(0)
    puc_percent = c(0)
    df = as.matrix(df)
    df = apply(df[,c("combinedFDR", "PUC")],2,function(x){as.numeric(as.vector(x))})
    #fdr_range = range(as.vector(df[,"combinedFDR"]))
    fdr_min = min(as.vector(df[,"combinedFDR"]))
    fdr_max = max(as.vector(df[,"combinedFDR"]))
    print(fdr_min)
    print(fdr_max)
    m_df = as.matrix(df[,c("combinedFDR", "PUC")])
    m_df = apply(m_df,2,function(x){as.numeric(as.vector(x))})

    for(r in sort(seq(as.numeric(fdr_min), as.numeric(fdr_max), length.out=100))){  ##length.out: desired length of the sequence
        combined_fdr = c(combined_fdr, r)

        #nrow(m_df[which(m_df[,"combinedFDR"] <= r & m_df[,"PUC"] == 1), , drop=FALSE])
        df_under_fdr = m_df[which(m_df[,"combinedFDR"] <= r), , drop=FALSE]
        DEN = nrow(df_under_fdr)

        goog_puc_in_df_under_fdr =  df_under_fdr[which(df_under_fdr[,"PUC"] == 1), , drop=FALSE]
        NUM = DEN - nrow(goog_puc_in_df_under_fdr)

        puc_percent = c(puc_percent, as.numeric(NUM*100/DEN))
    }

    plt = cbind(combined_fdr, puc_percent)

    pdf(paste(PUCoutfile, str, 'FDRvsPUC.pdf', sep='-'))
    plot(combined_fdr, puc_percent, type="o", ylab="PUC(%)", xlab="combined FDR", pch=10, cex=.2, ylim=c(0, 100) )
    #plot(plt, type="o", ylab="PUC(%)", xlab="combined FDR", pch=10, cex=.2 )

    dev.off()

    return(plt)

}


   # consistent pairs
   consist_elems = read.csv( argv$consistent, header=FALSE)
   consist_elems = as.vector(consist_elems[!duplicated(consist_elems[,1]),1] )	#delete duplicated array probes
   consist_elems = consist_elems[order(consist_elems)] # ascending sort
   consist_elems = grep("<==>", consist_elems, value=TRUE) # keep pairs

   # keep consistent pairs only
   data = read.csv( argv$file, header=TRUE, check.names=FALSE)
   subdata = data[data[,1] %in% consist_elems,]
   total_numb_input_files = length(grep("pairName", colnames(subdata)))
   subdata = remove_redundant_columns(subdata, total_numb_input_files)
   #print("Removed redundant id cols from corr file")

   outForPUC = subdata

   CorrAnalysColnames = colnames(outForPUC)
   # only keep the requested analysis number
   if(analysiscorr != "ALL"){
       CorrAnalysColnames = CorrAnalysColnames[grep(analysiscorr,CorrAnalysColnames)]
   }
   print(CorrAnalysColnames)
   outForPUC = outForPUC[,CorrAnalysColnames]
   #print(head(outForPUC))

   # calculate combined Pvalue for interest group
   PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
   PvalueColnames = PvalueColnames[grep(search_group,PvalueColnames)]
   #print(PvalueColnames)
   total_numb_input_files = length(PvalueColnames)
   interestedPvalueData = outForPUC[,PvalueColnames]
   interestedPvalueData = as.matrix(interestedPvalueData)
   interestedPvalueData = apply(interestedPvalueData,2,function(x){as.numeric(as.vector(x))})
   combinedPvalue = apply(interestedPvalueData,1
							,function(pvalues){
										pvalues = pvalues[!is.na(pvalues)]
										statistics = -2*log(prod(pvalues))
										degreeOfFreedom = 2*length(pvalues)
										combined = 1-pchisq(statistics,degreeOfFreedom)
									}
							)
   outForPUC = cbind(outForPUC,combinedPvalue)

   # calculate median coefficient for interest
   CoefficientColnames = colnames(outForPUC)[grep("Coefficient",colnames(outForPUC))]
   CoefficientColnames = CoefficientColnames[grep(search_group,CoefficientColnames)]
   #print(CoefficientColnames)
   interestedCoefficientData = outForPUC[,CoefficientColnames]
   interestedCoefficientData = as.matrix(interestedCoefficientData)
   interestedCoefficientData = apply(interestedCoefficientData,2,function(x){as.numeric(as.vector(x))})
   combinedCoefficient = apply(interestedCoefficientData,1, function(x){round(median(x, na.rm = TRUE), 3)})
   outForPUC = cbind(outForPUC,combinedCoefficient)

   result = outForPUC
   #calculate FDR for combined pvalue
   combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
   result = cbind(result,combinedFDR)
   #print(head(result))
   write.csv (result,paste( outputFile, search_group_santized, "tmp-out.csv", sep='-'))


   # calculate median FoldChange
   FoldChangeMetabolic =  ''
   if(noPUC){
           # no need to loop up the fold change
        } else {
   FoldChangeMetabolic =  read.csv(FoldChangeFile,header = TRUE,check.names=FALSE)

   if(argv$numbDataFromFCfile){
   		total_numb_input_files = length(grep("geneName", colnames(FoldChangeMetabolic)))
   }


   FoldChangeMetabolic = remove_redundant_columns(FoldChangeMetabolic, total_numb_input_files)
   #print("Removed redundant id cols from fc file")

	 #print(head(FoldChangeMetabolic))
   FoldChangeColnames = colnames(FoldChangeMetabolic)[grep(foldchVar,colnames(FoldChangeMetabolic))]
   # if an Analysis number is specified then select only that analysis number
   if(analysisfc != "ALL"){
       FoldChangeColnames = FoldChangeColnames[grep(analysisfc,FoldChangeColnames)]
   }
   #print(FoldChangeColnames)

   interestedFoldChangeData = FoldChangeMetabolic[,FoldChangeColnames]
   interestedFoldChangeData = as.matrix(interestedFoldChangeData)
   interestedFoldChangeData = apply(interestedFoldChangeData,2,function(x){as.numeric(as.vector(x))})
   combinedFoldChange = apply(interestedFoldChangeData,1, function(x){round(median(x, na.rm = TRUE), 3)})
   FoldChangeMetabolic$geneName = rownames(FoldChangeMetabolic)
   FoldChangeMetabolic = cbind(FoldChangeMetabolic,combinedFoldChange)
   }


pow <- function(x=10, y=6) {
   # function to print x raised to the power y
   result <- x^y
   return(result)
}


unlog <- function(FoldChangeData, base){
	if(logbase != 0){
	    if(logbase == 1) {
	          FoldChangeData = apply(FoldChangeData, 2, function(x) {exp(x)})    # using the base e
	    } else {
	          FoldChangeData = apply(FoldChangeData, 2, function(x) {pow(base, x)})
	    }
	}
	return(FoldChangeData)
}

###########################################################################################################
#calc PUC
###########################################################################################################

   Data = result
   head(result)
#--------------------------------------------------------------------------------
# calculate PUC
#--------------------------------------------------------------------------------
forPUC = function(FoldChangeMetabolic,noPUC){

	#change out format to partner1 partner2 for PUC
	row_names_Data = rownames(Data)
	pair = str_split( row_names_Data ,"<==>")
	pairs = t(as.data.frame(pair))

	colnames(pairs) = c("partner1","partner2")
	Data = apply(Data, 2, function(x) as.numeric(as.character(x))) # convert the chars to numeric
	rownames(pairs) = row_names_Data # remove this if you do not want row names
        outForPUC = cbind(pairs,Data)

        grep_cols_c = grep("pvalue", colnames(outForPUC), value=TRUE, fixed=TRUE) #ignore.case = TRUE ,
        grep_cols_c = append(grep_cols_c, grep("combined" , colnames(outForPUC), value=TRUE, fixed=TRUE) )

        g_grep_cols = c("partner1","partner2", grep_cols_c)

        outForPUC = outForPUC[,g_grep_cols]

        if(noPUC){
           # no need to loop up the fold change
        } else {
        # attach the foldChange information for each partner
        FoldChangeCol = grep("combined" , colnames(FoldChangeMetabolic), value=TRUE, fixed=TRUE)
	FoldMetab1_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner1"]), c("geneName", FoldChangeCol)]
	colnames(FoldMetab1_InPair) = c("partner1InFold","partner1_FoldChange")
	FoldMetab2_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner2"]),c("geneName", FoldChangeCol)]
	colnames(FoldMetab2_InPair) = c("partner2InFold","partner2_FoldChange")

    FoldMetab1_InPair = cbind(FoldMetab1_InPair[, 1, drop=F], unlog(FoldMetab1_InPair[, 2, drop=F], logbase))
	FoldMetab2_InPair = cbind(FoldMetab2_InPair[, 1, drop=F], unlog(FoldMetab2_InPair[, 2, drop=F], logbase))
	#print(head(FoldMetab1_InPair))
	#print(head(FoldMetab2_InPair))

	outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
	}

        # calculate correlation Direction For combined correlation coefficient of interest
        # at this point we only have the consistent pairs left, so the value of combined corr coeff is ok to use
        interestedCoefficientColnames = grep("Coefficient",colnames(outForPUC), value=TRUE, fixed=TRUE)
	print(interestedCoefficientColnames)
	interestedCorrelationData = outForPUC[,interestedCoefficientColnames, drop=FALSE]
	interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
	signOfInterestedCorrelationData = interestedCorrelationData/abs(interestedCorrelationData)
	rownames(signOfInterestedCorrelationData) = c()
	colnames(signOfInterestedCorrelationData) = paste(colnames(interestedCorrelationData),"correlationDirection",sep=".")
	matchedExpressionDirection = signOfInterestedCorrelationData


        if(noPUC){
           # no need to loop up the fold change direction
           outForPUC = cbind(outForPUC,signOfInterestedCorrelationData)
        } else {
	# calculate fold change direction for each partner
	FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))] # since this is using the combined fold change calculated above, you do not need the foldchVar variable
	FoldChangeData = outForPUC[,FoldChangeColnames]


	#print(head(FoldChangeData))

	FoldChangeDirection = NA

	#if (foldchthresh==0){
	#    FoldChangeDirection <- as.matrix(FoldChangeData)
	#    FoldChangeDirection[FoldChangeDirection<0] <- -1
	#    FoldChangeDirection[FoldChangeDirection>=0] <- 1
	#    FoldChangeDirection = as.data.frame(FoldChangeDirection)
	#    #print(head(FoldChangeDirection))
	#} else { # for unlog data
	    FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
	#}
	names(FoldChangeDirection) = c()
	colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"FoldChangeDirection",sep=".")

	# calculate if fold change direction are the same for the two partners
	IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
	names(IfFoldChangeDirectionMatch) = c()
	colnames(matchedExpressionDirection) = c()

	# use "matchedExpressionDirection" and "IfFoldChangeDirectionMatch" to calc PUC, i.e. if these two are the same PUC=1 (good)
	PUC = IfFoldChangeDirectionMatch * matchedExpressionDirection
	outForPUC = cbind(outForPUC,signOfInterestedCorrelationData,FoldChangeDirection,IfFoldChangeDirectionMatch,PUC)
        }

	return(outForPUC)

}

    PUCoutfile = ''
    if(noPUC){
           PUCoutfile = paste(outputFile, search_group_santized, "noPUC-output.csv",sep='-')
    } else {
         PUCoutfile = paste(outputFile, search_group_santized, "PUC-output.csv",sep='-')
    }

    if(file.exists(PUCoutfile))
    {
        print("Deleted old PUC output")
        file.remove(PUCoutfile)
    }

    result = forPUC(FoldChangeMetabolic,noPUC)
    data = result

    if(noPUC){
           # no need to calc stats
    } else {

    # calculate percentage of unexpected correlations in the entire file
    print("Total items")
    DEN = length(result[,"PUC"])
    print(DEN)

    print("Items not 1")
    NUM = DEN - length(result[which(result[,"PUC"]==1), "PUC"])
    print(NUM)

    PUC_Prop = as.numeric(NUM*100/DEN)
    print(PUC_Prop)
    PUC_Prop = paste(c("Percentage of Unexpected Correlation=" , PUC_Prop, "%") , collapse='')


    if(file.exists(PUCoutfile)){
      out = file(PUCoutfile ,'a')
      write.csv(result, file=out, row.names=FALSE)
      close(out)
    } else{
      write.csv(result, PUCoutfile ,row.names=FALSE)
    }

    if(file.exists(PUCoutfile)){
      out = file(PUCoutfile ,'a')
      write.csv(PUC_Prop, file=out, row.names=FALSE)
      close(out)
    }

    print("Done with PUC")

    # sort as per combined fdr
    sorted_result = result[order(result$"combinedFDR"), ]
    plt = calc_PUC_at_thresholds(PUCoutfile, sorted_result) # walk along fdr and calc puc at diff fdr
    write.csv(plt, paste(PUCoutfile, "all-edges.csv", sep='.'), row.names=FALSE)

    # if there are only pos or neg correlations, use the if condition to avoid error
    if(nrow(sorted_result[which(sorted_result$combinedCoefficient.correlationDirection == 1),])>0){
    pos_plt = calc_PUC_at_thresholds(PUCoutfile, sorted_result[which(sorted_result$combinedCoefficient.correlationDirection == 1),], "pos") # walk along fdr and calc puc at diff fdr
    write.csv(pos_plt, paste(PUCoutfile, "pos-edges.csv", sep='.'), row.names=FALSE)
    }
    if(nrow(sorted_result[which(sorted_result$combinedCoefficient.correlationDirection == -1),])>0){
    neg_plt = calc_PUC_at_thresholds(PUCoutfile, sorted_result[which(sorted_result$combinedCoefficient.correlationDirection == -1),], "neg") # walk along fdr and calc puc at diff fdr
    write.csv(neg_plt, paste(PUCoutfile, "neg-edges.csv", sep='.'), row.names=FALSE)
    }
    data = sorted_result
    }



###########################################################################################################
#generate network
###########################################################################################################


#----------------------------------------------------------------------
# calc stats of the generated network
#----------------------------------------------------------------------
calc_stats = function(inNet, PUC_Prop, correlThreshold=0){

        out = file(paste(networkFile,'-stats-log.txt',sep='') ,'a')

        combcoeffcol = as.vector(inNet[,"combinedCoefficient"])
        # not sure why it behaves this way!
        lowest_neg_corrl = min(combcoeffcol[combcoeffcol < 0])
        highest_neg_corrl = max(combcoeffcol[combcoeffcol < 0])

        lowest_pos_corrl = min(combcoeffcol[combcoeffcol > 0])
        highest_pos_corrl = max(combcoeffcol[combcoeffcol > 0])


        coefficientData = inNet[,grep("correlationDirection",colnames(inNet)),drop=FALSE]

		## artifically add row (with 0) incase it is only 1 row
		if(nrow(coefficientData) == 1){
			tmp_v = c("0")
			tmp_artif = data.frame(tmp_v)
			colnames(tmp_artif) = c("combinedCoefficient.correlationDirection")
			coefficientData = rbind(coefficientData, tmp_artif)
		}

		coefficientData = apply(coefficientData,2,function(x){as.numeric(as.vector(x))})

        res_pos = apply(coefficientData, 2, function(x) sum(x == 1))
        res_neg = apply(coefficientData, 2, function(x) sum(x == -1))
        res = cbind(res_pos, res_neg)
        #print(res)

        #print(paste(c("Total number of edges: ", nrow(inNet)), collapse=""))

        setpartner1 = unique(as.vector(inNet[,"partner1"]))
        setpartner2 = unique(as.vector(inNet[,"partner2"]))

        nodes = union(setpartner1, setpartner2)
        #print(paste(c("Number of unique nodes: ", length(nodes)), collapse=""))

        df_a = inNet[,c("partner1InFold","partner1_FoldChange"), drop=F]
        colnames(df_a) = c("partnerInFold","partner_FoldChange")
        rownames(df_a) <- NULL

        df_b = inNet[,c("partner2InFold","partner2_FoldChange"), drop=F]
        colnames(df_b) = c("partnerInFold","partner_FoldChange")
        rownames(df_b) <- NULL

        df_ab = rbind(df_a, df_b)
        uniq_df_ab = unique(df_ab)
        print("Unique nodes (partner, fc) df:")
        print(dim(uniq_df_ab))

        upreg = length(uniq_df_ab[as.numeric(uniq_df_ab[,"partner_FoldChange"]) > 1, "partner_FoldChange"])
        dnreg = length(uniq_df_ab[as.numeric(uniq_df_ab[,"partner_FoldChange"]) < 1, "partner_FoldChange"])

        potential_pos_corr_edges = choose(upreg, 2) + choose(dnreg, 2)
        potential_neg_corr_edges = upreg * dnreg

        edgesDistinctNodes = inNet[which(inNet[,"partner1"]!=inNet[,"partner2"]), ,drop=F]
        print(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse=""))

        write.csv(paste(c("Lowest pos corr:" , lowest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Highest pos corr:" , highest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Lowest neg corr:" , lowest_neg_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Highest neg corr:" , highest_neg_corrl), collapse='') , file=out, row.names=FALSE)

        write.csv(paste(c("Pos corr edges:" , toString(res_pos)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Neg corr edges:" , toString(res_neg)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Ratio of Pos to Neg corr edges:" , toString(res_pos/res_neg)), collapse='') , file=out, row.names=FALSE)

        write.csv(paste(c("Potential Pos corr edges:" , toString(potential_pos_corr_edges)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Potential Neg corr edges:" , toString(potential_neg_corr_edges)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Potential Ratio of Pos to Neg corr edges:" , toString(potential_pos_corr_edges/potential_neg_corr_edges)), collapse='') , file=out, row.names=FALSE)


        write.csv(paste(c("Total number of edges: ", nrow(inNet)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of unique nodes: ", length(nodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of upreg nodes: ", upreg), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of dwnreg nodes: ", dnreg), collapse='') , file=out, row.names=FALSE)

        write.csv(paste(c("Ratio of edges to nodes: ", nrow(inNet)/length(nodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("PUC with ip,fisher,fdr cuts: ", PUC_Prop), collapse='') , file=out, row.names=FALSE)

        close(out)

}


#----------------------------------------------------------------------
# for a given combinedPvalue and combinedfdr threshold, generate network
#----------------------------------------------------------------------
generateNetwork = function(){
    #generate network After Filter

	out = data

	out = as.matrix(out)
    #df = apply(df[,c("combinedFDR", "PUC")],2,function(x){as.numeric(as.vector(x))})

	# feb 10 2022: if it is a single dataset it becomes a vector so adding the drop=FALSE
	out = out[(as.numeric(out[,"combinedPvalue"])<combinedPvalueCutoff)==1,, drop=FALSE]
	out = out[(as.numeric(out[,"combinedFDR"])<combinedFDRCutoff)==1,, drop=FALSE]
	print(head(out))
	
	# find significant in individual pvalue: all of the pvalues for all of the datasets for each pair must be smaller than threshold
	# find pvalue data
  # if it is a single dataset it becomes a vector so adding the drop=FALSE
	pvalueData = out[,grep("pvalue",colnames(out)), drop=FALSE]
	pvalueData = pvalueData[,grep(search_group,colnames(pvalueData)), drop=FALSE]
	#print(head(pvalueData))


	pvalueData = as.matrix(pvalueData)
	#print(head(pvalueData))
	
	## artifically add row (with 0) incase it is only 1 row
	artif_flag = FALSE
	if(nrow(pvalueData) == 1){
			tmp_v = rep(NA, ncol(pvalueData))
			tmp_artif = data.frame(t(tmp_v))
			colnames(tmp_artif) = colnames(pvalueData)
			rownames(tmp_artif) = "DUMMY"
			pvalueData = rbind(tmp_artif, pvalueData) # this somehow converts it to numeric
			#print(head(pvalueData))
			artif_flag = TRUE
	} else{
		pvalueData = apply(pvalueData,2,function(x){as.numeric(as.vector(x))})
	}
	#print(head(pvalueData))
	
	
	# remove the artifical row
	if(artif_flag){
				pvalueData = pvalueData[-1,,drop=F]
	}

	# added this on july 17 2018, after I added this is in the create nets for diff corr
	if(is.null(nrow(pvalueData))){
		print("Qutting since no edges passed the fisher pval or fdr cuts.")
		q()
	}
	

	# calculate the largest pvalue among all datasets for each gene, this smallest pvalue must be smaller than threshold
	passIndevidualPvalue = apply(pvalueData,1,max,na.rm=T)<individualPvalueCutoff  # june 21 2020
	#print(head(passIndevidualPvalue))
	outNetwork = out[passIndevidualPvalue, , drop=FALSE]
	#print(head(outNetwork))

	#outNetwork$names<-rownames(outNetwork)
	#splt = str_split_fixed(outNetwork$names, "<==>", 2)
	#outNetwork = cbind(outNetwork, splt)

	write.csv(outNetwork,paste0(networkFile,"-prePUCcut.csv"), quote=FALSE)
	#print(head(outNetwork))

	 if(noPUC){
	    # do nothing
     } else {
		# find puc after pval, cp, fdr cuts
	 	DEN = length(outNetwork[,"PUC"])
		NUM = DEN - length(outNetwork[as.numeric(outNetwork[,"PUC"])==1, "PUC"])
		PUC_Prop = as.numeric(NUM*100/DEN)
		print("PUC after the ip, fisher, fd cuts:")
		print(PUC_Prop)

	    # keep PUC expected
	    outNetwork = outNetwork[as.numeric(outNetwork[,"PUC"])==1 ,,drop=F]
		outNetwork = outNetwork[!is.na(as.numeric(outNetwork[,"PUC"])),,drop=F] # remove the rows with 'NA' in PUC columns
		write.csv (outNetwork,networkFile, quote=FALSE)
		calc_stats(outNetwork, PUC_Prop)
     }

	#print(head(outNetwork))

    print("Done!")
}


    generateNetwork()
