#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# usage:
#cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/TransNet/analysis/corr/merged/ENR.ER.LR/p1_cp1_cfdr1/
#cut -d, -f 1,x-y merged_FolChMedian_merged-parallel-output.csv > merged_FolChMedian_merged-parallel-output-coeff-only.csv
#Rscript /data/rodriguesrr/scripts/R/check-consis-corrs-in-merged-file.R merged_FolChMedian_merged-parallel-output-coeff-only.csv 8

# adding additional cores increases the RAM requirements and job crashes
#swarm -g 200 -t 8 --time 1-00:00:00 --module R/3.6.0 --logdir swarmlogs -f consis-swarm-cmd.sh



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

libraries <- c("argparser", "hash", "psych", "gtools", "corpcor", "reshape", "doParallel", "readxl", "tools", "meta", "stringr")
uname = Sys.getenv(c("USER"))
pathRlibs = paste0("/data/",uname, "/R_LIBS")
print(pathRlibs)
ifelse(!dir.exists(pathRlibs), dir.create(pathRlibs), FALSE)

for(lib in libraries)
{
	checkInstallAndLoadRlibs(lib, pathRlibs)
}

foldchThreshold = 0
correlThreshold = 0
pvalThreshold = 1
combinedPvalueCutoff = 1
combinedFDRCutoff = 1



allowed_cores = args[2]
cl <<- makeCluster(allowed_cores, type="FORK")



#### trim the pvalues and FDR columns ####
merged_df = read.csv( args[1], header=TRUE, check.names=FALSE, row.names=1)
merged_df = merged_df[, grep(" Coefficient_Expt_",colnames(merged_df), value=T)]
head(merged_df)
total_numb_input_files = length( grep(" Coefficient_Expt_", 
                                                grep("Analys 1 ", colnames(merged_df), value=T),
                                                      value=T))
total_numb_input_files



#------------------------------------------------------------------
# this finds the analysis number
#------------------------------------------------------------------
find_analysis_number = function(x){
   #res = as.numeric(str_match(x, "Analys ([0-9]+) .*")[,2])
   res = str_match(x, "Analys ([0-9]+-?[0-9]*) .*")[,2]
   res
}




#------------------------------------------------------------------
# identify the rows where values are >/< than threshold
#------------------------------------------------------------------
calc_consistency = function(idx, s_df, Threshold, total_numb_input_files){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > Threshold)
			neg = sum(dfx < Threshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
			res
}






#------------------------------------------------------------------
# Parallelized: check which measurement ("gene") or pairs have the same trend across experiments
# identify the rows where values are >/< than threshold
#------------------------------------------------------------------
checkConsistency_across_expts_parallel = function(s_df, condition, total_numb_input_files, correlThreshold , pvalThreshold, foldchThreshold){


    list_rows_passing_consistency = list()

	resl=unique(unlist(lapply(colnames(s_df), find_analysis_number)))
	print(resl)

	consistent_rows_per_analysis = c()

    if(condition == "Coefficient" || condition == "DiffCorrs"){

		rows_passing_consistency = c()
		rows_passing_consistency = unlist(parLapply(cl, 1:nrow(s_df), calc_consistency, s_df, correlThreshold, total_numb_input_files))

		list_rows_passing_consistency$ConsisCorrel = rows_passing_consistency
		consistent_rows_per_analysis = append(consistent_rows_per_analysis, rows_passing_consistency)

		##### send this for per analysis #####
	    # this part is useful to calc. comb FDR on only the consistent genes
      #  CalcCombPvalFdrPerAnalysis(resl, consistent_rows_per_analysis, merged_df, "CORR")
      
      analyses = paste(c("Analys ", resl, " "), collapse='')
      o_f_name = gsub(' ' , '', paste("./per_analysis/",analyses, "-consis_genes.csv", sep=''))
      write(consistent_rows_per_analysis, o_f_name)  
    } 


	  #print(consistent_rows_per_analysis)
    #print(list_rows_passing_consistency)

    return(list_rows_passing_consistency)
}




#------------------------------------------------------------------
# check Consistency across ALL expts using the thresholds
#------------------------------------------------------------------
#result_dumper = list()  # to dump the consistent results in output files

cols_processed = 1
while(cols_processed < ncol(merged_df))  # loop over the merged dataset in steps of (number of input files)
{
   subset_df = merged_df[,cols_processed:(cols_processed+total_numb_input_files-1)]
   cols_processed = cols_processed + total_numb_input_files

   if(grepl("Coefficient", colnames(subset_df)[1])){
       #out_res = list()
       #out_res$name = colnames(subset_df) #paste(colnames(subset_df), collapse='<-->')

		   tmp_tmp = checkConsistency_across_expts_parallel(subset_df, "Coefficient", total_numb_input_files, correlThreshold , pvalThreshold, foldchThreshold)

       #out_res = append(out_res, tmp_tmp)
       #result_dumper[[length(result_dumper)+1]] <- out_res

   } 

}

print("Not dumping result_dumper into file.")

#stop cluster
stopCluster(cl)

print("Done")
