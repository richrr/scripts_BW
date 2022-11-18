# R version 4.1
args <- commandArgs(trailingOnly = TRUE)

cohort_w_cutpoint <- args[1]
print(cohort_w_cutpoint)

cohort_w_cutpoint_df = read.table(cohort_w_cutpoint, sep = "\t", header = TRUE, na.strings=c("", "NA", "NaN", "N_A"))
#head(cohort_w_cutpoint_df)
dim(cohort_w_cutpoint_df)

second_cohort <- args[2]
print(second_cohort)

second_cohort_df = read.table(second_cohort, sep = "\t", header = TRUE, na.strings=c("", "NA", "NaN", "N_A"))
#head(second_cohort_df)
dim(second_cohort_df)

lktnames = grep("LKT" ,colnames(second_cohort_df), value=T)
head(lktnames)

BIG_CMDS_LIST = c()
counter = 0
for(otu in lktnames){
	counter = counter + 1
  print(counter)
  print(otu)
  
  # if exists in first cohort
  if( !(otu %in% cohort_w_cutpoint_df[,"ID"]) ){
      print("Absent in First cohort")
      next
  } #else{
    #  print("Exists")
  #}
  

  info_in_first_cohort = cohort_w_cutpoint_df[which(cohort_w_cutpoint_df[,"ID"] == otu), ]
  print(info_in_first_cohort)
  cutpoint_value_in_first_cohort = sort(unlist(info_in_first_cohort[,2]))[1]  # take the "smallest" entry if more than one, e.g. for Roseburia_faecis which was twice 
  print(cutpoint_value_in_first_cohort)
  
  sdf = second_cohort_df[,c(otu, "dmfs_event", "dmfs_time")]
  print(head(sdf))
  tmp_lkt_file = paste0(getwd(), "/temp/cohort_two_per_lkt_infile/", otu, "_", cutpoint_value_in_first_cohort,".txt")
  write.table(sdf, tmp_lkt_file, row.names=F, sep="\t", quote=F)
  
  # calls
  cmd_to_run = paste("Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Multi-var-Cox_Regression.R ./", tmp_lkt_file, otu, cutpoint_value_in_first_cohort)
  print(cmd_to_run)
	BIG_CMDS_LIST = c(BIG_CMDS_LIST, cmd_to_run)
  system(cmd_to_run)
  
  #if(counter == 5)
  #{
  #    break
  #}
}

write(BIG_CMDS_LIST, "commands-ran.txt")
print("Done!")
