library(readxl)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

# needs (JAMS) excel file with sheets for abundance, metadata, and lkt names
# needs info about which "ColnameofGroups" to use for lefse

# usage #
# cd /data/rodriguesrr/FMT_pittsburg/PD1_CpG_network_for_DD/pd1-NM-pitts-and-other-cohorts/PD1_MelanomaCohorts_db22-selected/batch_correct_all-cohorts
# give full path as argument
# Rscript /data/rodriguesrr/scripts/R/prep_JAMS_for_lefse.R /data/rodriguesrr/FMT_pittsburg/PD1_CpG_network_for_DD/pd1-NM-pitts-and-other-cohorts/PD1_MelanomaCohorts_db22-selected/batch_correct_all-cohorts Study_Clin_Response



indir = paste0(args[1],"/")
GroupName = args[2]

#--- 
# list all the excel files
#--- 
all_files = list.files(indir, pattern=".xlsx", all.files=FALSE, full.names=F)
#all_files


BIGLIST = c()

for (f in all_files){
  print(f)
  #Sheets = excel_sheets(f)
  #print(Sheets)


  #--- 
  # create output directory with all needed files and for results
  #---
  odir = str_replace(f, ".xlsx", "")
  dir.create(odir)
  
  
  #--- 
  # metadata
  #--- 
  Metadata_df = read_excel(f, sheet = "Metadata", na = "N_A") # tibble
  #print(Metadata_df[, c("row.names", GroupName)])
  write.table(Metadata_df[, c("row.names", GroupName)], paste0(odir, "/map.txt"), row.names=F, quote=F, sep="\t")

  
  #--- 
  # LKT abundance
  #--- 
  LKT_PPM_df = read_excel(f, sheet = "LKT_PPM", na = "N_A") # tibble
  LKT_PPM_df = data.frame(LKT_PPM_df)
  colnames(LKT_PPM_df)[1] = c("#OTU ID")
  print(head(LKT_PPM_df))
  
  
  #--- 
  # lkt names
  #--- 
  LKT_featuretable_df = read_excel(f, sheet = "LKT_featuretable", na = "N_A") # tibble
  LKT_featuretable_df = data.frame(LKT_featuretable_df)
  #print(head(LKT_featuretable_df))
  
  taxonomy = paste(LKT_featuretable_df[, "Kingdom"], LKT_featuretable_df[, "Phylum"],  
                LKT_featuretable_df[, "Class"], LKT_featuretable_df[, "Order"],
                LKT_featuretable_df[, "Family"], LKT_featuretable_df[, "Genus"], 
                LKT_featuretable_df[, "Species"], sep="; ")
  #print(head(taxonomy))
  
  name_info = cbind(LKT_featuretable_df[, "row.names"], taxonomy)
  colnames(name_info) = c("#OTU ID", "taxonomy")
  print(head(name_info))
  #print(dim(name_info))
  
  
  #--- 
  # otu table for prep for lefse command
  #--- 
  otuTable = merge(LKT_PPM_df, name_info, by="#OTU ID", all.x=T)
  print(head(otuTable))
  print(dim(LKT_PPM_df))
  print(dim(otuTable))
  
  ofname = paste0(odir, "/", odir, ".tmp.tsv")
  write.table(otuTable, ofname, row.names=F, quote=F, sep="\t")
  
  outfname = paste0(odir, "/", odir, ".tsv")
  comand = paste("echo '# Constructed from biom file' >", outfname)
  print(comand)
  system(comand)
  
  comand2 = paste("cat", ofname, ">>", outfname)
  print(comand2)
  system(comand2)
  
  
  #--- 
  # commands for running lefse pipeline
  #--- 
  cmd1 = paste("cd", paste0(indir, odir))
  cmd2 = paste("Rscript /data/rodriguesrr/scripts/R/prep_qiime2_for_lefse.R", paste0(odir, ".tsv"), "map.txt 1")
  cmd3 = paste("chmod +x", paste0(odir, ".tsv.lefse.sh"))
  cmd4 = paste0("./", paste0(odir, ".tsv.lefse.sh"))
  
  BIGLIST = c(BIGLIST, cmd1, cmd2, cmd3, cmd4)
  
}

print(BIGLIST)
write(BIGLIST, "launch.cmds.for.lefse.sh", sep="\n")
