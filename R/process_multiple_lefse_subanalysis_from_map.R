
#ml R/3.6.0
#cd /data/MicrobiomeCore/analysis/GTAM45_20220714
#Rscript /data/rodriguesrr/scripts/R/process_multiple_lefse_subanalysis_from_map.R /data/MicrobiomeCore/analysis/GTAM45_20220714/map-files/map.xlsx /data/MicrobiomeCore/analysis/GTAM45_20220714/map-files Group
#bash launch_lefse_analysis.sh

# loop in R to create a map file for each sheet/comparison

library(readxl)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

infile = args[1] # map excel file with sheets representing subset of samples for analysis
mapdir = args[2] # dir to create the maps in
beta_sign_col = args[3] # column name to be used for lefse.

indir = getwd()


cmdlist = c()

shits = excel_sheets(infile)
for(f in shits){
  
  print(f)
  df = data.frame(read_excel(infile, sheet=f))
  #print(head(df))
  #print(dim(df))
  
  # [0 based] col number to be used for lefse
  colIndxforLefse = grep(paste0("^",beta_sign_col,"$"), colnames(df)) - 1   # ^Group$ to avoid other columns with "Group" in the name
  print(colIndxforLefse)

  omapfile = paste0(mapdir, "/", f, ".txt")
  #write.table(df, omapfile, row.names=F, quote=F, sep="\t")

  cmd1=paste0("cd ", indir)

  # get full path to otu_table_w_tax.txt in exported folder
  OtuwTax = list.files( path = paste0(indir, "/", f, "/"), pattern = "otu_table_w_tax.txt$", full.names = T, recursive = T)
  #print(OtuwTax)
  
  timestamp = str_replace_all( str_extract(OtuwTax, "/exported_\\d+/"), "[^0-9]", "")  #[^[:digit:]]
  #print(timestamp)
  
  expdir = paste0(indir, "/", f, "/exported_", timestamp, "/")
  q2visdir = paste0(indir, "/", f, "/q2_visualizations_", timestamp, "/")
  if( dir.exists(expdir) && dir.exists(q2visdir) ){
      #dir.create(file.path(q2visdir,"custom"))
      #dir.create(file.path(paste0(q2visdir,"/custom"), "lefse"))
  } else{
      q()
  }
  
  
  # prep_qiime2_for_lefse.R
  cmd2=paste("Rscript /data/rodriguesrr/scripts/R/prep_qiime2_for_lefse.R", OtuwTax, omapfile, colIndxforLefse)
  cmd3 = paste("chmod +x", paste0(OtuwTax, ".lefse.sh"))
  cmd4 = paste0("bash ", OtuwTax, ".lefse.sh")

  # copy to the folder
  cmd5=paste0("mkdir -p ", q2visdir,"/custom/lefse")
  cmd6=paste0("cp ", expdir, "/*.lefse*pdf ", paste0(q2visdir,"/custom/lefse/"))
  cmd7=paste0("cp ", expdir, "/*.lefse.results.tsv ", paste0(q2visdir,"/custom/lefse/"))
  
  # zip results
  cmd8=paste("zip -r", paste0(f, ".zip"), q2visdir)

  cmdlist = c(cmdlist, cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, "\n")

}

write(cmdlist, "launch_lefse_analysis.sh")

print("Review and Launch: bash launch_lefse_analysis.sh")
