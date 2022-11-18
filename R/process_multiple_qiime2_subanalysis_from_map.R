
#ml R/3.6.0
#cd /data/MicrobiomeCore/analysis/GTAM45_20220714
#Rscript /data/rodriguesrr/scripts/R/process_multiple_qiime2_subanalysis_from_map.R /data/MicrobiomeCore/analysis/GTAM45_20220714/map-files/map.xlsx /data/MicrobiomeCore/analysis/GTAM45_20220714/map-files /data/MicrobiomeCore/runs/M03213_20220801/q2d2_out 7500 Group
#bash launch_qiime2_analysis.sh

# loop in R to create a map file for each sheet/comparison

library(readxl)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

infile = args[1] # map excel file with sheets representing subset of samples for analysis
mapdir = args[2] # dir to create the maps in
sourcedir = args[3] # q2d2_out/ from runs folder
fixed_rarefaction_threshold = args[4]
beta_sign_col = args[5]

indir = getwd()  # e.g. cd /data/MicrobiomeCore/analysis/GTAM42_20220318/


cmdlist = c()

shits = excel_sheets(infile)
for(f in shits){
  
  print(f)
  df = data.frame(read_excel(infile, sheet=f))
  print(head(df))
  print(dim(df))

  omapfile = paste0(mapdir, "/", f, ".txt")
  write.table(df, omapfile, row.names=F, quote=F, sep="\t")

  cmd1=paste0("cd ", indir)

  # q2d2_extract_subset.pl
  workdir = paste0(indir, "/", f)
  cmd2=paste0("perl /data/MicrobiomeCore/scripts/q2d2/q2d2_extract_subset.pl -s ", sourcedir," -d ", workdir, " -m ", omapfile, " -p KEEP,yes")

  # q2d2_cd.pl [fixed threshold?]
  cmd3=paste0("perl /data/MicrobiomeCore/scripts/q2d2/q2d2_cd.pl -w ", workdir, " -m ", omapfile, " -d ", fixed_rarefaction_threshold, " -b ", beta_sign_col)

  cmdlist = c(cmdlist, cmd1, cmd2, cmd3, "\n")

}

write(cmdlist, "launch_qiime2_analysis.sh")

print("Review and Launch: bash launch_qiime2_analysis.sh")
