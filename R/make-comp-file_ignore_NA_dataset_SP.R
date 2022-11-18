args = commandArgs(trailingOnly=TRUE)

setwd(file.path("./"))

corrdir = args[1]
phenodir = args[2]
microdir = args[3]

corr_analy = args[4]
pheno_comp = args[5]
micro_comp = args[6]

basefolder = paste(corrdir, corr_analy, "p1_cp1_cfdr1.ignoreNaDatasets", sep="/")
merged_corr_file = paste(basefolder, "merged_FolChMedian_merged-parallel-output.csv", sep="/")

per_analys_folder = paste(basefolder, "per_analysis", sep="/")
consis_corr_file = paste(per_analys_folder, "Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv", sep="/")
fc_comp_file = paste0(per_analys_folder, "/", pheno_comp, ".AND.", micro_comp, "-medFC-comp-file.csv")

out_net_dir = paste(per_analys_folder, "net", sep="/")
dir.create(file.path(out_net_dir))



read_merged_file = function(dir, analysis){
  infile = list.files(paste(dir, analysis, sep="/"),
                      pattern = "FolChMedian_merged-parallel-output.csv$",
                      recursive = T,
                      full.names = T)
  print(infile)
  df = read.csv(infile, check.names=F, header=T, row.names=1)
  #print(head(df))
  return(df)
}


get_Median_FC = function(indf){
  FoldChangeColnames = grep("FolChMedian", colnames(indf), value=T)
  interestedFoldChangeData = indf[,FoldChangeColnames]
  interestedFoldChangeData = as.matrix(interestedFoldChangeData)
  interestedFoldChangeData = apply(interestedFoldChangeData,2,function(x){as.numeric(as.vector(x))})
  combinedFoldChange = apply(interestedFoldChangeData, 1, function(x){round(median(x, na.rm = TRUE), 3)})
  indf$geneName = rownames(indf)
  indf = cbind(indf,combinedFoldChange)
  outdf = indf[, c("geneName", "combinedFoldChange")]
  #print(head(outdf))
  return(outdf)
}


phenodf = read_merged_file(phenodir, pheno_comp)
phenofc = get_Median_FC(phenodf)

microdf = read_merged_file(microdir, micro_comp)
microfc = get_Median_FC(microdf)

outfc = rbind(phenofc, microfc)
colnames(outfc) = c("geneName", "Analys 1 FolChMedian ABX vs No_ABX")
head(outfc)

write.csv(outfc, fc_comp_file, quote=F, row.names=F)


setwd(file.path(out_net_dir))

create_net_cmd = paste0("Rscript /data/rodriguesrr/scripts/R/create-network_bw.R --file ", merged_corr_file, " --group GROUP --consistent ", consis_corr_file," --foldchange ", fc_comp_file," --indivPvalCutoff 1 --combPvalCutoff 0.15 --analysisfc 'Analys\ 1\ ' --analysiscorr 'Analys\ 1\ ' --combFDRCutoff 1 --numbDataFromFCfile --output ./net_")

print(create_net_cmd)

system(create_net_cmd)
