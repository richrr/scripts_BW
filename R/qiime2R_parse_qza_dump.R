
##### 
# this takes a */core_metrics_*/beta_sig_*/*_bsig.qzv file
# like /data/MicrobiomeCore/analysis/GTAM45_20220714/D15/core_metrics_202208082212/beta_sig_Group_20220808221611/unweighted_unifrac_anosim_bsig.qzv 
# and creates a filename */core_metrics_*/*_pcoa_results.qza to be processed

# this allows adding beta diversity pvalues into the pdf with plot

# potential downside: while anosim and permanova will both call this script,
# but each copy of the pdf plot will have either anosim or permanova
# the csv file from export_q2vis_qzv_to_pdf will have both anosim and permanova in one file if needed

#####


 
# Aug 10 2022
# ml R/4.1
# R --version
# R version 4.1.3 (2022-03-10) -- "One Push-Up"
# .libPaths()    /spin1/home/linux/rodriguesrr/R/4.1/library/

# cd /data/MicrobiomeCore/analysis/GTAM45_20220714
# usage: Rscript /data/rodriguesrr/scripts/R/qiime2R_parse_qza_dump.R /data/MicrobiomeCore/analysis/GTAM45_20220714/D15/core_metrics_202208082212/beta_sig_Group_20220808221611/unweighted_unifrac_anosim_bsig.qzv /data/MicrobiomeCore/analysis/GTAM45_20220714/ Group pdfdump pdfdump//D15/core_metrics_202208082212/beta_sig_Group_20220808221611/unweighted_unifrac_anosim_bsig.qzv/index.html_df.csv

# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(qiime2R)

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(qiime2R)

library(stringr)

library(gridExtra)



infile = args[1]   # "/data/MicrobiomeCore/analysis/GTAM45_20220714/D15/core_metrics_202208082212/beta_sig_Group_20220808221611/unweighted_unifrac_anosim_bsig.qzv"
infile = str_replace(infile, "beta_sig_.*/", "")

multivarTest = ""
if(grepl("anosim", infile)){
    infile = str_replace(infile, "_anosim_bsig.qzv", "_pcoa_results.qza")
    multivarTest = "anosim"
} else if (grepl("permanova", infile)){
    infile = str_replace(infile, "_permanova_bsig.qzv", "_pcoa_results.qza")
    multivarTest = "permanova"
} 



pcoa_res_infile = infile  
pcoa_res_infile     # "/data/MicrobiomeCore/analysis/GTAM45_20220714/D15/core_metrics_202208082212/unweighted_unifrac_pcoa_results.qza"
pcoaRes<-read_qza(pcoa_res_infile)    


#### make sure current directory (args[2]) has "/" at the end
args[2]

removeCd = str_replace(pcoa_res_infile, args[2], "")  # "D10/core_metrics_202208082211/unweighted_unifrac_pcoa_results.qza"
DirStruc = str_split( removeCd , "/")[[1]]  

shannon<-read_qza( paste0( paste( head(DirStruc, -1), collapse="/"), "/shannon_vector.qza") )$data %>% rownames_to_column("SampleID")    # "D10/core_metrics_202208082211/shannon_vector.qza"

metadata<-read.table(paste0("map-files/", DirStruc[1], ".txt"), header=T)  # read_q2metadata  # "map-files/D10.txt"


categCol = args[3]

outpdfname = paste0(args[4], "/", str_replace_all(removeCd, "/", "_") , "_", multivarTest, ".pdf")


betasignif = read.csv(args[5], header=T, check.names=F)
#betasignif

# https://stackoverflow.com/questions/41164675/adding-a-table-of-values-below-the-graph-in-ggplot2
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 10,
                     padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(betasignif, rows=NULL, theme=tt)



plt1 <- pcoaRes$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes_string(x="PC1", y="PC2", size="shannon", color=categCol, shape=categCol)) +  # https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="Antibiotic Usage") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Group") + ggtitle(removeCd)
#  ggsave(outpdfname, height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

pdf(outpdfname)
#plt1
grid.arrange(plt1, tbl, nrow = 2)
dev.off()


q()
