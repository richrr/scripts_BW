

#### this only plots the pcoa but does not provide beta diversity pvalues
#### when finding the map file, code also needs to account for the fact that the input file might start with /data/....
 

# Aug 10 2022
# ml R/4.1
# R --version
# R version 4.1.3 (2022-03-10) -- "One Push-Up"
# .libPaths()    /spin1/home/linux/rodriguesrr/R/4.1/library/

# cd /data/MicrobiomeCore/analysis/GTAM45_20220714
# usage: Rscript /data/rodriguesrr/scripts/R/qiime2R_parse_qza.R D15/core_metrics_202208082212/jaccard_pcoa_results.qza Group pdfdump

# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(qiime2R)

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(qiime2R)

library(stringr)


pcoa_res_infile = args[1]

pcoaRes<-read_qza(pcoa_res_infile)    #"D10/core_metrics_202208082211/unweighted_unifrac_pcoa_results.qza"

DirStruc = str_split(pcoa_res_infile, "/")[[1]]

shannon<-read_qza( paste0( paste( head(DirStruc, -1), collapse="/"), "/shannon_vector.qza") )$data %>% rownames_to_column("SampleID")    # "D10/core_metrics_202208082211/shannon_vector.qza"

metadata<-read.table(paste0("map-files/", DirStruc[1], ".txt"), header=T)  # read_q2metadata  # "map-files/D10.txt"


categCol = args[2]

outpdfname = paste0(args[3], "/", str_replace_all(pcoa_res_infile, "/", "_") , ".pdf")


pcoaRes$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes_string(x="PC1", y="PC2", size="shannon", color=categCol, shape=categCol)) +  # https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="Antibiotic Usage") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Group")
  ggsave(outpdfname, height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches




q()

# http://joey711.github.io/phyloseq/import-data.html
# rarefied table?
# taxa bar plot; # 
# alpha diversity metrics # http://joey711.github.io/phyloseq/plot_richness-examples.html
# ordination http://joey711.github.io/phyloseq/plot_ordination-examples.html
#library(phyloseq)


#physeq<-qza_to_phyloseq(
#    features="../D10/table.qza",
#    tree="../D10/rooted-tree.qza",
#    "../D10/taxonomy.qza",
#    metadata = "../map-files/D10.txt"
#    )

#physeq





#### do not use ####
heatmap = function(){


      SVs<-read_qza("../D10/table.qza")$data

      taxonomy<-read_qza("../D10/taxonomy.qza")$data

      SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent


      SVsToPlot<-  
        data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
        rownames_to_column("Feature.ID") %>%
        arrange(desc(MeanAbundance)) %>%
        #top_n(100, MeanAbundance) %>%
        pull(Feature.ID) #extract only the names from the table
        

        

      SVs %>%
        as.data.frame() %>%
        rownames_to_column("Feature.ID") %>%
        gather(-Feature.ID, key="SampleID", value="Abundance") %>%
        mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
        group_by(SampleID, Feature.ID) %>%
        summarize(Abundance=sum(Abundance)) %>%
        left_join(metadata) %>%
        mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
        left_join(taxonomy) %>%
        #mutate(Feature=paste(Feature.ID, Taxon)) %>%
        #mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
        #ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
        
        ggplot(aes(x=SampleID, y=Feature.ID, fill=NormAbundance)) +
        geom_tile() +
        facet_grid(~`Group`, scales="free_x") +
        theme_q2r() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        scale_fill_viridis_c(name="log10(% Abundance)")
        ggsave("heatmap.pdf", height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches


}
