library(igraph)
args = commandArgs(trailingOnly=TRUE)

# this file makes a different node file for every node for which you want to run bibc
# this is useful when you want to edit the node file of all phenos to run bibc for each pheno
# the ideal input file for this code can be from the output node files -w-header.csv created from make-bibc-files-cmds-all-pairs.R
# since they already have the other node groups set to the correct type

#ml R/4.4.1 
#cd /data/rodriguesrr/PD1_Milano/Sep-2024/Longitudinal_Multiomics_Project_LMP/analys/comp_corr/DCR_Y_merge_time/p1_cp1_cfdr1/per_analysis/net/bibc
#grep REALWORLD LKT__and__REALWORLD/LKT__and__REALWORLD-nodes.csv | cut -d, -f1 > REALWORLD-nodes-of-interest.txt
#Rscript /data/rodriguesrr/scripts/R/make-bibc-cmds-per-node_bw.R /data/rodriguesrr/PD1_Milano/Sep-2024/Longitudinal_Multiomics_Project_LMP/analys/comp_corr/DCR_Y_merge_time/p1_cp1_cfdr1/per_analysis/net/net_FolChMedian_FCAnalys_1_CorrAnalys_1___DCR_indiv-pval_0.3_comb-pval_0.05_comb-fdr_0.25_.csv LKT__and__REALWORLD/LKT__and__REALWORLD-nodes.csv-w-header.csv edges REALWORLD-nodes-of-interest.txt LKT


### this code is probably different than the ../python/make_per_pheno_cmds.py ###


edge_file = args[1] 
attrib_file = args[2] # three column file containing "name", "Type" and "Class" (otu or gene). Comes from the output node files -w-header.csv created from make-bibc-files-cmds-all-pairs.R
outstr = args[3]
node_file = args[4] # this file contains the node for which you want to make new node files
other_type = args[5] # to track which is the other group type

nodes <- as.vector(unlist(read.csv(node_file, header = FALSE, stringsAsFactors = FALSE)))
print(nodes)

attributesgraph = NULL
if(grepl(".csv$", attrib_file))
  {
    attributesgraph <- read.csv(attrib_file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    attributesgraph <- read.delim(attrib_file, header = TRUE, stringsAsFactors = FALSE, sep='\t')
  }
head(attributesgraph)
dim(attributesgraph)


cmds = c("ml python", paste0("python /data/rodriguesrr/scripts/python/NN_bibc/import_network_data.py --input ", edge_file, " --source partner1 --target partner2 --outputfileprefix ", outstr))


for( node in nodes){
    print(node)
    tmp = attributesgraph
    sanitize_node = gsub("[[:space:]]",'-',node)
    sanitize_node = gsub("[[:punct:]]",'-',sanitize_node)
    
    dname = paste(sanitize_node, other_type, sep="__and__")
    
    dir.create(dname)
    outfile = paste0(dname, "/", dname, "-node.csv")
    
    for (r in 1:nrow(tmp)){
        # change all items from nodes file except the node in iteration to "ignore"
        if(tmp[r, "name"] %in% nodes){
            if(tmp[r, "name"] != node){
            	tmp[r, "Class"] = "ignore"
        	}
        }   
    }
    
    #write.csv(tmp, outfile, row.names=F, quote=F)
    write.table(tmp[ , c(1,3)], outfile, quote=F, row.names=F, col.names=F, sep=",")  # the file that bibc code needs
    print(paste0("printed to ", outfile))

    cmd1 = paste0("python /data/rodriguesrr/scripts/python/NN_bibc/BiBC_degree_calculator.py edges.pickle --bibc_groups node_types --bibc_calc_type bibc --node_map ", outfile, " --node_groups gene micro --log") 

    cmds = c(cmds, cmd1)

}


cmds = c(cmds, "#run after all the above have completed; # convert results from score to perc rank, and then later you can merge w nodes.csv", "Rscript /data/rodriguesrr/scripts/R/parse-multiple-nolan-bibc-results.R edges_degree_BiBC_gene_micro.tsv all-bibc-results.csv")
write(cmds, paste0(outstr, other_type, "-run-per-node-bibc-cmds.txt"), sep = "\n")

print("Done!")
