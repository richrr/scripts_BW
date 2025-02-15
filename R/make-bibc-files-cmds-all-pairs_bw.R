#library(igraph)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

# ml R/4.4.1 
# cd /data/rodriguesrr/PD1_Milano/Sep-2024/Longitudinal_Multiomics_Project_LMP/analys/comp_corr/DCR_Y_merge_time/p1_cp1_cfdr1/per_analysis/net/bibc
# Rscript /data/rodriguesrr/scripts/R/make-bibc-files-cmds-all-pairs_bw.R /data/rodriguesrr/PD1_Milano/Sep-2024/Longitudinal_Multiomics_Project_LMP/analys/comp_corr/DCR_Y_merge_time/p1_cp1_cfdr1/per_analysis/net/net_FolChMedian_FCAnalys_1_CorrAnalys_1___DCR_indiv-pval_0.3_comb-pval_0.05_comb-fdr_0.25_.csv /data/rodriguesrr/PD1_Milano/Sep-2024/Longitudinal_Multiomics_Project_LMP/analys/comp_corr/mw_sp_comp-output-node-info.txt edges

edge_file = args[1] # full path
attrib_file = args[2] # two column file containing "name" and "Class" (otu or gene)
outstr = args[3] # this arg is used to prefix the "pickle" file created from the edges as well as to name the bibc cmds file

nodesDf <- read.delim(attrib_file, header = TRUE, stringsAsFactors = FALSE, sep='\t')
#head(nodesDf)
dim(nodesDf)

pairs_seen_so_far = c()
cmds = c("ml python", paste0("python /data/rodriguesrr/scripts/python/NN_bibc/import_network_data.py --input ", edge_file, " --source partner1 --target partner2 --outputfileprefix ", outstr))

# optional, only there in case you want to double check or trouble shoot later
results_file_builder_string = "#paste"

for(col in colnames(nodesDf)[2]){  # colnames(nodesDf)[-1]
  print(col)
  subdf = nodesDf[,c(colnames(nodesDf)[1],col)]
  #print(head(subdf))
  print(dim(subdf))

  grps = unique(as.vector(subdf[,2]))
  grps = grps[grps != "Ignore"]
  print(grps)

  pairs = t(combn(grps, 2))
  print(pairs)

  for(r in 1:nrow(pairs)){
    fr = paste(pairs[r,1], pairs[r,2], sep='__and__')
    rr = paste(pairs[r,2], pairs[r,1], sep='__and__')
    if(!(fr %in% pairs_seen_so_far) && !(rr %in% pairs_seen_so_far)){
        print(fr)
        pairs_seen_so_far = c(pairs_seen_so_far, fr, rr)

        nodeotu = subdf[which(subdf[, col] == pairs[r,1]), ]
        nodegene = subdf[which(subdf[, col] == pairs[r,2]), ]
        nodeignore = subdf[which(subdf[, col] != pairs[r,1] & subdf[, col] != pairs[r,2]), ]

        print(nrow(nodeotu))
        print(nrow(nodegene))
        print(nrow(nodeignore))
        print(sum(nrow(nodeotu), nrow(nodegene), nrow(nodeignore)))

        dir.create(fr)

        nodeotu$Class = "micro"
        nodegene$Class = "gene"
        nodeignore$Class = "ignore"

        BIGDF = rbind(nodeotu, nodegene, nodeignore)
        print(head(BIGDF))
        colnames(BIGDF) = c("name", "Type", "Class")

        onodefile = paste0(fr,'/',fr,'-nodes.csv')
        write.csv(BIGDF, paste0(onodefile, "-w-header.csv"), quote=F, row.names=F)
        write.table(BIGDF[ , c(1,3)], onodefile, quote=F, row.names=F, col.names=F, sep=",")  # the file that bibc code needs

        cmd1 = paste0("python /data/rodriguesrr/scripts/python/NN_bibc/BiBC_degree_calculator.py edges.pickle --bibc_groups node_types --bibc_calc_type bibc --node_map ", onodefile, " --node_groups gene micro --log") 

        cmds = c(cmds, cmd1)
        
        #print(cmds)
        
        results_file_builder_string = paste(results_file_builder_string, paste0(onodefile, "_edges_degree_BiBC_gene_micro.tsv"))

    }
  }

}

results_file_builder_string = paste(results_file_builder_string, "> results.tsv")

cmds = c(cmds, results_file_builder_string, "#run after all the above have completed; # convert results from score to perc rank, and then later you can merge w nodes.csv", "Rscript /data/rodriguesrr/scripts/R/parse-multiple-nolan-bibc-results.R edges_degree_BiBC_gene_micro.tsv all-bibc-results.csv")
write(cmds, paste0(outstr,"-run-bibc-cmds.txt"), sep = "\n")

print("Done!")
