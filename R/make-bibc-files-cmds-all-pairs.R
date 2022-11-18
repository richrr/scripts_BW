library(igraph)


args = commandArgs(trailingOnly=TRUE)

# ml R/3.6.0
# cd /data/rodriguesrr/FMT_pittsburg/TransNet/analys/network/lmfc/v2-mcode-cluster-8/w-interpro-only-connected-to-lkt/expanded/bibc/
# Rscript /data/rodriguesrr/scripts/R/make-bibc-files-cmds-all-pairs.R clust8-w-interpro-w-expanded-edges.csv node-w-group-info.txt clust8-w-ipr-expnd

edge_file = args[1] # full path
attrib_file = args[2] # two column file containing "name" and "Class" (otu or gene)
outstr = args[3]

nodesDf <- read.delim(attrib_file, header = TRUE, stringsAsFactors = FALSE, sep='\t')
#head(nodesDf)
dim(nodesDf)

pairs_seen_so_far = c()
cmds = c()

for(col in colnames(nodesDf)[-1]){
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

        nodeotu$Class = "otu"
        nodegene$Class = "gene"
        nodeignore$Class = "ignore"

        BIGDF = rbind(nodeotu, nodegene, nodeignore)
        print(BIGDF)
        colnames(BIGDF) = c("name", "Type", "Class")

        onodefile = paste0(fr,'/',fr,'-nodes.csv')
        write.csv(BIGDF, onodefile, quote=F, row.names=F)

        cmds = c(cmds, paste("/local/cluster/R-3.6.0/bin/Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/prep-for-betw-central.R", edge_file, onodefile, onodefile))

        cmds = c(cmds, paste0('SGE_Batch -c "/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/betweennessadapter node_type microbe degs ', onodefile,'.graphml > bibc-', fr, '-results.txt" -m 50G -F 50G -q transkingdom -M rodrrich@oregonstate.edu -r ', paste0("log-", fr)))

        #print(cmds)

    }
  }

}

cmds = c(cmds, "#run after all the above have completed; get the rel, rank, and perc rank files", "#Rscript ~/Morgun_Lab/richrr/scripts/R/parse-multiple-bibc-results.R -results.txt all-bibc-results.csv")
write(cmds, paste0(outstr,"-","bibc-cmds.txt"), sep = "\n")

print("Done!")
