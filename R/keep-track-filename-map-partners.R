args = commandArgs(trailingOnly=TRUE)


infile = args[1]
shotgunMapfile = args[2]
qiime2Mapfile = args[3]
genephenoMapfile = args[4]

# usage:
# cd /data/rodriguesrr/Marie_Vetizou/MV023M/analysis/transnet-FF-FM/all-w-shotgun/corr/net/v1/m-s-p/separateFishers/2.5percFDR.ExpressGr20.FM.degs/pass0/summary
# Rscript /data/rodriguesrr/scripts/R/keep-track-filename-map-partners.R corrs.FM.deg.csv /data/rodriguesrr/Marie_Vetizou/MV023M/data/JAMS_output_JM/PD1_MarieMouse_Inulin/Tables/FeatTable.Map.tsv /data/rodriguesrr/Marie_Vetizou/MV023M/analysis/16S/q2d2/exported_201908211423/asv_biom_taxonomy.tsv /data/rodriguesrr/Marie_Vetizou/MV023M/analysis/transnet/host-labels-three-coumn.txt


# /data/rodriguesrr/Marie_Vetizou/MV023M/analysis/16S/q2d2/blastnames/qiime2_asv_top_blast.csv
#blastMapfile = args[5]
# top blast hit
# dff = read.delim(blastMapfile, header=T, check.names=F, sep='\t')
# head(dff)



# get the file
df = read.csv(infile, header=T, check.names=F, stringsAsFactors = FALSE)
#head(df)
dim(df)



# shotgun map
shotgunnamespace = read.delim(shotgunMapfile, header=T, check.names=F, sep='\t')
#head(shotgunnamespace)

# get the default qiime2 taxa
dfq = read.delim(qiime2Mapfile, header=T, check.names=F, sep='\t')
# add dummy namespace
dfq$Type = "ASV"
# arrange as per the other namespace
dfq = dfq[, c("Type", "#OTUID", "taxonomy")]
colnames(dfq) = c("Type", "ID", "Description")
#head(dfq)

namespace = rbind(shotgunnamespace, dfq)

# get the gene pheno map (this file was created to fit the above file formats)
dfg = read.delim(genephenoMapfile, header=T, check.names=F, sep='\t')
namespace = rbind(namespace, dfg)
rownames(namespace) = namespace[,"ID"]

# full namespace map file
write.table(namespace, "namespace.tsv", quote=F, row.names=F, sep='\t')



# map partners
p1 = cbind( df[,"partner1", drop=F] , namespace[as.vector(df[,"partner1"]), c("Type", "Description")] )
colnames(p1) = paste0("p1.", colnames(p1))
#head(p1, 25)

p2 = cbind( df[,"partner2", drop=F] , namespace[as.vector(df[,"partner2"]), c("Type", "Description")] )
colnames(p2) = paste0("p2.", colnames(p2))
#head(p2, 25)

p = cbind(p1, p2)
dim(p)
#write.table(p, "partners.tsv", quote=F, row.names=F, sep='\t')



# write output
outdf = cbind(df, p)
write.table(outdf, paste0(infile, ".w.map.namespace.tsv"), quote=F, row.names=F, sep='\t')



print("Done!")
