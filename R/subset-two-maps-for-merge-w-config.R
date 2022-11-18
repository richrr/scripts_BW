args = commandArgs(trailingOnly=TRUE)
configFile = args[1]
tissue1_tissue2 = args[2]
list_tiss1 = args[3]
list_tiss2 = args[4]
outdir = args[5]
expt = args[6]

source(configFile)

# sinteractive --cpus-per-task=4 --mem=75G --time=10:00:00
# ml R/3.6.0
# cd /data/rodriguesrr/Stephanie_Prescott/12KM.Females/data/microbe-rnaseq_pheno_wk15.WD
# Rscript /data/rodriguesrr/scripts/R/subset-two-maps-for-merge-w-config.R config-files/12K-stool-rnaseq_pheno-config.txt stool-fisher5percpancreas /data/rodriguesrr/Stephanie_Prescott/12KM.Females/analysis/microbe/comp/merged/all-tissues-consis_fc-microbes.txt /data/rodriguesrr/Stephanie_Prescott/12KM.Females/analysis/rnaseq/comp/merged.WD.pancreas/p1_cp1_cfdr0.15/per_analysis/Analys1-consis.csv-comb-pval-output.csv.-p.1.combpval.0.05.combfdr.1.cutoff.csv.consis_genes.csv /data/rodriguesrr/Stephanie_Prescott/12KM.Females/analysis/microbe-rnaseq/corr/ 12K




map1 = read.delim(mapf1, header=T, check.names=F, as.is=T, colClasses = c("character"))
map2 = read.delim(mapf2, header=T, check.names=F, as.is=T, colClasses = c("character"))

smap1 = map1[which(map1[,grp1Col] == grp1), ]
smap2 = map2[which(map2[,grp2Col] == grp2), ]

#smap1
#smap2

# if any of the "group" columns are present in the other map file it will create .x and .y, but not otherwise
# e.g. map1 has Group and map2 has Treatment and Group, will result in Group.x and Group.y
# and you would have to account for either "Group" or "Group.x" or "Group.y"
# however, prefixing them before merge makes colnames unique and you only have to work with the prefixed colnames
colnames(smap1) = paste0(colnames(smap1), "_OmicsX")
colnames(smap2) = paste0(colnames(smap2), "_OmicsY")

merge1Col = paste0(merge1Col, "_OmicsX")
merge2Col = paste0(merge2Col, "_OmicsY")

grp1Col = paste0(grp1Col, "_OmicsX")
grp2Col = paste0(grp2Col, "_OmicsY")

mapSampId1 = paste0(mapSampId1, "_OmicsX")
mapSampId2 = paste0(mapSampId2, "_OmicsY")


# not doing merge (all=True) in maps because here we are using mouse ids for merging and the microbe/pheno have multiple timepoints for the same mice
# so it will create a problem. it may work if we create a new (column in the) map file with appropriate sample names to be used 
# like for combining pheno and rnaseq (where we only took the same time point and only added the NA in rnaseq for the samples that had pheno but no rnaseq)
# but it will be too much to make new map files for 5 tissues for microbes and 2 rnaseq
# so keep only common samples
# combined map file
smerged = merge(smap1, smap2, by.x= merge1Col, by.y= merge2Col)


# prefix the numerical animal ids so R doesn't treat them as index later
MID = paste0("Mouse_",as.vector(unlist(smerged[, merge1Col])))


NewGROUP = paste(smerged[, grp1Col] , smerged[, grp2Col], sep=".AND.")
smerged = cbind(MID, smerged)
smerged = cbind(smerged, NewGROUP)
#smerged


gname = paste0(grp1, ".AND.", grp2)
mapOutDir = paste0(tissue1_tissue2, "/", expt, "/map-file/")
dir.create(mapOutDir, recursive=TRUE)
mapfname = paste0(mapOutDir,gname, ".txt")
write.table(smerged, mapfname,sep="\t", quote=F, row.names=F)



data1 = read.csv(dataf1, check.names=F, header=T, row.names=1)
data1 = data1[, as.vector(unlist(smerged[, mapSampId1]))]
colnames(data1) = paste0("Mouse_",as.vector(unlist(smerged[, merge1Col])))
#head(data1)


data2 = read.csv(dataf2, check.names=F, header=T, row.names=1)
data2 = data2[, as.vector(unlist(smerged[, mapSampId2]))]
colnames(data2) = paste0("Mouse_",as.vector(unlist(smerged[, merge1Col])))
#head(data2)


# combined data file
rdata = rbind(data1, data2)
DID = rownames(rdata)
rdata = cbind(DID, rdata)

dataOutDir = paste0(tissue1_tissue2, "/", expt, "/data-file/")
dir.create(dataOutDir, recursive=TRUE)
datafname = paste0(dataOutDir,gname, ".csv")
write.csv(rdata, datafname, quote=F, row.names=F)


# make an analysis file:
H = c("ColumID","Group_A","Group_B","Detail","AnalysisType")
D = c("NewGROUP", gname, "", "", "correlation")

analysOutDir = paste0(tissue1_tissue2, "/", expt, "/analys_file/")
dir.create(analysOutDir, recursive=TRUE)
analyfname = paste0(analysOutDir,gname, ".tsv")
write(H, file = analyfname, sep = "\t", ncolumns = length(H))
write(D, file = analyfname, sep = "\t", ncolumns = length(H), append = T)



# build the command for perform analysis.

outputfname = paste0(outdir, "/", expt, "/",tissue1_tissue2,"/", gname, "/sp_")
cmd = paste("Rscript /data/rodriguesrr/scripts/R/perform-analyses_bw.R",
                      datafname,
                      "--lists", list_tiss1, list_tiss2,
                      "--mapFile",
                      mapfname,
                      "--mapFileSampleIdColName MID --AnalysToDoList", analyfname,
                      "--comparMethod NA --correlMethod sp --dataFileSymbolColName DID --pairedInfoColumn NA --output", outputfname, "--stringtorelativize NA --multicore --cores 4\n")
print(cmd)
#system(cmd)
cat(cmd,file="cmd-list.swarm",append=TRUE)

q()
