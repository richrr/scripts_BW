args = commandArgs(trailingOnly=TRUE)

map1 = read.delim(args[1], header=T, check.names=F, as.is=T, colClasses = c("character"))
map2 = read.delim(args[2], header=T, check.names=F, as.is=T, colClasses = c("character"))

grp1 = args[3]
grp2 = args[4]


smap1 = map1[which(map1[,"Group"] == grp1), ]
smap2 = map2[which(map2[,"Group"] == grp2), ]


#head(map1)
#head(map2)
#smap1
#smap2

# combined map file
smerged = merge(smap1, smap2, by.x= "Animal_Number", by.y= "AnimalNumber")
MID = paste0("Mouse_",as.vector(unlist(smerged[, "Animal_Number"])))
GROUP = paste(smerged[, "Group.x"] , smerged[, "Group.y"], sep=".AND.")
smerged = cbind(MID, smerged)
smerged = cbind(smerged, GROUP)

gname = paste0(grp1, ".AND.", grp2)
mapfname = paste0("map-file/",gname, ".txt")
write.table(smerged, mapfname,sep="\t", quote=F, row.names=F)



data1 = read.csv(args[5], check.names=F, header=T, row.names=1)
data1 = data1[, as.vector(unlist(smerged[, "ID"]))]
colnames(data1) = paste0("Mouse_",as.vector(unlist(smerged[, "Animal_Number"])))
head(data1)


data2 = read.csv(args[6], check.names=F, header=T, row.names=1)
data2 = data2[, as.vector(unlist(smerged[, "SampleID"]))]
colnames(data2) = paste0("Mouse_",as.vector(unlist(smerged[, "Animal_Number"])))
head(data2)

# combined data file
rdata = rbind(data1, data2)
DID = rownames(rdata)
rdata = cbind(DID, rdata)

datafname = paste0("data-file/",gname, ".csv")
write.csv(rdata, datafname, quote=F, row.names=F)



tissue = args[7]



# make an analysis file:
H = c("ColumID","Group_A","Group_B","Detail","AnalysisType")
D = c("GROUP", gname, "", "", "correlation")
analyfname = paste0("analys_file/",gname, ".tsv")
write(H, file = analyfname, sep = "\t", ncolumns = length(H))
write(D, file = analyfname, sep = "\t", ncolumns = length(H), append = T)


# build the command for perform analysis.

# was used for stool only.
#listmicr = paste0("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/",tissue,"/shortlisted-microbes.txt")

# for all other tissues
listmicr = "/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/comp/merged/all-tissues-shortlisted-microbes.txt"

outputfname = paste0("/data/rodriguesrr/Stephanie_Prescott/12KLMN/analys/microbe/","corr/",tissue,"/", gname, "/sp_")
cmd = paste("Rscript /data/rodriguesrr/scripts/R/perform-analyses_bw.R",
                      datafname,
                      "--lists", listmicr,
                      "/data/rodriguesrr/Stephanie_Prescott/06-08-2020/data/12KLMN-pheno-list.txt --mapFile",
                      mapfname,
                      "--mapFileSampleIdColName MID --AnalysToDoList", analyfname,
                      "--comparMethod NA --correlMethod sp --dataFileSymbolColName DID --pairedInfoColumn NA --output", outputfname, "--stringtorelativize NA")
print(cmd)
system(cmd)


q()
