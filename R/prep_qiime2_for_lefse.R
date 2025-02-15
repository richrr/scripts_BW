args = commandArgs(trailingOnly=TRUE)

# otu table with taxonomy is args[1]
# metadata file is args[2]
# comma separated list of (0-based) column numbers to use a class and subclass. only 2 elements allowed.
# map file with this format sampleID,class,subclass should be entered as 1,2 for args[3]

data = read.delim(args[1], header=2, skip=1, as.is = T, check.names = F, row.names=1, stringsAsFactors=F)
#head(data)
#dim(data)

# move the taxonomy column to the start
newdata = data[,c(ncol(data),1:(ncol(data)-1))]
colnames(newdata)[1] = 'ID'

# replace ';' with '|' as lefse likes it
#newdata$ID = gsub('D_.__', '', newdata$ID)
newdata$ID = gsub('; ', '|', newdata$ID)

# change rownames to include taxonomy and the exact OTU ID; delete the taxonomy (ID) column
newdata$ID = paste(newdata$ID , row.names(newdata), sep='_')
rownames(newdata) = newdata$ID
newdata = newdata[, -1]
#head(newdata)



map = read.delim(args[2], header=T, as.is = T, check.names = F, row.names=1, stringsAsFactors=F)
#head(map)

mapcols = as.numeric(as.vector(strsplit(args[3], ",")[[1]]))
mapredux = map[, c( mapcols ), drop=F]
SampleID = rownames(mapredux)
mapredux = cbind(mapredux, SampleID)
#head(mapredux)


tnewdata = t(newdata)
#head(tnewdata)

dim(tnewdata)
dim(mapredux)

fdf = t(merge(mapredux, tnewdata, by="row.names"))
fdf = fdf[-1,]
head(fdf)

write.table(fdf, paste0(args[1], ".for.lefse.tsv"), sep="\t", col.names=F, quote=F)




# build shell command file
#cmd = "#!/bin/bash\nexport TMPDIR=/lscratch/\\$SLURM_JOB_ID\n\nmodule load lefse\n"
cmd = "\n#run 'unset R_LIBS' with quotes on the terminal \nml R/4.3.0\nml python/3.10\nmodule load lefse/1.0.8\n"
if(length(mapcols) == 1){
    cmd = paste0(cmd, "format_input.py ", args[1], ".for.lefse.tsv" , " ", args[1], ".lefse.in", " -c 1 -u 2 -o 1000000")
} else if(length(mapcols) == 2){
    cmd = paste0(cmd, "format_input.py ", args[1], ".for.lefse.tsv" , " ", args[1], ".lefse.in", " -c 1 -s 2 -u 3 -o 1000000")
} else {
    print("Maximum of 2 columns allowed from metadata.")
    quit()
}

##### created one header file that will be used to write to the lefse output  /data/rodriguesrr/scripts/lefse_header_file.tsv

cmd = paste0(
cmd,"\n",
"run_lefse.py ", args[1], ".lefse.in", " ", args[1], ".lefse.res.tsv", "\n",
"plot_res.py ",  args[1], ".lefse.res.tsv", " ", args[1], ".lefse.pdf --dpi 600 --format pdf" , "\n",
"plot_cladogram.py ",  args[1], ".lefse.res.tsv", " ", args[1], ".lefse.cladogram.pdf --dpi 600 --format pdf" , "\n",
"cat /data/rodriguesrr/scripts/lefse_header_file.tsv ", args[1], ".lefse.res.tsv", " > ", args[1], ".lefse.results.tsv"
)
#print(cmd)


shellfile = paste0(args[1], ".lefse.sh")
fileConn<-file(shellfile)
writeLines(cmd, fileConn)
close(fileConn)

print("Manually run these commands on biowulf")

system(paste0("cat ",shellfile))

#execcmd = paste0("sbatch --cpus-per-task=2 --time=2:00:00 --mail-type=END ", shellfile)
#print(execcmd)
#system(execcmd)
