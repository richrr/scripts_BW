#!/usr/bin/env Rscript

# doesn't need JAMS installed.

# sinteractive --mem=150g -N 10 --time=10:00:00
# ml R/3.6.0
# Rscript /data/rodriguesrr/scripts/R/kraken_classify_sequences_from_file.R -f rep-seqs_100.fasta -p rep-sequences-100


suppressPackageStartupMessages(library(RCurl))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(benchmarkme))
#####################################
# Define System-specific Functions ##
#####################################
if ((.Platform$OS.type) != "unix"){
    stop("JAMS only works on UNIX. Install Linux and try again.")
}
#Get slurm job ID
slurmjobid <- as.character(Sys.getenv("SLURM_JOB_ID"))

#Decide which kind of system you are on.
if (nchar(slurmjobid) < 3){
   print("You are not on a Slurm Workload Cluster")
   #Define appropriate functions for non-slurm system
   detectBatchCPUs <- function() {
        ncores <- detectCores()
        if (is.na(ncores)) {
            stop("Could not determine how many CPUs you have. Aborting.")
        }
        return(ncores)
    }

    detectAvailRAM <- function(){
        totmembytes <- as.numeric(get_ram())

        return(totmembytes)
    }

} else {
    print(paste("You are on a Slurm Workload Manager Cluster under jobID", slurmjobid))
    #Define appropriate functions for slurm system
    detectBatchCPUs <- function() {
        slurmjobid <- as.character(Sys.getenv("SLURM_JOB_ID"))
        ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

        if (is.na(ncores)) {
            #Try plan B
            sacctraw <- system2("sacct", args = c("-j", slurmjobid, "-X", "-P"), stdout = TRUE)
            jobinforaw <- sacctraw[2]
            jobinfoheaders <- sacctraw[1]
            jobinfo <- unlist(strsplit(jobinforaw, split="\\|"))
            names(jobinfo) <- unlist(strsplit(jobinfoheaders, split="\\|"))
            ncores <- as.integer(jobinfo["AllocCPUS"])
            print(jobinfo)
            if (is.na(ncores)) {
                stop("Could not determine how many CPUs you have. Aborting.")
            }
        }

        return(ncores)
    }

    detectAvailRAM <- function(){
        mempercpu <- as.integer(Sys.getenv("SLURM_MEM_PER_CPU"))
        mempernode <- as.integer(Sys.getenv("SLURM_MEM_PER_NODE"))
        cpuspertask <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

        if (!(is.na(mempernode))){
            totmem <- mempernode
        } else {
            totmem <- mempercpu * cpuspertask
        }

        totmembytes <- totmem * 1000000

        return(totmembytes)
    }
}

############################
## Define other functions ##
############################
filetype <- function(path){
    f <- file(path)
    ext <- summary(f)$class
    close.connection(f)
    ext
}

# get path of running script
getScriptPath <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
    } else {
        return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
}



#Define defaults
defopt <- list()
defopt$outdir <- getwd()
defopt$fasta <- "sequences.fasta"
defopt$threads <- detectBatchCPUs()
defopt$dbdir <- "/data/MicrobiomeCore/JAMSdb/JAMSdbJan2021Full_96Gbk2db/krakendb"  #/data/Trinchieri_lab/JAMSdb/JAMSdbApr2020_96Gbk2db/krakendb
defopt$prefix <- "MyResults"
defopt$krakenconfidencescore <- 0

option_list <- list(
    make_option(c("-f", "--fasta"), action="store", default=defopt$fasta, type="character", help="Input multifasta file"),

    make_option(c("-d", "--dbdir"), default=NULL, action="store", type="character", help = str_c("path to JAMS database directory")),

    make_option(c("-k", "--krakenconfidencescore"), default=defopt$krakenconfidencescore, type = "numeric", action = "store", help = str_c("confidence score threshold for kraken2 classificaiton (default: ", defopt$krakenconfidencescore, "). Only change this if you really know what you are doing.")),

    make_option(c("-p", "--prefix"), default=defopt$prefix, action="store", type="character", help = str_c("output prefix (default: ", defopt$prefix, ")")),

    make_option(c("-x", "--threads"), default=defopt$threads, action="store", type="numeric", help = str_c("number of threads (default: ", defopt$threads, ")"))

)

# parse the options
args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(OptionParser(option_list = option_list), args)
opt <- merge.list(opt, defopt)
#Cheap trick
opt$analysis <- "isolate"
opt$testdependencies <- FALSE
opt$contigs <- NULL
opt$host <- "none"

#####################
## Set environment ##
#####################

#Get Script path
opt$bindir <- getScriptPath()

#Fix path relativity
fixrelpath <- function(JAMSpath = NULL){
    require(R.utils)
    if (!(isAbsolutePath(JAMSpath))){
        fixedpath <- getAbsolutePath(JAMSpath)
    } else {
        fixedpath <- JAMSpath
    }
    #Chomp a "/" from the end of paths
    fixedpath <- gsub("/$", "", fixedpath)

    return(fixedpath)
}
for (pathtofix in c("fasta", "outdir", "dbdir")){
    if (!is.null(opt[[pathtofix]])){
        opt[[pathtofix]] <- fixrelpath(opt[[pathtofix]])
    }
}

# give help if needed input option not provided
if (is.null(opt$dbdir)) {
    print("You must supply the JAMS database directory.")
    parse_args(OptionParser(option_list=option_list), c("-h"))
    q()
} else {
    #If db directory does exist, does it have the appropriate structure
    #Stop if unreasonable settings
    if (!(dir.exists(opt$dbdir))){
        flog.info("JAMS db directory not found. You must supply one. Aborting now.")
        q()
    } else {
        #DB path was supplied. Now validate the structure.
        #suppressPackageStartupMessages(library(JAMS))
        #detect resources
        opt$totmembytes <- detectAvailRAM()
        flog.info(paste("You have ~", round((opt$totmembytes)/1000000000, 1), "Gigabytes of RAM memory"))
        flog.info(paste("You have", opt$threads, "CPUs to run on"))
        #opt <- check_resources(opt = opt)
        opt$abort = FALSE
    }
}

#Abort if all is not good.
if (opt$abort == TRUE){
    flog.info("Aborting now.")
    q()
}


#library(JAMS)
#save.image()

opt$workdir <- opt$outdir
opt$workingkrakendb <- defopt$dbdir  # opt$krakendb
minlength <- 31
fastafile <- opt$fasta

renameheaders = FALSE
#Make sure youre in the right directory
setwd(opt$workdir)


filter_sequence_by_length <- function(sequence = NULL, minlength = 0, maxlength = Inf){
    sequence_lengths <- lapply(1:length(sequence), function(x){ length(sequence[[x]]) })
    seqsIwant <- which(sequence_lengths > minlength & sequence_lengths < maxlength)
    sequence <- sequence[seqsIwant]

    return(sequence)
}


kraken_classify_taxonomy <- function(opt = NULL, fastafile = NULL, confidence = 0.00001){

    #Default is kraken2
    #Count kmers
    flog.info(paste("Counting k-mers of query sequences with kraken2 and using confidence score", as.character(confidence)))
    krakenargs <- c("--db", opt$workingkrakendb, "--threads", opt$threads, "--confidence", confidence, fastafile)
    print(krakenargs)
    kraken2taxid <- system2("kraken2", args = krakenargs, stdout = TRUE, stderr = FALSE)
    kraken2taxid <- strsplit(kraken2taxid, split = "[\t]", fixed = FALSE)
    k2out <- plyr::ldply(kraken2taxid, rbind)
    colnames(k2out) <- c("ClassFlag", "Sequence", "Taxid", "Length", "kmers")
    k2out[] <- lapply(k2out, as.character)
    JAMStaxtablefile <- file.path(opt$workingkrakendb, "JAMS_taxtable.tsv")
    if (file.exists(JAMStaxtablefile)){
        JAMStaxtable <- read.table(file = JAMStaxtablefile, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = NULL, header = TRUE)
    } else {
        #Fall back on generic taxonomy table and warn user
        flog.info("JAMS taxonomy table not found. Falling back on generic JAMS taxtable.")
        data(JAMStaxtable)
    }
    JAMStaxtable[] <- lapply(JAMStaxtable, as.character)
    krakendf <- k2out[, c("Sequence", "Taxid")]
    krakendf <- left_join(krakendf, JAMStaxtable)
    #Fill in taxids which are NOT in the database with unclassifieds
    unclassdf <- JAMStaxtable[which(JAMStaxtable$Taxid == 0), 2:ncol(JAMStaxtable)]
    krakendf[which(is.na(krakendf$Domain == TRUE)), colnames(unclassdf)] <- unname(unclassdf)

    return(krakendf)
}


get_sequence_stats<-function(sequences=NULL){
    sequencenames<-names(sequences)
    sequencelengths<-getLength(sequences)
    sequenceGCs<-sapply(1:length(sequences), function(x){ GC(sequences[[x]]) })

    seqstats<-data.frame(Sequence=sequencenames, Length=sequencelengths, GC=sequenceGCs)

    return(seqstats)
}


library(seqinr)
library(plyr)
mycontigs <- read.fasta(file = fastafile, seqtype = "DNA", forceDNAtolower = FALSE)
mycontigs <- filter_sequence_by_length(sequence = mycontigs, minlength = minlength)
#Write a temporary file to the system so that kraken can classify the contigs
write.fasta(sequences = mycontigs, names = names(mycontigs), nbchar = 80, file.out = "tempseqstoclass.fa")

#Kraken classify contigs
krakendf <- kraken_classify_taxonomy(opt = opt, fastafile = "tempseqstoclass.fa", confidence = opt$krakenconfidencescore)
file.remove("tempseqstoclass.fa")

#Filter out any vertebrate DNA if exists
#opt$contigsdata <- subset(krakendf, Kingdom != "k__Metazoa")
seqstats <- get_sequence_stats(sequences = mycontigs)
krakendf <- left_join(krakendf, seqstats, by = "Sequence")

outfn <- paste(opt$prefix, "tsv", sep = ".")

write.table(krakendf, file=outfn, col.names=TRUE, row.names=FALSE, sep="\t")

lktdf = krakendf[,c("Sequence", "LKT")]
colnames(lktdf) = c("#OTUID",	"taxonomy")
write.table(lktdf, file=paste0(opt$prefix,"-taxa.txt"), col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
