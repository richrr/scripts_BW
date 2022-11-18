#!/usr/bin/env Rscript

## needs JAMS installed 

# from John McCulloch
# usage: ./kraken_classify_sequences_from_file -f rep-sequences-200.fasta -p rep-sequences-200 


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
defopt$dbdir <- "/data/Trinchieri_lab/JAMSdb/JAMSdbApr2020_96Gbk2db/krakendb"
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
        suppressPackageStartupMessages(library(JAMS))
        #detect resources
        opt$totmembytes <- detectAvailRAM()
        flog.info(paste("You have ~", round((opt$totmembytes)/1000000000, 1), "Gigabytes of RAM memory"))
        flog.info(paste("You have", opt$threads, "CPUs to run on"))
        opt <- check_resources(opt = opt)
    }
}

#Abort if all is not good.
if (opt$abort == TRUE){
    flog.info("Aborting now.")
    q()
}


library(JAMS)
#save.image()

opt$workdir <- opt$outdir
opt$workingkrakendb <- opt$krakendb
minlength <- 31
fastafile <- opt$fasta

renameheaders = FALSE
#Make sure youre in the right directory
setwd(opt$workdir)
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

outfn <- paste(opt$prefix, "xlsx", sep = ".")

write.xlsx(krakendf, file=outfn, col.names=TRUE, row.names=FALSE, borders="all", colWidths="auto")
