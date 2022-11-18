# Analysis
# R version 4.0.2
args <- commandArgs(trailingOnly = TRUE)
library('rmarkdown')

# only get output
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R ./ famous-infile.txt > log-famous-test.txt

# To put the output and error in the same file (assuming sh/bash)
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R ./ famous-infile.txt > log-famous.txt 2>&1

# swarm -g 10 -t 2 --time 4:00:00 --logdir swarmlogs --partition ccr,norm,quick --module R/4.0 -f swarm.sh
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R ./ OS-infile.txt
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R ./ PFS-infile.txt

# wget https://hackage.haskell.org/package/pandoc-2.11.1.1/pandoc-2.11.1.1.tar.gz
# tar xvzf pandoc-2.11.1.1.tar.gz
## /home/rodriguesrr/pandoc-2.11.1.1
# cd pandoc-2.11.1.1
# HOMEBREW_NO_AUTO_UPDATE=1 brew install haskell-stack
# stack setup
# stack install

##############
# OR did not work
# cabal update
# cabal install cabal-install
# cabal --version
# cabal install pandoc
###############


source("/data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_functions.R")
#source("/data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/ev-cutpoints/EvaluateCutpoints_multiple-analyses_functions.R")

# Step 0. Please set a working directory.

mainDir <- (args[1])
setwd(mainDir)

# Creates a directory for all future analyses.

if ( dir.exists("analysis-results")) {
  print("Directory already exists")
} else {
  dir.create("analysis-results")
}

# Step 1. Choose a file. Please change the name below.

file <- (args[2])

minimumSamps = 0
# force minimum number of samples in both categs for cutpoint
if(length(args) > 2){
  minimumSamps = as.numeric(args[3])
}

# Step 2. Read table. Please choose a type of separator: "", ",", ";", "\t".

table <- read.table(file, sep = "\t", header = TRUE, na.strings=c("", "NA", "NaN", "N_A"))
df <- data.frame(table)
#df = apply(df, 2, function(x){as.numeric(as.vector(x))})
#head(df)


# Step 3. Get date and prepare a folder for analysis.

currentDate <- gsub(':',"_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
setwd("analysis-results")
dir.create(currentDate)

# Step 5. Print all column names <- to simplify copy-pasting :)

colNames <- colnames(df)
colsCommaSeparated <- dput(colNames)

# Step 6. Choose if cutpoints should be counted or you would like to set your own cutpoint.
# Paste the value of the cutpoint to produce plots and analyses for this cutpoint.
# If setCutpoint variable is set to NULL, best cutpoint will be counted for each method:
# cutp, maxstat, ROC01, Youden for 2 groups and rolr for 3 groups
# If setCutpoint variable has one set value, histogram, ROC plot, Kaplan Meier plot and Standarized log-rank statistic plots will be produced.
# If setCutpoint variable has two set values, low vs high, low vs medium and medium vs high statistics will be performed and  plots (histogram, Kaplan-Meier) will be produced.

setCutpoint <- NULL
# Example 2: setCutpoint <- 5
# Example 3: setCutpoint <- c(5,6)

# Step 7. Choose biomarkers, time, event variables (see Step 5).

time <- c("dmfs_time")
event <- c("dmfs_event")
biomarkers <- setdiff(colsCommaSeparated, c(time, event))
biomarkers 

# Step 8. Creating directories tree, preparing data & performing analyses.

counter = 0
setwd(currentDate)
for (b in biomarkers) {
  counter = counter+1
  print(paste(counter, b))
  dir.create(b)
  timeLen <- length(time)
  eventLen <- length(event)
  if ((timeLen == 0) || (eventLen == 0)) {
    print("Please correct time or event variables")
  } else if ((timeLen == 1) && (eventLen == 1)) {
    print("One time and event variable - results will be produced in biomarker folder")
    mainFunction(setCutpoint, df, time, event, b, minimumSamps)
  } else if ((timeLen > 1) | (eventLen > 1)) {
    print("Please provide one time and one event variable")
  }

}
