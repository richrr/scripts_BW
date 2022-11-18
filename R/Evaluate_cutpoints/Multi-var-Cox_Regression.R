
## only get output
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Multi-var-Cox_Regression.R ./ infile.txt Age,BMI 66.02,36.34 

# ------------------------------------------------------------------------------
# Libraries

#library("OptimalCutpoints")
library("ggplot2")
library("survival")
library("survMisc")
library("scales")
library("survminer")
#library("maxstat")
library("plotly")
#library("rolr")



# ------------------------------------------------------------------------------

# Get biomarker counts

biomarkerCounts <- function(biomarker) {
  result = as.data.frame(table(biomarker))
  Percent = (result$Freq*100)/length(biomarker)
  length.results = cbind(result, Percent)
  return(length.results)
}


# Get cox model with coefficients and result table

resultTableSurvival <- function(df, biomarker, survival, event) {
  # https://win-vector.com/2018/09/01/r-tip-how-to-pass-a-formula-to-lm/
  form = as.formula (paste("Surv(survival, event)", paste(biomarker, collapse = " + "), sep = " ~ "))  
  res.cox <- coxph(form, data = df)
  model <- summary( res.cox )
  return(summary( res.cox ))
}


CoxMethod <- function(df, time, event, biomarker, cutpoints) {
    
  vector.biomarker <- df[, biomarker, drop=F ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  
  print(head(vector.biomarker))
  print(cutpoints)

  TALLDF = ''
  for(i in 1:ncol(vector.biomarker)){
    tmp = ifelse( vector.biomarker[,i] < cutpoints[i], 0, 1 )
    TALLDF = cbind(TALLDF, tmp)
  }

  TALLDF = TALLDF[,-1,drop=F]
  TALLDF = apply(TALLDF, 2, function(x){as.numeric(x)})
  colnames(TALLDF) = biomarker
  #print(head(TALLDF))

  
  category.df <- data.frame( df[, c(time, event)], TALLDF)
  print(category.df)
  result.table <- resultTableSurvival(category.df, biomarker, vector.survival, vector.event)
  print(result.table)
  
  sink("Results.csv")
  print(result.table)
  sink()
  
  print("Done!")
  q()
    
}





# ------------------------------------------------------------------------------
# Functions needed for analysis

# Function performed for default analyses.

defaultAnalyses <- function(df, time, event, biomarker, cutpoints) {
  setwd(paste(biomarker, collapse = "__"))
  CoxMethod(df, time, event, biomarker, cutpoints)
  setwd('..')
}


mainFunction <- function(df, time, event, biomarker, cutpoints) {
  defaultAnalyses(df, time, event, biomarker, cutpoints)
}





#----------------------------------------------

# Analysis
# R version 4.0.2
args <- commandArgs(trailingOnly = TRUE)
library('rmarkdown')



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


# Step 7. Choose biomarkers, time, event variables (see Step 5).

time <- c("dmfs_time")
event <- c("dmfs_event")
#biomarkers <- c("ESR1.205225_at", "ESR1.211233_x_at", "ESR1.211234_x_at", "ESR1.211235_s_at", "ESR1.211627_x_at", 
#                "ESR1.215551_at", "ESR1.215552_s_at", "ESR1.217163_at", "ESR1.217190_x_at", 
#                "PGR.208305_at", "ERBB2.210930_s_at", "ERBB2.216836_s_at")
biomarkers <- unlist(strsplit(args[3], ",")) #setdiff(colsCommaSeparated, c(time, event))
biomarkers #= c("ESR1")
cutpoints <- as.numeric(unlist(strsplit(args[4], ",")))
cutpoints

# Step 8. Creating directories tree, preparing data & performing analyses.

setwd(currentDate)
#for (b in biomarkers) {
#  print(b)
  b = biomarkers
  dir.create(paste(b, collapse = "__"))
  timeLen <- length(time)
  eventLen <- length(event)
  if ((timeLen == 0) || (eventLen == 0)) {
    print("Please correct time or event variables")
  } else if ((timeLen == 1) && (eventLen == 1)) {
    print("One time and event variable - results will be produced in biomarker folder")
    mainFunction(df, time, event, b, cutpoints)
  } else if ((timeLen > 1) | (eventLen > 1)) {
    print("Please provide one time and one event variable")
  }
#}
#----------------------------------------------
