# Analysis

# Step 0. Please set a working directory.

mainDir <- ("/Users/magda/Documents/ev-cutpoints-multiple-analyses/")
setwd(mainDir)

# Creates a directory for all future analyses.

if ( dir.exists("analysis-results")) {
  print("Directory already exists")
} else {
  dir.create("analysis-results")
}

# Step 1. Choose a file. Please change the name below.

file <- ("~/Documents/ev-cutpoints-multiple-analyses/breastcancer_GSE2034.txt")

# Step 2. Read table. Please choose a type of separator: "", ",", ";", "\t".

table <- read.table(file, sep = "\t", header = TRUE)
df <- data.frame(table)

# Step 3. Get date and prepare a folder for analysis.

currentDate <- format(Sys.time(), "%a %b %d %Y %X")
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
biomarkers <- c("ESR1.205225_at", "ESR1.211233_x_at", "ESR1.211234_x_at", "ESR1.211235_s_at", "ESR1.211627_x_at", 
                "ESR1.215551_at", "ESR1.215552_s_at", "ESR1.217163_at", "ESR1.217190_x_at", 
                "PGR.208305_at", "ERBB2.210930_s_at", "ERBB2.216836_s_at")

# Step 8. Creating directories tree, preparing data & performing analyses.

setwd(currentDate)
for (b in biomarkers) {
  dir.create(b)
  timeLen <- length(time)
  eventLen <- length(event)
  if ((timeLen == 0) || (eventLen == 0)) {
    print("Please correct time or event variables")
  } else if ((timeLen == 1) && (eventLen == 1)) {
    print("One time and event variable - results will be produced in biomarker folder")
    mainFunction(setCutpoint, df, time, event, b)
  } else if ((timeLen > 1) | (eventLen > 1)) {
    print("Please provide one time and one event variable")
  }
}