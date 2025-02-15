# ml R/4.0.5
## only get output
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_Regression_Plot_Kaplan_Meyer_for_categories_multiple-analyses.R ./ infile.txt > log-famous-test.txt

#ml R/4.0.5
#cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/final-AE-related-analysis/landmark-300-OS
#Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_Regression_Plot_Kaplan_Meyer_for_categories_multiple-analyses.R ./PFS/G_1_4/ infile.txt 


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

library(forcats)

source("/data/rodriguesrr/scripts/R/Evaluate_cutpoints/util-functions-Evaluate_cutpoints.R")

# ------------------------------------------------------------------------------

# Get biomarker counts

biomarkerCounts <- function(biomarker) {
  result = as.data.frame(table(biomarker))
  Percent = (result$Freq*100)/length(biomarker)
  length.results = cbind(result, Percent)
  write.csv(length.results, file="Input-stats.csv", row.names=F, quote=F)
  return(length.results)
}



# Get cox model with coefficients and result table

resultTableSurvival <- function(biomarker, survival, event, df, x, biomarkername) {
  
  input.table <-  biomarkerCounts(biomarker)
  print(input.table)
  
  

  # flipping the level order in coxph to match the order of levels and results with the KM plot
  #------------------------------------
  # by default:
      ## the levels are based on alphabetical sorting
      ## the first level is treated as "reference".
  # reorder the levels:
      ##  so the the last item becomes reference
  # df$x = factor(df$x)
  # original_levels_in_categ = levels(df$x)
  # number_of_levels_in_categ = length(levels(df$x))
  # last_categ_after_sort_is_ref = levels(df$x)[number_of_levels_in_categ]
  # df$x = relevel(df$x, ref = last_categ_after_sort_is_ref)
  # new_levels_in_categ = levels(df$x)
  #------------------------------------
  df$x = forcats::fct_rev(df$x) # reverses the order of levels and basically ends up doing the same as above
  
  res.cox <- do.call(
    coxph,
    list(formula = Surv(survival, event) ~ x, data = df)
  )
  

  # generate the stats and plots to check for coxph assumptions (method 1)
  # the solid line is a smoothing spline fit to the plot, with the dashed lines representing a +/- 2-standard-error band around the fit.
  # Note that, systematic departures from a horizontal line are indicative of non-proportional hazards, since proportional hazards assumes that estimates β1,β2,β3
  # do not vary much over time. 
  # From the graphical inspection, if there is no pattern with time, the assumption of proportional hazards appears to be supported.
  plot_scaled_Schoenfeld_residuals_with_time(res.cox, biomarkername)

  model <- summary( res.cox )
  #print(summary( res.cox ))
  
  out_string =  format_results(model)
  return(out_string)
 

}



CoxAndKaplanMeyerMethod <- function(df, time, event, biomarker) {
    
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  
  category <- vector.biomarker
  x <- vector.biomarker
  category.df <- data.frame( df, category, x )

  restext <- resultTableSurvival(vector.biomarker, vector.survival, vector.event, category.df, x, biomarker)
  print(restext)
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  # plots to check for coxph assumptions (method 2)
  # parallel lines means proportional hazards assumption is met.
  # intersecting lines means variable violates proportional hazards assumption.
  # log(-log(S)) method, using Kaplan-Meier estimator
  plot_log_neg_log_S_vs_log_time_using_Kaplan_Meier_estimator(category.df, category, biomarker)

  # http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
  ggsurvfit <- ggsurvplot(fit, pval = F, risk.table = TRUE, risk.table.col = "strata", title = biomarker, ggtheme=theme_survminer(base_family = "Helvetica")) 
  # center align labels
  ggsurvfit$plot <- ggsurvfit$plot + theme(plot.title = element_text(hjust = 0.5)) 
  ggsurvfit$table <- ggsurvfit$table + theme(plot.title = element_text(hjust = 0.5)) 
  # add annotation in plot
  ggsurvfit$plot <- ggsurvfit$plot+ 
              ggplot2::annotate("text", x = Inf, y = Inf, vjust = 1, hjust = 1,
                                label = restext, size = 2) # x = 0, y = 0.2

  #ggsave(paste(biomarker, "KaplanMeierPlot.png"), plot = print(ggsurvfit), dpi=600)
  #ggsave(paste(biomarker, "KaplanMeierPlot.eps"), plot = print(ggsurvfit), dpi=600)
  #ggsave(paste(biomarker, "KaplanMeierPlot_.pdf"), plot = print(ggsurvfit), dpi=600, onefile=F) # else produces an additional first blank page
  
  pdf(paste0(biomarker, "_KaplanMeierPlot.pdf"), onefile=F)
  plot = print(ggsurvfit)
  dev.off()
}






# ------------------------------------------------------------------------------
# Functions needed for analysis

# Function performed for default analyses.

defaultAnalyses <- function(df, time, event, biomarker) {
  setwd(biomarker)
  CoxAndKaplanMeyerMethod(df, time, event, biomarker)
  setwd('..')
}


mainFunction <- function(df, time, event, biomarker) {
  defaultAnalyses(df, time, event, biomarker)
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
biomarkers <- setdiff(colsCommaSeparated, c(time, event))
biomarkers #= c("ESR1")

# Step 8. Creating directories tree, preparing data & performing analyses.

setwd(currentDate)
for (b in biomarkers) {
  print(b)
  dir.create(b)
  timeLen <- length(time)
  eventLen <- length(event)
  if ((timeLen == 0) || (eventLen == 0)) {
    print("Please correct time or event variables")
  } else if ((timeLen == 1) && (eventLen == 1)) {
    print("One time and event variable - results will be produced in biomarker folder")
    mainFunction(df, time, event, b)
  } else if ((timeLen > 1) | (eventLen > 1)) {
    print("Please provide one time and one event variable")
  }
}
#----------------------------------------------
