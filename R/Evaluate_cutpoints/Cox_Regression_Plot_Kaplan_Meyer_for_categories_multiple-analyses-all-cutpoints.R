
# cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/ev-cutpoints/june-15-2021/pitts/lkt-pfs
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_Regression_Plot_Kaplan_Meyer_for_categories_multiple-analyses-all-cutpoints.R ./ infile.txt
#swarm -g 20 -t 4 --time 10-00:00:00 --module R/4.0.5 --partition=norm -f cmd.swarm --logdir log_swarm
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

biomarkerCountsCat <- function(biomarker) {
  result = as.data.frame(table(biomarker))
  Percent = (result$Freq*100)/length(biomarker)
  length.results = cbind(result, Percent)
  return(length.results)
}

biomarkerCounts <- function(biomarker, cutpoint) {
  length.total <- length(biomarker)
  length.lower <- length(biomarker[biomarker < cutpoint ])
  length.upper <- length.total - length.lower      # low is < cutpoint and everything else (including cutpoint) in high
  length.lower.ratio <- length.lower/length.total
  length.upper.ratio <- length.upper/length.total
    
  #length.lower.result <- paste(length.lower, "(", percent(length.lower.ratio), ")")
  #length.upper.result <- paste(length.upper, "(", percent(length.upper.ratio), ")")
  #length.results <- c(length.total, length.lower.result, length.upper.result)
  return(min(length.lower.ratio, length.upper.ratio))
}



# Get cox model with coefficients and result table

resultTableSurvival <- function(biomarker, survival, event, df, x) {
  
  length.results <-  biomarkerCountsCat(x)
  print(length.results)
  
  length.results <-  lapply(length.results, function(x) as.character(x))
  #print(length.results)
  
  length.results <-  do.call(c,length.results)
  #print(length.results)

  res.cox <- do.call(
    coxph,
    list(formula = Surv(survival, event) ~ x, data = df)
  )
  
  model <- summary( res.cox )
  #print(summary( res.cox ))
  coef <- model$coefficients
  HR <- signif( coef[ 2 ], digits = 3 )
  names(HR) = "HR"
  q <- 1-(1-95/100)/2
  z <- qnorm( q )
  HR.lower <- signif( exp( coef[1] - z * coef[3] ), digits = 3 )
  HR.upper <- signif( exp( coef[1] + z * coef[3] ), digits = 3 )
  CI <- paste(" ( ", HR.lower, "-", HR.upper, " ) " )
  names(CI) = "CI"
  p <- signif( model$sctest[ "pvalue" ], digits = 3 )
  
  result.table <- c(length.results, HR, CI, p)
  result.table <- data.frame(result.table)
  Estimate = rownames(result.table)
  result.table = cbind(Estimate, result.table)
  colnames(result.table) <- c("Estimate", "Result")
  #print(result.table)
  return(result.table)
}


CoxAndKaplanMeyerMethod <- function(df, time, event, biomarker) {
    
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  for(cutpoint in sort(unique(vector.biomarker))){
  
      if(biomarkerCounts(vector.biomarker, cutpoint) < 0.15){  # 15% percent minimum samples in both group
        next
      }
      
      category <- ifelse( vector.biomarker < cutpoint, "low", "high" )
      x <- ifelse( vector.biomarker < cutpoint, 0, 1 )
      category.df <- data.frame( df, category, x )

      result.table <- resultTableSurvival(vector.biomarker, vector.survival, vector.event, category.df, x)
      write.csv(result.table, file=paste0("Cutpoint-", cutpoint, "-Results.csv"), row.names=F, quote=F)

      restext = result.table[c("HR", "pvalue"), ]
      restext$Estimate[2] = 'p'
      restext = paste(restext$Estimate, restext$Result , sep=' = ', collapse='\n')
      restext = paste(paste("Cutpoint", cutpoint , sep=' = '), restext, sep='\n')
      print(restext)

    
      fit <- do.call(
        survfit,
        list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
      )
      
      # http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
      ggsurvfit <- ggsurvplot(fit, pval = F, risk.table = TRUE, risk.table.col = "strata", title = biomarker, ggtheme=theme_survminer(base_family = "Helvetica")) 
      # center align labels
      ggsurvfit$plot <- ggsurvfit$plot + theme(plot.title = element_text(hjust = 0.5)) 
      ggsurvfit$table <- ggsurvfit$table + theme(plot.title = element_text(hjust = 0.5)) 
      # add annotation in plot
      ggsurvfit$plot <- ggsurvfit$plot+ 
                  ggplot2::annotate("text", x = Inf, y = Inf, vjust = 1, hjust = 1,
                                    label = restext, size = 4) # x = 0, y = 0.2

      #ggsave(paste(biomarker, "Cutpoint", cutpoint,"KaplanMeierPlot.png",sep="_"), plot = print(ggsurvfit), dpi=600)
      ggsave(paste(biomarker, "Cutpoint", cutpoint,"KaplanMeierPlot.eps",sep="_"), plot = print(ggsurvfit), dpi=600)
      ggsave(paste(biomarker, "Cutpoint", cutpoint,"KaplanMeierPlot_.pdf",sep="_"), plot = print(ggsurvfit), dpi=600, onefile=F) # else produces an additional first blank page
      

    } # end loop over cutpoint

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
