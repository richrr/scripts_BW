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

resultTableSurvival <- function(biomarker, survival, event, df, x) {
  
  length.results <-  biomarkerCounts(biomarker)
  print(length.results)
  
  
  if(FALSE){
  length.results <-  lapply(length.results, function(x) as.character(x))
  #print(length.results)
  
  length.results <-  do.call(c,length.results)
  #print(length.results)
  }


  #---- flipping the level order in coxph to match the order of levels and results with the KM plot ----#
  # by default:
      ## the levels are based on alphabetical sorting
      ## the first level is treated as "reference".
  # reorder the levels:
      ##  so the the last item becomes reference
  df$x = factor(df$x)
  original_levels_in_categ = levels(df$x)
  number_of_levels_in_categ = length(levels(df$x))
  last_categ_after_sort_is_ref = levels(df$x)[number_of_levels_in_categ]
  df$x = relevel(df$x, ref = last_categ_after_sort_is_ref)
  new_levels_in_categ = levels(df$x)


  res.cox <- do.call(
    coxph,
    list(formula = Surv(survival, event) ~ x, data = df)
  )
  
#### to do:
## generate the stats and plots to check for coxph assumptions

  model <- summary( res.cox )
  #print(summary( res.cox ))
  #coef <- model$coefficients

  
  
  if(FALSE){  # number_of_levels_in_categ == 2
      HR <- signif( coef[ 2 ], digits = 3 )
      names(HR) = "HR"
      q <- 1-(1-95/100)/2
      z <- qnorm( q )
      HR.lower <- signif( exp( coef[1] - z * coef[3] ), digits = 3 )
      HR.upper <- signif( exp( coef[1] + z * coef[3] ), digits = 3 )
      CI <- paste0("(", HR.lower, " - ", HR.upper, ")" )
      names(CI) = "CI"
      p <- signif( model$sctest[ "pvalue" ], digits = 3 ) # Score (logrank) test
      
      print(paste0(length.results[1], ": n=", length.results[3]))
      print(paste0(length.results[2], ": n=", length.results[4]))
      print(paste0(HR, " ", CI, ", P=", p))
      
      result.table <- c(length.results, HR, CI, p)
      result.table <- data.frame(result.table)
      Estimate = rownames(result.table)
      result.table = cbind(Estimate, result.table)
      colnames(result.table) <- c("Estimate", "Result")
      #print(result.table)
      return(result.table)
  }


  #------ a generic block so it works for any # of levels ------#

print(model$conf.int)
print("Here")
print(model$coefficients)
print("There")
  resout <- signif( model$conf.int[, c(1,3,4), drop=F] , digits = 3 )
  pval = signif( model$coefficients[, "Pr(>|z|)", drop=F], digits = 3 )
  print(dim(resout))
  print(dim(pval))
  resout = cbind(resout, pval)
  colnames(resout) = c("HR", "HR.lower.95", "HR.upper.95", "Pr(>|z|).predictor.var.has.relation.w.response.var.in.the.model") 

  print(resout)
  print("There again")
  HRinfo = resout[, "HR", drop=F]
  HR_string = paste0( colnames(HRinfo), ":\n",  
                                          paste(rownames(HRinfo), HRinfo[,1], sep="=", collapse='\n')
                    )
  

  Score_logrank_test_pval <- signif( model$sctest[ "pvalue" ], digits = 3 )
  Likelihood_ratio_test_pval <- signif( model$logtest[ "pvalue" ], digits = 3 )
  pval_string = paste0("P-value:\nScore_(logrank)_test_p=", Score_logrank_test_pval, "\nLikelihood_ratio_test_p=", Likelihood_ratio_test_pval)

  write.csv(resout, file="Results-stats.csv", row.names=F, quote=F)
  cat(pval_string, file = "Results-stats.csv", append=T, sep = '\n')

  # return so that we can get HR and the approp pval
  out_string = paste( HR_string, pval_string, sep= "\n")
  return(out_string)
  

}


CoxAndKaplanMeyerMethod <- function(df, time, event, biomarker) {
    
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  
  category <- vector.biomarker
  x <- vector.biomarker
  category.df <- data.frame( df, category, x )

  restext <- resultTableSurvival(vector.biomarker, vector.survival, vector.event, category.df, x)
  
  
  #restext = result.table[c("HR", "pvalue"), ]
  #restext$Estimate[2] = 'p'
  #restext = paste(restext$Estimate, restext$Result , sep=' = ', collapse='\n')
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
                                label = restext, size = 2) # x = 0, y = 0.2

  #ggsave(paste(biomarker, "KaplanMeierPlot.png"), plot = print(ggsurvfit), dpi=600)
  #ggsave(paste(biomarker, "KaplanMeierPlot.eps"), plot = print(ggsurvfit), dpi=600)
  #ggsave(paste(biomarker, "KaplanMeierPlot_.pdf"), plot = print(ggsurvfit), dpi=600, onefile=F) # else produces an additional first blank page
  
  pdf(paste(biomarker, "KaplanMeierPlot_.pdf"), onefile=F)
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
