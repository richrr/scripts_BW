## specific version where landmark is done with "OS" and for deciding AE status but the KM and Cox is done using PFS




# (1) Landmark analysis (30, 60, 90, 180, 270, 300, 365)
    #- for say landmark 365 OS days
          # keep the patients with OS >= 365 (magical number/landmark)
          # those patients that have developed AE before 1 yr keep the actual time to AE and AE indicator
          # TA (time to AE) >= 365 is considered as no AE
          # do KM/Cox-regress with PFS

# (2) Time-dependent covariate




# ml R/4.0.5
# cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/landmark-analysis/june-25-2021
# Rscript /data/rodriguesrr/scripts/R/OS-landmark-PFS-KM.Cox-landmark-analysis.R ./ infile-for-specific-landmark-analysis.txt

# R version 4.0.5
args <- commandArgs(trailingOnly = TRUE)
library("survival")
library("survminer")
library("dplyr") 
library("gtsummary")
library("tidyverse")

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
# the row names are lost
table <- read.table(file, sep = "\t", header = TRUE, na.strings=c("", "NA", "NaN", "N_A"), row.names=1)
indf <- data.frame(table)
rownames_indf = rownames(indf)
indf = apply(indf, 2, function(x){as.numeric(as.vector(x))})
indf <- data.frame(indf)
rownames(indf) = rownames_indf
head(indf)


# Step 3. Get date and prepare a folder for analysis.

currentDate <- gsub(':',"_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
setwd("analysis-results")
dir.create(currentDate)
setwd(currentDate)



##----------------------------------------
# Landmark analysis
##----------------------------------------
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#Part_2:_Landmark_Analysis_and_Time_Dependent_Covariates
# it seems this tutorial keeps the AE as it is post landmark date (since in her case all AE happened before landmark date)
# Important: if the patient did not have AE by landmark date, they should be considered no AE (regardless whether they develop AE or not later)

landmark_analysis = function(file, BMT, lnd, substract){
  #print(head(BMT))
  print(lnd)
  
  # -------
  # step 1: keep only patients who survived more than this time
  # -------
  lm_dat <- 
  BMT %>% 
  filter(Overall_Survival_days >= lnd) 
  
  
  print(dim(BMT))
  print(dim(lm_dat))
  #print(lm_dat)
  
  
  # -------
  # step 2: only consider AE if happened before landmark time, everything else is no AE, 
  # -------
  # patients (who developed AE or not) after landmark are considered no AE.
  lm_dat$tmp_landmarkDeltaA  =  lm_dat[, "AE_indicator"]
  pats_to_be_considered_noAE = rownames(lm_dat[which(lm_dat[,"time_to_AE"] >=lnd), ])
  print(lm_dat[pats_to_be_considered_noAE,])
  lm_dat[pats_to_be_considered_noAE, "tmp_landmarkDeltaA"] = 0
  
  # same as below: (returns TRUE or FALSE)
  # patients who did not develop AE until landmark time are considered no AE (regardless if they develop it later or not)
  lm_dat$landmarkDeltaA = {{lm_dat[,"AE_indicator"] == 1} & {lm_dat[,"time_to_AE"] < lnd}}
  
  print(lm_dat)
  
  lm_dat = lm_dat[,c("PFS_days","Progression_Event","time_to_AE","landmarkDeltaA")]
  colnames(lm_dat) = c("T1", "delta1", "TA", "landmarkDeltaA")
    
  #print(lm_dat)
  
  

  # substract landmark from T1 (OS or PFS time)
  if(substract){
      lm_dat <- 
      lm_dat %>% 
      mutate(
        lm_T1 = T1 - lnd
        )
  } else{  # what GT was asking (only makes a difference in the plot, not the HR and pvalue)
      lm_dat <- 
      lm_dat %>% 
      mutate(
        lm_T1 = T1
        )
  }
  
  print(lm_dat)


  # -------
  # Kaplan Meier
  # -------
  lm_fit <- survfit(Surv(lm_T1, delta1) ~ landmarkDeltaA, data = lm_dat)
  print(lm_fit)
  pdf(paste(file,lnd,"substract-landmark-days",substract,"OS-landmark-PFS.KM.pdf",sep="_"), onefile=F)
  p = ggsurvplot(lm_fit, data = lm_dat,             # data used to fit survival curves.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = F )
  print(p)
  dev.off()
  

  # -------
  # Cox regression
  # -------  
  res.cox = coxph(
          Surv(lm_T1, delta1) ~ landmarkDeltaA, 
          data = lm_dat
          ) #%>% 
          #gtsummary::tbl_regression(exp = TRUE)
    
  
  
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
  
  Landmark = lnd
  names(Landmark) = "Landmark"
  result.table <- c(Landmark, HR, CI, p)
  result.table <- data.frame(result.table)
  Estimate = rownames(result.table)
  result.table = cbind(Estimate, result.table)
  colnames(result.table) <- c("Estimate", "Result")
  print(result.table)
  write.csv(result.table, file=paste(file,lnd,"substract-landmark-days",substract,"OS-landmark-PFS.Cox-regress.csv",sep="_"), row.names=F, quote=F)

}



for(lnd in c(60, 80, 90, 120, 150, 180, 210, 240, 270, 300, 365)){    
    landmark_analysis(file, indf, lnd, TRUE)
    landmark_analysis(file, indf, lnd, FALSE)
}


print("Done!")



q()
