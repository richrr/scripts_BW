# Analysis
# R version 4.0.5
args <- commandArgs(trailingOnly = TRUE)
#library('rmarkdown')
library("survival")
library("survminer")


# cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021

## swarm -g 10 -t 2 --time 4:00:00 --logdir swarmlogs --partition ccr,norm,quick --module R/4.0.5 -f swarm.sh
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint-w-wo-logscale.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021/early_63_samps/default/fact-os/ infile.txt
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint-w-wo-logscale.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021/early_63_samps/default/fact-pfs/ infile.txt
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint-w-wo-logscale.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021/early_63_samps/default/lkt-pfs/ infile.txt

# cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/houston
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint-w-wo-logscale.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/houston/lkt-pfs/ infile.txt


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
indf <- data.frame(table)
indf = apply(indf, 2, function(x){as.numeric(as.vector(x))})
indf <- data.frame(indf)
head(indf)



# Step 3. Get date and prepare a folder for analysis.

currentDate <- gsub(':',"_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
currentDate <- gsub(' ',"_",currentDate)
setwd("analysis-results")
#dir.create(currentDate)

# Step 5. Print all column names <- to simplify copy-pasting :)

colNames <- colnames(indf)
colsCommaSeparated <- dput(colNames)

# Step 6. Choose biomarkers, time, event variables (see Step 5).

time <- c("dmfs_time")
event <- c("dmfs_event")
biomarkers <- setdiff(colsCommaSeparated, c(time, event))
biomarkers 

# Step 8. Creating directories tree, preparing data & performing analyses.

doEverything = function(df){
      print(typeof(df))
      print(class(df))
      df <- data.frame(df)
      print(head(df))
      
      # keeps giving the error if the fit formula is written inline inside the ggsurvplot
      #Error in stats::get_all_vars(.formula, data = data) : 
      #'data' must be a data.frame, environment, or list
      #Calls: doEverything ... <Anonymous> -> .pvalue -> .extract.survfit -> <Anonymous>

      # https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#Kaplan-Meier_plot_-_ggsurvplot
      fit = survfit(Surv(dmfs_time, dmfs_event) ~ 1, data = df)
      pdf("Kaplan-Meier_Survival-Curve.pdf")
      km = ggsurvplot(
          fit, data = df, 
          xlab = "Time", 
          ylab = "Survival probability")
      print(km)
      dev.off()


      # http://www.sthda.com/english/wiki/cox-proportional-hazards-model
      univ_formulas <- sapply(biomarkers, function(x) {as.formula(paste('Surv(dmfs_time, dmfs_event)~', x))})

      univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df)})

      # Extract data 
      univ_results <- lapply(univ_models,
                             function(x){ 
                                
                                # the fitted cox model
                                y <- x
                                

                                x <- summary(x)
                                varname = rownames(x$coef)
                                print(varname)
                                
                                # add a plot using cox regression fit model
                                pdf(paste0("Cox-Regress_Survival-Curve_",varname,".pdf"))
                                g <- ggsurvplot(
                                    fit = survfit(y, data = df),  #https://github.com/kassambara/survminer/issues/67
                                    xlab = "Time", 
                                    ylab = "Survival probability",
                                    title = varname)
                                print(g)
                                dev.off()


                                p.value<-signif(x$wald["pvalue"], digits=3)
                                wald.test<-signif(x$wald["test"], digits=2)
                                beta<-signif(x$coef[1], digits=3);#coeficient beta
                                HR <-signif(x$coef[2], digits=5);#exp(beta)
                                HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                                HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                                HR <- paste0(HR, " (", 
                                             HR.confint.lower, "-", HR.confint.upper, ")")
                                res<-c(beta, HR, wald.test, p.value)
                                names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                              "p.value")
                                return(res)
                                #return(exp(cbind(coef(x),confint(x))))
                               })
      res <- t(as.data.frame(univ_results, check.names = FALSE))
      res = as.data.frame(res)
      res
      write.csv(res, paste0(currentDate,"-out.csv"))
}

# unlog both (abundance and time); log one; log another; log both

if ( dir.exists("unlog")) {
  print("Directory already exists")
} else {
  dir.create("unlog")
}
setwd("unlog")
doEverything(indf)
setwd("../")


# https://stackoverflow.com/questions/25150907/log-transforming-predictor-variables-in-survival-analysis
# log variable's abundance
if ( dir.exists("logVar")) {
  print("Directory already exists")
} else {
  dir.create("logVar")
}
setwd("logVar")

varlogdf = indf[ , -which(names(indf) %in% c("dmfs_time", "dmfs_event"))]
varlogdf = log(varlogdf + 1 , 2)
varlogdf = cbind(varlogdf, indf[, c("dmfs_time", "dmfs_event")])
head(varlogdf)
doEverything(varlogdf)
setwd("../")

print("Done with and without log transform of the predictor variable")

#q()

# Following may not be needed

# log time
if ( dir.exists("logTime")) {
  print("Directory already exists")
} else {
  dir.create("logTime")
}
setwd("logTime")

timelogdf = indf[ , -which(names(indf) %in% c("dmfs_time"))]
dmfs_time = log(indf[, "dmfs_time"] + 1 , 2)
timelogdf = cbind(timelogdf, dmfs_time)
head(timelogdf)
doEverything(timelogdf)
setwd("../")



# log both
if ( dir.exists("log")) {
  print("Directory already exists")
} else {
  dir.create("log")
}
setwd("log")

vartimelogdf = indf[ , -which(names(indf) %in% c("dmfs_event"))]
vartimelogdf = log(vartimelogdf + 1 , 2)
dmfs_event = indf[, "dmfs_event"]
vartimelogdf = cbind(vartimelogdf, dmfs_event)
head(vartimelogdf)
doEverything(vartimelogdf)
setwd("../")



print("Done")


q()
