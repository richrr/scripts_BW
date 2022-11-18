# Analysis
# R version 4.0.5
args <- commandArgs(trailingOnly = TRUE)
#library('rmarkdown')
library("survival")
library("survminer")


# cd /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021

# swarm -g 10 -t 2 --time 4:00:00 --logdir swarmlogs --partition ccr,norm,quick --module R/4.0.5 -f swarm.sh
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021/early_63_samps/default/fact-os/ infile.txt
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021/early_63_samps/default/fact-pfs/ infile.txt
# Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_regress_categ-or-continuous-vars-no-cutpoint.R /data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/Cox_regression_continuous/may-30-2021/early_63_samps/default/lkt-pfs/ infile.txt


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

#minimumSamps = 0
# force minimum number of samples in both categs for cutpoint
#if(length(args) > 2){
#  minimumSamps = as.numeric(args[3])
#}

# Step 2. Read table. Please choose a type of separator: "", ",", ";", "\t".

table <- read.table(file, sep = "\t", header = TRUE, na.strings=c("", "NA", "NaN", "N_A"))
df <- data.frame(table)
df = apply(df, 2, function(x){as.numeric(as.vector(x))})
df <- data.frame(df)
head(df)


# Step 3. Get date and prepare a folder for analysis.

currentDate <- gsub(':',"_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
setwd("analysis-results")
#dir.create(currentDate)

# Step 5. Print all column names <- to simplify copy-pasting :)

colNames <- colnames(df)
colsCommaSeparated <- dput(colNames)

# Step 6. Choose biomarkers, time, event variables (see Step 5).

time <- c("dmfs_time")
event <- c("dmfs_event")
biomarkers <- setdiff(colsCommaSeparated, c(time, event))
biomarkers 

# Step 8. Creating directories tree, preparing data & performing analyses.

# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#Kaplan-Meier_plot_-_ggsurvplot
pdf("Kaplan-Meier_Survival-Curve.pdf")
ggsurvplot(
    fit = survfit(Surv(dmfs_time, dmfs_event) ~ 1, data = df), 
    xlab = "Time", 
    ylab = "Survival probability")
dev.off()


# test #
#res.cox = coxph(Surv(dmfs_time, dmfs_event) ~ Age, data = df)
#pdf("Cox_Survival-Curve-test-age.pdf")
#ggsurvplot(survfit(res.cox, data = df))
#dev.off()



# http://www.sthda.com/english/wiki/cox-proportional-hazards-model
univ_formulas <- sapply(biomarkers, function(x) {    
                        f = as.formula(paste('Surv(dmfs_time, dmfs_event)~', x))

                        ## add a plot (KM)
                        #pdf(paste0("Kaplan-Meier_Survival-Curve_",x,".pdf"))
                        #g <- ggsurvplot(
                        #    fit = surv_fit(f, data = df),  #https://github.com/kassambara/survminer/issues/283
                        #    xlab = "Time", 
                        #    ylab = "Survival probability",
                        #    title = x)
                        #print(g)
                        #dev.off()
                        
                        f
})

        
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

print("Done")


q()
