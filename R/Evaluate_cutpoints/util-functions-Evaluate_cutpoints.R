


# https://www.sthda.com/english/wiki/cox-model-assumptions#testing-proportional-hazards-assumption
plot_scaled_Schoenfeld_residuals_with_time = function(res.cox, biomarker){
  
  info_string = "#If the test is not statistically significant for each of the covariates, and the global test is also not statistically significant, we can assume the proportional hazards."
  test.ph <- cox.zph(res.cox)
  print(test.ph)
  write.csv(test.ph$table, file="Test_proportional-hazards_assumption-stats.csv", row.names=T, quote=F)
  cat(info_string, file = "Test_proportional-hazards_assumption-stats.csv", append=T, sep = '\n')

  ssplot = survminer::ggcoxzph(test.ph)
  pdf(paste0(biomarker, "_scaled_Schoenfeld_residuals_with_time.pdf"), onefile=F)
  plot = print(ssplot)
  dev.off()

}


#------ a generic block so it works for any # of levels ------#
# formatting the results from the model
format_results = function(model){
  #print(model$conf.int)
  #print(model$coefficients)
  resout <- signif( model$conf.int[, c(1,3,4), drop=F] , digits = 3 )
  pval = signif( model$coefficients[, "Pr(>|z|)", drop=F], digits = 3 )
  resout = cbind(resout, pval)
  colnames(resout) = c("HR", "HR.lower.95", "HR.upper.95", "Pr(>|z|).predictor.var.has.relation.w.response.var.in.the.model") 
  print(resout)

  HRinfo = resout[, "HR", drop=F]
  HR_string = paste0( colnames(HRinfo), ":\n",  
                                          paste(rownames(HRinfo), HRinfo[,1], sep="=", collapse='\n')
                    )
  

  Score_logrank_test_pval <- signif( model$sctest[ "pvalue" ], digits = 3 )
  Likelihood_ratio_test_pval <- signif( model$logtest[ "pvalue" ], digits = 3 )
  pval_string = paste0("P-value:\nScore_(logrank)_test_p=", Score_logrank_test_pval, "\nLikelihood_ratio_test_p=", Likelihood_ratio_test_pval)
  print(pval_string)

  write.csv(resout, file="Results-stats.csv", row.names=T, quote=F)
  cat(pval_string, file = "Results-stats.csv", append=T, sep = '\n')

  # return so that we can get HR and the approp pval
  out_string = paste( HR_string, pval_string, sep= "\n")
  return(out_string)
  
}


# https://rpubs.com/kaz_yos/kleinbaum-ph-assumption (ignore since the method is no longer supported after rms v4.2)
# https://stackoverflow.com/questions/22422687/r-survival-package-plotting-log-logsurvival-against-logtime
plot_log_neg_log_S_vs_log_time_using_Kaplan_Meier_estimator =  function(catdf, category, biomarker){

  SurvObj <- Surv( catdf$dmfs_time, catdf$dmfs_event )

  levelsinCateg = levels(factor(category))
  numbOfLevels = length(levelsinCateg)
  
  ## Treatment variable (x)
  pdf(paste0(biomarker, "_log_neg_log_S_vs_log_time.pdf"), onefile=F)
  plot( survfit(SurvObj ~ category), fun="cloglog" , ylab = "-log(log(survival))" , xlab = "log time", col = seq(1:numbOfLevels))
  legend("topleft", legend = levelsinCateg, text.col=seq_along(levelsinCateg), lty=1, col = 1:numbOfLevels ) #, pch=1)
  dev.off()

 }
