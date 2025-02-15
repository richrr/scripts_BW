# ------------------------------------------------------------------------------
# Libraries

library("OptimalCutpoints")
library("ggplot2")
library("survival")
library("survMisc")
library("scales")
library("survminer")
library("maxstat")
library("plotly")
library("rolr")

# ------------------------------------------------------------------------------
# Functions for particular methods of cutpoint determination

# Youden Method

youdenMethod <- function(df, time, event, biomarker) {
  
  method <- "Youden"
  
  # Count optimal cutpoint and generate a table 
  
  optimal.cutpoint.Youden <- optimal.cutpoints(X = biomarker, status = event, tag.healthy = 0, methods = method, data = df, pop.prev = NULL, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
  Youdenptable <- summary(optimal.cutpoint.Youden)$p.table$Global
  
  summaryDf <- data.frame(Youdenptable[[1]])
  write.csv(summaryDf, file="YoudenTableOutput.csv")
  
  
  # Generate histogram and ROC plots
  
  cutpoint <- Youdenptable[[method]][[1]][[1]]
  
  histogramPlotYouden <- ggplot(df, aes(x=df[, biomarker])) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))

  histogramYoudenWithLabs <- histogramPlotYouden + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramYouden.png"))

  png(paste(biomarker, "ROCYouden.png"))
  youdenPlot <- plot(optimal.cutpoint.Youden, which = 1)
  dev.off()
  
}


# ROC Method

ROCMethod <- function(df, time, event, biomarker) {
  
  method <- "ROC01"
  
  # Count optimal cutpoint and generate a table 
  
  optimal.cutpoint.roc <- optimal.cutpoints(X = biomarker, status = event, tag.healthy = 0, methods = method, data = df, pop.prev = NULL, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
  Rocptable <- summary(optimal.cutpoint.roc)$p.table$Global
  
  summaryDf <- data.frame(Rocptable[[1]])
  write.csv(summaryDf, file="ROCTableOutput.csv")
  
  
  # Generate histogram and ROC plots
  
  cutpoint <- Rocptable[[method]][[1]][[1]]
  
  histogramPlotROC <- ggplot(df, aes(x=df[, biomarker])) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))

  histogramROCWithLabs <- histogramPlotROC + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramROC.png"))

  png(paste(biomarker, "ROCPlot.png"))
  rocPlot <- plot(optimal.cutpoint.roc, which = 1)
  dev.off()
  
}



# check minimumSamps is met in each categories
NumbSamplesCheck = function(cutpdf, sdf, minimumSamps){
    cutpointRet = ''  # the row number of value to be used as cutpoint
    for(r in 1:nrow(cutpdf)){
        cutpnow = cutpdf[[ r, 1 ]]
        #print(cutpnow)

        x <- ifelse( sdf[,3] < cutpnow, 0, 1 )
        grpA = sum(x == 0)
        grpB = sum(x == 1)
        if(min(grpA, grpB) >= minimumSamps)
          {
            cutpointRet = r
            print(paste(cutpnow, grpA, grpB))
            break
          }
    }
    return(cutpointRet)
}


# Get biomarker counts

biomarkerCounts <- function(biomarker, cutpoint) {
  length.total <- length(biomarker)
  length.lower <- length(biomarker[biomarker < cutpoint ])
  length.upper <- length.total - length.lower      # low is < cutpoint and everything else (including cutpoint) in high
  length.lower.ratio <- length.lower/length.total
  length.upper.ratio <- length.upper/length.total
  length.lower.result <- paste(length.lower, "(", percent(length.lower.ratio), ")")
  length.upper.result <- paste(length.upper, "(", percent(length.upper.ratio), ")")
  length.results <- c(length.total, length.lower.result, length.upper.result)
  return(length.results)
}

# Get biomarker counts for three groups

biomarkerCountsThree <- function(biomarker, low.cutoff, high.cutoff) {
  length.total <- length(biomarker)
  print(low.cutoff)
  print(length.total)
  length.lower <- length(biomarker[biomarker < low.cutoff])
  length.upper <- length(biomarker[biomarker > high.cutoff])
  length.medium <- length.total - length.lower - length.upper
  length.lower.ratio <- length.lower / length.total
  length.upper.ratio <- length.upper / length.total
  length.medium.ratio <- length.medium / length.total
  length.lower.result <- paste(length.lower, "(", percent(length.lower.ratio), ")")
  length.upper.result <- paste(length.upper, "(", percent(length.upper.ratio), ")")
  length.medium.result <- paste(length.medium, "(", percent(length.medium.ratio), ")")
  length.results <- c(length.total, length.lower.result, length.medium.result, length.upper.result)
  return(length.results)
}

# Get cox model with coefficients and result table

resultTableSurvival <- function(biomarker, cutpoint, survival, event, df, x) {
  
  length.results <- biomarkerCounts(biomarker, cutpoint)
  
  res.cox <- do.call(
    coxph,
    list(formula = Surv(survival, event) ~ x, data = df)
  )
  
  model <- summary( res.cox )
  #print(summary( res.cox ))
  coef <- model$coefficients
  HR <- signif( coef[ 2 ], digits = 3 )
  q <- 1-(1-95/100)/2
  z <- qnorm( q )
  HR.lower <- signif( exp( coef[1] - z * coef[3] ), digits = 3 )
  HR.upper <- signif( exp( coef[1] + z * coef[3] ), digits = 3 )
  CI <- paste(" ( ", HR.lower, "-", HR.upper, " ) " )
  p <- signif( model$sctest[ "pvalue" ], digits = 3 )
  
  result.table.col.names <- c( "Cutpoint", "Biomarker < Cutpoint", "Biomarker >= Cutpoint", "HR", "CI", "P-value" )
  result.table.row.names <- c( cutpoint, as.character(length.results[2]), as.character(length.results[3]), HR, CI, p )
  result.table <- data.frame(result.table.col.names, result.table.row.names)
  colnames(result.table) <- c("Estimate", "Result")
  return(result.table)
}

# CutP Method

cutPMethod <- function(df, time, event, biomarker, minimumSamps) {
  
  method <- "cutp"
  
  # NA are removed
  sdf = df[,c(time, event, biomarker)]
  sdf <- na.omit(sdf)
  print(paste0("Non-NA: ", nrow(sdf)))
  
  vector.biomarker <- sdf[, biomarker ]
  vector.survival <- sdf[, time ]
  vector.event <- sdf[, event ]
  
  cph1 <- do.call(
    coxph,
    list( formula = Surv( vector.survival, vector.event ) ~ vector.biomarker, data = sdf)
  )
  
  
  allCutpointsDf <- data.frame(cutp(cph1))
  if(ncol(allCutpointsDf) != 4) {print("Cutp did not return output. Moving to next"); return(1)}
  colnames(allCutpointsDf) <- c(biomarker, "U", "Q", "pvalue")
  write.csv(allCutpointsDf, file="allCutpointsOutput.csv")
  biomarkerPValueDf <- data.frame( allCutpointsDf[biomarker], allCutpointsDf$pvalue )
  colnames(biomarkerPValueDf) <- c("biomarker", "pvalue")
  biomarkerPValueDfSorted <- biomarkerPValueDf[order(biomarkerPValueDf$biomarker),]
  colnames(biomarkerPValueDfSorted) <- c(biomarker, "pvalue")
  write.csv(biomarkerPValueDfSorted, file="outputSortedbyBiomarkerValue.csv")
  
  named.biomarkerPValueDfSorted <- data.frame(biomarkerPValueDfSorted)
  colnames(named.biomarkerPValueDfSorted) <- c("biomarker", "pvalue")
  bind.cutpointsDf <- cbind(named.biomarkerPValueDfSorted$pvalue, named.biomarkerPValueDfSorted$pvalue)
  rownames(bind.cutpointsDf) <- named.biomarkerPValueDfSorted$biomarker
  rev.cutpointsDf <- t(bind.cutpointsDf)
  colnames(rev.cutpointsDf) <- rownames(bind.cutpointsDf)
  rownames(rev.cutpointsDf) <- c( ""," ")
  
  y.axisSettings <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    ticks = ""
  )
  x.axisSettings <- list(
    title = biomarker
  )
  
  heatmap <- plot_ly(
    z = data.matrix(rev.cutpointsDf),
    x = colnames(rev.cutpointsDf),
    y = rownames(rev.cutpointsDf),
    type = "heatmap",
    colorbar = list(
      xanchor = "left",
      yanchor = "top",
      ypad = 0,
      xpad = 0,
      lenmode = "pixels",
      len= 150,
      nticks=3
    )) %>% layout(
      yaxis = y.axisSettings,
      xaxis = x.axisSettings
    )
  
  htmlwidgets::saveWidget(heatmap, file = "heatmap.html")
  cutpointRet = NumbSamplesCheck(allCutpointsDf, sdf, minimumSamps)
  if(cutpointRet == '') {print("Cannot find cutpoint passing minimum number of samples criteria. Moving to next"); return(1)}
  cutpoint <- signif(allCutpointsDf[[ cutpointRet, 1 ]], digits = 4 )
  p.value <- signif(allCutpointsDf[[ cutpointRet, 4 ]], digits = 4 )
  
  category <- ifelse( vector.biomarker < cutpoint, "low", "high" )
  x <- ifelse( vector.biomarker < cutpoint, 0, 1 )
  category.df <- data.frame( sdf, category, x )

  result.table.cutp <- resultTableSurvival(vector.biomarker, cutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="CutPResults.csv")
  
  rownames(result.table.cutp) = result.table.cutp$Estimate
  restext = result.table.cutp[c("Cutpoint", "HR", "P-value"), ]
  restext$Estimate[3] = 'p'
  restext = paste(restext$Estimate, restext$Result , sep=' = ', collapse='\n')
  print(restext)


  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  
  histogramPlotCutp <- ggplot(sdf, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramCutpWithLabs <- histogramPlotCutp + labs(x=paste(biomarker))
  ggsave(paste0(biomarker, "_histogramCutP.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste0(biomarker, "_KaplanMeierPlotCutp_no_risk_table.png"))
  
  
  # http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
  #ggsurvfit <- ggsurvplot(fit, pval = TRUE, risk.table = TRUE)  #surv.col = c( "#2c3e50", "Red" ),   

  ggsurvfit <- ggsurvplot(fit, pval = F, risk.table = TRUE, risk.table.col = "strata", title = biomarker, ggtheme=theme_survminer(base_family = "Helvetica")) 
  # center align labels
  ggsurvfit$plot <- ggsurvfit$plot + theme(plot.title = element_text(hjust = 0.5)) 
  ggsurvfit$table <- ggsurvfit$table + theme(plot.title = element_text(hjust = 0.5)) 
  # add annotation in plot
  ggsurvfit$plot <- ggsurvfit$plot+ 
              ggplot2::annotate("text", x = Inf, y = Inf, vjust = 1, hjust = 1,
                                label = restext, size = 4) # x = 0, y = 0.2

  #ggsave(paste(biomarker, "KaplanMeierPlot.png"), plot = print(ggsurvfit), dpi=600)
  #ggsave(paste(biomarker, "KaplanMeierPlot.eps"), plot = print(ggsurvfit), dpi=600)
  #ggsave(paste(biomarker, "KaplanMeierPlot_.pdf"), plot = print(ggsurvfit), dpi=600, onefile=F) # else produces an additional first blank page
  
  pdf(paste0(biomarker, "_KaplanMeierPlot_.pdf"))
  print(ggsurvfit, newpage=FALSE)
  dev.off()
  
}


# Maxstat Method

maxstatMethod <- function(df, time, event, biomarker, minimumSamps) {
  
  method <- "maxstat"
  
  # NA are removed
  sdf = df[,c(time, event, biomarker)]
  sdf <- na.omit(sdf)
  print(paste0("Non-NA: ", nrow(sdf)))
  
  vector.biomarker <- sdf[, biomarker ]
  vector.survival <- sdf[, time ]
  vector.event <- sdf[, event ]
  

  mod <- maxstat.test( 
    Surv( vector.survival, vector.event ) ~ vector.biomarker, 
    data=sdf, smethod="LogRank", pmethod="Lau92", iscores=TRUE, minprop = minimumSamps/length(vector.biomarker)
  )
  
  cutpoint <- signif(mod$estimate[[ 1 ]], digits = 4)
  p.value <- mod$p.value
  
  modstats <- mod$stats
  modcuts <-mod$cuts
  dataxy <- data.frame(modcuts, modstats)
  
  maxstatplot <- ggplot( dataxy, aes( x=modcuts, y=modstats ) ) +
    geom_line(colour= "#2c3e50") +
    geom_point(colour= "#2c3e50") +
    geom_vline(aes(xintercept=cutpoint)) +
    labs(x=biomarker,y="Standarized log-rank statistics")
  ggsave(paste(biomarker, "MaxstatPlot.png"))
  
  category <- ifelse( vector.biomarker < cutpoint, "low", "high" )
  x <- ifelse(vector.biomarker < cutpoint, 0, 1)
  category.df <- data.frame( sdf, category, x )
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, cutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="MaxStatResults.csv")


  rownames(result.table.cutp) = result.table.cutp$Estimate
  restext = result.table.cutp[c("Cutpoint", "HR", "P-value"), ]
  restext$Estimate[3] = 'p'
  restext = paste(restext$Estimate, restext$Result , sep=' = ', collapse='\n')
  print(restext)


  histogramPlotMaxstat <- ggplot(sdf, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramMaxstatWithLabs <- histogramPlotMaxstat + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramMaxstat.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotMaxstat_no_risk_table.png"))
  
  # http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
  #ggsurvfit <- ggsurvplot(fit, pval = TRUE, risk.table = TRUE)  #surv.col = c( "#2c3e50", "Red" ),   
  
  ggsurvfit <- ggsurvplot(fit, pval = F, risk.table = TRUE, risk.table.col = "strata", title = biomarker, ggtheme=theme_survminer(base_family = "Helvetica")) 
  # center align labels
  ggsurvfit$plot <- ggsurvfit$plot + theme(plot.title = element_text(hjust = 0.5)) 
  ggsurvfit$table <- ggsurvfit$table + theme(plot.title = element_text(hjust = 0.5)) 
  # add annotation in plot
  ggsurvfit$plot <- ggsurvfit$plot+ 
              ggplot2::annotate("text", x = Inf, y = Inf, vjust = 1, hjust = 1,
                                label = restext, size = 4) # x = 0, y = 0.2

  ggsave(paste(biomarker, "KaplanMeierPlot.png"), plot = print(ggsurvfit), dpi=600)
  ggsave(paste(biomarker, "KaplanMeierPlot.eps"), plot = print(ggsurvfit), dpi=600)
  ggsave(paste(biomarker, "KaplanMeierPlot_.pdf"), plot = print(ggsurvfit), dpi=600, onefile=F) # else produces an additional first blank page

}

# Rolr Method

rolrMethod <- function(setCutpoint, df, time, event, biomarker) {
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  func <- function(df) {
    df <- df[, c(biomarker, time, event)]
    df <- na.omit(df)
    #res <- rhier( times = vector.survival, status = vector.event, x = vector.biomarker, ns = 15, alt = 'increase' )
    res <- rhier( times = df[, time ], status = df[, event ], x = df[, biomarker ], ns = 15, alt = 'increase' ) # RR
  }
  
  rolr <- func( df )
  
  if (is.null(setCutpoint)) {
    low.cutoff = rolr[[ 1 ]][[ 1 ]]
    high.cutoff = rolr[[ 1 ]][[ 2 ]]
  } else {
    low.cutoff = setCutpoint[1]
    high.cutoff = setCutpoint[2]
  }
  
  
  result.table.col.names <- c("Cutpoint 1", "Cutpoint 2", "Low", "Medium", "High")
  length.results <- biomarkerCountsThree(vector.biomarker, low.cutoff, high.cutoff)
  result.table.row.names <- c(low.cutoff, high.cutoff, as.character(length.results[2]), as.character(length.results[3]), as.character(length.results[4]))
  result.table.rolr <- data.frame(result.table.col.names, result.table.row.names)
  colnames(result.table.rolr) <- c("Estimate", "Result")
  write.csv(result.table.rolr, file="GroupCounts.csv")
  
  category <- ifelse(vector.biomarker < low.cutoff, "low", ifelse(vector.biomarker > high.cutoff, "high", "medium" ))
  category.df <- data.frame(vector.biomarker, vector.survival, vector.event, category )
  
  category.df.lowmed <- data.frame( subset( category.df, subset= category == "medium" | category == "low" ))
  category.df.lowhigh <- data.frame( subset( category.df, subset= category == "high" | category == "low" ))
  category.df.medhigh <- data.frame( subset( category.df, subset= category == "high" | category == "medium" ))
  
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  if (low.cutoff != high.cutoff) {
    
    lowhigh <- do.call(
      survfit,
      list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df.lowhigh)
    )
    
    lowmedium <- do.call(
      survfit,
      list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df.lowmed)
    )
    
    mediumhigh <- do.call(
      survfit,
      list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df.medhigh)
    )
    
    cat.df <- data.frame( vector.biomarker, vector.survival, vector.event )
    
    categorize.low <- subset( cat.df, subset = vector.biomarker < low.cutoff | biomarker > high.cutoff )
    category.low = ifelse( categorize.low$vector.biomarker < low.cutoff, "low", "high" )
    cat.low <- data.frame( categorize.low, category.low )
    
    lowhigh.surv <- do.call(
      coxph,
      list( formula = Surv( vector.survival, vector.event ) ~ category.low, data = cat.low)
    )
    
    categorize.med <- subset( cat.df, subset = vector.biomarker < high.cutoff )
    category.med = ifelse( categorize.med$vector.biomarker < low.cutoff, "medium", "high" )
    cat.med <- data.frame(categorize.med, category.med)
    
    lowmedium.surv <- do.call(
      coxph,
      list( formula = Surv( vector.survival, vector.event ) ~ category.med, data = cat.med)
    )
    
    categorize.high <- subset( cat.df, subset = vector.biomarker > low.cutoff )
    category.high = ifelse( categorize.high$vector.biomarker < high.cutoff, "medium", "high" )
    cat.high <- data.frame( categorize.high, category.high )
    
    mediumhigh.surv <- do.call(
      coxph,
      list( formula = Surv( vector.survival, vector.event ) ~ category.high, data = cat.high)
    )

    groupEstimates <- function (group, summaryModel) {
      model <- summary( summaryModel )
      coef <- model$coefficients
      HR <- signif(coef[ 2 ], digits = 3 )
      q <- 1-(1-95/100)/2
      z <- qnorm(q)
      HR.lower <- signif( exp( coef[ 1 ] - z * coef[ 3 ]), digits = 3 )
      HR.upper <- signif( exp( coef[ 1 ] + z * coef[ 3 ]), digits = 3 )
      CI <- paste( " (", HR.lower, "-", HR.upper, ") " )
      p <- signif( model$sctest[ "pvalue" ], digits = 3 )
      estimates <- c(group,HR,CI, p)
      return(estimates)
    }

    estimates.total.table <- rbind(
      groupEstimates("Low vs High", lowhigh.surv),
      groupEstimates("Low vs Medium", lowmedium.surv),
      groupEstimates("Medium vs High", mediumhigh.surv)
    )
    colnames(estimates.total.table) <- c("Group", "HR", "CI", "p-value")
    write.csv(estimates.total.table, file="EstimatesTable.csv")
    
    lowmediumfit.plot <- ggsurvplot(lowmedium, surv.col = c( "Red", "Blue" ))
    ggsave(paste(biomarker, "LowVsMedium.png"))
    lowhighfit.plot <- ggsurvplot(lowhigh, surv.col = c( "#2c3e50", "Red" ))
    ggsave(paste(biomarker, "LowVsHigh.png"))
    mediumhighfit.plot <- ggsurvplot(lowhigh, surv.col = c("#2c3e50", "Red"))
    ggsave(paste(biomarker, "MediumVsHigh.png"))
    lowmediumhighfit.plot <- ggsurvplot(fit, surv.col = c("#2c3e50", "Red", "Blue"))
    ggsave(paste(biomarker, "LowVsMediumVsHigh.png"))
    
    histogramPlotRolr <- ggplot(df, aes(x=vector.biomarker)) +
      geom_histogram(fill="#2c3e50") +
      geom_vline(aes(xintercept=low.cutoff)) +
      geom_vline(aes(xintercept=high.cutoff))

    histogramRolrWithLabs <- histogramPlotRolr + labs(x=paste(biomarker))
    ggsave(paste(biomarker, "histogramRolr.png"))
    
    
  } else {
    print("Please select two distinct cutpoints.")
  }
}

# Method for 1 selected cutpoint

oneCutpointAdapt <- function(setCutpoint, df, time, event, biomarker) {
  
 vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  cph1 <- do.call(
    coxph,
    list( formula = Surv( vector.survival, vector.event ) ~ vector.biomarker, data = df)
  )
  
  
  allCutpointsDf <- data.frame(cutp(cph1))
  colnames(allCutpointsDf) <- c(biomarker, "U", "Q", "pvalue")

  category <- ifelse( vector.biomarker < setCutpoint, "low", "high" )
  x <- ifelse( vector.biomarker < setCutpoint, 0, 1 )
  category.df <- data.frame( df, category, x )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, setCutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="CutPResults.csv")
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  histogramPlotCutp <- ggplot(df, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=setCutpoint))
  
  histogramCutpWithLabs <- histogramPlotCutp + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramCutP.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotCutp.png"))
  
  mod <- maxstat.test( 
    Surv( vector.survival, vector.event ) ~ vector.biomarker, 
    data=df, smethod="LogRank", pmethod="Lau92", iscores=TRUE
  )
  
  modstats <- mod$stats
  modcuts <-mod$cuts
  dataxy <- data.frame(modcuts, modstats)
  
  maxstatplot <- ggplot( dataxy, aes( x=modcuts, y=modstats ) ) +
    geom_line(colour= "#2c3e50") +
    geom_point(colour= "#2c3e50") +
    geom_vline(aes(xintercept=setCutpoint)) +
    labs(x=biomarker,y="Standarized log-rank statistics")
  ggsave(paste(biomarker, "MaxstatPlot.png"))
  
}

# ------------------------------------------------------------------------------
# Functions needed for analysis

# Function performed for default analyses.

defaultAnalyses <- function(df, time, event, biomarker, minimumSamps) {
  #numberOfGroups = c("2","3")
  numberOfGroups = c("2") # for the LKT
  setwd(biomarker)
  for (n in numberOfGroups) {
    dir.create(n)
    setwd(n)
    if (n == 2) {
      #methods = c("cutp", "maxstat") #c("Youden", "ROC01", "cutp", "maxstat")
      methods = c("cutp") #c("Youden", "ROC01", "cutp", "maxstat")
      for (m in methods) {
        print(m)
        dir.create(m)
        setwd(m)
        if (m == "Youden") {
          youdenMethod(df, time, event, biomarker)
        } else if (m == "ROC01") {
          ROCMethod(df, time, event, biomarker)
        } else if (m == "cutp") {
          cutPMethod(df, time, event, biomarker, minimumSamps)
        } else if (m == "maxstat") {
          maxstatMethod(df, time, event, biomarker, minimumSamps)
        }
          setwd('..')
        }
    }
    if (n == 3) {
      print("rolr")
      dir.create("rolr")
      setwd("rolr")
      rolrMethod(setCutpoint, df, time, event, biomarker)
      setwd('..')
    }
    setwd('..')
  }
  setwd('..')
}

# Function performed if one cutpoint is set.

adaptCutpoint2groups <- function(setCutpoint, df, time, event, biomarker) {
  setwd(biomarker)
  dirname <- paste(as.character(setCutpoint), "_2groups")
  dir.create(dirname)
  setwd(dirname)
  tryCatch({
    oneCutpointAdapt(setCutpoint, df, time, event, biomarker)
  }, error = function(e) {
    print("Please select cutpoints between minimum and maximum biomarker value")
  })
  setwd('../..')
}

# Function performed if two cutpoints are set.

adaptCutpoint3groups <- function(setCutpoint, df, time, event, biomarker) {
  setwd(biomarker)
  dirname <- paste(as.character(setCutpoint[1]), as.character(setCutpoint[2]), "_3groups")
  dir.create(dirname)
  setwd(dirname)

  tryCatch({
    rolrMethod(setCutpoint, df, time, event, biomarker)
  }, error = function(e) {
    print("Please select cutpoints between minimum and maximum biomarker value")
  })
  
  setwd('../..')
}

# Choosing number of groups & performing analysis based on the provided value

mainFunction <- function(setCutpoint, df, time, event, biomarker, minimumSamps) {
  if (is.null(setCutpoint)) {
    defaultAnalyses(df, time, event, biomarker, minimumSamps)
  } else if (length(setCutpoint) == 1) {
    adaptCutpoint2groups(setCutpoint, df, time, event, biomarker)
  } else if (length(setCutpoint) == 2) {
    adaptCutpoint3groups(setCutpoint, df, time, event, biomarker)
  } else {
    print("Please select a correct cutpoint value")
  }
}
