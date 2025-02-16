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
library("R.utils")

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

# Get biomarker counts

biomarkerCounts <- function(biomarker, cutpoint) {
  length.total <- length(biomarker)
  length.lower <- length(biomarker[biomarker < cutpoint ])
  length.upper <- length(biomarker[biomarker > cutpoint ])
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
  coef <- model$coefficients
  HR <- signif( coef[ 2 ], digits = 3 )
  q <- 1-(1-95/100)/2
  z <- qnorm( q )
  HR.lower <- signif( exp( coef[1] - z * coef[3] ), digits = 3 )
  HR.upper <- signif( exp( coef[1] + z * coef[3] ), digits = 3 )
  CI <- paste(" ( ", HR.lower, "-", HR.upper, " ) " )
  p <- signif( model$sctest[ "pvalue" ], digits = 3 )
  
  result.table.col.names <- c( "Cutpoint", "Biomarker < Cutpoint", "Biomarker > Cutpoint", "HR", "CI", "P-value" )
  result.table.row.names <- c( cutpoint, as.character(length.results[2]), as.character(length.results[3]), HR, CI, p )
  result.table <- data.frame(result.table.col.names, result.table.row.names)
  colnames(result.table) <- c("Estimate", "Result")
  return(result.table)
}

# CutP Method

cutPMethod <- function(df, time, event, biomarker) {
  
  method <- "cutp"
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  cph1 <- do.call(
    coxph,
    list( formula = Surv( vector.survival, vector.event ) ~ vector.biomarker, data = df)
  )
  
  
  allCutpointsDf <- data.frame(cutp(cph1))
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
  
  cutpoint <- signif(allCutpointsDf[[ 1, 1 ]], digits = 4 )
  p.value <- signif(allCutpointsDf[[ 1, 4 ]], digits = 4 )
  
  category <- ifelse( vector.biomarker < cutpoint, "low", "high" )
  x <- ifelse( vector.biomarker < cutpoint, 0, 1 )
  category.df <- data.frame( df, category, x )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, cutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="CutPResults.csv")
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  histogramPlotCutp <- ggplot(df, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramCutpWithLabs <- histogramPlotCutp + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramCutP.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotCutp.png"))
}


# Maxstat Method

maxstatMethod <- function(df, time, event, biomarker) {
  
  method <- "maxstat"
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  mod <- maxstat.test( 
    Surv( vector.survival, vector.event ) ~ vector.biomarker, 
    data=df, smethod="LogRank", pmethod="Lau92", iscores=TRUE
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
  category.df <- data.frame( df, category, x )
  
  fit <- do.call(
    survfit,
    list( formula = Surv( vector.survival, vector.event ) ~ category, data = category.df )
  )
  
  result.table.cutp <- resultTableSurvival(vector.biomarker, cutpoint, vector.survival, vector.event, category.df, x)
  write.csv(result.table.cutp, file="MaxStatResults.csv")
  
  histogramPlotMaxstat <- ggplot(df, aes(x=vector.biomarker)) +
    geom_histogram(fill="#2c3e50") +
    geom_vline(aes(xintercept=cutpoint))
  
  histogramMaxstatWithLabs <- histogramPlotMaxstat + labs(x=paste(biomarker))
  ggsave(paste(biomarker, "histogramMaxstat.png"))
  
  ggsurvfit <- ggsurvplot(fit, surv.col = c( "#2c3e50", "Red" ))
  ggsave(paste(biomarker, "KaplanMeierPlotMaxstat.png"))
  
}

# Rolr Method

rolrMethod <- function(setCutpoint, df, time, event, biomarker) {
  
  vector.biomarker <- df[, biomarker ]
  vector.survival <- df[, time ]
  vector.event <- df[, event ]
  
  func <- function(df) {
    df <- df[, c(biomarker, time, event)]
    df <- na.omit(df)
    res <- rhier(times = df[, time ], status = df[, event ], x = df[, biomarker ], ns = 15, alt = 'increase')
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

defaultAnalyses <- function(df, time, event, biomarker) {
  numberOfGroups = c("2","3")
  setwd(biomarker)
  for (n in numberOfGroups) {
    dir.create(n)
    setwd(n)
    if (n == 2) {
      methods = c("Youden", "ROC01", "cutp", "maxstat")
      for (m in methods) {
        dir.create(m)
        setwd(m)
        if (m == "Youden") {
          youdenMethod(df, time, event, biomarker)
        } else if (m == "ROC01") {
          ROCMethod(df, time, event, biomarker)
        } else if (m == "cutp") {
          cutPMethod(df, time, event, biomarker)
        } else if (m == "maxstat") {
          maxstatMethod(df, time, event, biomarker)
        }
        setwd('..')
      }
    }
    if (n == 3) {
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

mainFunction <- function(cutpoints, df, time, event, biomarker) {
  if (is.null(cutpoints)) {
    defaultAnalyses(df, time, event, biomarker)
  } else if (length(cutpoints) == 1) {
    adaptCutpoint2groups(cutpoints, df, time, event, biomarker)
  } else if (length(cutpoints) == 2) {
    adaptCutpoint3groups(cutpoints, df, time, event, biomarker)
  } else {
    print("Please select a correct cutpoint value")
  }
}



# ------------------------------------------------------------------------------
# Analysis



evaluateCutpoints <- function (biomarkerList = biomarkers, time, event, setCutpoint, setCutpoint2, df) {
  
  mainDir <- ("/data/rodriguesrr/FMT_pittsburg/Pittsburg_Cohorts_paper/ev-cutpoints/test")
  setwd(mainDir)
  
  if ( dir.exists("analysis-results")) {
    print("Directory already exists")
  } else {
    dir.create("analysis-results")
  }
  
  # https://raw.githubusercontent.com/budczies/CutoffFinder/master/inst/exampledata/breastcancer_GSE2034.txt
  file <- ("/data/rodriguesrr/scripts/R/Evaluate_cutpoints/breastcancer_GSE2034.txt")
  
  table <- read.table(file, sep = "\t", header = TRUE)
  df <- data.frame(table)
  
  currentDate <- format(Sys.time(), "%a %b %d %Y %X")
  setwd("analysis-results")
  dir.create(currentDate)
  
  colNames <- colnames(df)
  colsCommaSeparated <- dput(colNames)
  
  # Step 6. Choose if cutpoints should be counted or you would like to set your own cutpoint.
  # Paste the value of the cutpoint to produce plots and analyses for this cutpoint.
  # If setCutpoint variable is set to NULL, best cutpoint will be counted for each method:
  # cutp, maxstat, ROC01, Youden for 2 groups and rolr for 3 groups
  # If setCutpoint variable includes a list of values that is equal in length to biomarker list, histogram, ROC plot, Kaplan Meier plot and Standarized log-rank statistic plots will be produced.
  # If setCutpoint variable includes a list of values that is equal in length to biomarker list, low vs high, low vs medium and medium vs high statistics will be performed and  plots (histogram, Kaplan-Meier) will be produced.
  
  #setCutpoint <- NULL
  #setCutpoint2 <- NULL
  setCutpoint <- c(12.52, 5.33, 5.12, 4.83, 4.2, 7, 5.95, 4.08, 4.78, 2.54, 2.66, 11)
  setCutpoint2 <- c(14, 6.66, 7.52, 6.38, 5.73, 8.35, 6.68, 5.65, 6.61, 6.42, 4.9, 14.06)
  
  # Step 7. Choose biomarkers, time, event variables (look Step 5).
  
  time <- c("dmfs_time")
  event <- c("dmfs_event")
  biomarkers <- c("ESR1.205225_at", "ESR1.211233_x_at", "ESR1.211234_x_at", "ESR1.211235_s_at", "ESR1.211627_x_at", 
                  "ESR1.215551_at", "ESR1.215552_s_at", "ESR1.217163_at", "ESR1.217190_x_at", "PGR.208305_at", "ERBB2.210930_s_at", "ERBB2.216836_s_at")
  
  print("A")
  print(biomarkers)
  print("B")
  print(biomarkerList)
  
  # Step 8. Creating directories tree, preparing data & performing analyses.
  
  setwd(currentDate)
  print("C")
  
  prepareCutpointsDf <- function(biomarker, biomarkerList, cutpoint = NULL, cutpoint2 = NULL) {
    
    if (is.null(cutpoint) && is.null(cutpoint2)) {
      return (NULL)
    } else {
      cutpointsFunc(biomarker, biomarkerList, cutpoint, cutpoint2)
    }
  }
  
  cutpointsFunc <- function(biomarker, biomarkerList, cutpoint = NULL, cutpoint2 = NULL) {
    
    if (is.null(cutpoint2)) {
      cutpointDf <- data.frame(rbind(cutpoint))
    } else {
      cutpointDf <- data.frame(rbind(cutpoint, cutpoint2))
    }
    colnames(cutpointDf) <- biomarkerList
    rownames(cutpointDf) <- NULL
    
    result <- cutpointDf[[biomarker]]
    
  }
  
  if(is.null(setCutpoint) || is.null(setCutpoint2) || (length(setCutpoint) = length(biomarkerList)) || (length(setCutpoint2) = length(biomarkerList))) {
    for (b in biomarkerList) {
    print("C")
      dir.create(b)
      timeLen <- length(time)
      eventLen <- length(event)
      if ((timeLen == 0) || (eventLen == 0)) {
        print("Please correct time or event variables")
      } else if ((timeLen == 1) && (eventLen == 1)) {
        print("One time and event variable - results will be produced in biomarker folder")
        cutpointResult <- R.utils::doCall(prepareCutpointsDf, args=list(biomarker = b, biomarkerList=biomarkerList, cutpoint = setCutpoint, cutpoint2 = setCutpoint2))
        mainFunction(cutpoints = cutpointResult, df, time, event, b)
      } else if ((timeLen > 1) | (eventLen > 1)) {
        print("Please provide one time and one event variable")
      }
    }
  } else {
    print("Please select correct biomarkers, time, event and cutpoint values")
  }
  
  
}

evaluateCutpoints(
  biomarkers,
  time,
  event,
  setCutpoint,
  setCutpoint2,
  df
)
