# basic count data QC

#` Run the count data QC pipeline
count_data_QC <- function(df, conditions) {
  # create a directory for the outputs
  wd <- getwd()
  dir.create('count_qc', showWarnings=FALSE)
  setwd('./count_qc')
  
  # run the QC analysis
  plot_log_dists(df, 
                 conditions, 
                 'with zero rows', 
                 save=TRUE);
  df <- remove_zero_rows(df);
  df <- df + 1
  plot_log_dists(df, 
                 conditions, 
                 'without zero rows', 
                 save=TRUE);
  # check_nbinom(df);
  
  # back to the previous directory
  setwd(wd)
  return(df)
}

#` Plot the distribution of log counts for visual QC
plot_log_dists <- function(df, conditions, suffix='', save=FALSE) {
  get_package('reshape2')
  get_package('ggplot2')
  title=paste("distribution of log counts by sample", suffix)
  logcounts <- log(df)
  logcounts <- melt(logcounts)
  conditiontable <- data.frame(conditions=conditions, variable=names(df))
  logcounts <- merge(logcounts, conditiontable, by='variable')
  suppressWarnings({
    p <- ggplot(data=logcounts) +
      geom_line(aes(x=value, group=variable, colour=conditions), stat="density") +
      xlab("log expression count") +
      ggtitle(title)
    print(p)
    if (save) ggsave(plot=p, file=paste(gsub(title, pattern=' ', replacement='_'), ".pdf", sep=""))
  })
  return(p)
}

#` Check whether the data are negative-binomially distributed
check_nbinom <- function(df, plots=FALSE) {
  get_package('fitdistrplus')
  get_package('plyr')

  # fit the distribution to each sample
  results <- lapply(df, best_distribution, plots=plots)
  problems = 0
  for(result in results) {
    if (result != "nbinom") {
      problems ++
      warning(paste("Sample", names(df)[i], "might be", results[[i]], "distributed"))
    }
  }
  return(problems)
}

#` returns the distribution most likely to fit the data
best_distribution <- function(x, plots=TRUE) {
  distributions <- c("pois", "nbinom", "binom")
  
  suppressMessages({
    suppressWarnings({
      fitdata <- lapply(distributions, fit_distribution, x=x)
    })
  })
  
  distributions <- distributions[!is.na(fitdata)]
  fitdata <- remove_na_list(fitdata)
  aics <- lapply(fitdata, function(y) y$aic)
  aics <- remove_na_list(aics)
  distributions <- lapply(fitdata, function(y) y$distname)
  
  if (plots) {
    cdfcomp(fitdata, xlogscale=TRUE, legendtext=distributions, discrete=TRUE)
    denscomp(fitdata, xlim=c(1,100), legendtext=distributions)
  }
  
  return(distributions[which.min(aics)])
}

#` Remove NAs from a list
remove_na_list <- function(x) {
  return(x[!is.na(x)])
}

#` fits parameters for a distribution to the data
#`
#` return NA if fitting produced an error
fit_distribution <- function(distr, x, method="mle") {
  ret <- try(fitdist(data=x, distr=distr, method=method),
      silent=TRUE)
  if (class(ret) == "try-error") ret <- NA
  return(ret)
}

#` Remove any rows containing all-zeros
remove_zero_rows <- function(df) {
  df[apply(df, 1, function(x) !all(x==0)),]
}
