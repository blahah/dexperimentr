# Differential expression

#' Perform the differential expression workflow
#' 
#' 
#'
#' @param counts data frame of counts with samples as named columns 
#' and genes as named rows
#' @param conditions vector of condition names
#' @param annotation_file string path to a TSV containing annotation data to add
#' to the output data frame
#' @param method string; the differential expression method to use (only 'EBSeq' currently)
#' @param emrounds integer; number of rounds of expectation maximisation to use
#' @param prob_cutoff double; posterior probability to use as a cutoff for differential expression
#' @return list:
#' - final: a data frame containing counts, DE probability estimates,
#'          fold changes, and annotations
#' - results: the results object output by the DE test (useful for 
#'            diagnostic plots and QC)
#' - prob_cols: indices of columns containing posterior probabilities
#' - mean_cols: indices of columns containing mean expression counts
infer_DE <- function(counts, 
                     conditions, 
                     annotation_file,
                     method="EBSeq",
                     emrounds=5,
                     patterns=NULL,
                     named_patterns=list(),
                     prob_cutoff=0.95) {
  # create a directory for the outputs
  wd <- getwd()
  dir.create('de_data', showWarnings=FALSE)
  setwd('./de_data')
  
  # run the analysis
  if (method == "EBSeq") {
    de_data <- infer_EBSeq(counts, conditions, emrounds=emrounds, patterns=patterns)
  } else {
    stop("method not supported")
  }
  de_data <- add_means_and_errors(de_data, conditions)
  pattern_data <- output_pattern_sets(de_data, conditions, 
                                      named_patterns, prob_cutoff)
  de_data[['final']] <- pattern_data[['final']]
  de_data[['num_patterns']] <- pattern_data[['num_patterns']]
  de_data[['word_patterns']] <- replace_patterns(pattern_data[['num_patterns']],
                                                 named_patterns)
  de_data[['final']] <- merge_annotation(de_data[['final']], annotation_file)
  write_results(de_data)
  
  # back to the previous directoy
  setwd(wd)
  
  return(de_data)
}

#' Perform the differential expression experiment using EBSeq, automatically
#' selecting either the binary or multiway workflow based on the number of
#' unique conditions.
#'
#' @param counts data frame of counts with samples as named columns 
#' and genes as named rows
#' @param conditions vector of condition names
#' @param emrounds integer; number of rounds of expectation maximisation to use
#' @return list:
#' - final: a data frame containing counts, DE probability estimates,
#'          fold changes, and annotations
#' - results: the results object output by the DE test (useful for 
#'            diagnostic plots and QC)
#' - prob_cols: indices of columns containing posterior probabilities
#' - mean_cols: indices of columns containing mean expression counts
infer_EBSeq <- function(counts, conditions, emrounds=5, patterns=NULL) {
  get_package('EBSeq', bioconductor=TRUE)
  ncond = length(unique(conditions))
  if (ncond < 2) {
    error("There must be at least two different conditions to perform DE")
  } else if (ncond == 2) {
    # binary
    res <- infer_binary_EBSeq(counts, conditions, emrounds)
  } else {
    # multiway
    res <- infer_multiway_EBSeq(counts, conditions, emrounds, patterns)
  }
  return(res)
}

#' Perform the binary differential expression experiment using EBSeq
#'
#' @param counts data frame of counts with samples as named columns 
#' and genes as named rows
#' @param conditions vector of condition names
#' @param emrounds integer; number of rounds of expectation maximisation to use
#' @return list:
#' - final: a data frame containing counts, DE probability estimates,
#'          fold changes, and annotations
#' - results: the results object output by the DE test (useful for 
#'            diagnostic plots and QC)
#' - prob_cols: indices of columns containing posterior probabilities
#' - mean_cols: indices of columns containing mean expression counts
infer_binary_EBSeq <- function(counts, conditions, emrounds=5) {  
  counts <- data.matrix(counts)
  
  # normalization factors
  normfactors <- MedianNorm(counts)
  
  # run
  results <- EBTest(Data=counts,
                    Conditions=conditions,
                    sizeFactors=normfactors,
                    maxround=emrounds)
  
  # posterior probabilities
  pp <- GetPPMat(results)
  
  # calculate fold change
  fc <- PostFC(results)
  PlotPostVsRawFC(EBOut=results, FCOut=fc)

  # get normalised counts
  normcounts <- normalise_counts(counts, normfactors)

  # merge
  merged <- merge(normcounts, pp, by="row.names", all.x=T)
  merged <- merge(merged, fc$PostFC, 
                  by.x="Row.names", by.y="row.names", 
                  all.x=T)

  # tidy up merged data
  n <- ncol(merged)
  merged[,2:n] <- 
    data.frame(apply(merged[,2:n], 2, as.numeric)) # all cols are strings after merge

  names(merged)[1] <- "gene.id"
  names(merged)[ncol(merged)] <- "log2.FC"
  merged$log2.FC <- log2(merged$log2.FC)
  
  return(list(final=merged, 
              results=results,
              de_prob_cols='PPDE',
              ee_prob_col='PPEE'))
}

normalise_counts <- function(counts, normfactors) {
  round(t(t(counts) / normfactors))
}

#' Perform the multiway differential expression experiment using EBSeq
#'
#' @param counts data frame of counts with samples as named columns 
#' and genes as named rows
#' @param conditions vector of condition names
#' @param emrounds integer; number of rounds of expectation maximisation to use
#' @return list:
#' - final: a data frame containing counts, DE probability estimates,
#'          fold changes, and annotations
#' - results: the results object output by the DE test (useful for 
#'            diagnostic plots and QC)
#' - prob_cols: indices of columns containing posterior probabilities
#' - mean_cols: indices of columns containing mean expression counts
infer_multiway_EBSeq <- function(counts, conditions, emrounds=5, patterns=NULL) {
  counts <- data.matrix(counts)
  
  # conditions
  if (is.null(patterns)) {
    patterns = GetPatterns(conditions)
  }
  
  # normalization factors
  normfactors <- MedianNorm(counts)
  
  # run
  results <- EBMultiTest(counts,
                         NgVector=NULL,
                         Conditions=conditions,
                         AllParti=patterns,
                         sizeFactors=normfactors,
                         maxround=emrounds)
#                         Alpha=1.1,
#                         Beta=0.68,
                         # fixHyper=FALSE)
  
  # parse results
  pp <- GetMultiPP(results)

  # add fold-changes
  fc <- GetMultiFC(results)

  # normalise counts
  normcounts <- normalise_counts(counts, normfactors)
  
  # merge counts and DE
  merged <- merge(normcounts, pp$PP, by="row.names", all.x=T)
  merged <- merge(merged, fc$PostFCMat, 
                  by.x="Row.names", by.y="row.names", 
                  all.x=T)
  merged <- merge(merged, fc$Log2PostFCMat, 
                  by.x="Row.names", by.y="row.names", 
                  all.x=T, suffixes=c('', '.log2'))
  names(merged)[1] <- "gene.id"
  
  # store the probability column indices for pattern detection
  de_prob_cols <- rownames(patterns)
  # get the no-difference pattern
  ee_prob_col <- get_EE_pattern(patterns)
  # remove the ee col from the de cols
  de_prob_cols <- setdiff(de_prob_cols, ee_prob_col)
  
  return(list(final=merged, 
              results=results, 
              de_prob_cols=de_prob_cols,
              ee_prob_col=ee_prob_col))
}

#' Add means and standard errors to de_data by condition
#'
#' @param de_data list containing result data frame, DE test output and column indices
#' @param conditions vector of condition names
#' @return list, de_data with added mean and error columns for each condition
add_means_and_errors <- function(de_data, conditions) {
  final <- de_data[['final']]
  n_conds <- length(unique(conditions))
  n_cols <- ncol(final)
  mean_cols <- c()
  for (cond in unique(conditions)) {
    cols <- which(conditions == cond)
    # shift right by one to account for gene names
    cols <- cols + 1
    final[,paste(cond, 'mean', sep='.')] <- rowMeans(final[,cols])
    mean_cols <- c(mean_cols, ncol(final))
    final[,paste(cond, 'stderr', sep='.')] <- apply(final[,cols], 1, function(x) sd(x)/sqrt(length(x)))
  }
  de_data[['final']] <- final
  de_data[['mean_cols']] <- mean_cols
  return(de_data)
}

#' Create a data frame with all expression, DE and annotation data
merge_annotation <- function(de_data, annotation_file, by='gene.id') {
  annot <- read.csv(annotation_file, head=T, as.is=T)

  # annotation file must have a column with name gene.id
  # otherwise we assume the first column is the gene id  
  if (!(by %in% names(annot))) {
    names(annot)[1] <- by
  }

  output <- merge(de_data, annot, by=by, all.x=T)
  
  return(unique(output))
}

#' Write out a file containing all DE data collected
write_results <- function(de_data) {
  print("Writing out all results")
  final <- de_data[['final']]
  write.table(x=final,
              file="all_DE_data_with_patterns.csv",
              sep=",",
              row.names=F,
              col.names=T)
  patterns <- unique(final$pattern)
  # export
  for (pat in patterns) {
    patrows <- subset(final, pattern == pat)
    if (nrow(patrows)) {
      print(paste("Saving", nrow(patrows), " results for pattern:", pat))
      write.table(x=patrows,
                  file=paste(pat, ".csv", sep=""),
                  sep=",",
                  row.names=F,
                  col.names=T)
    } else {
      print(paste("No signficant results for pattern:", pat))
    }
  }
}

