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
                     emrounds=25,
                     prob_cutoff=0.95) {
  if (method == "EBSeq") {
    de_data <- infer_EBSeq(counts, conditions, emrounds=emrounds)
  } else {
    error("method not supported")
  }
  de_data <- add_means_and_errors(de_data, conditions)
  de_data[['final']] <- output_pattern_sets(de_data, conditions, prob_cutoff)
  de_data[['final']] <- merge_annotation(de_data[['final']], annotation_file)
  write_results(de_data)
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
infer_EBSeq <- function(counts, conditions, emrounds=25) {
  ncond = length(unique(conditions))
  if (ncond < 2) {
    error("There must be at least two different conditions to perform DE")
  } else if (ncond == 2) {
    # binary
    res <- infer_binary_EBSeq(counts, conditions, emrounds)
  } else {
    # multiway
    res <- infer_multiway_EBSeq(counts, conditions, emrounds)
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
infer_binary_EBSeq <- function(counts, conditions, emrounds=25) {
  get_package('EBSeq')
  
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
  # merge expression and probabilities
  final <- data.frame(gene.id=rownames(counts),
                      counts,
                      pp, 
                      log2FC=log2(fc$PostFC),
                      FC=fc$PostFC)
  
  # store the probability column indices for pattern detection
  prob_cols <- (ncol(counts)+1):(ncol(counts)+1+ncol(pp))
  
  return(list(final=final, 
              results=results,
              prob_cols=prob_cols))
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
infer_multiway_EBSeq <- function(counts, conditions, emrounds=25) {
  get_package('EBSeq')
  
  counts <- data.matrix(counts)
  
  # conditions
  patterns = GetPatterns(conditions)
  
  # normalization factors
  normfactors <- MedianNorm(counts)
  
  # run
  results <- EBMultiTest(counts,
                         NgVector=NULL,
                         Conditions=conditions,
                         AllParti=patterns,
                         sizeFactors=normfactors,
                         maxround=emrounds)
  
  # parse results
  pp <- GetMultiPP(results)
  
  # merge counts and DE
  final <- data.frame(gene.id=rownames(counts),
                      counts,
                      pp$PP)
  
  # store the probability column indices for pattern detection
  prob_cols <- (ncol(final)-ncol(pp$PP)):(ncol(final))
  
  # add fold-changes
  fc <- GetMultiFC(results)
  final <- cbind(final, fc$PostFCMat, fc$Log2PostFCMat)
  
  return(list(final=final, 
              results=results, 
              prob_cols=prob_cols))
}

#' Add means and standard errors to de_data by condition
#'
#' @param de_data list containing result data frame, DE test output and column indices
#' @param conditions vector of condition names
#' @return list, de_data with added mean and error columns for each condition
add_means_and_errors <- function(de_data, conditions) {
  final <- de_data[['final']]
  n_conds <- length(unique(conds))
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
  output <- merge(de_data, annot, by=by)
  
  return(unique(output))
}

#' Write out DE information for each DE gene expression pattern
output_pattern_sets <- function(de_data, conditions, prob_cutoff) {
  final <- de_data[['final']]
  mean_cols <- de_data[['mean_cols']]
  # add patterns
  final$pattern <- apply(final[,mean_cols],
                       1,
                       function(x) paste(pattern(x), collapse='_'))
  # select probable DE genes above cutoff
  n <- length(unique(conditions))
  if (n == 2) {
    prob_cols = which(names(final)=="PPDE")
    sig <- final[which(final[,prob_cols] >= prob_cutoff),]
  } else {
    alln <- length(conditions)
    prob_cols <- (alln + 2):(alln + n + 2)
    sig <- final[which(apply(final[,prob_cols], 1, function(x) any(x >= prob_cutoff))),]
  }
  patterns <- unique(sig$pattern)
  # export
  for (pat in patterns) {
    patsig <- sig[sig$patterns == pat,]
    if (nrow(patsig))
      write.table(x=patsig,
                  file=paste(paste(unique(conditions), collapse="_"), pat, ".csv", sep=""),
                  sep=",",
                  row.names=F,
                  col.names=T)
  }
  return(final)
}

#' Reduce a sequence of numbers to its pattern of changes.
#' 
pattern <- function(x, p=c(1), i=1, j=2) {
  if (x[i] == x[j]) {
    # no change
    n <- p[length(p)]
  } else if (x[i] > x[j]) {
    # decrease
    n <- p[length(p)] - 1
  } else {
    # increase
    n <- p[length(p)] + 1
  }
  # recurse along the pattern
  if (length(x) > j) {
    n <- pattern(x, n, i+1, j+1)
  }
  # add this to the pattern so far
  pattern <- c(p, n)
  # before the top-level return, adjust to so min is 1
  if (length(pattern) == length(x) && min(pattern) < 1) {
    pattern = pattern + abs(min(pattern)) + 1
  }
  return(pattern)
}

#' Write out a file containing all DE data collected
write_results <- function(de_data) {
  write.table(x=de_data[['final']],
              file="all_DE_data_with_patterns.csv",
              sep=",",
              row.names=F,
              col.names=T)
}