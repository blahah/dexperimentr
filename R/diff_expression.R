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
                     named_patterns=list(),
                     prob_cutoff=0.95) {
  # create a directory for the outputs
  wd <- getwd()
  dir.create('de_data', showWarnings=FALSE)
  setwd('./de_data')
  
  # run the analysis
  if (method == "EBSeq") {
    de_data <- infer_EBSeq(counts, conditions, emrounds=emrounds)
  } else {
    stop("method not supported")
  }
  de_data <- add_means_and_errors(de_data, conditions)
  de_data[['final']] <- output_pattern_sets(de_data, conditions, 
                                            named_patterns, prob_cutoff)
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
infer_EBSeq <- function(counts, conditions, emrounds=25) {
  get_package('EBSeq', bioconductor=TRUE)
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
  merged <- merge(merged, fc$PostFC, by.x="Row.names", by.y="row.names", all.x=T)

  # tidy up merged data
  n <- ncol(merged)
  merged[,2:n] <- 
    data.frame(apply(merged[,2:n], 2, as.numeric)) # all cols are strings after merge

  names(merged)[1] <- "gene.id"
  names(merged)[ncol(merged)] <- "log2.FC"
  merged$log2.FC <- log2(merged$log2.FC)
  
  # store the probability column indices for pattern detection
  prob_cols <- (ncol(counts)+1):(ncol(counts)+1+ncol(pp))
  
  return(list(final=merged, 
              results=results,
              prob_cols=prob_cols))
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
infer_multiway_EBSeq <- function(counts, conditions, emrounds=25) {
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

  # normalise counts
  normcounts <- normalise_counts(counts, normfactors)
  
  # merge counts and DE
  final <- data.frame(gene.id=rownames(counts),
                      normcounts,
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
  output <- merge(de_data, annot, by=by, all.x=T)
  
  return(unique(output))
}

#' Write out DE information for each DE gene expression pattern
output_pattern_sets <- function(de_data, conditions, 
                                named_patterns, prob_cutoff) {
  final <- de_data[['final']]
  mean_cols <- de_data[['mean_cols']]
  # add patterns
  final$pattern <- apply(final[,mean_cols],
                       1,
                       function(x) paste(pattern(x), collapse='_'))
  # genes below cutoff should have flat pattern
  final$pattern[which(final$PPDE < prob_cutoff)] <- 
    paste(rep(1, length(mean_cols)), collapse="_")
  # replace named patterns
  final$pattern <- sapply(final$pattern,
                           function(x) {
                             if (x %in% names(named_patterns)) {
                               return(named_patterns[[x]])
                             } else {
                               return(x)
                             }
                           })
  # select probable DE genes above cutoff
  n <- length(unique(conditions))
  if (n == 2) {
    sig <- final[which(final$PPDE >= prob_cutoff),]
  } else {
    alln <- length(conditions)
    prob_cols <- (alln + 2):(alln + n + 2)
    sig <- final[which(apply(final[,prob_cols], 1, function(x) any(x >= prob_cutoff))),]
  }
  if (!nrow(sig)) {
    stop("There are no rows with significantly differential expression")
  } else {
    print(paste("There were", nrow(sig), "rows (out of", nrow(final), 
      "tested) with significantly differential expression (>=", prob_cutoff, ")"))
    print(table(sig$pattern))
  }
  patterns <- unique(sig$pattern)
  # export
  for (pat in patterns) {
    patsig <- sig[sig$pattern == pat,]
    if (nrow(patsig)) {
      print(paste("Saving", nrow(patsig), "significant results for pattern:", pat))
      write.table(x=patsig,
                  file=paste(paste(unique(conditions), collapse="_"), pat, ".csv", sep=""),
                  sep=",",
                  row.names=F,
                  col.names=T)
    } else {
      print(paste("No signficant results for pattern:", pat))
    }
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
  print("Writing out all results")
  write.table(x=de_data[['final']],
              file="all_DE_data_with_patterns.csv",
              sep=",",
              row.names=F,
              col.names=T)
}