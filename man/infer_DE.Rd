\name{infer_DE}
\alias{infer_DE}
\title{Perform the differential expression workflow}
\usage{
infer_DE(counts, conditions, annotation_file, method = "EBSeq",
  plot = FALSE, emrounds = 25, prob_cutoff = 0.95)
}
\arguments{
  \item{counts}{data frame of counts with samples as named
  columns and genes as named rows}

  \item{conditions}{vector of condition names}

  \item{annotation_file}{string path to a TSV containing
  annotation data to add to the output data frame}

  \item{method}{string; the differential expression method
  to use (only 'EBSeq' currently)}

  \item{plot}{logical; if TRUE, diagnostic plots will be
  produced}

  \item{emrounds}{integer; number of rounds of expectation
  maximisation to use}
}
\value{
list: - final: a data frame containing counts, DE
probability estimates, fold changes, and annotations -
results: the results object output by the DE test (useful
for diagnostic plots and QC) - prob_cols: indices of
columns containing posterior probabilities - mean_cols:
indices of columns containing mean expression counts
}
\description{
Perform the differential expression workflow
}

