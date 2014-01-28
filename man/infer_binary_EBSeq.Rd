\name{infer_binary_EBSeq}
\alias{infer_binary_EBSeq}
\title{Perform the binary differential expression experiment using EBSeq}
\usage{
infer_binary_EBSeq(counts, conditions, emrounds = 25)
}
\arguments{
  \item{counts}{data frame of counts with samples as named
  columns and genes as named rows}

  \item{conditions}{vector of condition names}

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
Perform the binary differential expression experiment using
EBSeq
}

