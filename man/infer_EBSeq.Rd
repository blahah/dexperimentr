\name{infer_EBSeq}
\alias{infer_EBSeq}
\title{Perform the differential expression experiment using EBSeq, automatically
selecting either the binary or multiway workflow based on the number of
unique conditions.}
\usage{
infer_EBSeq(counts, conditions, emrounds = 25, plot = FALSE)
}
\description{
Return a list: - final: a data frame containing counts, DE
probability estimates, and fold changes - results: the
results object output by the DE test (useful for diagnostic
plots and QC)
}

