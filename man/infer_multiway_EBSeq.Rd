\name{infer_multiway_EBSeq}
\alias{infer_multiway_EBSeq}
\title{Perform the multiway differential expression experiment using EBSeq}
\usage{
infer_multiway_EBSeq(counts, conditions, emrounds = 25, plot = FALSE)
}
\description{
Return a list: - final: a data frame containing counts, DE
probability estimates, and fold changes - results: the
results object output by the DE test (useful for diagnostic
plots and QC)
}

