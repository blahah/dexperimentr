\name{infer_DE}
\alias{infer_DE}
\title{Perform the differential expression workflow using the specified method}
\usage{
infer_DE(counts, conditions, annotation_file, method = "EBSeq",
  plot = FALSE, emrounds = 25)
}
\description{
Return a list: - final: a data frame containing counts, DE
probability estimates, fold changes, and annotations -
results: the results object output by the DE test (useful
for diagnostic plots and QC)
}

