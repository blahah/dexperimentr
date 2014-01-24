# defining conditions

#` Generate the possible expression patterns for a binary experiment (A vs. B)
binary_conditions <- function(samples) {
  return(list(DE=samples,
              NDE=rep(1, length(samples))))
}

#` Generate all possible expression patterns for the supplied samples (All vs All)
all_conditions <- function(samples) {
  library(EBSeq)
  return(GetPatterns(samples))
}