# generate all combinatorial patterns and return
# them as a list suitable for use with baySeq
generate_patterns <- function(conditions) {
  get_package('EBSeq', bioconductor=TRUE)
  patterns <- GetPatterns(conditions)
  library(plyr)
  patterns <- alply(patterns, 1)
  patterns <- as.list(patterns)
  patterns <- lapply(patterns, function(x) {
    z <- sapply(x, function(y) { rep(y, 3) })
    dim(z) <- NULL
    return(z)
  })
}

# given a countdata object (CD) and a numeric
# vector specifying the rows to use as a subset,
# learn the patterns associated with those genes
# and return them as a list, ready to be used
# to specify the groups in a subsequent baySeq run
learn_patterns <- function(CD, subset, samplesize=1000,
                           cluster=NULL, consensus=TRUE,
                           posterior_cutoff=0.9, plot=TRUE) {
  get_package("baySeq", bioconductor=TRUE)
  CD <- getPriors.NB(CD,
                     samplesize=1000,
                     estimation="QL",
                     cl=cluster,
                     consensus=consensus)

  # don't specify a distribution for getLikelihoods
  CD <- getLikelihoods(CD,
                       pET='BIC',
                       cl=cluster,
                       subset=subset)

  # get best pattern for each row with posterior cutoff
  row_patterns <- get_row_patterns(CD, subset=subset, cutoff=posterior_cutoff)
  row_patterns <- row_patterns[row_patterns != 'none']

  list(cd=CD, patterns=CD@groups[names(CD@groups) %in% row_patterns])
}

# assign the most likely pattern to each row
# optionally return only a subset of the rows
get_row_patterns <- function(CD, subset=TRUE, cutoff=NULL) {
  row_patterns <- apply(CD@posteriors[subset,],
                    MARGIN=1,
                    FUN=function(x) {
                      if(all(is.na(x))) {
                        # not in the subset
                        return(NA)
                      }
                      posts <- exp(x)
                      if (!is.null(cutoff) &&
                          max(posts) < cutoff) {
                        # probability below cutoff
                        return("none")
                      }
                      maxidx <- which.max(exp(x))
                      names(x)[maxidx]
                  })
  return(row_patterns[!is.na(row_patterns)])
}

# get the posterior associated with the most likely pattern for each row
# optionally return only a subset of the rows
get_row_pattern_posts <- function(CD, subset=TRUE) {
  row_pattern_posts <- apply(CD@posteriors[subset,],
                        MARGIN=1,
                        FUN=function(x) {
                          if(all(is.na(x))) {
                            return(NA)
                          }
                          max(exp(x))
                        })
  return(row_pattern_posts[!is.na(row_pattern_posts)])
}