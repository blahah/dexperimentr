#' Load a package, installing it if it isn't already installed
#'
#' This function loads a package. If the package isn't installed,
#' it is automatically installed from CRAN (default) or BioconductoR
#' (if bioconductor=TRUE).
#' @examples
#' get_package('ggplot2')
#' get_package('EBSeq', bioconductor=TRUE)
get_package <- function(package, bioconductor=FALSE) {
  if (!package %in% installed.packages()) {
    if (bioconductor) {
      return(install.packages(package))
    } else {
      source("http://bioconductor.org/biocLite.R")
      return(biocLite(package))
    }
  }
  return(library(package, character.only=T))
}