#' Load a package, installing it if it isn't already installed
#'
#' This function loads a package. If the package isn't installed,
#' it is automatically installed from CRAN (default) or BioconductoR
#' (if bioconductor=TRUE).
#' @examples
#' get_package('ggplot2')
#' get_package('EBSeq', bioconductor=TRUE)
get_package <- function(package, bioconductor=FALSE) {
  inst <- !package %in% installed.packages()
  if (!inst) {
    tryCatch({
        library(package, character.only=T)
      }, 
      silent=TRUE,
      error=function(e){
        inst <- TRUE
      })
  }
  if (inst) {
    print("installing")
    if (bioconductor) {
      source("http://bioconductor.org/biocLite.R")
      return(biocLite(package))
    } else {
      install.packages(package)
    }
    return(library(package, character.only=T))
  }
}