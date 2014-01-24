\name{get_package}
\alias{get_package}
\title{Load a package, installing it if it isn't already installed}
\usage{
get_package(package, bioconductor = FALSE)
}
\description{
This function loads a package. If the package isn't
installed, it is automatically installed from CRAN
(default) or BioconductoR (if bioconductor=TRUE).
}
\examples{
get_package('gpplot2')
get_package('EBSeq', bioconductor=TRUE)
}

