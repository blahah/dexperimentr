# Create all intersections of GO enriched GO terms between patterns
# This script should be run from the directory containing all GO
# enrichment analysis result CSVs from dexperimentr

#' find the GO terms intersecting between
#' each possible pair of pattern enrichment sets
go_intersections <- function(ontology) {
  x <- list.files(pattern=ontology)
  pairs <- combn(x, 2)
  overlap <- data.frame()
  for (col in 1:ncol(pairs)) {
    lfile <- pairs[1,col]
    l <- read.csv(lfile, as.is=T)
    rfile <- pairs[2,col]
    r <- read.csv(rfile, as.is=T)
    i <- intersect(l$GO.ID, r$GO.ID)
    if(length(i) > 0) {
      this <- data.frame(ontology=ontology,
                         left=clean_filename(lfile),
                         right=clean_filename(rfile),
                         id=i,
                         term=l[l$GO.ID%in%i,'Term'])
      overlap <- rbind(overlap,
                       this)
    }
  }
  write.csv(x=overlap, file=paste(ontology, "intersections.csv", sep="_"))
}

#' Strip out file extension, unnecessary words and underscores
#' from a filename, returning the cleaned string
clean_filename <- function(x) {
  x <- gsub(x, pattern="_GO\\.csv", replacement="")
  x <- gsub(x, pattern=paste("_", ontology, "_", sep=""), replacement="")
  x <- gsub(x, pattern="_", replacement=" ")
  return(x)
}

# run for each ontology
go_intersections("BP")
go_intersections("CC")
go_intersections("MF")
