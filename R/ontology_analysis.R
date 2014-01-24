# GO analysis

#` Run GO-term enrichment analysis
#`
#` expressiondata must be a data frame with columns corresponding to
#` expression patterns and rows corresponding to genes or transcripts.
#` Both rows and columns must be named. Row names must match the format
#` used in the mappingfile.
#`
#` mappingsfile must be the path to a file containing gene/transcript name
#` -> GO ID mappings.
#`
#` ppcutoff is the cutoff for posterior probability of being associated
#` with each pattern, used to filter out non-significant rows.
#`
#` alpha is the p-value cutoff to use in the Fisher's exact test.
ontology_enrichment <- function(expressiondata,
                                mappingsfile,
                                ppcutoff=0.95,
                                alpha=0.05) {
  get_package(topGO)
  geneID2GO <- readMappings(file = mappingsfile)
  
  # define the gene cutoff function for posterior probabilities
  topDiffGenes <- function(pp) {
    return(pp >= ppcutoff)
  }
  
  # iterate through patterns performing GO analysis
  for (pattern in colnames(expressiondata)) {
    allgenes <- expressiondata[,pattern]
    names(allgenes) <- rownames(expressiondata)
    ontologies <- c('BP', 'CC', 'MF')
    
    # iterate through ontologies testing each separately
    for (ontology in ontologies) {
      GOdata <- new("topGOdata",
                    description = "Test", ontology = ontology,
                    allGenes = allgenes, geneSel = topDiffGenes,
                    nodeSize = 10, annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
      
      # run fisher test
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      
      # write out short-form results
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                         orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
      allRes <- allRes[allRes$classicFisher <= 0.05,]
      write.table(file=paste(ontology, pattern, "GO.csv", sep='_'), 
                  x=allRes,
                  row.names=F,
                  col.names=T,
                  sep=",")
      
      # print graph of signficant nodes
      printGraph(GOdata, resultFisher, 
                 firstSigNodes = 15, 
                 fn.prefix = paste(ontology, pattern, sep='_'), 
                 useInfo = "all", 
                 pdfSW = TRUE)
    }
  }
}