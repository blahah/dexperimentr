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
ontology_enrichment <- function(de_data,
                                mappingsfile,
                                conditions,
                                ppcutoff=0.95,
                                alpha=0.05) {
  go_results <- perform_GO_enrichment(de_data,
                        mappingsfile,
                        conditions,
                        ppcutoff,
                        alpha)
  plot_high_level_GO(go_results)
  plot_detailed_GO(de_data, go_results)
}

plot_high_level_GO <- function(go_results) {
  for (samplename in names(go_results)) {
    sample <- go_results[[samplename]]
    for (ontname in names(sample)) {
      ontology <- sample[[ontname]]
      plot_ontology(ontology)
    }
  }
}

plot_detailed_GO <- function(go_results,
                             maxsize = 50,
                             minsize = 2) {
  terms <- terms[which(terms < maxsize)]
  terms <- terms[which(terms > minsize)]
  length(terms)
  for(term in names(terms)) {
    tryCatch({
      print(term)
      plot_term(term)
    }, error=function(err) {
      print(paste('failed to plot', term, "due to error:", err))
    })
  }
}

plot_term <- function(term, save=TRUE) {
  d <- unique(alldata[alldata$Term == term,]) 
  d <- d[-is.na(d$id),c("fc", "Annotation", "Ontology", "id", 'up')]
  names(d) <- c('fc', 'gene', "ontology", "id", 'up')
  d$gene <- paste(d$gene, " (", d$id, ")", sep="") # add locus IDs to gene names
  d$gene <- gsub('(.{1,50})(\\s|$)', '\\1\n', d$gene) # split lines at 50 chars
  d$gene <- factor(d$gene, levels=d$gene[order(d$up, -d$fc)], ordered=TRUE)
  p <- ggplot(d, aes(x=gene, y=fc, fill=up)) +
    geom_bar(stat="identity") +
    scale_fill_manual(name="Enriched in", values = colours, labels = labels) +
    theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), # top, right, bottom, left
          legend.text.align = 0, # left align
          axis.text.y = element_text(hjust=1,
                                     vjust=0.8)) +
    theme(panel.border = theme_border()) +
    theme(axis.ticks.y=element_blank()) +
    theme(panel.grid.major.y = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    scale_y_continuous(name=expression('Log2 fold change (' * italic('O. sativa') * '/' *
                                         italic('E. glabrescens') * ')', sep=""),
                       expand = c(0, 0)) +
    scale_x_discrete(name="Gene",
                     expand=c(0,0)) +
    coord_flip()
  if (save) {
    term <- gsub(x=term, pattern="'", replacement="")
    term <- gsub(x=term, pattern=" ", replacement="_")
    ggsave(paste('GOplots/', d$ontology[1], '_', term, '.png', sep=""),
           height=7*length(d$fc)+80, width=350, units="mm")
    ggsave(paste('GOplots/', d$ontology[1], '_', term, '.pdf', sep=""),
           height=7*length(d$fc)+80, width=350, units="mm")
  }
  return(p)
}

plot_ontology <- function(ontology, save=TRUE) {
  d <- both.pc[both.pc$ontology == "BP",]
  d <- d[with(d, order(variable)),]
  p <- ggplot(data=d, aes(x=description, y=value, fill=variable)) +
    geom_bar(stat='identity') +
    theme_bw() +
    theme(legend.text.align = 0) + # left align
    theme(panel.border = theme_border()) +
    theme(axis.ticks.y=element_blank()) +
    scale_y_continuous(name="Percent of genes differentially expressed",
                       expand = c(0, 0)) +
    scale_x_discrete(name="GO annotation (count of genes)",
                     expand=c(0,0)) +
    scale_fill_manual(name="Enriched in", values=colours, labels=labels) +
    coord_flip()
  if(save) {
    ggsave(paste(ontology, '_GO_enrichment.png', sep=''))
    ggsave(paste(ontology, "_GO_enrichment.pdf", sep=''))
  }
  return(p)
}

.pt <- 1 / 0.352777778

len0_null <- function(x) {
  if (length(x) == 0)  NULL
  else                 x
}

theme_border <- function(
  type = c("left", "right", "bottom", "top", "none"),
  colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + opts( panel.border=theme_border(type=c("bottom","left")) ) + ...
  type <- match.arg(type, several.ok=TRUE)
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border", "element_blank", "element")
  )
}

element_grob.theme_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  type = NULL,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  if (is.null(type)) type = element$type
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if ("bottom" %in% type) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (length(type)==0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

perform_GO_enrichment <- function(de_data,
                                mappingsfile,
                                conditions,
                                ppcutoff=0.95,
                                alpha=0.05) {
  final <- de_data[['final']]
  if (length(unique(conditions)) == 2)
    prob_cols <- which(names(final)=='PPDE')
  else
    prob_cols <- de_data[['prob_cols']]
  
  get_package('topGO')
  geneID2GO <- readMappings(file = mappingsfile)
  
  # define the gene cutoff function for posterior probabilities
  topDiffGenes <- function(pp) {
    return(pp >= ppcutoff)
  }
  
  results = list()
  
  # iterate through patterns performing GO analysis
  for (pattern in colnames(final)[prob_cols]) {
    print(paste("GO enrichment testing for pattern", pattern))
    pattern_results <- list()
    
    allgenes <- final[,pattern]
    names(allgenes) <- final$gene.id
    ontologies <- c('BP', 'MF', 'CC')
    
    # iterate through ontologies testing each separately
    for (ontology in ontologies) {
      print(paste("Fisher testing GO enrichment for ontology", ontology, "in condition", pattern))
      GOdata <- new("topGOdata",
                    description = "Test", ontology = ontology,
                    allGenes = allgenes, geneSel = topDiffGenes,
                    nodeSize = 10, annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
      
      # run fisher test
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      
      # write out short-form results
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                         orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
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
      
      # save ontology
      pattern_results[[ontology]] <- allRes
    }
    # save pattern
    results[['pattern']] <- pattern_results
  }
  return(results)
}