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
                                longmappingsfile,
                                named_patterns=list(),
                                ppcutoff=0.95,
                                alpha=0.05) {
  wd <- getwd()
  dir.create('functional_analysis', showWarnings=FALSE)
  setwd('./functional_analysis')
  for (d in c('results', 'graphs', 'plots')) {
    dir.create(d, showWarnings=FALSE)
  }

  go_results <- perform_GO_enrichment(de_data,
                        mappingsfile,
                        conditions,
                        ppcutoff,
                        alpha)
  # replace named patterns
  print('replacing pattern names..')
  go_results$Pattern <- sapply(go_results$Pattern,
                          function(x) {
                            if (x %in% names(named_patterns)) {
                              return(named_patterns[[x]])
                            } else {
                              return(x)
                            }
                        })
  # preserve input order of pattern names
  go_results$Pattern <- factor(go_results$Pattern, levels=named_patterns, ordered=TRUE)
  # make plots
  plot_high_level_GO(go_results, de_data, longmappingsfile)
  plot_detailed_GO(de_data, go_results)
  setwd(wd)
}

#' Produce a horizontal stacked barplot showing the
#' proportion of genes in each enriched GO term that
#' fall into each expression pattern
plot_high_level_GO <- function(go_results, de_data, longmappingsfile) {
  print("Making high-level GO plots")
  long_go <- read.csv(longmappingsfile, sep="\t", as.is=TRUE)
  names(long_go)[2] <- "gene.id"
  de_data[['final']] <- merge(de_data[['final']], long_go, by="gene.id")
  final <- de_data[['final']]
  for (pattern in unique(go_results$Pattern)) {
    subset <- go_results[go_results$Pattern == pattern,]
    for (ontname in unique(subset$Ontology)) {
      ontology <- subset[subset$Ontology == ontname,]
      plot_ontology(ontname, ontology, de_data, pattern)
    }
  }
}

#' For each enriched GO term produce a plot showing
#' the expression and direction of each gene.
#' 
#' If there are two conditions, a horizontal barplot is produced,
#' with log2 fold change on the x-axis, genes on the y-axis, and
#' one gene per line.
#' 
#' If there are >2 conditions, a grid of barplots is produced,
#' with one barplot per gene, one bar per condition, and the y-axis
#' within each barplot being log expression count.
plot_detailed_GO <- function(go_results,
                             de_data,
                             maxsize = 50,
                             minsize = 2) {
  print("Making detailed GO plots")
  terms <- terms[which(terms < maxsize)]
  terms <- terms[which(terms > minsize)]
  length(terms)
  for(term in names(terms)) {
    tryCatch({
      print(term)
      plot_term(term, de_data)
    }, error=function(err) {
      print(paste('failed to plot', term, "due to error:", err))
    })
  }
}

#' Produce a plot showing the expression and direction of each gene
#' in the specified GO term.
#' 
#' If there are two conditions, a horizontal barplot is produced,
#' with log2 fold change on the x-axis, genes on the y-axis, and
#' one gene per line.
#' 
#' If there are >2 conditions, a grid of barplots is produced,
#' with one barplot per gene, one bar per condition, and the y-axis
#' within each barplot being log expression count.
plot_term <- function(term, save=TRUE) {
  print(paste("Plotting", term, "GO term"))
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
  print(p)
  return(p)
}


#' Produce a plot with all enriched GO terms and the prop of each pattern
plot_ontology <- function(ontname, go_results, de_data, pattern, save=TRUE) {
  print(paste("Plotting", ontname, "ontology"))
  # extract data from de_data
  final <- de_data[['final']]
  mean_cols <- de_data[['mean_cols']]
  # subset to this ontology
  d <- final[final$Ontology == ontname,]
  # subset to enriched go terms
  d <- d[d$go.id %in% go_results$GO.ID,]
  # we only want the pattern and GOid columns
  d <- d[,c('pattern', 'go.id')]
  # count each pattern per GOid
  d <- as.data.frame(t(table(d)))
  # convert to percentages, count totals
  d <- ddply(d, .(go.id), function(x) {
    goid <- x$go.id[1]
    total <- go_results$Annotated[go_results$GO.ID==goid][1]
    # TODO: fix these proportions - are we leaving out some counts because
    # they're not significant?
    # START HERE
    data.frame(x, prop = x$Freq / total, total)
  })
  # add descriptions
  d <- merge(d, 
             data.frame(go.id=go_results$GO.ID, 
                        description=go_results$Term), 
             by='go.id')
  # order
  d$pattern <- factor(d$pattern, levels=levels(go_results$Pattern), ordered=TRUE)
  d$description <- paste(d$description, "\n", d$go.id, " (", d$total, ")", sep="")
  d <- d[with(d, order(-total, pattern)),]
  d <- transform(d, description = reorder(description, -total))
  p <- ggplot(data=d, aes(x=description, y=prop, fill=pattern)) +
    geom_bar(stat='identity') +
    theme_bw() +
    theme(legend.text.align = 0) + # left align
    theme(panel.border = theme_border()) +
    theme(axis.ticks.y=element_blank()) +
    scale_y_continuous(name="Percent of genes differentially expressed",
                       expand = c(0, 0)) +
    scale_x_discrete(name="GO annotation (count of genes)",
                     expand=c(0, 0)) +
    scale_fill_brewer(name="Enriched in", type="div") +
    coord_flip()
  if(save) {
    height = (3*dim(d)[1])+60
    width = 300
    units = "mm"
    ggsave(paste("plots/", pattern, "_", ontname, '_GO_enrichment.png', sep=''),
           width=width,
           height=height,
           units=units)
    ggsave(paste("plots/", pattern, "_", ontname, "_GO_enrichment.pdf", sep=''),
           width=width,
           height=height,
           units=units)
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
  print("Performing GO enrichment tests")
  final <- de_data[['final']]
  if (length(unique(conditions)) == 2) {
    prob_cols <- which(names(final)=='PPDE')
  } else {
    prob_cols <- de_data[['prob_cols']]
  }

  get_package('topGO', bioconductor=TRUE)
  geneID2GO <- readMappings(file = mappingsfile)

  results = data.frame()
    
  # iterate through patterns performing GO analysis
  for (pattern in unique(final$pattern)) {
    print(paste("GO enrichment testing for pattern", pattern))
    
    allgenes <- 1:dim(final)[1]
    names(allgenes) <- final$gene.id
    topDiffGenes <- function(row) {
      return(final$pattern[row] == pattern)
    }
    ontologies <- c('BP', 'MF', 'CC')
    
    # iterate through ontologies testing each separately
    for (ontology in ontologies) {
      print(paste("Fisher testing GO enrichment for ontology", ontology, "in condition:", pattern))
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
      allRes <- allRes[allRes$classicFisher <= alpha,]
      write.table(file=paste("results/", ontology, pattern, "GO.csv", sep='_'), 
                  x=allRes,
                  row.names=F,
                  col.names=T,
                  sep=",")
      
      # print graph of signficant nodes
      printGraph(GOdata, resultFisher, 
                 firstSigNodes = 15, 
                 fn.prefix = paste("graphs/", ontology, pattern, sep='_'), 
                 useInfo = "all", 
                 pdfSW = TRUE)
      
      # save ontology and pattern
      allRes$Ontology <- ontology
      allRes$Pattern <- pattern
      results <- rbind(results, allRes)
    }
  }
  return(results)
}