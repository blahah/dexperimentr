
diff_expression_QC <- function(de_data, conditions) {
  results <- de_data['results']
  if (length(unique(conditions)) <= 2) {
    # binary
    
  } else {
    # multiway

  }
}

binary_diagnostic_plots <- function(results) {
  QQP(results.2way)
  DenNHist(results.2way)
}

multiway_diagnostic_plots <- function(results) {
  # basic diagnostic plots
  PlotPattern(pp.patterns)
  QQP(results)
  DenNHist(results)
  PlotPostVsRawFC(EBOut=results, FCOut=fc)
}

multiway_pattern_plots <- function(de_data) {
  # write out pattern plot
  library(ggplot2)
  library(reshape2)
  patterns <- as.data.frame(patterns)
  patterns$pattern <- rownames(patterns)
  patterns.melted <- melt(patterns, id='pattern')
  p <- ggplot(patterns.melted, aes(pattern, variable, fill=factor(value))) +
    geom_tile() +
    scale_fill_brewer(palette="Paired", guide=F) +
    xlab('Pattern') +
    ylab('Condition') +
    labs(fill='') +
    coord_flip()
  ggsave('expression_patterns.pdf', p, width=6, height=6)
  
  # Count the genes in each pattern (PP >= 0.95)
  final.df <- as.data.frame(final)
  names(final.df)[1] <- 'Echi1'
  final.df$rgi <- rownames(final)
  final.df <- final.df[-which(is.na(allgenes)),]
  countsig <- function(f) {
    length(which(f >= 0.95))
  }
  
  pcs <- apply(final.df[9:13],2,countsig)
  pattern.counts <- data.frame(pattern=names(pcs), count=pcs)
  q<- ggplot(pattern.counts, aes(pattern, count)) +
    geom_bar(stat='identity') +
    geom_text(aes(y=count-100, label=count, colour=T)) +
    scale_colour_brewer(palette="Paired", guide=F) +
    ggtitle('Count of genes with a posterior probability >= 0.95 of\n following each expression pattern')
  ggsave('gene_count_by_expression_pattern.pdf', q)
}