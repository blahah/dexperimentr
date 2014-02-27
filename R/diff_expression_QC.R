
diff_expression_QC <- function(de_data, conditions) {
  # create a directory for the outputs
  wd <- getwd()
  dir.create('de_qc', showWarnings=FALSE)
  setwd('./de_qc')

  # basic universal diagnostics
  basic_diagnostics(de_data, conditions)
  
  # run the QC analysis
  results <- de_data[['results']]
  if (length(unique(conditions)) <= 2) {
    # binary
    binary_diagnostic_plots(results)
  } else {
    # multiway
    multiway_diagnostic_plots(results)
  }
  
  setwd(wd)
}

basic_diagnostics <- function(de_data, conditions) {
  final <- de_data[['final']]
  # plot (normalised) count sanity diagnostics
  counts <- final[,2:(1+length(conditions))]
  plot_correlation_matrix(counts)
  plot_pca(counts, final$pattern)
  plot_count_collision(counts)
  # plot probability diagnostics
  if (length(unique(conditions)) == 2) {
    probs <- final[,'PPDE']
  } else {
    probs <- final[,de_data[['prob_cols']]]
  }
  plot_prob_dist(probs)
}

binary_diagnostic_plots <- function(results) {
  
  FC = PostFC(results)
  pdf('posterior_FC_vs_raw_FC.pdf')
  PlotPostVsRawFC(results,FC)
  dev.off()
  
  pdf('quantile_plots.pdf')  
  par(mfrow=c(1,1))
  QQP(results, GeneLevel=T)
  dev.off()
  
  pdf('hyperparameter_histograms_vs_prob_density.pdf')  
  par(mfrow=c(1,1))
  DenNHist(results, GeneLevel=T)
  dev.off()
  
  par(mfrow=c(1,1))
  plot_convergence(results)
}

multiway_diagnostic_plots <- function(results) {
  # basic diagnostic plots
  PlotPattern(pp.patterns)
  QQP(results)
  DenNHist(results)
  PlotPostVsRawFC(EBOut=results, FCOut=fc)
  plot_convergence(results)
}

multiway_pattern_plots <- function(de_data) {
  # write out pattern plot
  get_package('ggplot2')
  get_package('reshape2')
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

plot_convergence <- function(results) {
  get_package('ggplot2')
  get_package('reshape2')
  d <- data.frame(alpha=results$Alpha,
                  beta=results$Beta,
                  iteration=1:(length(results$Alpha)))
  d <- melt(d, id='iteration')
  p <- qplot(data=d, x=iteration, 
              y=value, colour=variable, 
              geom="line", ymin=0) + 
          theme_bw() +
          ggtitle('Convergence of hyperparameters')
  ggsave(plot=p, filename="hyperparameter_convergence.pdf")
}

plot_pca <- function(counts, pattern) {
  get_package('ggplot2')
  get_package('gridExtra')
  pca <- prcomp(counts, scale=T)
  scores <- as.data.frame(pca$x[,1:3])
  scores$pattern <- pattern
  pdf('pca.pdf')
  pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=pattern) + 
              theme(legend.position="none", plot.margin = unit(c(2, 1, 1, 1), "cm"))
  pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=pattern) + 
              theme(legend.direction = "horizontal", 
                    legend.position = c(0.1, 1.05),
                    plot.margin = unit(c(2, 1, 1, 1), "cm"))
  pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=pattern) + 
              theme(legend.position="none",
                    plot.margin = unit(c(2, 1, 1, 1), "cm"))
  print(grid.arrange(pc1.2, pc1.3, pc2.3, ncol=3, nrow=1))
  dev.off()
}

plot_correlation_matrix <- function(counts) {
  get_package('ggplot2')
  get_package('reshape2')
  c <- melt(cor(counts))
  p <- ggplot(data=c, aes_string(x=names(c)[1], y=names(c)[2], fill="value")) +
   geom_tile() +
   xlab('') +
   ylab('')
  ggsave(plot=p, filename='correlation_matrix.pdf')
}

plot_count_collision <- function(counts) {
  get_package('ggplot2')
  get_package('GGally')
  p <- ggpairs(counts)
  pdf('all_vs_all_counts_scatter.pdf')
  print(p)
  dev.off()
}

plot_prob_dist <- function(probs) {
  get_package('ggplot2')
  get_package('reshape2')
  d <- melt(probs)
  if (ncol(d) == 1) {
    p <- ggplot(data=d, aes(x=value)) + geom_density()
  } else {
    p <- ggplot(data=d, aes(x=value, colour=variable)) + geom_density()
  }
  print(p)
  ggsave(plot=p, filename='prob_dist.pdf')
}