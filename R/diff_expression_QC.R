
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
    multiway_pattern_plots(de_data)
  }

  setwd(wd)
}

basic_diagnostics <- function(de_data, conditions) {
  final <- de_data[['final']]
  # plot (normalised) count sanity diagnostics
  counts <- de_data$results@data
  try(plot_correlation_matrix(counts))
  try(plot_pca(counts, final$pattern))
  try(plot_count_collision(counts))
  # plot probability diagnostics
  if (length(unique(conditions)) == 2) {
    probs <- final[,'PPDE']
  } else {
    prob_cols <- c(de_data[['ee_prob_col']], de_data[['de_prob_cols']])
    print(prob_cols)
    print(summary(final))
    probs <- final[,prob_cols]
  }
  try(plot_prob_dist(probs))
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
  # TODO: plot all pattern shapes along with how many
  # genes follow each pattern

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

multiway_pattern_plots <- function(de_data) {
  # write out pattern plot
  get_package('ggplot2')
  get_package('reshape2')

  num_patterns <- de_data[['num_patterns']]
  num_patterns <- sapply(num_patterns, function(x) strsplit(x, '_', fixed=T))

  len.cond <- sapply(num_patterns, function(x) length(x) > 1)
  num_patterns <- num_patterns[len.cond]
  num_patterns <- sapply(num_patterns, as.numeric)

  expanded_np <- t(as.data.frame(num_patterns))
  print(expanded_np)
  names(expanded_np) <- unique(de_data[['conditions']])
  if (is.null(names(de_data[['word_patterns']]))) {
    names(de_data[['word_patterns']]) <- de_data[['word_patterns']]
  }
  patterns <- data.frame(pattern=de_data[['word_patterns']],
                         original=names(de_data[['word_patterns']]),
                         expanded_np)
  patterns.melted <- melt(patterns, id=c("pattern", "original"))
  patterns.melted$pattern <- clean_strings(patterns.melted$pattern)
  patterns.melted$variable <- underscore_to_space(patterns.melted$variable)
  patterns.melted$value <- as.numeric(patterns.melted$value)
  p <- ggplot(patterns.melted, aes(variable, factor(value))) +
          geom_line(aes(colour=original, group=original)) +
          facet_grid(pattern~., scales="free_y") +
          xlab('Condition') +
          ylab('Level') +
          guides(colour=FALSE) +
          scale_x_discrete(expand=c(0.1,0.1)) +
          theme_bw(base_size=9)
  print(p)
  ggsave('expression_patterns.pdf', p,
          width=length(names(expanded_np))*0.8,
          height=length(unique(patterns$pattern))*0.8)

  # Plot number of genes in each pattern
  pattern.counts <- melt(table(de_data[['final']]$pattern))
  names(pattern.counts) <- c('pattern', 'count')
  pattern.counts <- pattern.counts[-which(pattern.counts$pattern == "no significant pattern")]
  pattern.counts$pattern <- clean_strings(pattern.counts$pattern)
  q <- ggplot(pattern.counts, aes(reorder(pattern, -count), count)) +
    geom_bar(stat='identity') +
    # geom_text(aes(y=count+(max(count)/80), label=count, colour=T)) +
    # ggtitle('Count of genes following each expression pattern') +
    theme_bw(base_size=9) +
    # theme(axis.text.x = element_text(angle=90, hjust=1,
    #                                  vjust=0.5)) +
    xlab('pattern') + ylab('number of genes')
  print(q)
  ggsave('gene_count_by_expression_pattern.pdf', q,
         width=length(pattern.counts$pattern)*0.8,
         height=2)
}

plot_convergence <- function(results) {
  get_package('ggplot2')
  get_package('reshape2')
  d <- data.frame(alpha=results$Alpha[,1],
                  beta=results$Beta[,1],
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
  names(c)[3] <- 'correlation'
  p <- ggplot(data=c,
              aes_string(x=names(c)[1], y=names(c)[2], fill="correlation")) +
   geom_tile() +
   xlab('') +
   ylab('')
  ggsave(plot=p, filename='correlation_matrix.pdf')
}

plot_count_collision <- function(counts) {
  get_package('ggplot2')
  get_package('GGally')
  p <- ggpairs(counts)
  pdf('all_vs_all_counts_scatter.pdf', width=10, height=10)
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
    p <- ggplot(data=d, aes(x=value, colour=variable)) +
            geom_density() +
            scale_y_log10() +
            ylab("log density") +
            ylim(1e-02, 1e+04)
  }
  print(p)
  ggsave(plot=p, filename='prob_dist.pdf')
}

underscore_to_space <- function(x) {
  return(sapply(x, function(y) {
    a <- gsub(as.character(y),
            pattern="_",
            replacement=" ")
    return(a)
  }))
}

wrap_long_strings <- function(x, length=25) {
  sapply(x, function(y) {
    z <- gsub('(.{1,12})(\\s)', '\\1\n', y)
    print(z)
  })
}

clean_strings <- function(x) {
  x <- underscore_to_space(x)
  x <- wrap_long_strings(x)
  return(x)
}
