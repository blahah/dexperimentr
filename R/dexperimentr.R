### Two-way rice vs. echi
binary_DE_workflow <- function(expression_file,
                               conditions,
                               annotation_file,
                               GO_mappings_file,
                               long_mappings_file,
                               named_patterns,
                               emrounds=10) {
  ## 0. Import and prepare the data

  # load the expression counts file
  expression_data <- read.csv(expressionfile, as.is=T, sep="\t")
  
  # convert first column to row names if necessary
  if (length(conds) < ncol(expression_data)) {
    # name rows with gene IDs
    rownames(expression_data) <- expression_data[,1]
    
    # convert to matrix, discard old gene id column
    expression_data <- as.data.frame(expression_data[,-1])
  } else {
    expression_data <- as.data.frame(expression_data)
  }
  
  ## 1. Count-data quality control
  expression_data <- count_data_QC(expression_data, conditions)
  
  ## 2. Infer differential expression
  res <- infer_DE(counts=expression_data, 
                  conditions=conditions, 
                  annotation_file=annotation_file,
                  emrounds=emrounds,
                  named_patterns=named_patterns,
                  prob_cutoff=0.99)
  
  ## 3. DE quality control
  diff_expression_QC(res, conds)
  
  ## 4. Gene set enrichment analysis
  go_res <- ontology_enrichment(res,
                      GO_mappings_file, 
                      conditions,
                      long_mappings_file,
                      named_patterns,
                      ppcutoff=0.99,
                      alpha=0.01)

  return(list(de_res=res, go_res=go_res))
}

### Three-way rice-echi-maize
multiway_DE_workflow <- function() {
  ## 0. Import and prepare the data
  
  # load the expression counts file
  expression_data <- read.csv('/data2/rnaseq/echi_rice/rsem/merged_expression_file.txt', as.is=T, sep='\t')
  
  # name rows with gene IDs
  rownames <- expression_data[,1]
  
  # convert to matrix, discard old gene id column
  expression_data <- as.matrix(expression_data[,-1])
  
  ## 1. Define the conditions
  conditions <- binary_conditions(c(1, 1, 1, 2, 2, 2, 3, 3))
  
  ## 2. Count-data quality control
  count_data_QC(expression_data)
  
  ## 3. Perform differential expression analysis
  diff_exp_data <- infer_multiway_DE(expression_data)
  
  ## 4. Differential expression quality control
  multiway_DE_QC(diff_exp_data)
  
  ## 5. Gene set enrichment analysis
  ontology_enrichment(expressiondata = diff_exp_data,
                    mappingsfile = "/data/genomes/os/annotation/Rice_topGO.txt",
                    ppcutoff = 0.95,
                    alpha = 0.02)
}