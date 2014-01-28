### Two-way rice vs. echi
binary_DE_workflow <- function(expression_file,
                               conditions,
                               annotation_file,
                               GO_mappings_file) {
  ## 0. Import and prepare the data

  # load the expression counts file
  expression_data <- read.table(expression_file)
  
  # convert first column to row names if necessary
  if (conditions < ncol(expression_data)) {
    # name rows with gene IDs
    rownames <- expression_data[,1]
    
    # convert to matrix, discard old gene id column
    expression_data <- as.data.frame(expression_data[,-1])
  } else {
    expression_data <- as.data.frame(expression_data)
  }
  
  ## 1. Count-data quality control
  expression_data <- count_data_QC(expression_data)
  
  ## 2. Infer differential expression
  res <- infer_DE(counts=expression_data, 
                  conditions=conds, 
                  annotation_file=annotation_file,
                  emrounds=25)
  
  ## 3. DE quality control
  diff_expression_QC(res, conds)
  
  ## 4. Gene set enrichment analysis
  ontology_enrichment(expressiondata = diff_exp_data,
                    mappingsfile = GO_mappings_file,
                    ppcutoff = 0.95,
                    alpha = 0.02)
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