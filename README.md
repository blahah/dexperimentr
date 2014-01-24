dexperimentr
============

dexperimentr is a collection of tools for automating best-practise RNAseq differential expression workflows.


## Usage

To run the entire workflow, you need the following:

1. your count data as a dataframe, samples as named columns and genes/transcripts as named rows
2. the path to your annotation file, which is a TSV with the first column being gene/transcript identifiers, and subsequent columns containing information about those genes that you want to include in the output
3. Your long-form GO term file, which is a TSV with columns: go.id, id, term, ontology, definition (see below)

Then you just run the pipeline:

```R
library(dexperimentr)
run_DE_workflow(counts, annotation, gene_ontology)
```

An HTML report, `report.html`, will be generated with a log of the process and diagnostic plots and guidance to help you interpret your data.

In addition, the following directories will be created and populated:

- `count_qc/`: Basic quality control information about your count data
- `de_data/`: Expression data, differential expression probabilities and expression patterns
- `de_qc/`: Quality control information about the differential expression experiment
- `functional_analysis/`: Gene Ontology functional category enrichment analysis data and plots