rmarkdown::render(
  input       = '/home/thermo1/nextflow_pipelines/fragpipe_msstats_pipeline_dev/report/run_msstats.Rmd',
  output_file = 'msstats_report_ludovic.html',
  output_dir  = './results',
  params = list(
    input_file           = '/mnt/d/MSF/nextflow_test/MSF_02_3reps/results/diann/msstats.csv',
    protein_map_file     = '/mnt/d/MSF/nextflow_test/MSF_02_3reps/results/diann/protein_map.csv',
    comparisons_file     = '/mnt/d/MSF/nextflow_test/MSF_02_3reps/design/comparisons.csv',
    highlight_genes_file = '/mnt/d/MSF/nextflow_test/MSF_02_3reps/design/genes_to_highlight.txt',
    fasta_file           = '/home/thermo1/database_fasta/2026-03-03-decoys-contam-uniprotkb_proteome_UP000005640_AND_revi_2026_03_03.fasta.fas',
    scripts_dir          = '/home/thermo1/nextflow_pipelines/fragpipe_msstats_pipeline_dev/report/scripts',
    imputation_method    = 'ludovic',
    pvalue_cutoff        = 0.05,
    log2fc_cutoff        = 1,
    report_title         = 'MS Poject Analysis',
    report_author        = 'Jacek Szymanski'
  )
)


