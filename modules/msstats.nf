/*
========================================================================================
    MSstats Report Process Module
========================================================================================
    Runs an RMarkdown-based MSstats analysis on DIA-NN output (msstats.csv).
----------------------------------------------------------------------------------------
*/

process MSSTATS_REPORT {
    tag "MSstats Analysis"
    label 'process_medium'
    publishDir "${params.outdir}/msstats",      mode: 'copy', pattern: '*.{html,xlsx,csv.gz}'
    publishDir "${params.outdir}/msstats", mode: 'copy', pattern: '*_files/**'
    publishDir "${params.outdir}/msstats/logs",  mode: 'copy', pattern: '*.log'

    container 'bioconductor/bioconductor_docker:devel'

    input:
    path msstats_file          // msstats.csv from DIANN_TO_MSSTATS
    path rmd_template          // run_msstats.Rmd
    path scripts               // scripts/ directory (R scripts + comparison_template.Rmd)
    path logo_path
    path comparisons_file      // comparisons.csv (condition_A, condition_B)
    path highlight_genes_file  // genes_to_highlight.txt (optional, can be empty)
    path fasta
    path protein_map           // protein_map.csv from DIANN_TO_MSSTATS (can be empty)

    output:
    path "msstats_report.html",      emit: report
    path "*_files/**",               emit: report_files, optional: true
    path "*.xlsx",                   emit: excel, optional: true
    path "*.csv.gz",                 emit: tables, optional: true
    path "*.log",                    emit: logs,  optional: true

    script:
    """
    Rscript -e "rmarkdown::render(
        input       = '${rmd_template}',
        output_file = 'msstats_report.html',
        output_dir  = '\$PWD',
        params = list(
            input_file           = '\$PWD/${msstats_file}',
            comparisons_file     = '\$PWD/${comparisons_file}',
            highlight_genes_file = '\$PWD/${highlight_genes_file}',
            scripts_dir          = '\$PWD/${scripts}',
            log_dir              = '\$PWD',
            logo_path            = '\$PWD/${logo_path}',
            fasta_file           = '\$PWD/${fasta}',
            protein_map_file     = '\$PWD/${protein_map}',
            imputation_method    = '${params.imputation_method}',
            pvalue_cutoff        = ${params.fdr_threshold},
            log2fc_cutoff        = ${params.log2fc_threshold},
            report_title         = '${params.report_title}',
            report_author        = '${params.report_author}'
        )
    )"
    """
}
