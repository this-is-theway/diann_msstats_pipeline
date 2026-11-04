/*
========================================================================================
    DIA-NN to R Analysis Pipeline
========================================================================================
    A Nextflow pipeline for proteomics data analysis using native DIA-NN
    followed by statistical analysis and HTML report generation in R
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    INCLUDE MODULES
========================================================================================
*/

include { RAW_TO_MZML       } from './modules/diann'
include { DIANN             } from './modules/diann'
include { DIANN_TO_MSSTATS  } from './modules/diann'
include { MSSTATS_REPORT    } from './modules/msstats'



/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    log.info """\

    ===========================================================
                  D I A - N N  +  M S S T A T S
    -----------------------------------------------------------
              P R O T E O M I C S  A N A L Y S I S 
    ===========================================================
    input                : ${params.input}
    fasta                : ${params.fasta}
    experimental design  : ${params.experimental_design}
    comparisons          : ${params.comparisons_file}
    outdir               : ${params.outdir}
    threads              : ${params.diann_threads}
    ram (GB)             : ${params.diann_ram}
    MBR                  : ${params.diann_mbr}
    ===========================================================
    """
    .stripIndent()

        // -------------------------------------------------------------------
    // Validate required parameters
    // -------------------------------------------------------------------
    if (!params.input && !params.diann_report && !params.msstats_csv) {
        error "Please specify input MS files with --input (glob pattern)"
    }
    if (!params.fasta) {
        error "Please specify FASTA database with --fasta"
    }
    if (!params.experimental_design) {
        error "Please specify experimental design CSV with --experimental_design"
    }
    if (!params.rmd_template) {
        error "Please specify an RMarkdown template with --rmd_template"
    }
    if (!params.scripts) {
        error "Please specify R scripts directory with --scripts"
    }
    if (!params.comparisons_file) {
        error "Please specify comparisons CSV with --comparisons_file"
    }

    // -------------------------------------------------------------------
    // Create input channels
    // -------------------------------------------------------------------

    // Experimental design, Rmd template, and R scripts as separate channels
    design_ch       = channel.fromPath(params.experimental_design, checkIfExists: true)
    rmd_ch          = channel.fromPath(params.rmd_template, checkIfExists: true)
    scripts_ch      = channel.fromPath(params.scripts, checkIfExists: true)
    logo_ch          = channel.fromPath(params.logo_path, checkIfExists: true)
    comparisons_ch  = channel.fromPath(params.comparisons_file, checkIfExists: true)

    // Highlight genes file (optional — create empty file if not provided)
    highlight_ch = params.highlight_genes_file
        ? channel.fromPath(params.highlight_genes_file, checkIfExists: true)
        : channel.of(file("NO_HIGHLIGHT")).map { file -> 
            def f = file("${workDir}/empty_highlights.txt")
            f.text = ""
            return f
        }
    fasta_ch    = channel.fromPath(params.fasta, checkIfExists: true)

    
    if (params.msstats_csv) {
        // ---- MSstats-only: reuse a pre-computed msstats.csv ----
        log.info "MSstats-only mode: using ${params.msstats_csv} (skipping DIA-NN)"
        msstats_ch     = channel.fromPath(params.msstats_csv, checkIfExists: true)
        if (params.protein_map) {
            protein_map_ch = channel.fromPath(params.protein_map, checkIfExists: true)
        } else {
            // empty placeholder: annotation then relies on the FASTA only
            def empty_map = file("${workDir}/empty_protein_map.csv")
            empty_map.text = ""
            protein_map_ch = channel.value(empty_map)
        }
    }
    else {

                if (params.diann_report) {
                    // ---- reuse an existing DIA-NN report.parquet ----
                    log.info "Conversion-only mode: using ${params.diann_report} (skipping DIA-NN search)"
                    report_ch = channel.fromPath(params.diann_report, checkIfExists: true)
                }
                else {
                    // input files are collected and converted to mzML if .raw format
                    input_ch = channel.fromPath(params.input, type:'any')

                    if (params.convert_raw) {
                        RAW_TO_MZML(input_ch)
                        ms_files_ch = RAW_TO_MZML.out.mzml
                            .mix(input_ch.filter { !it.name.endsWith('.raw') })
                            .collect()
                    } else {
                        // DIA-NN reads .raw natively when the Thermo runtime is available
                        ms_files_ch = input_ch.collect()
                    }

                    // -------------------------------------------------------------------
                    // STEP 1: Run DIA-NN (library-free, predicted library from FASTA)
                    // -------------------------------------------------------------------
                    DIANN(
                        ms_files_ch,
                        fasta_ch
                    )

                    report_ch = DIANN.out.report
                }

                // -------------------------------------------------------------------
                // STEP 2: report.parquet -> msstats.csv
                // -------------------------------------------------------------------
                DIANN_TO_MSSTATS(
                    report_ch,
                    design_ch
                )

                msstats_ch     = DIANN_TO_MSSTATS.out.msstats_file
                protein_map_ch = DIANN_TO_MSSTATS.out.protein_map
    }

    // -------------------------------------------------------------------
    // STEP 3: MSstats analysis & HTML report
    // -------------------------------------------------------------------
    MSSTATS_REPORT(
        msstats_ch,
        rmd_ch,
        scripts_ch,
        logo_ch,
        comparisons_ch,
        highlight_ch,
        fasta_ch,
        protein_map_ch
    )

    // -------------------------------------------------------------------
    // Completion handler
    // -------------------------------------------------------------------
    workflow.onComplete = {
        println "Pipeline completed at : ${workflow.complete}"
        println "Execution status      : ${workflow.success ? 'Success' : 'Failed'}"
        println "Results directory     : ${params.outdir}"
    }
}
