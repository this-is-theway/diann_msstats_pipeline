/*
========================================================================================
    DIA-NN Process Module
========================================================================================
    Native DIA-NN library-free search (in silico predicted spectral library from FASTA)
    and conversion of report.parquet into an MSstats input CSV.
    --export-quant adds the fragment intensities (Fr.N.Quantity) to report.parquet,
    which MSstats uses as features (DIA-NN developers' recommended route).
----------------------------------------------------------------------------------------
*/



process RAW_TO_MZML {
    tag "${raw_file.simpleName}"
    label 'process_medium'

    input:
    path raw_file

    output:
    path "*.mzML", emit: mzml

    when:
    raw_file.name.endsWith('.raw')

    script:
    """
    thermorawfileparser \
        -i=${raw_file} \
        -o=. \
        -f=2
    """
}



process DIANN {
    tag "DIA-NN Analysis"
    label 'process_high'
    publishDir "${params.outdir}/diann", mode: 'copy'

    input:
    path ms_files
    path fasta

    output:
    path "report.parquet",          emit: report
    path "report*.tsv",             emit: matrices, optional: true
    path "report-lib.parquet",      emit: library,  optional: true
    path "*.log.txt",               emit: logs,     optional: true
    path "target.fasta",            emit: fasta

    script:
    def files   = (ms_files instanceof List ? ms_files : [ms_files])
    def f_args  = files.collect { "--f ${it}" }.join(' ')
    def var_mod = params.diann_oxidation ? "--var-mod UniMod:35,15.994915,M" : ""
    def mbr     = params.diann_mbr ? "--reanalyse" : ""
    """
    # DIA-NN generates its own decoys: drop FragPipe-style decoys (rev_/DECOY_)
    # so they are not searched as targets
    awk '/^>/{keep = !(\$0 ~ /^>(rev_|DECOY_)/)} keep' ${fasta} > target.fasta
    echo "FASTA entries: \$(grep -c '^>' ${fasta}) -> \$(grep -c '^>' target.fasta) after decoy removal"

    ${params.diann_bin} \\
        ${f_args} \\
        --fasta target.fasta \\
        --fasta-search \\
        --predictor \\
        --gen-spec-lib \\
        --out report.parquet \\
        --out-lib report-lib.parquet \\
        --threads ${task.cpus} \\
        --verbose 1 \\
        --qvalue ${params.diann_qvalue} \\
        --matrices \\
        --cut 'K*,R*' \\
        --missed-cleavages ${params.diann_missed_cleavages} \\
        --min-pep-len ${params.diann_min_pep_len} \\
        --max-pep-len ${params.diann_max_pep_len} \\
        --min-pr-charge ${params.diann_min_pr_charge} \\
        --max-pr-charge ${params.diann_max_pr_charge} \\
        --min-pr-mz ${params.diann_min_pr_mz} \\
        --max-pr-mz ${params.diann_max_pr_mz} \\
        --min-fr-mz ${params.diann_min_fr_mz} \\
        --max-fr-mz ${params.diann_max_fr_mz} \\
        --met-excision \\
        --unimod4 \\
        --var-mods ${params.diann_max_var_mods} \\
        ${var_mod} \\
        ${mbr} \\
        --rt-profiling \\
        --export-quant \\
        ${params.diann_extra_args}
    """
}



process DIANN_TO_MSSTATS {
    tag "report.parquet -> msstats.csv"
    label 'process_medium'
    publishDir "${params.outdir}/diann", mode: 'copy'

    input:
    path report               // DIA-NN report.parquet (with --export-quant fragment columns)
    path experimental_design  // CSV: file, experiment, bioreplicate, data_type

    output:
    path "msstats.csv",         emit: msstats_file
    path "protein_map.csv",     emit: protein_map
    path "MSstats_import.log",  emit: log, optional: true

    script:
    // Fragment-level conversion with MSstats' DIANNtoMSstatsFormat()
    """
    diann_to_msstats.R \\
        --report=${report} \\
        --design=${experimental_design} \\
        --out=msstats.csv \\
        --protein-map=protein_map.csv \\
        --mbr=${params.diann_mbr ? 'TRUE' : 'FALSE'} \\
        --qvalue=${params.diann_qvalue} \\
        --pg-qvalue=${params.diann_pg_qvalue} \\
        --global-qvalue=${params.diann_global_qvalue} \\
        --log=MSstats_import.log
    """
}
