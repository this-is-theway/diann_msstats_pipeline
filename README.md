# DIA-NN + MSstats proteomics pipeline

A Nextflow pipeline that takes raw mass-spectrometry files all the way to an
interactive HTML report of differentially abundant proteins:

```
 .raw / .mzML / .d files
        │
        ▼
 RAW_TO_MZML ──► DIANN ──► DIANN_TO_MSSTATS ──► MSSTATS_REPORT
 (.raw only)    search     report.parquet        MSstats statistics,
                           → msstats.csv         annotation, HTML report
```

---

## Contents

1. [Requirements](#1-requirements)
2. [Quick start](#2-quick-start)
3. [Preparing the input files](#3-preparing-the-input-files)
4. [What each process does](#4-what-each-process-does)
5. [Run modes (skipping steps)](#5-run-modes-skipping-steps)
6. [Parameters](#6-parameters)
7. [Output files](#7-output-files)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Requirements

- **Nextflow** ≥ 23.04
- **Singularity / Apptainer** with the container images below (or the tools
  installed locally)

| Process            | Container (set in `nextflow.config`)                                 | Needs                         |
|--------------------|----------------------------------------------------------------------|-------------------------------|
| `RAW_TO_MZML`      | `${container_dir}/thermorawfileparser/thermorawfileparser_1.4.5.sif`  | ThermoRawFileParser           |
| `DIANN`            | `${container_dir}/diann/diann_2.sif`                                  | DIA-NN 2.x                    |
| `DIANN_TO_MSSTATS` | `${container_dir}/msstats/msstats.sif`                                | R, MSstatsAlt, arrow, data.table |
| `MSSTATS_REPORT`   | `${container_dir}/msstats/msstats.sif`                                | R, MSstatsAlt, data.table, rmarkdown, readthedown, tidyverse, Biostrings, PCAtools, plotly, crosstalk, DT, openxlsx |

> Check the image names in `nextflow.config` (`process { withName: ... }`)
> match the files in your `container_dir`, and enable Singularity
> (`-with-singularity`, or `singularity.enabled = true` in your Nextflow config).

---

## 2. Quick start

```bash
nextflow run main.nf -profile local -with-singularity \
    --input               '/data/project/raw/*.raw' \
    --fasta               /data/fasta/human_uniprot.fasta \
    --experimental_design /data/project/design/experimental_design.csv \
    --comparisons_file    /data/project/design/comparisons.csv \
    --rmd_template        report/run_msstats.Rmd \
    --scripts             report/scripts \
    --logo_path           report/img/logo.png \
    --outdir              /data/project/results \
    --report_title        'MS Project Analysis' \
    --report_author       'Your Name'
```

- Put the `--input` glob in **single quotes**, so Nextflow expands it rather than the shell.
- Add `-resume` to continue an interrupted run without repeating finished steps.
- Instead of a long command line you can put the parameters in a file and run
  `nextflow run main.nf -profile local -params-file params.yaml`:

```yaml
# params.yaml
input:               '/data/project/raw/*.raw'
fasta:               '/data/fasta/human_uniprot.fasta'
experimental_design: '/data/project/design/experimental_design.csv'
comparisons_file:    '/data/project/design/comparisons.csv'
rmd_template:        'report/run_msstats.Rmd'
scripts:             'report/scripts'
logo_path:           'report/img/logo.png'
outdir:              '/data/project/results'
report_title:        'MS Project Analysis'
report_author:       'Your Name'
diann_threads:       16
diann_ram:           64
```

---

## 3. Preparing the input files

### 3.1 MS data files (`--input`)

DIA data in one of these formats:

| Format          | Notes                                                                        |
|-----------------|------------------------------------------------------------------------------|
| `.raw` (Thermo) | Converted to `.mzML` first (`RAW_TO_MZML`). Set `--convert_raw false` if your DIA-NN build reads `.raw` natively. |
| `.mzML`         | Used directly.                                                               |
| `.d` (Bruker)   | Used directly (timsTOF).                                                     |

### 3.2 Protein database (`--fasta`)

A **UniProt FASTA** (headers like `>sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2`).

- DIA-NN generates its own decoys. Entries starting with `rev_` or `DECOY_`
  (e.g. a FragPipe "decoys-contam" FASTA) are **removed automatically** before the search,
  so you can reuse the same FASTA.
- Contaminants (`contam_sp|...`) are kept.
- The same FASTA is used for annotation in the report (protein name, gene, organism),
  so UniProt-style headers give the best annotation.

### 3.3 Experimental design (`--experimental_design`)

CSV with one row per MS file:

```csv
file,experiment,bioreplicate,data_type
Control_01.raw,Control,1,DIA
Control_02.raw,Control,2,DIA
Control_03.raw,Control,3,DIA
Treated_01.raw,Treated,1,DIA
Treated_02.raw,Treated,2,DIA
Treated_03.raw,Treated,3,DIA
```

| Column         | Meaning                                                                          |
|----------------|----------------------------------------------------------------------------------|
| `file`         | MS file name. Matched to DIA-NN runs **without the extension**, so `Control_01.raw`, `Control_01.mzML` and `Control_01` are all equivalent. Directories are ignored. |
| `experiment`   | Condition / group name. Becomes MSstats `Condition` and is what you compare in `comparisons.csv`. |
| `bioreplicate` | Replicate number within the condition. Becomes MSstats `BioReplicate` as `<experiment>_<bioreplicate>` (e.g. `Control_1`). If empty, the run name is used. |
| `data_type`    | Informational only (e.g. `DIA`).                                                 |

Every file must appear only once. Files missing from the design are dropped (with a warning).

### 3.4 Comparisons (`--comparisons_file`)

Which conditions to compare. One row per comparison; `condition_A` is compared **against** `condition_B`
(positive log2FC = higher in `condition_A`):

```csv
condition_A,condition_B
Treated,Control
```

Names must match the `experiment` column of the design exactly.

### 3.5 Optional files

| File                       | Parameter                | Format                                              |
|----------------------------|--------------------------|-----------------------------------------------------|
| Proteins of interest       | `--highlight_genes_file` | Plain text, one gene symbol per line (e.g. `ALB`). Entry names (`ALBU_HUMAN`) or accessions (`P02768`) also work. These proteins get `Highlighted = Yes` in the report tables and Excel files. |
| Report logo                | `--logo_path`            | PNG image. The bundled `report/img/logo.png` can be used. |

---

## 4. What each process does

### `RAW_TO_MZML`
Converts Thermo `.raw` files to `.mzML` with **ThermoRawFileParser**. Runs once per
`.raw` file and is skipped for `.mzML` / `.d` input (or with `--convert_raw false`).

### `DIANN`
Runs **DIA-NN** in *library-free* mode on all files together:

1. Removes `rev_` / `DECOY_` entries from the FASTA (`target.fasta`).
2. Digests the FASTA *in silico* (trypsin, `K*,R*`) and predicts a spectral library
   with DIA-NN's deep-learning predictor.
3. Searches all runs against that library, with match-between-runs (`--reanalyse`) if `--diann_mbr true`.
4. Writes the main precursor-level report `report.parquet`, plus protein/gene matrices
   (`report.pg_matrix.tsv`, …) and the generated library. `--export-quant` adds the
   fragment intensities (`Fr.0.Quantity` … `Fr.11.Quantity`) to `report.parquet`; MSstats needs them.

Fixed modification: carbamidomethyl C (`--unimod4`). Variable: Met oxidation (`--diann_oxidation`),
N-terminal Met excision (`--met-excision`).

### `DIANN_TO_MSSTATS`
Converts `report.parquet` to the MSstats input table `msstats.csv`
(script `bin/diann_to_msstats.R`, run in the MSstats container). It follows the route
recommended by the DIA-NN developers, since MSstats does not read the DIA-NN 2.x format directly:

1. Joins the fragment intensities `Fr.0.Quantity` … `Fr.N.Quantity` with `;` into
   `Fragment.Quant.Corrected`, the column MSstats' DIA-NN importer expects.
2. Builds the annotation (`Run`, `Condition`, `BioReplicate`) from the experimental design.
3. Runs MSstats' `DIANNtoMSstatsFormat()`. This function:
   - removes decoys and applies the q-value filters: `Q.Value ≤ diann_qvalue`; with MBR
     `Lib.Q.Value ≤ diann_global_qvalue` and `Lib.PG.Q.Value ≤ diann_pg_qvalue`; without MBR the
     `Global.*` equivalents
   - keeps unique peptides only and removes features with too few measurements
   - uses `Protein.Names` (UniProt entry name, e.g. `ALBU_HUMAN`) as `ProteinName`. Each
     fragment of each precursor becomes one MSstats feature.
4. Writes the result to `msstats.csv` and the MSstats log to `MSstats_import.log`.
5. Also writes `protein_map.csv`: for each `ProteinName`, the DIA-NN protein group
   (UniProt accessions) and gene names. The report uses it as an annotation fallback.

### `MSSTATS_REPORT`
Renders `report/run_msstats.Rmd` into `msstats_report.html`:

1. **Import** – `msstats.csv` is read as-is (it is already filtered by `DIANNtoMSstatsFormat()`).
2. **Preprocessing** – log2 transformation, median normalisation (`equalizeMedians`), feature selection as
   recommended by the DIA-NN developers (`featureSubset = "all"`, `min_feature_count = 2`,
   `remove_uninformative_feature_outlier = FALSE`),
   Tukey median polish summarisation to protein level, model-based imputation of missing values.
3. **Quality control** – intensity box plots and PCA, before and after normalisation.
4. **Differential abundance** – MSstats linear model for each comparison in `comparisons.csv`;
   Benjamini–Hochberg adjusted p-values.
   - `--imputation_method ludovic` (default): each comparison is processed separately using only its two
     conditions, with its own PCA.
   - `--imputation_method AFT`: all samples are processed together and every comparison is tested in one model.
5. **Annotation** – each protein is matched to the FASTA to add `Protein_ID` (accession), `Protein_entry`,
   `Protein_name`, `Gene_name`, `Organism` and `Length`. Lookup order:
   entry name (`ALBU_HUMAN`) → accession → DIA-NN `protein_map.csv` → gene/accession reported by DIA-NN.
   For protein groups (`A_HUMAN;B_HUMAN`) the first member is used.
6. **Highlighting** – proteins from `--highlight_genes_file` get `Highlighted = Yes` (matched on gene name,
   entry name or accession, case-insensitive); all others `No`.
7. **Results** – interactive volcano and MA plots linked to a searchable table, and Excel exports.
   Filter the table's `Highlighted` column to `Yes` to see only your proteins of interest.

---

## 5. Run modes (skipping steps)

| Mode | How to start it | Runs |
|------|-----------------|------|
| Full pipeline | `--input` | `RAW_TO_MZML` → `DIANN` → `DIANN_TO_MSSTATS` → `MSSTATS_REPORT` |
| From a DIA-NN report | `--diann_report results/diann/report.parquet` | `DIANN_TO_MSSTATS` → `MSSTATS_REPORT` |
| Report only | `--msstats_csv results/diann/msstats.csv` (optionally `--protein_map results/diann/protein_map.csv`) | `MSSTATS_REPORT` |

Use the shorter modes to change statistics, comparisons or filters without repeating the
DIA-NN search. `--fasta` and `--experimental_design` are still required in every mode.
`--diann_report` needs a `report.parquet` from DIA-NN 2.x (not the older `report.tsv`).

---

## 6. Parameters

Set any parameter on the command line (`--name value`) or in a `-params-file`.
Defaults are in `nextflow.config`.

### Input / output

| Parameter              | Default     | Description |
|------------------------|-------------|-------------|
| `input`                | –           | MS files (glob in quotes), e.g. `'/data/raw/*.raw'`. **Required** for the full pipeline. |
| `fasta`                | –           | Protein FASTA database. **Required.** |
| `experimental_design`  | –           | Design CSV (section 3.3). **Required.** |
| `outdir`               | `results`   | Output directory. |

### DIA-NN search

| Parameter                | Default   | Description |
|--------------------------|-----------|-------------|
| `diann_bin`              | `diann`   | DIA-NN executable name or full path. |
| `diann_threads`          | `8`       | CPU threads for DIA-NN. |
| `diann_ram`              | `20`      | Memory for DIA-NN (GB). |
| `convert_raw`            | `true`    | Convert `.raw` to `.mzML` before DIA-NN. |
| `diann_mbr`              | `true`    | Match-between-runs (`--reanalyse`): a second pass that uses an empirical library built from all runs. Usually increases data completeness. |
| `diann_qvalue`           | `0.01`    | Precursor FDR (1 %) used in the search and in the `msstats.csv` filter. |
| `diann_missed_cleavages` | `1`       | Maximum trypsin missed cleavages. |
| `diann_min_pep_len` / `diann_max_pep_len` | `7` / `30` | Peptide length range. |
| `diann_min_pr_charge` / `diann_max_pr_charge` | `1` / `4` | Precursor charge range. |
| `diann_min_pr_mz` / `diann_max_pr_mz` | `300` / `1800` | Precursor m/z range (match your acquisition method). |
| `diann_min_fr_mz` / `diann_max_fr_mz` | `200` / `1800` | Fragment m/z range. |
| `diann_max_var_mods`     | `1`       | Maximum variable modifications per peptide. |
| `diann_oxidation`        | `true`    | Variable Met oxidation. |
| `diann_extra_args`       | `''`      | Any other DIA-NN options, passed as-is, e.g. `'--var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1'`. |

### Conversion to MSstats (`msstats.csv`)

| Parameter                | Default              | Description |
|--------------------------|----------------------|-------------|
| `diann_pg_qvalue`        | `0.01`               | Experiment-wide protein-group FDR (`Lib.PG.Q.Value` with MBR, `Global.PG.Q.Value` without). |
| `diann_global_qvalue`    | `0.01`               | Experiment-wide precursor FDR (`Lib.Q.Value` with MBR, `Global.Q.Value` without). |
| `diann_report`           | `''`                 | Existing `report.parquet`: skip the search (section 5). It must come from DIA-NN run with `--export-quant`. |

### Statistics and report

| Parameter              | Default        | Description |
|------------------------|----------------|-------------|
| `comparisons_file`     | –              | Comparisons CSV (section 3.4). **Required.** |
| `rmd_template`         | –              | `report/run_msstats.Rmd`. **Required.** |
| `scripts`              | –              | `report/scripts` directory. **Required.** |
| `logo_path`            | –              | Logo image for the report, e.g. `report/img/logo.png`. **Required.** |
| `highlight_genes_file` | `''`           | Optional list of proteins of interest (section 3.5); adds the `Highlighted` (Yes/No) column. |
| `imputation_method`    | `ludovic`      | `ludovic` (each comparison processed separately) or `AFT` (all samples together). See section 4. |
| `fdr_threshold`        | `0.05`         | Adjusted p-value cut-off for significance. |
| `log2fc_threshold`     | `1`            | \|log2FC\| cut-off for significance: `1` = 2-fold, `0.585` = 1.5-fold, `2` = 4-fold. |
| `report_title`         | `Proteomics Analysis Report` | Report title. |
| `report_author`        | `Author Name`  | Author shown in the report. |
| `msstats_csv`          | `''`           | Existing `msstats.csv`: run only the report (section 5). |
| `protein_map`          | `''`           | `protein_map.csv` to use with `msstats_csv`. |

### Execution

| Parameter       | Default | Description |
|-----------------|---------|-------------|
| `container_dir` | `/mnt/CF/MSF/nextflow_pipelines/containers/` | Folder containing the `.sif` images. |

Profiles: `-profile local` (this machine) or `-profile slurm` (each process submitted as a SLURM job).

### SLURM (`-profile slurm`)

| Parameter               | Default  | Description |
|-------------------------|----------|-------------|
| `slurm_queue`           | `''`     | Partition (`sbatch --partition`). Empty = cluster default. |
| `slurm_account`         | `''`     | Account to charge (`sbatch --account`). Empty = your default account. |
| `slurm_queue_size`      | `20`     | Maximum number of pipeline jobs queued or running at once. |
| `slurm_submit_rate`     | `10 sec` | Job submission rate limit. |
| `slurm_cluster_options` | `''`     | Any extra `sbatch` options, e.g. `'--qos=normal --exclude=node01'`. |

```bash
nextflow run main.nf -profile slurm -with-singularity -params-file params.yaml \
    --slurm_queue compute --slurm_account my_lab
```

Run Nextflow itself on a login node inside `screen`/`tmux` (or submit it as its own small job),
because it must keep running while the jobs execute.

Process resources come from the labels in `nextflow.config`: `process_low` (2 CPU / 4 GB),
`process_medium` (12 CPU / 48 GB; `RAW_TO_MZML`, `MSSTATS_REPORT`) and `process_high`
(`diann_threads` / `diann_ram`; `DIANN`). Failed tasks are retried twice.

---

## 7. Output files

```
results/
├── diann/
│   ├── report.parquet          # DIA-NN precursor-level report (main result)
│   ├── report.pg_matrix.tsv    # DIA-NN protein-group matrix (+ other *.tsv matrices / stats)
│   ├── report-lib.parquet      # predicted spectral library
│   ├── target.fasta            # FASTA used for the search (decoys removed)
│   ├── msstats.csv             # MSstats input (fragment level)
│   ├── MSstats_import.log      # DIANNtoMSstatsFormat log
│   └── protein_map.csv         # ProteinName → accession / gene (DIA-NN)
├── msstats/
│   ├── msstats_report.html     # interactive report – open in a browser
│   ├── msstats_report_files/   # report resources (keep next to the .html)
│   ├── *.xlsx                  # comparison results, protein/feature-level data
│   ├── pepetide_level_data.csv.gz  # feature-level data instead of .xlsx when > 1 million rows
│   └── logs/                   # MSstats log files
└── pipeline_info/              # Nextflow execution report, timeline, trace, DAG
```

Keep `msstats_report.html` and `msstats_report_files/` together when moving or sharing the report.

**Main result columns** (per comparison):

| Column | Meaning |
|--------|---------|
| `Protein` | Protein identifier from DIA-NN (entry name; protein groups joined by `;`) |
| `Protein_ID`, `Protein_entry`, `Protein_name`, `Gene_name`, `Organism`, `Length` | Annotation from the FASTA |
| `Highlighted` | `Yes` for proteins from `--highlight_genes_file`, otherwise `No` |
| `log2FC` | log2 fold change, `condition_A` vs `condition_B` |
| `SE`, `Tvalue`, `DF` | Model statistics |
| `pvalue`, `adj.pvalue` | Raw and Benjamini–Hochberg adjusted p-values |
| `Significance` | `Up` / `Down` / `NS` using `fdr_threshold` and `log2fc_threshold` |

---

## 8. Troubleshooting

| Problem | What to check |
|---------|---------------|
| `no Fr.N.Quantity columns in report.parquet` | The report was made without `--export-quant`. Re-run DIA-NN (the pipeline adds it), or add `--export-quant` to your own DIA-NN run. |
| `no rows left after matching DIA-NN runs to the experimental design` | The `file` column must match the MS file names (extension does not matter). |
| `Condition X not found in available conditions` | Names in `comparisons.csv` must match `experiment` in the design exactly (case, spaces). |
| Many proteins without `Gene_name` | The FASTA is not in UniProt format. The DIA-NN gene names from `protein_map.csv` are used where available. |
| Container not found | Image names in `nextflow.config` vs files in `container_dir`; Singularity enabled. |
| DIA-NN out of memory | Increase `--diann_ram`, or reduce `--diann_threads`. |
| Want to change statistics only | Re-run with `--msstats_csv results/diann/msstats.csv --protein_map results/diann/protein_map.csv`, or simply use `-resume`. |

### Citations

- **DIA-NN** – Demichev V, et al. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nature Methods*, 17, 41–44.
- **MSstats** – Choi M, et al. (2014). MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments. *Bioinformatics*, 30(17), 2524–2526.
- **ThermoRawFileParser** – Hulstaert N, et al. (2020). ThermoRawFileParser: modular, scalable, and cross-platform RAW file conversion. *J. Proteome Res.*, 19(1), 537–542.
