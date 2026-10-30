#!/usr/bin/env python3
"""
Convert a DIA-NN report.parquet into an MSstats input CSV.

Output (msstats.csv) uses the standard 10-column MSstats layout, the same one
FragPipe writes, so it can be read with FragPipetoMSstatsFormat() or passed to
MSstats directly:

    ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
    IsotopeLabelType, Condition, BioReplicate, Run, Intensity

A second file (protein_map.csv) links every ProteinName to the DIA-NN protein
group / accession(s) / gene names, used as a fallback in the R annotation step.

Filtering follows MSstats' DIANNtoMSstatsFormat():
    Q.Value         <= --qvalue
    PG.Q.Value      <= --pg-qvalue
    with --mbr:  Lib.Q.Value / Lib.PG.Q.Value <= --global-qvalue
    otherwise:   Global.Q.Value / Global.PG.Q.Value <= --global-qvalue
(each filter is applied only if the column is present in the report).
"""

import argparse
import os
import sys

import pandas as pd
import pyarrow.parquet as pq

MS_EXTENSIONS = (".raw", ".mzml", ".mzxml", ".d", ".wiff", ".wiff2", ".dia")


def run_name(path):
    """File name without directory and MS extension, as DIA-NN reports Run."""
    name = os.path.basename(str(path).strip().rstrip("/\\"))
    lower = name.lower()
    for ext in MS_EXTENSIONS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def read_design(path):
    design = pd.read_csv(path, dtype=str).fillna("")
    design.columns = [c.strip().lower() for c in design.columns]
    for col in ("file", "experiment"):
        if col not in design.columns:
            sys.exit(f"ERROR: experimental design is missing column '{col}'")
    if "bioreplicate" not in design.columns:
        design["bioreplicate"] = ""

    design["Run"] = design["file"].map(run_name)
    design["Condition"] = design["experiment"].str.strip()
    rep = design["bioreplicate"].str.strip()
    # BioReplicate must be unique per sample (otherwise MSstats treats the
    # design as repeated measures and the QC pivot collapses samples).
    design["BioReplicate"] = [
        f"{c}_{r}" if r else run for c, r, run in zip(design["Condition"], rep, design["Run"])
    ]
    if design["Run"].duplicated().any():
        dups = design.loc[design["Run"].duplicated(), "Run"].tolist()
        sys.exit(f"ERROR: duplicated runs in experimental design: {dups}")
    return design[["Run", "Condition", "BioReplicate"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, help="DIA-NN report.parquet")
    ap.add_argument("--design", required=True, help="CSV: file, experiment, bioreplicate[, data_type]")
    ap.add_argument("--out", default="msstats.csv")
    ap.add_argument("--protein-map", default="protein_map.csv")
    ap.add_argument("--protein-column", default="Protein.Names",
                    help="DIA-NN column used as MSstats ProteinName (default: Protein.Names, "
                         "i.e. UniProt entry names; falls back to Protein.Group when empty)")
    ap.add_argument("--intensity-column", default="Precursor.Quantity",
                    help="Precursor.Quantity (default, MSstats normalises itself) or Precursor.Normalised")
    ap.add_argument("--qvalue", type=float, default=0.01)
    ap.add_argument("--pg-qvalue", type=float, default=0.01)
    ap.add_argument("--global-qvalue", type=float, default=0.01)
    ap.add_argument("--mbr", action="store_true", help="DIA-NN was run with --reanalyse (MBR)")
    ap.add_argument("--proteotypic-only", action="store_true", help="keep only proteotypic precursors")
    args = ap.parse_args()

    available = set(pq.read_schema(args.report).names)

    required = ["Run", "Protein.Group", "Modified.Sequence", "Precursor.Charge",
                args.intensity_column, args.protein_column]
    missing = [c for c in required if c not in available]
    if missing:
        sys.exit(f"ERROR: columns missing from {args.report}: {missing}")

    q_filters = [("Q.Value", args.qvalue), ("PG.Q.Value", args.pg_qvalue)]
    if args.mbr:
        q_filters += [("Lib.Q.Value", args.global_qvalue), ("Lib.PG.Q.Value", args.global_qvalue)]
    else:
        q_filters += [("Global.Q.Value", args.global_qvalue), ("Global.PG.Q.Value", args.global_qvalue)]
    q_filters = [(c, t) for c, t in q_filters if c in available]

    optional = [c for c in ("Protein.Ids", "Protein.Names", "Genes", "Decoy", "Proteotypic") if c in available]
    columns = list(dict.fromkeys(required + optional + [c for c, _ in q_filters]))

    df = pq.read_table(args.report, columns=columns).to_pandas()
    print(f"Read {len(df):,} precursor rows from {args.report}")

    keep = pd.Series(True, index=df.index)
    if "Decoy" in df.columns:
        keep &= df["Decoy"].fillna(0).astype(int) == 0
    for col, thr in q_filters:
        keep &= df[col] <= thr
        print(f"  filter {col} <= {thr}")
    if args.proteotypic_only and "Proteotypic" in df.columns:
        keep &= df["Proteotypic"].fillna(0).astype(int) == 1
    keep &= df[args.intensity_column].notna() & (df[args.intensity_column] > 0)
    df = df.loc[keep].copy()
    print(f"{len(df):,} rows after filtering")

    # DIA-NN 2.x reports Run without extension; strip it anyway in case of paths
    df["Run"] = df["Run"].map(run_name)

    protein = df[args.protein_column].astype("string").fillna("").str.strip()
    no_name = protein == ""
    if no_name.any():
        print(f"  {no_name.sum():,} rows without {args.protein_column}; using Protein.Group instead")
        protein = protein.where(~no_name, df["Protein.Group"].astype("string"))
    df["ProteinName"] = protein

    design = read_design(args.design)
    report_runs = set(df["Run"].unique())
    design_runs = set(design["Run"])
    if report_runs - design_runs:
        print(f"WARNING: runs not in experimental design (dropped): {sorted(report_runs - design_runs)}",
              file=sys.stderr)
    if design_runs - report_runs:
        print(f"WARNING: design runs absent from DIA-NN report: {sorted(design_runs - report_runs)}",
              file=sys.stderr)
    df = df.merge(design, on="Run", how="inner")
    if df.empty:
        sys.exit("ERROR: no rows left after matching DIA-NN runs to the experimental design "
                 "(check the 'file' column against DIA-NN Run names)")

    out = pd.DataFrame({
        "ProteinName": df["ProteinName"],
        "PeptideSequence": df["Modified.Sequence"],
        "PrecursorCharge": df["Precursor.Charge"].astype(int),
        "FragmentIon": pd.NA,
        "ProductCharge": pd.NA,
        "IsotopeLabelType": "L",
        "Condition": df["Condition"],
        "BioReplicate": df["BioReplicate"],
        "Run": df["Run"],
        "Intensity": df[args.intensity_column],
    })
    out.to_csv(args.out, index=False, na_rep="NA")
    print(f"Wrote {args.out}: {len(out):,} rows, {out['ProteinName'].nunique():,} proteins, "
          f"{out['Run'].nunique()} runs")

    map_cols = {"ProteinName": "ProteinName", "Protein.Group": "Protein_group"}
    if "Protein.Ids" in df.columns:
        map_cols["Protein.Ids"] = "Protein_ids"
    if "Genes" in df.columns:
        map_cols["Genes"] = "Genes_diann"
    pmap = df[list(map_cols)].rename(columns=map_cols).drop_duplicates("ProteinName")
    pmap.to_csv(args.protein_map, index=False)
    print(f"Wrote {args.protein_map}: {len(pmap):,} proteins")


if __name__ == "__main__":
    main()
