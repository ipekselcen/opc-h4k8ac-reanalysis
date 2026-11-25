#!/usr/bin/env python3
"""
Fetch SRA run metadata for every GEO series in this project and write
DRAFT sample sheets for you to review.

Outputs
  config/metadata/<GSE>_runs.tsv   raw run table from SRA (keep for provenance)
  config/samples_rnaseq.draft.tsv  draft RNA-seq sheet (review, then rename)
  config/samples_chipseq.draft.tsv draft ChIP-seq sheet (review, then rename)

Labels (age, treatment, input/ChIP) are GUESSED from sample titles.
Anything the script can't place is marked CHECK. Nothing runs until you
review the drafts and save them as samples_rnaseq.tsv / samples_chipseq.tsv.

Usage:  python scripts/fetch_metadata.py
"""
import re
import sys
from pathlib import Path

import pandas as pd
from pysradb.sraweb import SRAweb

# GEO series -> (assay, cell type). Cell type comes from the study design,
# not from the titles, so it's set here once.
DATASETS = {
    "GSE263808": ("rnaseq", "O4"),       # O4+ OPCs, P5/P60, DMSO/KATi/SUV420i
    "GSE66047":  ("rnaseq", "PDGFRa"),   # PDGFRa-GFP sorted OPCs, P2 (SuperSeries: also other assays)
    "GSE122446": ("rnaseq", "PDGFRa"),   # adult PDGFRa+ OPCs
    "GSE240522": ("chipseq", "PDGFRa"),  # H4K8ac ChIP + inputs, P5/P60
}
KEEP_STRATEGY = {"rnaseq": "RNA-Seq", "chipseq": "ChIP-Seq"}

OUT = Path("config")
META = OUT / "metadata"


def find_gsm_column(df):
    """pysradb column names vary by version; find the column holding GSM IDs."""
    for col in df.columns:
        if df[col].astype(str).str.fullmatch(r"GSM\d+").any():
            return col
    sys.exit("Could not find a GSM column in the SRA metadata.")


def fetch_runs(db, gse):
    """Return one row per GSM with its runs joined by ';'."""
    srps = db.gse_to_srp(gse)["study_accession"].dropna().unique()
    meta = pd.concat([db.sra_metadata(s, detailed=True) for s in srps])
    meta.to_csv(META / f"{gse}_runs.tsv", sep="\t", index=False)

    gsm_col = find_gsm_column(meta)
    title_col = next(c for c in ["experiment_title", "sample_title"] if c in meta)
    return (
        meta.groupby(gsm_col)
        .agg(
            runs=("run_accession", lambda r: ";".join(sorted(set(r)))),
            layout=("library_layout", "first"),
            strategy=("library_strategy", "first"),
            title=(title_col, "first"),
        )
        .reset_index()
        .rename(columns={gsm_col: "gsm"})
    )


def guess_age(title):
    m = re.search(r"\bP(\d+)\b", title, re.I)
    if m:
        return f"P{m.group(1)}"
    if re.search(r"neonat", title, re.I):
        return "CHECK_neonatal"
    if re.search(r"adult", title, re.I):
        return "CHECK_adult"
    return "CHECK"


def guess_treatment(title):
    if re.search(r"DMSO|vehicle", title, re.I):
        return "DMSO"
    if re.search(r"KAT", title, re.I):
        return "KATi"
    if re.search(r"SUV", title, re.I):
        return "SUV420i"
    return "untreated"


def guess_rep(title):
    m = re.search(r"rep(?:licate)?\s*(\d+)", title, re.I)
    return m.group(1) if m else ""


def build_rnaseq(df):
    df = df.copy()
    df["age"] = df["title"].map(guess_age)
    df["treatment"] = df["title"].map(guess_treatment)
    df["include"] = "CHECK"
    # replicate = order within each group, unless the title states one
    df["replicate"] = df["title"].map(guess_rep)
    grp = ["cell_type", "age", "treatment"]
    auto = df.groupby(grp).cumcount() + 1
    df["replicate"] = df["replicate"].where(df["replicate"] != "", auto.astype(str))
    df["sample"] = df[grp].agg("_".join, axis=1) + "_rep" + df["replicate"]
    cols = ["sample", "include", "dataset", "gsm", "runs", "layout",
            "cell_type", "age", "treatment", "replicate", "title"]
    return df[cols]


def build_chipseq(df):
    df = df.copy()
    df["age"] = df["title"].map(guess_age)
    df["mark"] = df["title"].str.contains("input", case=False).map(
        {True: "input", False: "H4K8ac"})
    df["replicate"] = df["title"].map(guess_rep)
    df["include"] = "yes"
    df["sample"] = (df["cell_type"] + "_" + df["age"] + "_" + df["mark"]
                    + "_rep" + df["replicate"])
    # each ChIP sample points to the input from the same age and replicate
    inputs = df[df["mark"] == "input"].set_index(["age", "replicate"])["sample"]
    df["control"] = [
        inputs.get((a, r), "CHECK") if m != "input" else ""
        for a, r, m in zip(df["age"], df["replicate"], df["mark"])
    ]
    cols = ["sample", "include", "dataset", "gsm", "runs", "layout",
            "cell_type", "age", "mark", "replicate", "control", "title"]
    return df[cols]


def main():
    META.mkdir(parents=True, exist_ok=True)
    db = SRAweb()
    tables = {"rnaseq": [], "chipseq": []}

    for gse, (assay, cell_type) in DATASETS.items():
        print(f"Fetching {gse} ...")
        runs = fetch_runs(db, gse)
        runs = runs[runs["strategy"] == KEEP_STRATEGY[assay]]
        runs.insert(0, "dataset", gse)
        runs.insert(1, "cell_type", cell_type)
        print(f"  {len(runs)} {KEEP_STRATEGY[assay]} samples")
        tables[assay].append(runs)

    rna = build_rnaseq(pd.concat(tables["rnaseq"]))
    chip = build_chipseq(pd.concat(tables["chipseq"]))
    rna.to_csv(OUT / "samples_rnaseq.draft.tsv", sep="\t", index=False)
    chip.to_csv(OUT / "samples_chipseq.draft.tsv", sep="\t", index=False)

    n_check = (rna.astype(str).apply(lambda c: c.str.contains("CHECK")).any(axis=1).sum()
               + chip.astype(str).apply(lambda c: c.str.contains("CHECK")).any(axis=1).sum())
    print(f"\nWrote draft sheets. {n_check} rows contain CHECK and need your review.")


if __name__ == "__main__":
    main()
