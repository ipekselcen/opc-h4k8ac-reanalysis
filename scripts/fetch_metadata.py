#!/usr/bin/env python3
"""
Fetch SRA run metadata for every GEO series in this project and write
the sample sheets the workflow reads.

Outputs
  config/metadata/<GSE>_runs.tsv   raw run table from SRA (kept for provenance)
  config/samples_rnaseq.tsv        one row per RNA-seq sample
  config/samples_chipseq.tsv       one row per ChIP-seq sample (ChIP + input)

Labels come from explicit rules written after reading every GEO title by
hand (see label_rnaseq). A title that matches no rule stops the script, so
a new or renamed sample can never be labeled silently.

Usage:  python scripts/fetch_metadata.py
"""
import re
import sys
from pathlib import Path

import pandas as pd
from pysradb.sraweb import SRAweb

# GEO series -> assay
DATASETS = {
    "GSE263808": "rnaseq",   # O4+ OPCs, P5/P60, untreated/DMSO/KATi/SUV420i (Dansu, Selcen et al. 2024)
    "GSE66047":  "rnaseq",   # SuperSeries (Moyon et al. 2016): P2 OPCs, P18 OLs, Dnmt1 sub-study
    "GSE122446": "rnaseq",   # Tet1 study: P60 OPCs and P60 pre-myelinating OLs
    "GSE240522": "chipseq",  # H4K8ac ChIP + inputs, PDGFRa+ OPCs, P5/P60
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


def clean_title(title):
    """'GSM123: DMSO O4+ aOPC rep1; Mus musculus; RNA-Seq' -> 'DMSO O4+ aOPC rep1'"""
    return re.sub(r"^GSM\d+:\s*|;\s*Mus musculus;.*$", "", title).strip()


def fetch_runs(db, gse):
    """Return one row per GSM, with its runs joined by ';'."""
    srps = db.gse_to_srp(gse)["study_accession"].dropna().unique()
    meta = pd.concat([db.sra_metadata(s, detailed=True) for s in srps])
    meta.to_csv(META / f"{gse}_runs.tsv", sep="\t", index=False)

    meta["bases"] = pd.to_numeric(meta["run_total_bases"], errors="coerce")
    meta["spots"] = pd.to_numeric(meta["run_total_spots"], errors="coerce")
    gsm_col = find_gsm_column(meta)
    df = (
        meta.groupby(gsm_col)
        .agg(
            runs=("run_accession", lambda r: ";".join(sorted(set(r)))),
            layout=("library_layout", "first"),
            strategy=("library_strategy", "first"),
            instrument=("instrument_model", "first"),
            bases=("bases", "sum"),
            spots=("spots", "sum"),
            title=("experiment_title", "first"),
        )
        .reset_index()
        .rename(columns={gsm_col: "gsm"})
    )
    # read length per mate: total bases / spots, halved for paired-end
    mates = df["layout"].eq("PAIRED").map({True: 2, False: 1})
    df["read_length"] = (df["bases"] / df["spots"] / mates).round().astype("Int64")
    df["million_reads"] = (df["spots"] / 1e6).round(1)
    df["title"] = df["title"].map(clean_title)
    return df.drop(columns=["bases", "spots"])


def replicate(title):
    m = re.search(r"(?:rep|_n)\s*(\d+)", title)
    if not m:
        raise ValueError(f"No replicate number in title: {title}")
    return m.group(1)


def label_rnaseq(dataset, title):
    """Labels for one RNA-seq sample, from its GEO title.

    include = yes : quantified with salmon. Each DE analysis then picks its
                    samples by cell_type / age / treatment in config.yaml.
    include = no  : never downloaded.
    """
    lab = dict(genotype="WT", treatment="untreated", include="yes", note="")

    if dataset == "GSE263808":
        stage = re.search(r"\b([an])OPC\b", title).group(1)
        lab.update(cell_type="O4", age={"n": "P5", "a": "P60"}[stage])
        drug = title.split()[0]                       # Untreated / DMSO / KATi / SUV420i
        lab["treatment"] = "untreated" if drug == "Untreated" else drug
        return lab

    if dataset == "GSE66047":
        if title.startswith("P2 brain oligodendrocyte progenitor"):
            return {**lab, "cell_type": "PDGFRa", "age": "P2"}
        if title.startswith("P18 brain oligodendrocytes"):
            return {**lab, "cell_type": "OL", "age": "P18",
                    "note": "Plp1-GFP oligodendrocytes; lineage reference, not used in DE"}
        if "Dnmt1" in title:
            ko = "Olig1cre/+" in title
            return {**lab, "cell_type": "PDGFRa", "age": "P5",
                    "genotype": "Dnmt1cKO" if ko else "Dnmt1flox",
                    "include": "no", "note": "Dnmt1 sub-study; not part of this reanalysis"}

    if dataset == "GSE122446":
        if "oligodendrocyte progenitor" in title:
            return {**lab, "cell_type": "PDGFRa", "age": "P60"}
        if "pre-myelinating" in title:
            return {**lab, "cell_type": "preOL", "age": "P60",
                    "note": "Plp1-GFP pre-myelinating OLs; lineage reference, not used in DE"}

    raise ValueError(f"No labeling rule for {dataset}: '{title}'")


def build_rnaseq(df):
    labels = pd.DataFrame(
        [label_rnaseq(d, t) for d, t in zip(df["dataset"], df["title"])],
        index=df.index,
    )
    df = pd.concat([df, labels], axis=1)
    df["replicate"] = df["title"].map(replicate)
    geno = df["genotype"].map(lambda g: "" if g == "WT" else f"_{g}")
    df["sample"] = (df["cell_type"] + "_" + df["age"] + "_" + df["treatment"]
                    + geno + "_rep" + df["replicate"])
    cols = ["sample", "include", "dataset", "gsm", "runs", "layout", "read_length",
            "million_reads", "instrument", "cell_type", "age", "treatment",
            "genotype", "replicate", "title", "note"]
    return df[cols]


def build_chipseq(df):
    df = df.copy()
    df["cell_type"] = "PDGFRa"
    df["age"] = df["title"].str.extract(r"\((P\d+)\)", expand=False)
    df["mark"] = df["title"].str.contains("input", case=False).map(
        {True: "input", False: "H4K8ac"})
    df["replicate"] = df["title"].str.extract(r"rep\s*(\d+)", expand=False)
    if df[["age", "replicate"]].isna().any().any():
        raise ValueError("ChIP-seq title without age or replicate:\n"
                         + df.loc[df[["age", "replicate"]].isna().any(axis=1), "title"].to_string())
    df["include"] = "yes"
    df["sample"] = (df["cell_type"] + "_" + df["age"] + "_" + df["mark"]
                    + "_rep" + df["replicate"])
    # each ChIP sample uses the input from the same age and replicate
    inputs = df[df["mark"] == "input"].set_index(["age", "replicate"])["sample"]
    df["control"] = [
        "" if m == "input" else inputs[(a, r)]
        for a, r, m in zip(df["age"], df["replicate"], df["mark"])
    ]
    cols = ["sample", "include", "dataset", "gsm", "runs", "layout", "read_length",
            "million_reads", "instrument", "cell_type", "age", "mark",
            "replicate", "control", "title"]
    return df[cols]


def check_unique(df, name):
    dup = df["sample"][df["sample"].duplicated(keep=False)]
    if len(dup):
        raise ValueError(f"Duplicate sample names in {name}:\n{dup.to_string()}")


def main():
    META.mkdir(parents=True, exist_ok=True)
    db = SRAweb()
    tables = {"rnaseq": [], "chipseq": []}

    for gse, assay in DATASETS.items():
        print(f"Fetching {gse} ...")
        runs = fetch_runs(db, gse)
        runs = runs[runs["strategy"] == KEEP_STRATEGY[assay]].drop(columns="strategy")
        runs.insert(0, "dataset", gse)
        print(f"  {len(runs)} {KEEP_STRATEGY[assay]} samples")
        tables[assay].append(runs)

    rna = build_rnaseq(pd.concat(tables["rnaseq"], ignore_index=True))
    chip = build_chipseq(pd.concat(tables["chipseq"], ignore_index=True))
    check_unique(rna, "RNA-seq")
    check_unique(chip, "ChIP-seq")
    rna.to_csv(OUT / "samples_rnaseq.tsv", sep="\t", index=False)
    chip.to_csv(OUT / "samples_chipseq.tsv", sep="\t", index=False)

    used = rna[rna["include"] == "yes"]
    print(f"\nRNA-seq: {len(used)} of {len(rna)} samples included")
    print(used.groupby(["cell_type", "age", "treatment"]).size().to_string())
    print(f"\nChIP-seq: {len(chip)} samples")
    print(chip.groupby(["age", "mark"]).size().to_string())


if __name__ == "__main__":
    main()
