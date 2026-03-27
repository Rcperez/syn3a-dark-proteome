#!/usr/bin/env python3
"""
01_fetch_fasta.py
-----------------
Pull all uncharacterized / hypothetical protein sequences from JCVI-syn3A
(GenBank accession CP016816.2) and export:
  - syn3a_unknown_genes.csv     metadata table with is_strict flag
  - syn3a_unknown_all.fasta     all 132 sequences (strict + putative)
  - syn3a_unknown_strict.fasta  strictly uncharacterized only
  - syn3a_unknown_putative.fasta putative only
  - syn3a_all_cds_ordered.csv   full genome CDS in genomic order
                                (used for neighborhood context in later steps)

Usage:
    python 01_fetch_fasta.py --email your@email.com [--outdir ./data]
"""

import argparse
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO

ACCESSION = "CP016816.2"
STRAIN    = "JCVI-syn3A"

UNKNOWN_KEYWORDS = [
    "hypothetical protein",
    "unknown function",
    "uncharacterized",
    "domain of unknown function",
    "duf",
    "putative",
]

STRICT_KEYWORDS = [
    "hypothetical protein",
    "unknown function",
    "uncharacterized",
    "domain of unknown function",
    "duf",
]


def matches(qualifiers: dict, keywords: list) -> bool:
    product = qualifiers.get("product", [""])[0].lower()
    note    = qualifiers.get("note",    [""])[0].lower()
    return any(kw in (product + " " + note) for kw in keywords)


def fetch_record(accession: str) -> SeqIO.SeqRecord:
    print(f"Fetching {accession} from NCBI ...")
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gbwithparts",
        retmode="text",
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()
    time.sleep(0.5)
    print(f"  Loaded: {record.id} ({len(record.features)} features)")
    return record


def parse_unknown(record: SeqIO.SeqRecord) -> pd.DataFrame:
    rows = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        if not matches(feat.qualifiers, UNKNOWN_KEYWORDS):
            continue
        q           = feat.qualifiers
        protein_seq = q.get("translation", ["N/A"])[0]
        rows.append({
            "locus_tag":      q.get("locus_tag",  [""])[0],
            "gene_name":      q.get("gene",        [""])[0],
            "product":        q.get("product",     [""])[0],
            "note":           q.get("note",        [""])[0],
            "protein_id":     q.get("protein_id",  [""])[0],
            "start_bp":       int(feat.location.start),
            "end_bp":         int(feat.location.end),
            "strand":         feat.location.strand,
            "protein_seq":    protein_seq,
            "protein_len_aa": len(protein_seq) if protein_seq != "N/A" else None,
            "is_strict":      matches(q, STRICT_KEYWORDS),
        })
    return pd.DataFrame(rows)


def parse_all_cds(record: SeqIO.SeqRecord) -> pd.DataFrame:
    rows = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        q = feat.qualifiers
        rows.append({
            "locus_tag": q.get("locus_tag", [""])[0],
            "gene_name": q.get("gene",       [""])[0],
            "product":   q.get("product",    [""])[0],
            "start_bp":  int(feat.location.start),
            "end_bp":    int(feat.location.end),
            "strand":    feat.location.strand,
        })
    return pd.DataFrame(rows).sort_values("start_bp").reset_index(drop=True)


def write_fasta(df: pd.DataFrame, path: Path, strain: str) -> None:
    written = 0
    with open(path, "w") as fh:
        for _, row in df.iterrows():
            if row["protein_seq"] == "N/A":
                continue
            tier   = "strict" if row["is_strict"] else "putative"
            header = (
                f">{row['locus_tag']} [{strain}] "
                f"tier={tier} protein_id={row['protein_id']} "
                f"{row['product']}"
            )
            seq     = row["protein_seq"]
            wrapped = "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))
            fh.write(f"{header}\n{wrapped}\n")
            written += 1
    print(f"  Saved: {path} ({written} sequences)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--email",  required=True, help="Email for NCBI Entrez")
    parser.add_argument("--outdir", default="data", help="Output directory")
    args = parser.parse_args()

    Entrez.email = args.email
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    record     = fetch_record(ACCESSION)
    df         = parse_unknown(record)
    df_strict  = df[df["is_strict"]].copy()
    df_putative = df[~df["is_strict"]].copy()
    df_all_cds = parse_all_cds(record)

    print(f"\nUnknown-function CDS: {len(df)} total")
    print(f"  Strictly uncharacterized: {len(df_strict)}")
    print(f"  Putative only:            {len(df_putative)}")
    print(f"\nFull CDS (genomic order): {len(df_all_cds)}")

    # CSV exports
    df.to_csv(outdir / "syn3a_unknown_genes.csv", index=False)
    print(f"\n  Saved: {outdir}/syn3a_unknown_genes.csv")

    df_all_cds.to_csv(outdir / "syn3a_all_cds_ordered.csv", index=False)
    print(f"  Saved: {outdir}/syn3a_all_cds_ordered.csv")

    # FASTA exports
    write_fasta(df,           outdir / "syn3a_unknown_all.fasta",      STRAIN)
    write_fasta(df_strict,    outdir / "syn3a_unknown_strict.fasta",    STRAIN)
    write_fasta(df_putative,  outdir / "syn3a_unknown_putative.fasta",  STRAIN)

    print("\nDone.")


if __name__ == "__main__":
    main()
