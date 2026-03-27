#!/usr/bin/env python3
"""
02_fetch_prost.py
-----------------
Download and parse PROST homology results for JCVI-syn3A directly
from the PROST GitHub CDN. No web server or API key required.

Output:
  data/prost_lookup.json   per-protein dict keyed by locus tag

Usage:
    python 02_fetch_prost.py [--outdir ./data]
"""

import argparse
import gzip
import json
import re
from pathlib import Path

import requests

PROST_URL = "https://raw.githubusercontent.com/MesihK/minweb/master/jsonwp/PROST.json.gz"


def fetch_prost(url: str) -> dict:
    print(f"Fetching PROST data from GitHub CDN ...")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    data = json.loads(gzip.decompress(resp.content))
    print(f"  Downloaded: {len(resp.content)/1e6:.1f} MB")
    print(f"  Keys: {list(data.keys())[:5]}")
    return data


def clean_md(text: str) -> str:
    text = re.sub(r'\[.*?\]\(.*?\)', '', str(text))
    text = re.sub(r'[#*_`]', '', text)
    return text.strip()


def build_lookup(data: dict) -> dict:
    table_rows = data['table:results']['rows']
    row_by_jcvi = {}
    for row in table_rows:
        jcvi_num = str(row[1]).zfill(4)
        row_by_jcvi[jcvi_num] = row

    def extract_literature(page_data):
        for key, val in page_data.items():
            if "s3:" in key and "Literature" in str(val):
                return clean_md(val)
        return ""

    def extract_structural(page_data):
        for key, val in page_data.items():
            if key == "md:a1":
                return clean_md(val)
        return ""

    lookup = {}
    for key in data.keys():
        if not key.startswith("page:"):
            continue
        jcvi_num  = key.replace("page:", "")
        locus_tag = f"JCVISYN3A_{jcvi_num}"
        row       = row_by_jcvi.get(jcvi_num)
        if row is None:
            continue
        page_data = data[key]
        lookup[locus_tag] = {
            "ncbi_id":             str(row[0]).split("@")[-1],
            "prost_function":      row[2],
            "classification":      row[3],
            "best_homolog":        str(row[4]).split("@")[-1],
            "homolog_function":    row[5],
            "fatcat_p_score":      row[6],
            "seq_identity":        row[7],
            "homolog_source":      row[8],
            "n_prost_homologs":    row[9],
            "n_blast_homologs":    row[10],
            "n_foldseek_homologs": row[11],
            "literature":          extract_literature(page_data),
            "structural_homolog":  extract_structural(page_data),
        }

    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="data", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data   = fetch_prost(PROST_URL)
    lookup = build_lookup(data)

    out_path = outdir / "prost_lookup.json"
    with open(out_path, "w") as f:
        json.dump(lookup, f, indent=2)

    print(f"\nPROST lookup built: {len(lookup)} proteins")
    print(f"Saved: {out_path}")

    # Quick sample
    sample_tag = "JCVISYN3A_0005"
    if sample_tag in lookup:
        p = lookup[sample_tag]
        print(f"\nSample — {sample_tag}:")
        print(f"  Function : {p['prost_function']}")
        print(f"  Homolog  : {p['best_homolog']} — {p['homolog_function']}")
        print(f"  FATCAT p : {p['fatcat_p_score']}")


if __name__ == "__main__":
    main()
