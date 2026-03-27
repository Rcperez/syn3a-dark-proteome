#!/usr/bin/env python3
"""
04_extract_structured.py
------------------------
Use Claude Sonnet (Anthropic API) to extract structured functional
annotations from BioReason-Pro reasoning traces.

BioReason-Pro generates verbose reasoning that consistently exceeds
token limits before reaching structured output fields. This script
reads the raw reasoning traces and extracts clean JSON annotations.

Requires:
  - ANTHROPIC_API_KEY environment variable
  - data/bioreason_combined_checkpoint.json
  - data/prost_lookup.json
  - data/syn3a_unknown_genes.csv

Output:
  data/syn3a_master_annotations.csv    full annotation table
  data/syn3a_master_annotations.json   same, JSON format
  data/parsed_results.json             raw extracted JSON per protein

Usage:
    export ANTHROPIC_API_KEY=sk-ant-...
    python 04_extract_structured.py [--outdir ./data]
                                    [--model claude-sonnet-4-6]
"""

import argparse
import json
import os
import time
from pathlib import Path

import anthropic
import pandas as pd

FUNCTIONAL_CATEGORIES = [
    "Membrane transport", "Lipoprotein/membrane", "Proteolysis/peptidase",
    "RNA modification", "DNA metabolism", "Redox/oxidoreductase",
    "Hydrolase/phosphatase", "Kinase/signaling", "Methyltransferase",
    "Glycosyl transferase", "Acetyltransferase", "Transcriptional regulator",
    "Protein secretion", "Adaptor/scaffold", "Cytoskeletal/division",
    "Nucleic acid binding", "Membrane scaffold", "ECF transporter", "Unknown",
]


def build_extraction_prompt(locus_tag: str, length_aa: int,
                             interpro: str, prost_function: str,
                             prost_homolog: str, prost_fatcat: str,
                             bioreason_trace: str) -> str:
    categories = ", ".join(FUNCTIONAL_CATEGORIES)
    return f"""Extract a structured functional annotation from this BioReason-Pro reasoning trace for a JCVI-syn3A minimal cell protein.

Protein: {locus_tag} ({length_aa} aa)
InterPro summary: {interpro[:300]}
PROST best homolog: {prost_homolog} (FATCAT p={prost_fatcat})
PROST function: {prost_function}

BioReason-Pro reasoning trace:
{bioreason_trace[:3000]}

Return ONLY a valid JSON object with exactly these keys:
{{
  "molecular_function": "one concise line describing the molecular activity",
  "biological_process": "one concise line describing the cellular process",
  "functional_category": "exactly one of: {categories}",
  "confidence": "exactly one of: high, medium, low",
  "rationale": "2-3 sentences synthesizing the key evidence from InterPro, PROST homology, genomic neighborhood, and literature"
}}

Return only the JSON object. No preamble, no markdown fences, no other text."""


def extract_one(client: anthropic.Anthropic, prompt: str,
                model: str, retries: int = 3) -> dict:
    for attempt in range(retries):
        try:
            response = client.messages.create(
                model=model,
                max_tokens=512,
                messages=[{"role": "user", "content": prompt}]
            )
            text = response.content[0].text.strip()
            text = text.replace("```json", "").replace("```", "").strip()
            return json.loads(text)
        except json.JSONDecodeError as e:
            if attempt < retries - 1:
                time.sleep(2)
                continue
            raise e
    return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="data")
    parser.add_argument("--model",  default="claude-sonnet-4-6")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    client = anthropic.Anthropic()

    # Load inputs
    with open(outdir / "bioreason_combined_checkpoint.json") as f:
        br_traces = json.load(f)
    with open(outdir / "prost_lookup.json") as f:
        prost_lookup = json.load(f)

    df_meta = pd.read_csv(outdir / "syn3a_unknown_genes.csv")
    meta    = {row["locus_tag"]: row for _, row in df_meta.iterrows()}

    # Interpro results (optional — gracefully handle missing)
    interpro_path = outdir / "interpro_results.json"
    if interpro_path.exists():
        with open(interpro_path) as f:
            interpro_results = json.load(f)
    else:
        interpro_results = {}
        print("Warning: interpro_results.json not found — using empty dict")

    # Checkpoint
    parsed_path = outdir / "parsed_results.json"
    if parsed_path.exists():
        with open(parsed_path) as f:
            parsed = json.load(f)
        print(f"Resuming: {len(parsed)} already extracted")
    else:
        parsed = {}

    tags = list(br_traces.keys())
    errors = []

    for i, tag in enumerate(tags):
        if tag in parsed and "error" not in parsed[tag]:
            continue

        row    = meta.get(tag, {})
        prost  = prost_lookup.get(tag, {})
        length = int(row.get("protein_len_aa", 0)) if row else 0

        prompt = build_extraction_prompt(
            locus_tag      = tag,
            length_aa      = length,
            interpro       = interpro_results.get(tag, "No InterPro domains found"),
            prost_function = prost.get("prost_function", "N/A"),
            prost_homolog  = prost.get("homolog_function", "N/A"),
            prost_fatcat   = prost.get("fatcat_p_score", "N/A"),
            bioreason_trace= br_traces[tag],
        )

        print(f"[{i+1}/{len(tags)}] {tag} ...")
        try:
            result = extract_one(client, prompt, args.model)
            parsed[tag] = result
            print(f"  -> [{result.get('confidence', '?')}] "
                  f"{result.get('functional_category', '?')} | "
                  f"{result.get('molecular_function', '')[:80]}")
        except Exception as e:
            parsed[tag] = {"error": str(e)}
            errors.append(tag)
            print(f"  -> ERROR: {e}")

        if (i + 1) % 20 == 0:
            with open(parsed_path, "w") as f:
                json.dump(parsed, f)
            print(f"  Checkpoint saved ({i+1} done)")

    with open(parsed_path, "w") as f:
        json.dump(parsed, f, indent=2)

    if errors:
        print(f"\nErrors on {len(errors)} proteins: {errors}")

    # Assemble master CSV
    rows = []
    for tag, row in meta.items():
        prost  = prost_lookup.get(tag, {})
        result = parsed.get(tag, {})
        rows.append({
            "locus_tag":             tag,
            "tier":                  "strict" if row.get("is_strict") else "putative",
            "length_aa":             row.get("protein_len_aa", ""),
            "genbank_annotation":    row.get("product", ""),
            "interpro_domains":      interpro_results.get(tag, ""),
            "prost_function":        prost.get("prost_function", ""),
            "prost_classification":  prost.get("classification", ""),
            "prost_best_homolog":    prost.get("best_homolog", ""),
            "prost_homolog_function":prost.get("homolog_function", ""),
            "prost_fatcat_p":        prost.get("fatcat_p_score", ""),
            "prost_seq_identity":    prost.get("seq_identity", ""),
            "prost_literature":      prost.get("literature", ""),
            "bioreason_combined":    br_traces.get(tag, ""),
            "molecular_function":    result.get("molecular_function", ""),
            "biological_process":    result.get("biological_process", ""),
            "functional_category":   result.get("functional_category", ""),
            "confidence":            result.get("confidence", ""),
            "rationale":             result.get("rationale", ""),
        })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(outdir / "syn3a_master_annotations.csv", index=False)

    with open(outdir / "syn3a_master_annotations.json", "w") as f:
        json.dump(
            {r["locus_tag"]: {k: v for k, v in r.items() if k != "locus_tag"}
             for r in rows},
            f, indent=2
        )

    print(f"\nDone. {len(df_out)} proteins in master table.")
    print(f"Saved: {outdir}/syn3a_master_annotations.csv")
    print(f"Saved: {outdir}/syn3a_master_annotations.json")
    print("\n=== Functional category distribution ===")
    print(df_out["functional_category"].value_counts().to_string())
    print("\n=== Confidence distribution ===")
    print(df_out["confidence"].value_counts().to_string())


if __name__ == "__main__":
    main()
