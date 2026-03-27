#!/usr/bin/env python3
"""
03_bioreason_annotate.py
------------------------
Run BioReason-Pro RL on all uncharacterized syn3A proteins using a
combined prompt integrating InterPro domains, PROST homology, genomic
neighborhood context, and prior literature annotations.

Requires:
  - NVIDIA A100 GPU (T4 will work but is very slow due to CPU offload)
  - BioReason-Pro cloned: git clone https://github.com/bowang-lab/BioReason-Pro
  - HuggingFace access to wanglab/bioreason-pro-rl (sign ESM3 license)
  - data/syn3a_unknown_all.fasta
  - data/interpro_results.json
  - data/prost_lookup.json
  - data/syn3a_all_cds_ordered.csv

Output:
  data/bioreason_combined_checkpoint.json   raw reasoning traces per protein

Usage:
    python 03_bioreason_annotate.py --repo-root /path/to/BioReason-Pro
                                    [--outdir ./data]
                                    [--model wanglab/bioreason-pro-rl]
                                    [--max-new-tokens 512]
"""

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

import pandas as pd
import torch
from Bio import SeqIO
from transformers import AutoModelForCausalLM, AutoTokenizer


ORGANISM = "Mycoplasma mycoides"
BR_MODEL  = "wanglab/bioreason-pro-rl"


def load_gogpt(repo_root: str):
    """Load GOGPTPredictor via importlib to avoid stale module cache."""
    for key in list(sys.modules.keys()):
        if "gogpt" in key or "bioreason" in key:
            del sys.modules[key]

    gogpt_init = os.path.join(repo_root, "gogpt", "src", "gogpt", "__init__.py")
    gogpt_src  = os.path.join(repo_root, "gogpt", "src", "gogpt")
    spec = importlib.util.spec_from_file_location(
        "gogpt", gogpt_init,
        submodule_search_locations=[gogpt_src]
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gogpt"] = mod
    spec.loader.exec_module(mod)
    return mod.GOGPTPredictor


def load_bioreason(model_name: str):
    print(f"Loading BioReason-Pro from {model_name} ...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model     = AutoModelForCausalLM.from_pretrained(
        model_name,
        dtype=torch.bfloat16,
        device_map="auto",
    )
    model.eval()
    device = next(model.parameters()).device
    print(f"  Loaded — device: {device}")
    return tokenizer, model


def get_neighborhood(locus_tag: str, all_cds_df: pd.DataFrame,
                     locus_to_idx: dict, n: int = 3) -> str:
    idx = locus_to_idx.get(locus_tag)
    if idx is None:
        return "Neighborhood not available"
    lines = []
    for i in range(max(0, idx - n), idx):
        row = all_cds_df.iloc[i]
        lines.append(f"  [{i-idx}] {row['locus_tag']} — {row['product']}")
    row = all_cds_df.iloc[idx]
    lines.append(f"  [0] {row['locus_tag']} — {row['product']} <- TARGET")
    for i in range(idx + 1, min(len(all_cds_df), idx + n + 1)):
        row = all_cds_df.iloc[i]
        lines.append(f"  [+{i-idx}] {row['locus_tag']} — {row['product']}")
    return "\n".join(lines)


def build_prompt(prot: dict, interpro_results: dict,
                 prost_lookup: dict, all_cds_df: pd.DataFrame,
                 locus_to_idx: dict) -> str:
    tag       = prot["locus_tag"]
    interpro  = interpro_results.get(tag, "No InterPro domains found")
    prost     = prost_lookup.get(tag, {})
    neighborhood = get_neighborhood(tag, all_cds_df, locus_to_idx)

    if prost:
        prost_str = (
            f"PROST function: {prost.get('prost_function', 'N/A')}\n"
            f"Classification: {prost.get('classification', 'N/A')}\n"
            f"Best homolog: {prost.get('best_homolog', 'N/A')} — "
            f"{prost.get('homolog_function', 'N/A')}\n"
            f"FATCAT p-value: {prost.get('fatcat_p_score', 'N/A')}  |  "
            f"Sequence identity: {prost.get('seq_identity', 'N/A')}\n"
            f"Homologs: PROST={prost.get('n_prost_homologs', 0)} / "
            f"BLAST={prost.get('n_blast_homologs', 0)} / "
            f"Foldseek={prost.get('n_foldseek_homologs', 0)}"
        )
        lit_str    = prost.get("literature", "")
        struct_str = prost.get("structural_homolog", "")
    else:
        prost_str  = "No PROST data available"
        lit_str    = ""
        struct_str = ""

    return (
        f"Protein: {tag}\n"
        f"Organism: JCVI-syn3A (synthetic minimal cell, Mycoplasma mycoides, 493 genes)\n"
        f"Length: {prot['length']} amino acids\n"
        f"GenBank annotation: {prot['description'].split('protein_id=')[0].strip()}\n\n"
        f"=== Genomic neighborhood (±3 genes) ===\n{neighborhood}\n\n"
        f"=== InterPro domain architecture ===\n{interpro}\n\n"
        f"=== PROST homology analysis ===\n{prost_str}\n\n"
        f"=== Structural homolog ===\n{struct_str}\n\n"
        f"=== Prior literature annotations ===\n{lit_str}\n\n"
        f"Every gene in syn3A is essential or quasi-essential. "
        f"Integrating all evidence, provide a concise structured annotation:\n\n"
        f"Molecular function: [one line]\n"
        f"Biological process: [one line]\n"
        f"Functional category: [one of: Membrane transport, Lipoprotein/membrane, "
        f"Proteolysis/peptidase, RNA modification, DNA metabolism, Redox/oxidoreductase, "
        f"Hydrolase/phosphatase, Kinase/signaling, Methyltransferase, Glycosyl transferase, "
        f"Acetyltransferase, Transcriptional regulator, Protein secretion, "
        f"Adaptor/scaffold, Cytoskeletal/division, Nucleic acid binding, "
        f"Membrane scaffold, ECF transporter, Unknown]\n"
        f"Confidence: [high/medium/low]\n"
        f"Rationale: [2-3 sentences integrating all evidence sources]"
    )


def annotate(prompt: str, tokenizer, model, max_new_tokens: int = 512) -> str:
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
    with torch.no_grad():
        output = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=False,
            pad_token_id=tokenizer.eos_token_id,
        )
    new_tokens = output[0][inputs["input_ids"].shape[1]:]
    return tokenizer.decode(new_tokens, skip_special_tokens=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root",      required=True,
                        help="Path to cloned BioReason-Pro repo")
    parser.add_argument("--outdir",         default="data")
    parser.add_argument("--model",          default=BR_MODEL)
    parser.add_argument("--max-new-tokens", type=int, default=512)
    args = parser.parse_args()

    outdir   = Path(args.outdir)
    repo_root = args.repo_root

    # Fix sys.path
    sys.path.insert(0, repo_root)
    sys.path.insert(0, os.path.join(repo_root, "gogpt", "src"))

    # Load data
    fasta_path = outdir / "syn3a_unknown_all.fasta"
    records    = list(SeqIO.parse(fasta_path, "fasta"))

    def get_tier(r):
        for part in r.description.split():
            if part.startswith("tier="):
                return part.split("=")[1]
        return "unknown"

    proteins = [
        {
            "locus_tag":   r.id,
            "description": r.description,
            "tier":        get_tier(r),
            "sequence":    str(r.seq),
            "length":      len(r.seq),
        }
        for r in records
    ]
    print(f"Loaded {len(proteins)} proteins from {fasta_path}")

    with open(outdir / "interpro_results.json") as f:
        interpro_results = json.load(f)
    with open(outdir / "prost_lookup.json") as f:
        prost_lookup = json.load(f)

    all_cds_df  = pd.read_csv(outdir / "syn3a_all_cds_ordered.csv")
    locus_to_idx = {row["locus_tag"]: i for i, row in all_cds_df.iterrows()}

    # Load model
    tokenizer, model = load_bioreason(args.model)

    # Checkpoint
    checkpoint_path = outdir / "bioreason_combined_checkpoint.json"
    if checkpoint_path.exists():
        with open(checkpoint_path) as f:
            results = json.load(f)
        print(f"Resuming: {len(results)} already done")
    else:
        results = {}

    # Run
    for i, prot in enumerate(proteins):
        tag = prot["locus_tag"]
        if tag in results:
            continue
        print(f"[{i+1}/{len(proteins)}] {tag} ({prot['length']} aa) ...")
        try:
            prompt      = build_prompt(prot, interpro_results, prost_lookup,
                                       all_cds_df, locus_to_idx)
            results[tag] = annotate(prompt, tokenizer, model, args.max_new_tokens)
            print(f"  -> {results[tag][:120]}")
        except Exception as e:
            results[tag] = f"ERROR: {e}"
            print(f"  -> ERROR: {e}")

        if (i + 1) % 10 == 0:
            with open(checkpoint_path, "w") as f:
                json.dump(results, f)
            print(f"  Checkpoint saved ({i+1} done)")

    with open(checkpoint_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nDone. {len(results)} proteins annotated.")
    print(f"Saved: {checkpoint_path}")


if __name__ == "__main__":
    main()
