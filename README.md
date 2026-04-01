# Functional annotation of the JCVI-syn3A dark proteome

> Multimodal biological reasoning over 132 uncharacterized proteins in the world's smallest synthetic cell

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

JCVI-syn3A is the smallest self-replicating synthetic cell ever constructed — 493 genes, 543 kbp, derived from *Mycoplasma mycoides* — yet 132 of its protein-coding genes remain functionally uncharacterized in GenBank. This represents a fundamental barrier to complete whole-cell computational modeling: the most recent 4D whole-cell model of syn3A (Thornburg et al., Cell 2026) simulates all 493 gene products but must assign placeholder parameters to these uncharacterized proteins.

This repository contains the full computational pipeline and results for the first systematic application of multimodal biological reasoning (BioReason-Pro) to the syn3A dark proteome, combined with ESMFold structure prediction and Foldseek structural search for the most intractable cases.

## Results summary

| Metric | Value |
|--------|-------|
| Proteins annotated | 131/132 (99.2%) |
| Truly intractable | 1/132 (JCVISYN3A_0416) |
| Functional categories identified | 20 |
| Largest dark category | Adaptor/scaffold (18 proteins, 14%) |
| Membrane-associated overall | ~58% |
| Resolved by structural search (Tier 1) | 7/12 previously unknown |
| Confidence: high / medium / low | 9 / 50 / 73 |

## Key findings

- **Adaptor/scaffold proteins are the largest dark category** (18 proteins, 14%) — not previously characterized as such in the literature. These non-catalytic structural organizers are invisible to homology-based methods.
- **58% of the dark proteome is membrane-associated**, suggesting the missing functional knowledge concerns membrane organization rather than core metabolism.
- **RNA modification is overrepresented** (11 proteins, 8%) for a 493-gene genome, consistent with the importance of translational fine-tuning in minimal cells.
- **Kinase/signaling proteins** identified in a cell previously considered stimulus-free: `_0264` (serine/threonine kinase, FATCAT p=1.3e-8), `_0805` (Fic/Fido adenylyltransferase), `_0906` (phosphotransferase).
- **`_0353`** — DivIVA-class cytoskeletal organizer (FATCAT p=8.73e-14, RMSD 1.25 Å).
- **`_0433`** — CutC-like TIM-barrel hydrolase with metal cofactor binding (Foldseek alnTM=0.965 against *Pasteurella multocida*).
- **`_0138`** — Strong structural hit (alnTM=0.812) against *Streptococcus mutans* at only 10.6% sequence identity — completely invisible to BLAST.
- **`_0873`** — 66 aa, DUF951, structural hit (alnTM=0.830) against *Bacillus subtilis* ribosomal scaffold.
- **`_0416`** — The single remaining intractable protein: putative toxin or endoribonuclease with no sequence or structural homologs.

## Pipeline

The annotation pipeline runs in two tiers:

### Tier 1 — Sequence-based annotation (132 proteins)
```
NCBI GenBank CP016816.2
        ↓
Sequence retrieval (Biopython/Entrez)
        ↓
InterPro domain annotation (EBI InterProScan 5)
        ↓
PROST remote homology (GitHub CDN pull)
        ↓
Genomic neighborhood context (±3 flanking genes)
        ↓
BioReason-Pro RL inference (A100 GPU, combined prompt)
        ↓
Claude Sonnet structured extraction (Anthropic API)
        ↓
120/132 proteins annotated
```

### Tier 2 — Structural follow-up (12 remaining unknowns)
```
ESMFold structure prediction (facebook/esmfold_v1, A100)
        ↓
Foldseek structural search (PDB + AlphaFold Swiss-Prot)
        ↓
BioReason-Pro rerun with structural context
        ↓
Claude Sonnet structured extraction
        ↓
7/12 resolved → 131/132 total annotated
```

## Repository structure
```
syn3a-dark-proteome/
├── data/
│   ├── syn3a_master_annotations.csv        # Tier 1 annotations (132 proteins, 18 columns)
│   ├── syn3a_master_annotations_v2.csv     # Tier 2 annotations (updated with structural hits)
│   ├── syn3a_master_annotations.json       # Same, JSON format
│   ├── syn3a_unknown_all.fasta             # 132 uncharacterized protein sequences
│   ├── interpro_all.tsv                    # Raw InterProScan 5 output
│   ├── prost_lookup.json                   # PROST homology data (all 132 proteins)
│   ├── syn3a_all_cds_ordered.csv           # Full genome CDS in genomic order
│   ├── foldseek_results.json               # Foldseek top hits per protein (Tier 2)
│   ├── enhanced_parsed_results.json        # Tier 2 structured extraction results
│   ├── bioreason_enhanced_checkpoint.json  # Tier 2 BioReason-Pro reasoning traces
│   ├── esmfold_pdbs/                       # 12 ESMFold-predicted PDB structures
│   │   ├── JCVISYN3A_0138.pdb
│   │   ├── JCVISYN3A_0248.pdb
│   │   └── ... (12 total)
│   └── foldseek_results/                   # Per-protein Foldseek TSV files
│       ├── JCVISYN3A_0138_pdb.m8
│       ├── JCVISYN3A_0138_afsp.m8
│       └── ... (24 total, 12 proteins × 2 databases)
├── pipeline/
│   ├── 01_fetch_fasta.py                   # Pull sequences from NCBI GenBank
│   ├── 02_fetch_prost.py                   # Download PROST results from GitHub CDN
│   ├── 03_bioreason_annotate.py            # BioReason-Pro inference (requires A100)
│   └── 04_extract_structured.py            # Structured extraction via Anthropic API
├── results/
│   ├── fig1_functional_landscape.pdf/png   # Functional category distribution
│   ├── fig2_confidence.pdf/png             # Annotation confidence distribution
│   └── fig3_unknowns.pdf/png              # Tier 2 unknown proteins table
├── README.md
├── requirements.txt
└── LICENSE
```

## Data description

### `syn3a_master_annotations_v2.csv` (primary output)

| Column | Description |
|--------|-------------|
| `locus_tag` | JCVISYN3A_XXXX identifier |
| `tier` | strict or putative |
| `length_aa` | Protein length in amino acids |
| `genbank_annotation` | Original GenBank product field |
| `interpro_domains` | InterPro domain hits |
| `prost_function` | PROST-assigned function |
| `prost_best_homolog` | Best structural homolog (UniProt ID) |
| `prost_homolog_function` | Homolog functional description |
| `prost_fatcat_p` | FATCAT structural alignment p-value |
| `prost_seq_identity` | Sequence identity to best homolog |
| `prost_literature` | 5-study literature annotation history |
| `molecular_function` | Extracted molecular function |
| `biological_process` | Extracted biological process |
| `functional_category` | Assigned functional category (20 categories) |
| `confidence` | Annotation confidence (high/medium/low) |
| `rationale` | Evidence synthesis (2-3 sentences) |
| `esmfold_plddt` | ESMFold mean pLDDT score (Tier 2 only) |
| `foldseek_best_hit` | Best Foldseek structural hit |
| `foldseek_alntmscore` | Alignment TM-score |
| `foldseek_lddt` | LDDT score |
| `foldseek_resolved` | Whether structural search changed the annotation |
| `structural_insight` | What ESMFold/Foldseek adds |

## Methods summary

### Tier 1
1. **Sequence retrieval** — CDS features from CP016816.2 fetched via Biopython/Entrez. 132 proteins with unknown function retained.
2. **InterPro annotation** — EBI InterProScan 5 batch submission. 122/132 proteins received hits.
3. **PROST homology** — Results pulled directly from GitHub CDN (`raw.githubusercontent.com/MesihK/minweb/master/jsonwp/PROST.json.gz`). Provides remote homolog detection down to ~16% sequence identity.
4. **Genomic neighborhood** — ±3 flanking genes extracted from full ordered CDS list.
5. **BioReason-Pro inference** — Combined prompt integrating all four evidence streams. Inference on NVIDIA A100-SXM4-80GB in bfloat16.
6. **Structured extraction** — Claude Sonnet (Anthropic API) parses BioReason-Pro reasoning traces into structured JSON fields.

### Tier 2 (12 remaining unknowns)
7. **ESMFold structure prediction** — `facebook/esmfold_v1` via HuggingFace transformers. All 12 structures predicted and cached.
8. **Foldseek structural search** — `easy-search` against PDB and AlphaFold Swiss-Prot databases. 7/12 proteins received significant hits (alnTM ≥ 0.4).
9. **Enhanced BioReason-Pro** — Rerun with structural context added to prompt. 7/12 previously unknown proteins reclassified.

### Key negative results
- **GO-GPT** returned only root-level GO terms for all 132 proteins — no discriminatory signal. Expected for highly divergent mycoplasma proteins underrepresented in UniProt training data.
- **`pip install -e .`** for BioReason-Pro fails on Colab due to `flash-attn`, `deepspeed`, `trl[vllm]` build dependencies. Use `sys.path` insertion instead (see pipeline scripts).
- **EBI InterPro REST API** is rate-limited from shared cloud IPs. Use web batch interface or local InterProScan for scale.

## Quickstart

### Use pre-computed results (no GPU required)
```python
import pandas as pd
df = pd.read_csv("data/syn3a_master_annotations_v2.csv")
print(df[["locus_tag", "functional_category", "confidence",
          "molecular_function"]].head(20))
print(df["functional_category"].value_counts())
```

### Run the full pipeline (A100 GPU required for steps 3-4)
```bash
git clone https://github.com/Rcperez/syn3a-dark-proteome
cd syn3a-dark-proteome
pip install -r requirements.txt
python pipeline/01_fetch_fasta.py --email your@email.com
python pipeline/02_fetch_prost.py
python pipeline/03_bioreason_annotate.py --repo-root /path/to/BioReason-Pro
export ANTHROPIC_API_KEY=sk-ant-...
python pipeline/04_extract_structured.py
```

## Motivation and context

The Thornburg et al. (2026) 4D whole-cell model of syn3A (Cell 189, 1-16) simulates the complete cell cycle including all genetic information processes, metabolism, growth, and division. Despite this achievement, the model must assign placeholder parameters to uncharacterized gene products. That paper notes that "the majority of [genes that go untranscribed in some simulated cells] are genes of unknown function" — highlighting that functional annotation gaps directly limit model completeness and predictive power. Our annotations of 131/132 previously uncharacterized proteins provide the functional context needed to incorporate these proteins into the next generation of syn3A whole-cell models.

## Citation

If you use this dataset or pipeline, please cite:
```bibtex
@article{perez2026syn3a,
  title={Functional annotation of the JCVI-syn3A dark proteome
         using multimodal biological reasoning},
  author={Perez, Rolando},
  journal={bioRxiv},
  year={2026},
  note={Blue Marble Space Institute of Science}
}
```

Please also cite the tools this work depends on:

- **BioReason-Pro**: Fallahpour et al. 2026, bioRxiv 2026.03.19.712954
- **PROST**: Kilinc et al. 2025, Methods Mol Biol 2867:153-168
- **InterProScan**: Paysan-Lafosse et al. 2023, Nucleic Acids Res 51:D418-D427
- **ESMFold**: Lin et al. 2023, Science 379:1123-1130
- **Foldseek**: van Kempen et al. 2024, Nature Biotechnology 42:243-246
- **JCVI-syn3A**: Breuer et al. 2019, eLife 8:e36842
- **4D whole-cell model**: Thornburg et al. 2026, Cell 189:1-16

## License

MIT. See [LICENSE](LICENSE).

## Contact

Rolando Perez
Blue Marble Space Institute of Science
GitHub: [@Rcperez](https://github.com/Rcperez)
