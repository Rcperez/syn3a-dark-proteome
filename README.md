# Functional annotation of the JCVI-syn3A dark proteome

> Multimodal biological reasoning over 132 uncharacterized proteins in the world's smallest synthetic cell

[![bioRxiv](https://img.shields.io/badge/bioRxiv-preprint-red)](https://biorxiv.org)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Rcperez/syn3a-dark-proteome/blob/main/notebooks/01_fetch_and_interpro.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

JCVI-syn3A is the smallest self-replicating synthetic cell ever created — 493 genes, 543 kbp, derived from *Mycoplasma mycoides*. Every retained gene is essential or quasi-essential. Yet 132 of its protein-coding genes remain functionally uncharacterized in GenBank, representing a fundamental gap in our understanding of minimal life and blocking complete whole-cell computational modeling.

This repository contains the full computational pipeline and results for the first application of a multimodal reasoning large language model ([BioReason-Pro](https://github.com/bowang-lab/BioReason-Pro)) to the syn3A dark proteome. For each uncharacterized protein, we integrate:

- **InterPro domain architecture** (EBI InterProScan 5, batch submission)
- **PROST remote homology** (protein language model embeddings, down to 16% sequence identity)
- **Genomic neighborhood context** (±3 flanking genes from CP016816.2)
- **Five-study literature annotation history** (Danchin & Fang 2016 through Kilinc et al. 2025)

BioReason-Pro synthesizes these evidence streams into mechanistic reasoning traces, which are parsed into structured functional annotations via Claude Sonnet (Anthropic API).

## Key findings

- 132 proteins annotated (130 strictly uncharacterized, 2 putative)
- 122/132 (92%) received InterPro domain hits
- 20 functional categories identified
- Adaptor/scaffold proteins are the largest dark category (18 proteins, 14%) — not previously characterized as such in the literature
- Membrane-associated proteins dominate overall (58%)
- RNA modification machinery is overrepresented (11 proteins, 8%) for a 493-gene genome
- Kinase/signaling proteins identified in a supposedly stimulus-free minimal cell (`_0264`, `_0495`, `_0906`)
- Novel mechanistic hypotheses for the most intractable proteins: `_0353` as a DivIVA-class cytoskeletal organizer, `_0805` as a Fic/Fido-domain adenylyltransferase
- 12 proteins remain unresolvable with current evidence

## Repository structure
```
syn3a-dark-proteome/
├── data/
│   ├── syn3a_master_annotations.csv     # Full annotation table (132 proteins, 18 columns)
│   ├── syn3a_master_annotations.json    # Full annotation table (JSON)
│   ├── syn3a_unknown_all.fasta          # 132 uncharacterized protein sequences
│   ├── interpro_all.tsv                 # Raw InterProScan 5 output (merged batches)
│   ├── prost_lookup.json                # PROST homology data
│   └── syn3a_all_cds_ordered.csv        # Full genomic CDS order (458 genes)
├── pipeline/
│   ├── 01_fetch_fasta.py                # Pull unknown proteins from NCBI GenBank
│   ├── 02_fetch_prost.py                # Download and parse PROST results
│   ├── 03_bioreason_annotate.py         # BioReason-Pro inference (requires A100)
│   └── 04_extract_structured.py         # Structured extraction via Anthropic API
├── notebooks/
│   ├── 01_fetch_and_interpro.ipynb      # Colab: FASTA pull + InterPro batch submission
│   ├── 02_bioreason_inference.ipynb     # Colab: BioReason-Pro annotation (A100 GPU)
│   └── 03_analysis.ipynb               # Colab: analysis and figures
├── results/
│   └── (figures — coming soon)
├── README.md
├── requirements.txt
└── LICENSE
```

## Quickstart

### Option A: Use pre-computed results (no GPU required)
```python
import pandas as pd
df = pd.read_csv("data/syn3a_master_annotations.csv")
print(df[["locus_tag", "functional_category", "confidence", "molecular_function"]].head())
```

### Option B: Run the full pipeline

| Step | Notebook | Runtime |
|------|----------|---------|
| 1. Fetch sequences + InterPro | notebooks/01 | CPU, ~6 hours (InterPro web batch) |
| 2. BioReason-Pro annotation | notebooks/02 | A100 GPU, ~2 hours |
| 3. Analysis + figures | notebooks/03 | CPU, ~10 minutes |

### Option C: Run locally
```bash
git clone https://github.com/Rcperez/syn3a-dark-proteome
cd syn3a-dark-proteome
pip install -r requirements.txt
python pipeline/01_fetch_fasta.py --email your@email.com
python pipeline/02_fetch_prost.py
python pipeline/03_bioreason_annotate.py   # requires A100 GPU
export ANTHROPIC_API_KEY=sk-ant-...
python pipeline/04_extract_structured.py
```

## Data description

### `syn3a_master_annotations.csv`

| Column | Description |
|--------|-------------|
| `locus_tag` | JCVISYN3A\_XXXX identifier |
| `tier` | strict (zero annotation) or putative |
| `length_aa` | Protein length in amino acids |
| `genbank_annotation` | Original GenBank product field |
| `interpro_domains` | InterPro domain hits (text summary) |
| `prost_function` | PROST-assigned function |
| `prost_classification` | PROST TIGRfam classification tier |
| `prost_best_homolog` | Best structural homolog (UniProt ID) |
| `prost_homolog_function` | Homolog functional description |
| `prost_fatcat_p` | FATCAT structural alignment p-value |
| `prost_seq_identity` | Sequence identity to best homolog |
| `prost_literature` | 5-study literature annotation history |
| `bioreason_combined` | Full BioReason-Pro reasoning trace |
| `molecular_function` | Extracted molecular function |
| `biological_process` | Extracted biological process |
| `functional_category` | Assigned functional category |
| `confidence` | Annotation confidence (high/medium/low) |
| `rationale` | Evidence synthesis (2-3 sentences) |

## Methods summary

1. **Sequence retrieval**: CDS features from CP016816.2 (JCVI-syn3A) fetched via Biopython/Entrez. Features matching "hypothetical protein", "uncharacterized", "unknown function", or "DUF" retained (132 proteins, 130 strictly uncharacterized).

2. **InterPro annotation**: Sequences submitted in two batches of 66 to EBI InterProScan 5 via web interface (REST API rate-limited from Colab IPs). Results merged and parsed into per-protein domain summaries.

3. **PROST homology**: Results pulled directly from the PROST GitHub CDN (`raw.githubusercontent.com/MesihK/minweb/master/jsonwp/PROST.json.gz`) — no web server required. Provides remote homolog detection down to ~16% sequence identity with FATCAT structural alignment p-values and aggregated literature annotations from five prior studies (Danchin & Fang 2016, Yang & Tsui 2018, Antczak et al. 2019, Zhang et al. 2021, Bianchi et al. 2022).

4. **Genomic neighborhood**: Full ordered CDS list from CP016816.2 used to extract ±3 flanking genes per uncharacterized protein.

5. **BioReason-Pro inference**: Each protein annotated using BioReason-Pro RL (Fallahpour et al. 2026, bioRxiv) with a combined prompt integrating InterPro + PROST + genomic neighborhood + literature history. Inference on NVIDIA A100 GPU.

6. **Structured extraction**: BioReason-Pro reasoning traces parsed into structured fields using Claude Sonnet (Anthropic API) — necessary because BioReason-Pro's verbose output consistently exceeds token limits before reaching structured fields.

## Known limitations

- GO-GPT returned only root-level GO terms for all 132 proteins (GO:0003674, GO:0008150, GO:0005575) — no discriminatory signal. This is expected for highly divergent mycoplasma proteins underrepresented in GO-GPT's UniProt training data.
- InterPro batch submission via EBI web interface limited to 100 sequences per job — at larger scale, local InterProScan installation is required.
- BioReason-Pro confidence scores reflect reasoning trace quality, not experimental validation. All annotations are computational predictions requiring experimental confirmation.
- 12 proteins remain unresolvable (functional_category = "Unknown") with current evidence.

## Citation

If you use this dataset or pipeline, please cite:
```bibtex
@article{perez2026syn3a,
  title={Functional annotation of the JCVI-syn3A dark proteome using multimodal biological reasoning},
  author={Perez, Rolando},
  journal={bioRxiv},
  year={2026},
  note={Blue Marble Space Institute of Science}
}
```

Please also cite the tools this work depends on:

- BioReason-Pro: Fallahpour et al. 2026, bioRxiv 2026.03.19.712954
- PROST: Kilinc et al. 2025, Methods Mol Biol 2867:153-168
- InterProScan: Paysan-Lafosse et al. 2023, Nucleic Acids Res 51:D418-D427
- JCVI-syn3A: Breuer et al. 2019, eLife 8:e36842

## License

MIT. See [LICENSE](LICENSE).

## Contact

Rolando Perez  
Blue Marble Space Institute of Science  
GitHub: [@Rcperez](https://github.com/Rcperez)
