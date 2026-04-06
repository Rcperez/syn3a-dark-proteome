# Functional annotation of the JCVI-syn3A dark proteome

> Multimodal biological reasoning over 132 uncharacterized proteins in the world's smallest synthetic cell

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

JCVI-syn3A is the smallest self-replicating synthetic cell ever constructed -- 493 genes, 543 kbp, derived from Mycoplasma mycoides -- yet 132 of its protein-coding genes remain functionally uncharacterized in GenBank. This represents a fundamental barrier to complete whole-cell computational modeling: the most recent 4D whole-cell model of syn3A (Thornburg et al., Cell 2026) simulates all 493 gene products but must assign placeholder parameters to these uncharacterized proteins.

This repository contains the full computational pipeline and results for the first systematic application of multimodal biological reasoning (BioReason-Pro RL) to the syn3A dark proteome, combined with ESMFold structure prediction and Foldseek structural search for the most intractable cases.

## Results summary

| Metric | Value |
|--------|-------|
| Proteins processed | 132/132 (100%) |
| Proteins with category assignments | 131/132 |
| Sole computationally intractable protein | JCVISYN3A_0416 |
| Functional categories identified | 20 |
| Largest dark category | Adaptor/scaffold (24 proteins, 18%) |
| Membrane-associated overall | ~58% |
| Resolved by structural search (alnTM >= 0.4) | 7/12 |
| Not resolved by structural search (foldseek_resolved=False) | 5/12 |
| Confidence: high / medium / low | 9 / 50 / 73 |

## Key findings

- **Adaptor/scaffold proteins are the largest dark category** (24 proteins, 18%) -- non-catalytic structural organizers invisible to homology-based methods.
- **58% of the dark proteome is membrane-associated**, suggesting the missing functional knowledge concerns membrane organization rather than core metabolism.
- **RNA modification is overrepresented** (12 proteins, 9%), consistent with the importance of translational fine-tuning in minimal cells.
- **Kinase/signaling proteins** identified in a cell previously considered stimulus-free: JCVISYN3A_0264 (serine/threonine kinase, FATCAT p=1.3e-8), JCVISYN3A_0805 (Fic/Fido adenylyltransferase).
- **JCVISYN3A_0353** -- DivIVA-class cytoskeletal organizer (FATCAT p=8.73e-14).
- **JCVISYN3A_0433** -- CutC-like TIM-barrel hydrolase (Foldseek alnTM=0.965 against Pasteurella multocida).
- **JCVISYN3A_0138** -- Strong structural hit (alnTM=0.812) against Streptococcus mutans at only 10.6% sequence identity -- invisible to BLAST.
- **5 proteins not resolved by structural search** (foldseek_resolved=False, alnTM < 0.4): 4 carry low-confidence BioReason-Pro hypotheses; JCVISYN3A_0416 alone has no category assignment.

## Pipeline

### Tier 1 -- Sequence-based annotation (132 proteins)

```
NCBI GenBank CP016816.2
        |
Sequence retrieval (Biopython/Entrez)
        |
InterPro domain annotation (EBI InterProScan 5)
        |
PROST remote homology (GitHub CDN pull)
        |
Genomic neighborhood context (+-3 flanking genes)
        |
BioReason-Pro RL inference (A100 GPU, combined prompt)
        |
Claude Sonnet structured extraction (Anthropic API)
        |
120/132 proteins annotated
```

### Tier 2 -- Structural follow-up (12 remaining unknowns)

```
ESMFold structure prediction (facebook/esmfold_v1, A100)
        |
Foldseek structural search (PDB + AlphaFold Swiss-Prot)
        |
BioReason-Pro rerun with structural context
        |
Claude Sonnet structured extraction
        |
7/12 resolved (foldseek_resolved=True, alnTM >= 0.4)
5/12 not resolved (foldseek_resolved=False, alnTM < 0.4)
        |
132/132 proteins processed; 131 with category assignments
```

## Repository structure

```
syn3a-dark-proteome/
├── data/
│   ├── syn3a_master_annotations_v2.csv     # Final annotations (132 proteins, 31 columns)
│   ├── syn3a_master_annotations.csv        # Tier 1 annotations (18 columns)
│   ├── syn3a_master_annotations.json
│   ├── syn3a_unknown_all.fasta
│   ├── interpro_all.tsv
│   ├── prost_lookup.json
│   ├── syn3a_all_cds_ordered.csv
│   ├── foldseek_results.json
│   ├── enhanced_parsed_results.json
│   ├── bioreason_enhanced_checkpoint.json
│   ├── esmfold_pdbs/                       # 12 ESMFold-predicted PDB structures
│   └── foldseek_results/                   # 24 Foldseek TSV files (12 x 2 databases)
├── pipeline/
│   ├── 01_fetch_fasta.py
│   ├── 02_fetch_prost.py
│   ├── 03_bioreason_annotate.py
│   └── 04_extract_structured.py
├── notebooks/
│   ├── Syn3A_Tier1_annotation_v2.ipynb     # Tier 1 pipeline + all figures
│   └── Syn3A_Tier2_structural_annotation.ipynb
├── results/
│   ├── fig1_functional_landscape.pdf/png
│   ├── fig2_confidence.pdf/png
│   ├── fig3_unknowns.pdf/png               # 5 unresolved proteins (foldseek_resolved=False)
│   ├── fig3_unresolved_proteins.csv        # Same as CSV with full rationale column
│   └── fig4_tier2_foldseek.pdf/png        # alnTM/LDDT/pident for 7 resolved proteins
├── README.md
├── requirements.txt
└── LICENSE
```

## Key columns in syn3a_master_annotations_v2.csv

| Column | Description |
|--------|-------------|
| locus_tag | JCVISYN3A_XXXX identifier |
| tier | strict or putative |
| length_aa | Protein length in amino acids |
| functional_category | Assigned functional category (20 categories) |
| confidence | high / medium / low |
| molecular_function | Extracted molecular function |
| biological_process | Extracted biological process |
| rationale | Evidence synthesis (2-3 sentences) |
| esmfold_plddt | ESMFold mean pLDDT (Tier 2 proteins only) |
| foldseek_alntmscore | Alignment TM-score (primary structural confidence metric) |
| foldseek_resolved | True if alnTM >= 0.4; False if no significant structural hit |
| structural_insight | What ESMFold/Foldseek adds to the annotation |

## Motivation and context

The Thornburg et al. (2026) 4D whole-cell model of syn3A (Cell 189, 1-16) simulates the complete cell cycle including all genetic information processes, metabolism, growth, and division. Despite this achievement, the model must assign placeholder parameters to uncharacterized gene products, and the majority of genes that go untranscribed in at least one simulated cell cycle are genes of unknown function. Our annotations -- with explicit foldseek_resolved flags distinguishing structurally confirmed annotations from low-confidence hypotheses -- provide the functional context needed to incorporate the syn3A dark proteome into the next generation of whole-cell models.

## Quickstart

```python
import pandas as pd
df = pd.read_csv("data/syn3a_master_annotations_v2.csv")

# Final category distribution
print(df["functional_category"].value_counts())

# Structurally confirmed (Tier 2, alnTM >= 0.4)
resolved = df[df["foldseek_resolved"] == True]

# Not resolved by structural search
unresolved = df[df["foldseek_resolved"] == False]

# Sole truly intractable protein
intractable = df[df["functional_category"] == "Unknown"]
```

## Citation

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

Please also cite: BioReason-Pro (Fallahpour et al. 2026), PROST (Kilinc et al. 2025), InterProScan (Paysan-Lafosse et al. 2023), ESMFold (Lin et al. 2023), Foldseek (van Kempen et al. 2024), JCVI-syn3A (Breuer et al. 2019), 4D whole-cell model (Thornburg et al. 2026, Cell 189:1-16).

## License

MIT. See [LICENSE](LICENSE).

## Contact

Rolando Perez | Blue Marble Space Institute of Science | [@Rcperez](https://github.com/Rcperez)
