# Safe Organism Substitute Pipeline for iGEM

In iGEM, teams are restricted to **Risk Group 1/2 organisms**. However, many interesting proteins are only characterized in unsafe species.

This repository contains a three-step computational pipeline designed to help iGEM teams **find safe substitute organisms** when their protein of interest is normally found in organisms not permitted by iGEM (Risk Group 2+).  

By starting with a **whitelist of safe organisms** and a **core protein of interest**, the pipeline identifies 
- homologous proteins
- checks for functional motifs
- evaluates structural similarity of active pockets.

---
## Environment Requirements

- Python >= 3.8  
- Dependencies:
  - requests  
  - biopython  
  - numpy  
  - scipy  
  - tqdm  
- External software:
  - [mmseqs2](https://github.com/soedinglab/MMseqs2) (for fast homology search)  
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/) (for multiple sequence alignment)  

---

## Input Preparation
All input files should be stored in the same directory  `data/`:
- **whitelist.csv**  
  A list of safe species (RG1). Must contain:  
  - `Taxon`: NCBI TaxID (required)  
  - `species`: Species name (preferred)  
- **core.fasta**  
  The reference protein sequence you care about (e.g., Chitin synthase from a non-whitelist fungus).  
- **motif.txt**  
  A list of required motifs (regular expressions), one per line. Example: 
```
QRRRW
DxE
```
- **pocket.txt**  
Residue indices (1-based, based on `core.fasta`) that define the functional pocket. Example:  
```
123
124
130
```
---
## Running the Pipeline

### Step 1. Download homologous proteins from UniProt
Fetch sequences from UniProt for all whitelist species that match the search keyword, and merge with your core protein.

```bash
python Step1.py --data_dir data --searchword "Chitin synthase"
```
**Output (data/)**:
- meta.csv : metadata for downloaded proteins
- all.fasta : combined whitelist homologs + core protein

### Step 2. Homology and motif screening
Run mmseqs2 (core vs all) and MAFFT to align sequences, check for motif conservation, and filter out weak hits.
```bash
python Step2.py --data data --outdir result
```
**Output (result/)**:
- msa.fa : multiple sequence alignment
- core_vs_all.csv : list of homologs with alignment stats and motif QC
**Key fields**:
- hit_core(0/1): significant hit to core protein
- motifs_ok(0/1): motifs conserved
- fail_tag: mark if the sequence fails homology/motif QC

### Step 3. Pocket mapping and structural QC

Download AlphaFold structures, map the pocket residues, check pocket confidence (pLDDT), and calculate 3D similarity scores.
```bash
python Step3.py \
  --core-fasta data/core.fasta \
  --msa result/msa.fa \
  --pocket data/pocket.txt \
  --motif-file data/motif.txt \
  --outdir outputs
```
**Output**:
- pocket_residues.json : mapping of pocket positions across all sequences
- qc_plddt.csv : pocket quality (mean pLDDT, %low confidence residues)
- pocket_scores.csv : RMSD, shape overlap, sidechain conservation, overall pocket score
- final_result.csv : integrated result (only proteins that pass Step2 + Step3)
---
## Result Interpretation
- **core_vs_all.csv**: Filters candidates based on sequence homology and motif presence.
- **qc_plddt.csv**: Checks structural reliability of the pocket (fail if >30% residues have pLDDT < 70).
- **pocket_scores.csv**: Quantifies structural similarity to the unsafe core protein:
    - RMSD: pocket rigid-body deviation
    - ShapeOverlap: 3D pocket overlap (0–1)
    - SidechainChem: chemical similarity of pocket residues
    - Pocket3D_Score: weighted overall score (higher = more similar)
- **final_result.csv**: The safe candidate proteins (from whitelist organisms) that are both functionally conserved and structurally similar to the unsafe protein of interest.
    - These are the recommended substitutes for wet-lab experiments.

## Scoring Details
### 1. Quality Control (QC) by pLDDT
For each protein pocket region (mapped from the core sequence):
- **Mean pocket pLDDT**:  
  The average pLDDT score of pocket residues.

- **Low-confidence percentage** (`pocket_low70%`):  
  Fraction of pocket residues with pLDDT < 70.

- **Pass/fail rule**:  
  QC passes if  
  \[
  \text{low70\%} \leq 30\%
  \]

Proteins failing QC are excluded from scoring.

---

### 2. Weight Assignment for Pocket Residues

Each aligned residue pair (core ↔ target) receives a weight:

\[
w = \left(\frac{\min(\text{pLDDT}_\text{core}, \text{pLDDT}_\text{target})}{100}\right)^{\gamma} \times M
\]

- **γ (PLDDT_GAMMA)** = 1.5  
  (emphasizes high-confidence sites)

- **Anchor multiplier (M)** = 2.0  
  (applied if residue lies within the anchor motif inside the pocket)

- **Quality metric**: mean of all weights for a protein.

---

### 3. Geometry-Based Scoring

1. **Weighted Kabsch alignment**:  
   Aligns pocket residues using weights `w`.

   - Output: weighted RMSD (wRMSD).  
   - Normalized RMSD:  
     \[
     \text{norm\_wRMSD} = \frac{1}{1 + (wRMSD / 3.0)^2}
     \]
     (maps to [0,1] range, robust to outliers)

2. **Weighted shape overlap (ShapeOverlap)**:  
   - Symmetric coverage between two residue sets within a radius (default = 3.0 Å).  
   - Each side weighted by `w`.

---

### 4. Chemistry-Based Scoring

- **Weighted chemical conservation (SidechainChem)**:  
  For each aligned pair:

  - Score = 1.0 if identical AA  
  - Score = 0.7 if same chemical category (e.g., both hydrophobic)  
  - Score = 0.5 if conservative substitution (e.g., D↔E, K↔R)  
  - Score = 0.0 otherwise  

  Weighted average using `w`.

---

### 5. Final Composite Pocket Score

The final `Pocket3D_Score` is a weighted sum:

\[
\text{Pocket3D\_Score} =
0.45 \times \text{ShapeOverlap} +
0.30 \times \text{norm\_wRMSD} +
0.20 \times \text{SidechainChem} +
0.05 \times \text{Quality}
\]

- **Weights** reflect priority:  
  - Geometry overlap (45%)  
  - RMSD robustness (30%)  
  - Chemistry conservation (20%)  
  - Overall quality (5%)

---

### 6. Output Metrics

For each protein (`protein_id`):

- **RMSD**: Weighted RMSD after alignment  
- **ShapeOverlap**: Weighted geometric overlap  
- **SidechainChem**: Weighted chemical conservation  
- **Pocket3D_Score**: Final composite score