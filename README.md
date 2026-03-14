# Genome wide eccDNA hotspot propensity estimation

This repository contains the Jupyter notebook used to construct genome wide eccDNA hotspot propensity targets, annotate genomic bins with locus features, and train predictive models on hg38.

The main notebook is organized as a complete analysis workflow, from raw eccDNA interval tables to model evaluation and interpretation.

---

## 1. Notebook organization

The notebook is organized into the following main sections:

### 1. Configuration
This section defines:
- random seed
- bin size (`BIN_SIZES`)
- sub-bin resolution used for bigWig summaries (`SUBBIN_BP`)
- target parameters such as pseudocount (`EPS`) and number of shuffles (`N_SHUFFLES`)
- chromosome split used for evaluation
- paths to all required input files
- output directory

### 2. Genome bin construction
This part builds a fixed genomic grid on hg38, typically using 25 kb bins, and excludes invalid bins overlapping:
- assembly gaps
- blacklist regions

### 3. eccDNA event loading and preprocessing
This section loads eccDNA intervals from public sources and harmonizes them into a common coordinate representation. Events are grouped by experiment using the `(PubMed ID, assay)` combination, which defines the replicate unit used throughout the analysis.

### 4. Target construction
For each genomic bin, the notebook computes:
- **observed replicate support**, i.e. how many independent experiments report at least one eccDNA interval in the bin
- **chance support**, estimated with a shuffle-based null model that repeatedly randomizes interval positions within the same chromosome while preserving interval lengths
- a continuous enrichment target derived from observed support relative to chance support

### 5. Feature annotation
Each valid genomic bin is annotated with locus-level features derived from:
- reference sequence
- gene annotation
- repeat annotation
- mappability
- centromere/telomere context
- optional additional genome tracks

These features are divided into:
- **primary interpretable covariates** used for modelling
- **technical or genome-complexity covariates** used mainly for filtering and bias diagnostics

### 6. Dataset assembly
The notebook merges targets and features into modelling tables, including:
- genome-wide datasets
- optional assay-specific datasets
- sequence-feature augmented versions used in downstream experiments

### 7. Model training and evaluation
This section trains regression models on the continuous enrichment target and evaluates:
- random split performance
- chromosome holdout performance
- ablation settings
- permutation feature importance

### 8. DeepCircle-like benchmark
This final part builds a protocol-matched comparison dataset and evaluates lightweight models under a 1 kb balanced classification setting inspired by DeepCircle.

---

## 2. Packages to install

The notebook requires Python and the following main packages:

- `numpy`
- `pandas`
- `tqdm`
- `pyranges`
- `pyBigWig`
- `pyfaidx`
- `scikit-learn`
- `matplotlib`
- `torch`

A minimal installation can be done with:

```bash
pip install numpy pandas tqdm pyranges pyBigWig pyfaidx scikit-learn matplotlib torch
