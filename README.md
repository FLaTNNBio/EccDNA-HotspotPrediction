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
```

## Required input data

The notebook is driven by a configuration cell where all input paths are defined.  
Inputs fall into three groups:

1. **Required core files** for genome binning and target construction  
2. **Required eccDNA event tables** for building observed replicate support  
3. **Optional annotation tracks** for feature extraction and diagnostics  

### 1. Required core files

These files are necessary to define the genomic grid and exclude invalid regions.

| Variable | Expected file | Description | Required |
|---|---|---|---|
| `CHROM_SIZES_PATH` | `hg38.chrom.sizes` | Chromosome sizes for hg38 | Yes |
| `GAPS_BED_PATH` | `hg38_gaps.bed` | Assembly gaps in BED format (0-based, half-open) | Yes |
| `BLACKLIST_BED_PATH` | `encBlacklist.bed` | ENCODE blacklist regions in BED format (0-based, half-open) | Yes |

These files are used to:
- build the 25 kb genomic grid
- restrict the analysis to standard chromosomes
- remove bins overlapping assembly gaps or blacklist regions

---

### 2. Required eccDNA event tables

These are the main biological input tables used to compute replicate support and the enrichment target.

| Variable | Expected file | Description | Required |
|---|---|---|---|
| `CIRCLE_PATH` | `Circle_Seq.tsv.gz` | eccDNA intervals from Circle-seq experiments | Yes |
| `ATAC_PATH` | `ATAC_seq.tsv.gz` | eccDNA-related intervals from ATAC-seq experiments | Yes |
| `WGS_PATH` | `WGS.tsv.gz` | eccDNA/ecDNA intervals from whole-genome sequencing experiments | Yes |

Each event table should be a TSV or TSV.GZ file and contain at least the following columns:

| Column | Meaning | Required |
|---|---|---|
| `chr_hg38` | Chromosome name in hg38 coordinates | Yes |
| `start_hg38` | Interval start coordinate | Yes |
| `end_hg38` | Interval end coordinate | Yes |
| `pubmed_id` | Study identifier used to define replicate provenance | Yes |
| `library_type` | Assay label (e.g. Circle-seq, ATAC-seq, WGS) | Yes |
| `Length` | Interval length | Optional |

These tables are used to:
- define replicate units as unique `(PubMed ID, assay)` pairs
- compute observed replicate support for each bin
- build the shuffle-based null model

---

### 3. Required reference sequence and gene annotation

These files are needed to compute primary genomic features.

| Variable | Expected file | Description | Required |
|---|---|---|---|
| `FASTA_HG38` | `hg38.fa` | hg38 reference FASTA used for sequence-derived features | Yes |
| `GTF_GENCODE` | `gencode.v49.annotation.gtf.gz` | GENCODE annotation used for gene and TSS features | Yes |

These files are used to compute:
- GC fraction
- CpG-related summaries
- strand skews
- gene coverage
- gene count
- transcription start site density
- sequence complexity descriptors

---

### 4. Optional annotation tracks

These files are not strictly required to run the core target construction, but they are used for feature extraction, diagnostics, and bias assessment.

| Variable | Expected file | Description | Required |
|---|---|---|---|
| `RMSK_TSV` | `rmsk.txt.gz` | RepeatMasker annotations | Optional |
| `SEGDUP_BED` | `hg38_genomicSuperDups.bed` | Segmental duplication intervals | Optional |
| `MAPPABILITY_BW` | `k100.Umap.MultiTrackMappability.bw` | Umap k=100 mappability track | Optional |
| `CENTROMERE_FILE` | `centromeres_hg38.txt.gz` | Centromere intervals for hg38 | Optional but strongly recommended |
| `REPTIME_BW` | replication timing bigWig | Optional replication timing signal | Optional |

These resources are used to compute:
- repeat content
- segmental duplication coverage
- mappability summaries
- centromere distance
- optional replication timing features

---

### 5. Minimal set to reproduce the main paper results

To reproduce the main results reported in the manuscript, the minimum required inputs are:

- `hg38.chrom.sizes`
- `hg38_gaps.bed`
- `encBlacklist.bed`
- `Circle_Seq.tsv.gz`
- `ATAC_seq.tsv.gz`
- `WGS.tsv.gz`
- `hg38.fa`
- `gencode.v49.annotation.gtf.gz`
- `centromeres_hg38.txt.gz`

The following files are strongly recommended for the full feature set and bias diagnostics:

- `rmsk.txt.gz`
- `hg38_genomicSuperDups.bed`
- `k100.Umap.MultiTrackMappability.bw`


## Running the main regression models

The main modelling results in the manuscript are obtained from the processed genome-wide dataset built in the earlier sections of the notebook.  
Once bin construction, target construction, and feature annotation have been completed, the notebook trains regression models directly on the final merged table.

### Main input dataset

The main input used for model training is the processed file:

- `out_bins_experiment/dataset_global_hg38_25000_seqfeat.csv.gz`

This table contains one row per genomic bin and includes:
- bin coordinates
- validity flags
- observed support and null-based target columns
- sequence-derived features
- gene and transcription start site features
- positional features
- optional technical or genome-complexity covariates

The primary target used for modelling is:

- `y_log2`

which corresponds to the continuous enrichment signal derived from observed replicate support relative to the shuffle-based null expectation.

---

### Which models are trained

The notebook trains and evaluates the following regression models:

- `HGB_reg` (Histogram Gradient Boosting Regressor)
- `MLP_reg` (Multilayer Perceptron regressor)
- `Ridge`

These models are used for:
- the main random split benchmark
- the primary versus full feature comparison
- the positional ablation experiment
- permutation feature importance analysis

---

### Main modelling settings

The standard regression evaluation uses:

- `VAL_FRAC = 0.20`
- `TEST_FRAC = 0.20`
- `N_REPEATS = 10`
- `Y_COL = "y_log2"`

The benchmark is stratified by deciles of the target and repeated across multiple random seeds.

---

### Feature configurations

The notebook evaluates two feature configurations:

#### Primary features
These correspond to the interpretable locus covariates used in the main paper results:
- sequence composition and CpG context
- gene and transcription start site summaries
- telomere and centromere distances
- sequence complexity descriptors

#### Full features
These extend the primary set with technical / genome-complexity covariates, including:
- repeat fractions
- mappability summaries
- segmental duplication coverage

The **Primary** configuration is the main analysis setting.  
The **Full** configuration is used as a diagnostic comparison.

---


## DeepCircle matched comparison

The notebook also includes a protocol-matched comparison against DeepCircle-like datasets.  
This benchmark is meant to complement the main genome-wide regression setting by evaluating lightweight models on fixed 1 kb windows under a balanced classification protocol.

### What this section does

The DeepCircle comparison reproduces a classification setting in which:

- positive examples are genomic windows associated with reported eccDNA events
- negative examples are matched genomic windows not overlapping eccDNA regions
- all windows are represented with engineered features rather than raw sequence
- models are trained and evaluated under a DeepCircle-like protocol

This analysis is implemented in the notebook section:

- **`# DeepCircle confronto`**

and is organized into the following steps:

1. **Dataset discovery and setup**  
   Load DeepCircle BED files and define paths, genome resources, and split logic.

2. **Window construction**  
   Convert positive events into fixed 1 kb windows and generate matched negative windows.

3. **Feature extraction**  
   Compute sequence-derived and genomic-context features for each 1 kb window.

4. **Train/test split export**  
   Save processed datasets and split indices.

5. **Model training and evaluation**  
   Train lightweight classifiers and compute the classification metrics reported in the manuscript.

---

### Files required for the DeepCircle comparison

The following files are required in addition to the main genome resources already described above.

#### A. DeepCircle BED files

The notebook expects DeepCircle-style BED inputs, typically organized under a directory such as:

- `DeepCircle/DNABERT/eccdna/db/.../*.bed`

These BED files provide the positive eccDNA regions used to construct the benchmark windows.

You should update the corresponding root path in the notebook so that it points to your local DeepCircle dataset directory.

#### B. Genome reference and masking files

These are required to construct valid windows and sequence-based features:

- `hg38.fa`
- `hg38.chrom.sizes`
- `hg38_gaps.bed`
- `encBlacklist.bed`

These files are used to:
- build and validate 1 kb windows
- remove invalid regions
- extract sequence-derived features

#### C. Optional annotation tracks for engineered features

If you want the full engineered feature set for the DeepCircle benchmark, the same annotation resources used in the main analysis should also be available:

- `gencode.v49.annotation.gtf.gz`
- `centromeres_hg38.txt.gz`
- `rmsk.txt.gz`
- `hg38_genomicSuperDups.bed`
- `k100.Umap.MultiTrackMappability.bw`

These are used only if the selected benchmark feature set includes genomic context and diagnostic covariates.

---

### Minimum data needed to run this benchmark

At minimum, you need:

- DeepCircle BED files
- `hg38.fa`
- `hg38.chrom.sizes`
- `hg38_gaps.bed`
- `encBlacklist.bed`

If you want to reproduce the same engineered-feature benchmark reported in the manuscript, you should also provide:

- `gencode.v49.annotation.gtf.gz`
- `centromeres_hg38.txt.gz`
- optional repeat, mappability, and segmental duplication tracks

---

### How to run it

After completing the main configuration and ensuring all paths are correct:

1. Open the notebook section:
   - **`# DeepCircle confronto`**

2. Update the input paths for:
   - DeepCircle BED files
   - genome FASTA
   - chromosome sizes
   - masking files
   - optional annotation tracks

3. Run the DeepCircle cells in order:
   - dataset setup
   - window construction
   - feature extraction
   - split export
   - model training and evaluation

No manual reordering is required if the cells are run sequentially.
