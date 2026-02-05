# circRNA-m6A Deep Learning Analysis

We present here TCIEm6A for “translatable circRNAs identification and editing via m6A modification”. For the first time, we explored different strategies for encoding sub-molecular geographic information of translation-enhancing elements on circRNAs, from which a prediction framework was built upon a hybrid deep learning architectures combining Transformer, recurrent neural network, and LSTM-attention. 
### System
- Linux / macOS / Windows (WSL recommended)
- Python ≥ 3.7
- R ≥ 4.0

---

### Python Dependencies
The deep learning framework is implemented based on PyTorch and related scientific libraries.

- torch ≥ 1.8
- numpy
- pandas
- scikit-learn
- shap
- tqdm
- matplotlib
- seaborn

Optional (for GPU acceleration):
- CUDA ≥ 11.0
- cuDNN

---

### R Dependencies
All genomic feature extraction and annotation steps are performed in R.

#### CRAN packages
- data.table
- dplyr
- tidyr
- ggplot2
- stringr
- reticulate
- rtracklayer

#### Bioconductor packages
- GenomicRanges
- IRanges
- Biostrings
- GenomicFeatures
- BSgenome
- AnnotationDbi

---

### External Tools
- TransDecoder (v5.5.0 or later)

TransDecoder is used to extract candidate ORFs from circRNA sequences.
## Step 1: Extract ORFs with TransDecoder

```bash
TransDecoder.LongOrfs -t circ_sequence.fa
```

Generates `longest_orfs.cds`.

## Step 2: Align Sequences with R

Load `circrna_positive.rds` and align ORF sequences. Match ORFs to m6A sites in `m6a_circ_olp.rds`.

```r
hits <- findOverlaps(m6a_sites, result_gr)
m6a_sites$ORF <- NA; m6a_sites$ORF_Count <- 0; m6a_sites$ORF_Length <- 0
m6a_sites$ORF[hits@from] <- mcols(result_gr)$ORF_Genomic[hits@to]
m6a_sites$ORF_Count[hits@from] <- mcols(result_gr)$ORF_Count[hits@to]
m6a_sites$ORF_Length[hits@from] <- mcols(result_gr)$ORF_Length[hits@to]
```

Determine:
- Whether the m6A lies inside ORFs.
- Distance between m6A and both ends of ORFs.

## Step 3: Check Overlap with IRES

```r
hits <- findOverlaps(m6a_sites, result_gr)
m6a_sites$IRES <- NA; m6a_sites$IRES_Count <- 0; m6a_sites$IRES_Length <- 0
m6a_sites$IRES[hits@from] <- mcols(result_gr)$IRES[hits@to]
m6a_sites$IRES_Count[hits@from] <- mcols(result_gr)$IRES_Count[hits@to]
m6a_sites$IRES_Length[hits@from] <- mcols(result_gr)$IRES_Length[hits@to]
```

## Step 4: Check Overlap with RBP

```r
rbp <- readRDS("RBP_site_gr_38.rds")
rbp_list <- c("HNRNPC","YTHDC1","YTHDF1","YTHDF2","YTHDF3","METTL3","METTL14",
              "METTL16","WTAP","ALKBH5","FTO")

```

## Step 5: Overlap with miRNA、phastcons、SplicingSite_Num Sites

```r
miRNA_ALL_gr <- readRDS("miRNA_ALL_gr.rds")
mcols(candidate_sites)$miRNA <- FALSE
overlaps <- findOverlaps(candidate_sites, miRNA_ALL_gr)
mcols(candidate_sites)$miRNA[queryHits(overlaps)] <- TRUE
hg38_splicing_sites_100bp_gr <- readRDS("hg38_splicing_sites_100bp_gr.rds")
mcols(candidate_sites)$SplicingSite_Num <- FALSE
overlaps <- findOverlaps(candidate_sites, hg38_splicing_sites_100bp_gr)
mcols(candidate_sites)$SplicingSite_Num[queryHits(overlaps)] <- TRUE

```

## Step 6: Feature Preprocessing (preprocess_m6A_features.R)
```
Rscript server script/veb.R example output.fa test fasta 01 data/output mode2
```

## Step 7: Deep Learning Training & SHAP Interpretation

```
python predict_m6acirc.py test002 /you/path/data/output
```

## Final Workflow Summary

- Train model via `predict_m6acirc.py`.
- Visualize SHAP feature importance via `shap_analysis.py`.
- Evaluate using metrics such as AUC, accuracy, and F1 score.


