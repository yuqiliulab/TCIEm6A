# circRNA-m6A Deep Learning Analysis

This project develops a deep learning model to investigate the association between m6A sites and functional elements (ORFs/IRES/RBP/miRNA) in circRNA. The workflow involves R-based preprocessing followed by training and interpretation using the SRAMP model.

## Step 1: Extract ORFs with TransDecoder

```bash
TransDecoder.LongOrfs -t circrna_positive.fasta
TransDecoder.LongOrfs -t circrna_negative.fasta
```

Generates `longest_orfs.cds`.

## Step 2: Align Sequences with R

Load `circrna_positive.rds` and align ORF sequences. Match ORFs to m6A sites in `m6a_circ_olp.rds`.

```r
circrna_positive <- readRDS("circrna_positive.rds")
orf_sequences <- readDNAStringSet("longest_orfs.cds")
circrna_negative <- readRDS("circrna_negative.rds")
orf_sequences <- readDNAStringSet("longest_orfs.cds")
overlaps <- findOverlaps(m6a_circ_olp, circrna_pos)
overlaps <- findOverlaps(m6a_circ_olp_neg, circrna_neg)
```

Determine:
- Whether the m6A lies inside ORFs.
- Distance between m6A and both ends of ORFs.

## Step 3: Check Overlap with IRES

```r
overlap_results_pos <- findOverlaps(m6a_circ_olp_updated, ires_granges)
overlap_results_neg <- findOverlaps(m6a_circ_olp_neg_updated, ires_granges)
```

## Step 4: Check Overlap with RBP

```r
rbp <- readRDS("RBP_site_gr_38.rds")
rbp_list <- c("HNRNPC","YTHDC1","YTHDF1","YTHDF2","METTL3","METTL14",
              "METTL16","WTAP","ALKBH5","FTO")
```

## Step 5: Overlap with miRNA Sites

```r
miRNA_ALL_gr <- readRDS("miRNA_ALL_gr.rds")
mcols(m6a_circ_olp_updated)$miRNA <- FALSE
overlaps <- findOverlaps(m6a_circ_olp_updated, miRNA_ALL_gr)
mcols(m6a_circ_olp_updated)$miRNA[queryHits(overlaps)] <- TRUE
```

## Step 6: Feature Preprocessing (preprocess_m6A_features.R)

Output includes:
- `score_matrices_*.npy`
- `phenotype_features_*.npy`
- `labels_*.npy`

## Step 7: Deep Learning Training & SHAP Interpretation

```
deepSRAMP/
├── config/
│   └── params.py
├── models/
│   └── sramp.py
├── utils/
│   ├── data_loader.py
│   └── metrics.py
├── explain/
│   └── shap_analysis.py
├── train.py
```

## Final Workflow Summary

- Train model via `train.py`.
- Visualize SHAP feature importance via `shap_analysis.py`.
- Analyze positive and negative samples separately.
- Evaluate using metrics such as AUC, accuracy, and F1 score.

## Requirements

Python ≥ 3.7
R ≥ 4.0
PyTorch, NumPy, pandas, SHAP, Bioconductor packages
