# TCIEm6A

We present  TCIEm6A, the first deep learning framework for “translatable circRNAs identification and editing via m6A modification” at base resolution. For the first time, we explored different strategies for encoding the sub-molecular topological genomic contexts of translation-enhancing elements on circRNA transcripts, building a prediction framework upon a hybrid architecture combining Transformer, recurrent neural network, and LSTM-attention.
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
## Step 1: Data Pre-processing
CircRNA sequences are first annotated and converted into model input format.

Example command:
```bash
Rscript server_script/web1.R circrna_positive.fa test_positive data/output mode1
```

This step generates processed sequence files for downstream model training.

## Step 2: Model Training



```python
python sramp2.py
```

This script:

1.loads positive and negative datasets

2.performs random negative sampling

3.trains the SRAMP-based deep learning model

4.saves the best model checkpoints

## Step 3: Evaluation

Model performance is evaluated automatically during training.

Evaluation metrics include:

AUROC

AUPRC

Accuracy

Precision

Recall

F1-score

MCC

Results are saved in:

train_runs_all/

including:

metrics_run*.csv
summary_runs.csv
sramp_run*_best.pth


