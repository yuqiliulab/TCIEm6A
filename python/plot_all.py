import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import os

# =========================================================
# 字体与风格（在上一版基础上再 +4，投稿级）
# =========================================================
plt.rcParams.update({
    "font.size": 22,          # 原 18 → 22
    "axes.titlesize": 26,     # 原 22 → 26
    "axes.labelsize": 22,     # 原 18 → 22
    "legend.fontsize": 17,    # 原 13 → 17
    "xtick.labelsize": 22,    # 原 18 → 22
    "ytick.labelsize": 22,    # 原 18 → 22
    "lines.linewidth": 3.0    # 稍微加粗，适合 PDF
})

# =========================================================
# 1. 准备所有模型的 y_true / y_score
# =========================================================

models = {}

# ---------- TCIEm6A-full ----------
y_true = np.load("train_runs_TCIE/run8_y_true.npy")
y_score = np.load("train_runs_TCIE/run8_y_score.npy")
models["TCIEm6A-full"] = (y_true, y_score)

# ---------- DeepSRAMP-full ----------
y_true = np.load("train_runs_all/run8_y_true.npy")
y_score = np.load("train_runs_all/run8_y_score.npy")
models["DeepSRAMP-full"] = (y_true, y_score)

# ---------- CNN-full ----------
y_true = np.load("cnn_runs_full/cnn_run8_y_true.npy")
y_score = np.load("cnn_runs_full/cnn_run8_y_score.npy")
models["CNN-full"] = (y_true, y_score)

# ---------- 5 个传统 ML ----------
ml_df = pd.read_csv("ml_runs_norm/run8_predictions.csv")

ml_models = ["GNB-full", "k-NN-full", "LR-full", "PLS-DA-full", "SVM-full"]

for m in ml_models:
    tmp = ml_df[ml_df["Model"].str.contains(m)]
    models[m] = (tmp["y_true"].values, tmp["y_prob"].values)

# =========================================================
# 2. 颜色（论文固定）
# =========================================================

colors = {
    "TCIEm6A-full": "#1f77b4",
    "DeepSRAMP-full": "#ff7f0e",
    "CNN-full": "#2ca02c",
    "GNB-full": "#d62728",
    "k-NN-full": "#9467bd",
    "LR-full": "#8c564b",
    "PLS-DA-full": "#e377c2",
    "SVM-full": "#7f7f7f",
}

# =========================================================
# 3. AUROC（已去掉 0.5 虚线）
# =========================================================

plt.figure(figsize=(8, 8))

for name, (y_true, y_score) in models.items():
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    plt.plot(
        fpr, tpr,
        color=colors[name],
        label=f"{name} (AUROC={roc_auc:.3f})"
    )

# ❌ 删除随机参考线
# plt.plot([0, 1], [0, 1], "--", color="grey", lw=1.5)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("AUROC")
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig("AUROC_run8_all_models.png", dpi=300)
plt.savefig("AUROC_run8_all_models.pdf")
plt.close()

# =========================================================
# 4. AUPRC
# =========================================================

plt.figure(figsize=(8, 8))

for name, (y_true, y_score) in models.items():
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    pr_auc = average_precision_score(y_true, y_score)
    plt.plot(
        recall, precision,
        color=colors[name],
        label=f"{name} (AUPRC={pr_auc:.3f})"
    )

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("AUPRC")
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig("AUPRC_run8_all_models.png", dpi=300)
plt.savefig("AUPRC_run8_all_models.pdf")
plt.close()

print("✅ Saved:")
print("  AUROC_run8_all_models.png / .pdf")
print("  AUPRC_run8_all_models.png / .pdf")
