
import os
import random
import numpy as np
import pandas as pd
from tqdm import tqdm

import torch
from torch import nn, optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, roc_auc_score, f1_score,
    precision_score, recall_score, matthews_corrcoef
)

# =============================
# 基础配置
# =============================
#POS_CSV = "m6a_pos_one_hot_minmax.csv"
#NEG_CSV = "m6a_neg_one_hot_minmax.csv"
POS_CSV = "shap_最开始minmax/run_data_round8/pos_data.csv"
NEG_CSV = "shap_最开始minmax/run_data_round8/neg_sampled_data.csv"
OUT_DIR = "cnn_runs_full"
os.makedirs(OUT_DIR, exist_ok=True)

# 超参数
N_REPEATS = 10
BATCH_SIZE = 64
EPOCHS = 8
LR = 1e-3
TEST_SIZE = 0.2
SEED = 42
SEQ_LEN = 41
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


# =============================
# One-hot 编码函数
# =============================
def one_hot_encode(seq):
    """将 RNA 序列编码为 one-hot，通道顺序为 [A, C, G, T/U]"""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    one_hot = np.zeros((len(seq), 4), dtype=np.float32)
    for i, base in enumerate(seq.upper()):
        if base in mapping:
            one_hot[i, mapping[base]] = 1.0
    return one_hot


# =============================
# Dataset
# =============================
class SeqFeatureDataset(Dataset):
    def __init__(self, df, seq_col="seq_count", label_col="label", L=41):
        self.seq_col = seq_col
        self.label_col = label_col
        self.L = L
        self.feature_cols = [c for c in df.columns if c not in [seq_col, label_col]]
        self.seq_list = df[seq_col].astype(str).tolist()
        self.features = df[self.feature_cols].fillna(0.0).to_numpy(dtype=np.float32)
        self.labels = df[label_col].astype(np.float32).to_numpy()

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        seq = self.seq_list[idx]
        seq = seq[:self.L].ljust(self.L, 'N')  # pad/truncate
        onehot = one_hot_encode(seq)           # (L,4)
        label = np.array([self.labels[idx]], dtype=np.float32)
        return (
            torch.tensor(onehot).float(),          # (L,4)
            torch.tensor(self.features[idx]).float(), # (feature_dim,)
            torch.tensor(label).float()            # (1,)
        )


# =============================
# CNN 模型
# =============================
class CNN_SRAMPMixed(nn.Module):
    def __init__(self, feature_dim, seq_channels=4, num_classes=1):
        """
        feature_dim: 数值基因组特征的维度
        seq_channels: one-hot 通道数 (A,C,G,T/U -> 4)
        """
        super().__init__()
        in_channels = seq_channels + feature_dim  # 拼接后的通道数
        self.conv = nn.Sequential(
            nn.Conv1d(in_channels, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(64, 128, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1)
        )
        self.fc = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, num_classes),
            nn.Sigmoid()
        )

    def forward(self, seq_onehot, features):
        """
        seq_onehot: (B, L, 4)
        features: (B, feature_dim)
        """
        B, L, _ = seq_onehot.shape
        if features.dim() == 2:
            features = features.unsqueeze(1).repeat(1, L, 1)  # (B, L, feature_dim)
        x = torch.cat([seq_onehot, features], dim=-1)         # (B, L, 4+feature_dim)
        x = x.permute(0, 2, 1)                                # (B, C, L)
        x = self.conv(x)
        x = x.view(x.size(0), -1)
        out = self.fc(x)
        return out


# =============================
# 训练 & 评估函数
# =============================
def train_one_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    for seq, feat, label in loader:
        seq, feat, label = seq.to(device), feat.to(device), label.to(device)
        optimizer.zero_grad()
        out = model(seq, feat)
        loss = criterion(out, label)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * seq.size(0)
    return total_loss / len(loader.dataset)


def eval_model(model, loader, device):
    model.eval()
    ys, probs = [], []
    with torch.no_grad():
        for seq, feat, label in loader:
            seq, feat = seq.to(device), feat.to(device)
            out = model(seq, feat)
            ys.extend(label.squeeze(-1).cpu().numpy())
            probs.extend(out.squeeze(-1).cpu().numpy())
    ys, probs = np.array(ys), np.array(probs)
    preds = (probs >= 0.5).astype(int)
    metrics = {
        "ACC": float(accuracy_score(ys, preds)),
        "AUC": float(roc_auc_score(ys, probs)) if len(np.unique(ys)) > 1 else float("nan"),
        "F1": float(f1_score(ys, preds, zero_division=0)),
        "Prec": float(precision_score(ys, preds, zero_division=0)),
        "Recall": float(recall_score(ys, preds, zero_division=0)),
        "MCC": float(matthews_corrcoef(ys, preds)) if len(np.unique(preds)) > 1 else float("nan"),
        "N": int(len(ys)),
        "Pos_ratio": float(ys.mean())
    }
    return metrics


# =============================
# 主流程
# =============================
def main():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    pos_df = pd.read_csv(POS_CSV)
    neg_df = pd.read_csv(NEG_CSV)
    print(f"正样本: {len(pos_df)}，负样本: {len(neg_df)}")

    results_all = []

    for run in range(N_REPEATS):
        print("=" * 50)
        print(f"Run {run+1}/{N_REPEATS}")

        # 随机抽取负样本与正样本平衡
        neg_sample = neg_df.sample(n=len(pos_df), replace=False, random_state=SEED+run).reset_index(drop=True)
        df = pd.concat([pos_df, neg_sample], axis=0).reset_index(drop=True)

        train_df, test_df = train_test_split(
            df, test_size=TEST_SIZE, stratify=df["label"], random_state=SEED+run
        )

        train_ds = SeqFeatureDataset(train_df, "seq_count", "label", L=SEQ_LEN)
        test_ds = SeqFeatureDataset(test_df, "seq_count", "label", L=SEQ_LEN)
        train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True, num_workers=0)
        test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=0)

        feature_dim = len(train_ds.feature_cols)
        print(f"Detected feature_dim = {feature_dim}")

        model = CNN_SRAMPMixed(feature_dim=feature_dim).to(DEVICE)
        criterion = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=LR)

        best_auc = -1
        best_metrics = None

        for epoch in range(EPOCHS):
            loss = train_one_epoch(model, train_loader, optimizer, criterion, DEVICE)
            val_metrics = eval_model(model, test_loader, DEVICE)
            print(f"Run {run+1}  Epoch {epoch+1}/{EPOCHS}  "
                  f"loss={loss:.4f}  valAUC={val_metrics['AUC']:.4f}  ACC={val_metrics['ACC']:.4f}")

            # 保存最优模型
            if val_metrics['AUC'] > best_auc:
                best_auc = val_metrics['AUC']
                best_metrics = val_metrics.copy()
                torch.save(model.state_dict(), os.path.join(OUT_DIR, f"cnn_run{run+1}_best.pth"))

        best_metrics['run'] = run + 1
        results_all.append(best_metrics)
        pd.DataFrame([best_metrics]).to_csv(os.path.join(OUT_DIR, f"metrics_run{run+1}.csv"), index=False)
        print(f"✅ Saved best model for run {run+1}")

    # 汇总结果
    dfres = pd.DataFrame(results_all)
    dfres.to_csv(os.path.join(OUT_DIR, "summary_runs.csv"), index=False)
    print("🎯 All done! Summary saved to:", os.path.join(OUT_DIR, "summary_runs.csv"))


if __name__ == "__main__":
    main()
