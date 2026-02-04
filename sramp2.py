# train_deepsramp_custom.py
import os
import random
import math
import numpy as np
import pandas as pd

import torch
from torch import nn, optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, roc_auc_score, f1_score,average_precision_score,
    precision_score, recall_score, matthews_corrcoef
)

# ------- 修改文件 -------
POS_CSV = "m6a_pos_one_hot_minmax.csv"
NEG_CSV = "m6a_neg_one_hot_minmax.csv"
#NEG_CSV = "train_runs_9_seq/neg_sample_run8.csv"
#POS_CSV = "m6a_pos_norm.csv"
#NEG_CSV = "m6a_neg_norm.csv"
OUT_DIR = "train_runs_all"   # 保存结果的文件夹
os.makedirs(OUT_DIR, exist_ok=True)

# ------- 超参数 -------
N_REPEATS = 10                # 总共抽负样本训练多少遍
BATCH_SIZE = 64
EPOCHS = 8
LR = 1e-3
TEST_SIZE = 0.2               # 划分测试集比例
SEED = 42
L = 41                       # 序列定长（模型输入序列长度）
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# ------- 取碱基到索引的映射（保留 N 和 gap） -------
base2idx = {'A':0, 'C':1, 'G':2, 'T':3, 'U':3, 'N':4, '-':5}

def seq_to_idx(seq, L):
    """把字符串序列转为长度为 L 的整数索引向量（居中截断或居中pad，pad用'N'）"""
    s = seq.strip().upper().replace('U', 'T')  # 把 U 当做 T 处理
    if len(s) >= L:
        # 取中间 L 个碱基
        start = (len(s) - L) // 2
        sub = s[start:start+L]
    else:
        pad_left = (L - len(s)) // 2
        pad_right = L - len(s) - pad_left
        sub = ('N' * pad_left) + s + ('N' * pad_right)
    # map to ints (如果出现未知字符，映为 'N')
    ids = [base2idx.get(ch, base2idx['N']) for ch in sub]
    return np.array(ids, dtype=np.int64)  # long for torch

class SeqGenomeDataset(Dataset):
    """
    Dataset 返回三元组 (seq_idx_tensor, emb_tensor, label_tensor)
    - seq_idx_tensor: LongTensor (L,)   — 输入到 model 中，model 内会 one-hot
    - emb_tensor: FloatTensor (L, emb_dim)  — 基因组特征重复到每个位置
    - label: FloatTensor (1,)
    """
    def __init__(self, df, seq_col="seq_count", label_col="label", L=41):
        # df: 包含序列列、label列、其余列为基因组特征
        self.seq_col = seq_col
        self.label_col = label_col
        self.L = L
        # genomic features 列：所有除 seq_col, label_col 的数值列
        cols = [c for c in df.columns if c not in [seq_col, label_col]]
        self.emb_cols = cols
        # convert DataFrame to arrays for speed
        self.seq_list = df[seq_col].astype(str).tolist()
        self.emb_array = df[cols].fillna(0.0).to_numpy(dtype=np.float32)
        self.labels = df[label_col].astype(np.float32).to_numpy()
        # emb_dim
        self.emb_dim = self.emb_array.shape[1]
    def __len__(self):
        return len(self.labels)
    def __getitem__(self, idx):
        seq = self.seq_list[idx]
        seq_idx = seq_to_idx(seq, self.L)                 # (L,)
        emb_vec = self.emb_array[idx]                     # (emb_dim,)
        emb_rep = np.tile(emb_vec.reshape(1, -1), (self.L, 1))  # (L, emb_dim)
        label = np.array([self.labels[idx]], dtype=np.float32)  # (1,)
        return torch.from_numpy(seq_idx).long(), torch.from_numpy(emb_rep).float(), torch.from_numpy(label).float()

from model import SRAMP
def eval_model_with_scores(model, loader, device):
    model.eval()
    ys, probs = [], []
    with torch.no_grad():
        for seq_idx, emb, label in loader:
            seq_idx = seq_idx.to(device)
            emb = emb.to(device)
            out = model(seq_idx, emb)

            if out.dim() == 2 and out.size(1) == 1:
                p = out.squeeze(1).cpu().numpy()
            else:
                p = out.cpu().numpy().ravel()

            ys.extend(label.squeeze(-1).numpy().tolist())
            probs.extend(p.tolist())

    ys = np.array(ys)
    probs = np.array(probs)
    return ys, probs

# --------- 训练与评估函数 ----------
def train_one_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0.0
    for seq_idx, emb, label in loader:
        seq_idx = seq_idx.to(device)             # (B, L) long
        emb = emb.to(device)                     # (B, L, emb_dim)
        label = label.to(device)                 # (B, 1)
        optimizer.zero_grad()
        out = model(seq_idx, emb)                # 假设 model 返回概率 (B,1)
        # 有些 SRAMP 可能返回 shape (B,1) 或 (B,)，确保形状匹配
        if out.dim() == 2 and out.size(1) == 1:
            pred = out
        elif out.dim() == 1:
            pred = out.unsqueeze(-1)
        else:
            pred = out  # try
        loss = criterion(pred, label)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * seq_idx.size(0)
    return total_loss / len(loader.dataset)

def eval_model(model, loader, device):
    model.eval()
    ys, probs = [], []
    with torch.no_grad():
        for seq_idx, emb, label in loader:
            seq_idx = seq_idx.to(device)
            emb = emb.to(device)
            out = model(seq_idx, emb)
            if out.dim() == 2 and out.size(1) == 1:
                p = out.squeeze(1).cpu().numpy()
            else:
                p = out.cpu().numpy().ravel()
            ys.extend(label.squeeze(-1).numpy().tolist())
            probs.extend(p.tolist())
    ys = np.array(ys)
    probs = np.array(probs)
    preds = (probs >= 0.5).astype(int)
    metrics = {
        "ACC": float(accuracy_score(ys, preds)),
        "AUC": float(roc_auc_score(ys, probs)) if len(np.unique(ys))>1 else float("nan"),
        "F1": float(f1_score(ys, preds, zero_division=0)),
        "Prec": float(precision_score(ys, preds, zero_division=0)),
        "Recall": float(recall_score(ys, preds, zero_division=0)),
        "MCC": float(matthews_corrcoef(ys, preds)) if len(np.unique(preds))>1 else float("nan"),
        "AUPRC": float(average_precision_score(ys, probs)),  # ⭐ 新增
        "N": int(len(ys)),
        "Pos_ratio": float(ys.mean())
    }
    return metrics

# --------- 主流程：10 次重复采样训练并记录结果 ---------
def main():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    # 读文件
    pos_df = pd.read_csv(POS_CSV)
    neg_df = pd.read_csv(NEG_CSV)
    # 检查列
    assert "seq_count" in pos_df.columns or "seq" in pos_df.columns, "找不到序列列，请确认列名为 seq_count 或 seq"
    seq_col = "seq_count" if "seq_count" in pos_df.columns else "seq"
    label_col = "label"
    # 统计正样本数
    n_pos = len(pos_df)
    print(f"正样本数: {n_pos}, 负样本总数: {len(neg_df)}")
    # ==== 检查列名和数量 ====
    print("\n[列检查]")
    print("pos_df 列数:", len(pos_df.columns))
    print("neg_df 列数:", len(neg_df.columns))
    print("pos_df 列名:", pos_df.columns.tolist())
    print("neg_df 列名:", neg_df.columns.tolist())

    # 检查是否有空列或未命名列
    empty_cols = [c for c in pos_df.columns if "Unnamed" in c or c.strip() == ""]
    if empty_cols:
        print("⚠️ 检测到未命名或空列:", empty_cols)

    # 检查特征列数量（不含 seq_col 和 label_col）
    feature_cols = [c for c in pos_df.columns if c not in [seq_col, label_col]]
    print(f"特征列数（不含 seq_count 和 label）: {len(feature_cols)}")

    # 显示前几行数据，确认是否有异常值
    print("\npos_df 前几行：")
    print(pos_df.head(2))
    print(pos_df.columns)


    results_all = []
    for run in range(N_REPEATS):
        print("="*40)
        print(f"Run {run+1}/{N_REPEATS}")

        # 1) 从 neg 中随机抽取 n_pos 条
        neg_sample = neg_df.sample(n=n_pos, replace=False, random_state=SEED+run).reset_index(drop=True)
        print(neg_sample.columns)
        neg_file = os.path.join(OUT_DIR, f"neg_sample_run{run + 1}.csv")
        neg_sample.to_csv(neg_file, index=False)
        print(f"Saved negative sample for run {run + 1} to {neg_file}")

        # 2) 合并
        train_df = pd.concat([pos_df, neg_sample], axis=0).reset_index(drop=True)
        # 3) 将 label 放在最后或确保存在
        assert label_col in train_df.columns, "label 列缺失"
        # 4) split 为训练/测试
        train_sub, test_sub = train_test_split(train_df, test_size=TEST_SIZE, stratify=train_df[label_col], random_state=SEED+run)
        # 5) Dataset & Dataloader
        train_ds = SeqGenomeDataset(train_sub, seq_col, label_col, L=L)
        test_ds = SeqGenomeDataset(test_sub, seq_col, label_col, L=L)
        train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True, num_workers=0)
        test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=0)
        for seq_idx, emb, label in train_loader:
            print("seq_idx shape:", seq_idx.shape)
            print("emb shape:", emb.shape)
            print("label shape:", label.shape)
            print("seq_idx[0]:", seq_idx[0])  # 打印第一条序列索引看看是否正常
            print("emb_dim:", emb.shape[-1])
            break
        # 6) model & optimizer & loss
        model = SRAMP(mode='full', halfseqlen=(L-1)//2, oh=True).to(DEVICE)  # halfseqlen 与 L 对应
        criterion = nn.BCELoss()   # 因为 SRAMP 默认在返回前用了 sigmoid
        optimizer = optim.Adam(model.parameters(), lr=LR)

        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.5, patience=2, verbose=True
        )
        best_val_auc = -1
        best_metrics = None
        # 7) 训练
        for epoch in range(EPOCHS):
            train_loss = train_one_epoch(model, train_loader, optimizer, criterion, DEVICE)
            val_metrics = eval_model(model, test_loader, DEVICE)
            print(f"Run {run+1} Epoch {epoch+1}/{EPOCHS}  loss={train_loss:.4f}  valAUC={val_metrics['AUC']:.4f} ACC={val_metrics['ACC']:.4f}")
            # 保存最好模型（按 AUC）
            if not math.isnan(val_metrics['AUC']) and val_metrics['AUC'] > best_val_auc:
                best_val_auc = val_metrics['AUC']
                best_metrics = val_metrics.copy()
                torch.save(model.state_dict(), os.path.join(OUT_DIR, f"sramp_run{run+1}_best.pth"))
                # ====== 如果是第 8 次 run，保存 y_true 和 y_score ======
                if run + 1 == 8:
                    print("Saving ROC/PR data for run 8 ...")

                    best_model = SRAMP(mode='full', halfseqlen=(L - 1) // 2, oh=True).to(DEVICE)
                    best_model.load_state_dict(
                        torch.load(os.path.join(OUT_DIR, "sramp_run8_best.pth"), map_location=DEVICE)
                    )

                    y_true, y_score = eval_model_with_scores(best_model, test_loader, DEVICE)

                    np.save(os.path.join(OUT_DIR, "run8_y_true.npy"), y_true)
                    np.save(os.path.join(OUT_DIR, "run8_y_score.npy"), y_score)

                    print("Saved run8_y_true.npy and run8_y_score.npy")

        # 8) 记录
        if best_metrics is None:
            best_metrics = eval_model(model, test_loader, DEVICE)
        best_metrics['run'] = run+1
        results_all.append(best_metrics)
        # 写每轮结果到 csv
        pd.DataFrame([best_metrics]).to_csv(os.path.join(OUT_DIR, f"metrics_run{run+1}.csv"), index=False)
        print("Saved best model/metrics for run", run+1)
    # 汇总所有 run 的结果
    dfres = pd.DataFrame(results_all)
    dfres.to_csv(os.path.join(OUT_DIR, "summary_runs.csv"), index=False)
    print("All done. Summary saved to", os.path.join(OUT_DIR, "summary_runs.csv"))

if __name__ == "__main__":
    main()
