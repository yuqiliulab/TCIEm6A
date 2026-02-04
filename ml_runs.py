import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, precision_score, recall_score, matthews_corrcoef
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVC
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

# -------------------- 配置 --------------------
NEG_CSV = "shap_最开始minmax/run_data_round8/neg_sample.csv"
POS_CSV = "m6a_pos_one_hot_minmax.csv"
OUT_DIR = "ml_runs_norm"
os.makedirs(OUT_DIR, exist_ok=True)

SEED = 42
TEST_SIZE = 0.2
L = 41  # 序列长度（如果不够会自动填充N）

# -------------------- One-hot 编码 --------------------
def one_hot_encode(seq, alphabet="ACGTU"):
    mapping = {ch: i for i, ch in enumerate(alphabet)}
    onehot = np.zeros((L, len(alphabet)), dtype=np.float32)
    seq = seq.upper()[:L].ljust(L, 'N')  # 截断/填充
    for i, base in enumerate(seq):
        if base in mapping:
            onehot[i, mapping[base]] = 1.0
    return onehot.flatten()  # 展平成1维向量

# -------------------- 模型列表 --------------------
models = {
    "E_GNB": lambda: ("GaussianNB", GaussianNB()),
    "F_KNN": lambda: ("KNN", KNeighborsClassifier(n_neighbors=5)),
    "H_Logit": lambda: ("LogisticRegression", LogisticRegression(max_iter=500, random_state=SEED)),
    "I_PLSDA": lambda: ("PLSDA", PLSRegression(n_components=2)),
    "K_SVM": lambda: ("SVM", SVC(probability=True, random_state=SEED)),
}

# -------------------- 数据读取 --------------------
pos_df = pd.read_csv(POS_CSV)
neg_df = pd.read_csv(NEG_CSV)
n_pos = len(pos_df)
neg_df = neg_df.sample(n=n_pos, random_state=SEED).reset_index(drop=True)

data = pd.concat([pos_df, neg_df], axis=0).reset_index(drop=True)
label_col = "label"
seq_col = "seq_count"
feature_cols = [c for c in data.columns if c not in [seq_col, label_col]]

# 序列转 One-hot
print("Encoding sequences...")
seq_features = np.array([one_hot_encode(seq) for seq in tqdm(data[seq_col])])

# 数值特征
num_features = data[feature_cols].fillna(0).to_numpy(dtype=np.float32)

# 拼接
X = np.concatenate([seq_features, num_features], axis=1)
y = data[label_col].astype(int).to_numpy()

# 划分训练/测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE, random_state=SEED, stratify=y)

# 标准化
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# -------------------- 训练与评估 --------------------
def evaluate_model(name, clf):
    clf.fit(X_train, y_train)

    # 预测概率
    if hasattr(clf, "predict_proba"):
        probs = clf.predict_proba(X_test)[:, 1]
    else:
        probs = clf.predict(X_test)  # 对于 PLSRegression 等模型
    probs = probs.ravel()  # ✅ 展平成1维
    preds = (probs >= 0.5).astype(int)

    # 保存模型预测概率 CSV
    df_prob = pd.DataFrame({
        "Model": name,
        "y_true": y_test.ravel(),
        "y_prob": probs
    })
    df_prob.to_csv(os.path.join(OUT_DIR, f"{name}_run8_probs.csv"), index=False)

    # 计算指标
    metrics = {
        "ACC": accuracy_score(y_test, preds),
        "AUC": roc_auc_score(y_test, probs) if len(np.unique(y_test)) > 1 else np.nan,
        "F1": f1_score(y_test, preds),
        "Precision": precision_score(y_test, preds),
        "Recall": recall_score(y_test, preds),
        "MCC": matthews_corrcoef(y_test, preds),
        "N_test": len(y_test)
    }
    out_dir = os.path.join(OUT_DIR, name)
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame([metrics]).to_csv(os.path.join(out_dir, f"{name}_metrics.csv"), index=False)

    print(f"[{name}] AUC={metrics['AUC']:.4f} ACC={metrics['ACC']:.4f} F1={metrics['F1']:.4f}")
    return metrics


# -------------------- 主程序 --------------------
results = []
for key, build_fn in models.items():
    name, clf = build_fn()
    metrics = evaluate_model(name, clf)
    results.append({"Model": name, **metrics})

# 汇总所有模型预测概率到一个 CSV
all_probs = []
for key, build_fn in models.items():
    name, _ = build_fn()
    tmp = pd.read_csv(os.path.join(OUT_DIR, f"{name}_run8_probs.csv"))
    all_probs.append(tmp)

df_all_probs = pd.concat(all_probs, axis=0).reset_index(drop=True)
df_all_probs.to_csv(os.path.join(OUT_DIR, "run8_predictions.csv"), index=False)
print("✅ All model predictions saved to run8_predictions.csv")

# 汇总
pd.DataFrame(results).to_csv(os.path.join(OUT_DIR, "summary_results.csv"), index=False)
print("✅ All models done! Summary saved to:", os.path.join(OUT_DIR, "summary_results.csv"))
