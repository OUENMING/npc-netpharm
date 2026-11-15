"""E120/01 · 队列 A（GSE132112）：PGE2–Treg 轴检验 + PTGS2 稳健性。预注册见 预注册_先于看数据.md。
运行：<python> -I 01_bulk.py <geo目录> <输出目录>
"""
import gzip, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, rankdata, norm, t as tdist
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

GEO = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
rng = np.random.default_rng(20261116)

TIT = ["TNFRSF9", "TNFRSF18", "TNFRSF4", "CTLA4", "CCR8"]
TREG_ID = ["FOXP3", "IL2RA"]
CD3 = ["CD3D", "CD3E", "CD2"]
LIN = {"髓系": "CD14 CD68 CSF1R LYZ CD163 MSR1 ITGAM FCGR3A VCAN S100A8 S100A9".split(),
       "上皮/肿瘤": "EPCAM KRT5 KRT19 KRT17 TP63".split(),
       "成纤维": "COL1A1 COL1A2 DCN LUM FAP".split()}
ALL = sorted(set(TIT + TREG_ID + CD3 + ["PTGS2"] + [g for v in LIN.values() for g in v]))

# ---- 载入（同 E119/02）
p2g = {}
with open(GEO / "GPL17585_data.txt", errors="ignore") as fh:
    cols = None
    for l in fh:
        if l.startswith("ID\t"):
            cols = l.rstrip("\n").split("\t"); iS = cols.index("gene_symbols"); continue
        if cols is None or l.startswith(("#", "!", "^")): continue
        f = l.rstrip("\n").split("\t")
        if len(f) <= iS or not f[0].startswith("PSR") or f[iS] == "---": continue
        s = {x.strip() for x in re.split(r"\s*///\s*|;|,", f[iS]) if x.strip()} & set(ALL)
        if s: p2g[f[0]] = s
sm_ = gzip.open(GEO / "GSE132112_series_matrix.txt.gz", "rt", errors="ignore").read().split("\n")
src = [l for l in sm_ if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
lab = {}
for s in src:
    m = re.match(r"(LA-NPC_\d+)-(.+)", s.strip('"')); lab[m.group(1)] = m.group(2)
sums, cnt = {}, {}
with gzip.open(GEO / "GSE132112_nor_data_probelevel.txt.gz", "rt", errors="ignore") as fh:
    h = fh.readline().rstrip("\n").split("\t")[1:]
    for l in fh:
        k, rest = l.split("\t", 1)
        if k in p2g:
            v = np.array(rest.rstrip("\n").split("\t"), dtype=float)
            for g in p2g[k]: sums[g] = sums.get(g, 0) + v; cnt[g] = cnt.get(g, 0) + 1
expr = pd.DataFrame({g: sums[g] / cnt[g] for g in sums}, index=h).T
y = pd.Series([0 if lab[s] in ("Complete Response", "Partial Response") else 1 for s in h], index=h)
print("n 非应答/应答", int(y.sum()), int((1 - y).sum()), "| 未测到:", [g for g in ALL if g not in expr.index])
z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1), axis=0)
def score(gs): return z.loc[[g for g in gs if g in z.index]].mean(axis=0)
S = pd.DataFrame({"PTGS2": z.loc["PTGS2"], "TIT": score(TIT), "CD3": score(CD3), "FOXP3": z.loc["FOXP3"]})
for k, g in LIN.items(): S[k] = score(g)

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; r = np.empty(n); r[o] = np.minimum(q, 1); return r
def auc(x):
    a = x[y == 1].values; b = x[y == 0].values; u = mannwhitneyu(a, b, alternative="two-sided")
    return u.statistic / (len(a) * len(b)), u.pvalue
def auc_ci(x, B=2000):
    a = x[y == 1].values; b = x[y == 0].values; r = []
    for _ in range(B):
        aa = rng.choice(a, len(a)); bb = rng.choice(b, len(b))
        r.append((aa[:, None] > bb[None, :]).mean() + 0.5 * (aa[:, None] == bb[None, :]).mean())
    return np.percentile(r, [2.5, 97.5])
def pcor_spearman(x, yv, covs):
    """偏相关：先对秩做线性回归取残差，再求相关（P 用 n-2-k 自由度 t 检验）。"""
    R = lambda v: rankdata(v)
    C = np.column_stack([np.ones(len(x))] + [R(c) for c in covs])
    res = lambda v: R(v) - C @ np.linalg.lstsq(C, R(v), rcond=None)[0]
    rx, ry = res(x), res(yv); r = np.corrcoef(rx, ry)[0, 1]
    dfree = len(x) - 2 - len(covs); tt = r * np.sqrt(dfree / (1 - r * r)); return r, 2 * tdist.sf(abs(tt), dfree)

# ================= 第 1 部分
print("\n=== 第 1 部分 ===")
r1, p1 = spearmanr(S.PTGS2, S.TIT); print(f"P1 Spearman(PTGS2,TIT): rho={r1:.3f} P={p1:.3g}")
pr, pp = pcor_spearman(S.PTGS2.values, S.TIT.values, [S.CD3.values, S.FOXP3.values])
print(f"P1 偏相关(控制 CD3、FOXP3): rho={pr:.3f} P={pp:.3g}")
A2, p2 = auc(S.TIT); lo, hi = auc_ci(S.TIT)
print(f"P2 TIT→非应答 AUC={A2:.3f} (95%CI {lo:.3f}-{hi:.3f}) P={p2:.3g}")
q_main = bh([p1, p2]); print("主检验 BH(2): q(P1)=%.4f q(P2)=%.4f" % tuple(q_main))
rho_cd3 = spearmanr(S.PTGS2, S.CD3); rho_f = spearmanr(S.PTGS2, S.FOXP3)
print(f"参考：PTGS2 与 CD3 rho={rho_cd3[0]:.3f} P={rho_cd3[1]:.3g}；与 FOXP3 rho={rho_f[0]:.3f} P={rho_f[1]:.3g}；TIT 与 CD3 rho={spearmanr(S.TIT,S.CD3)[0]:.3f}")
print(f"参考：CD3 对非应答 AUC={auc(S.CD3)[0]:.3f} P={auc(S.CD3)[1]:.3g}；FOXP3 AUC={auc(S.FOXP3)[0]:.3f} P={auc(S.FOXP3)[1]:.3g}")

sec = []
for g in TIT + TREG_ID:
    if g not in z.index: sec.append((g, np.nan, np.nan)); continue
    a, p = auc(z.loc[g]); sec.append((g, a, p))
SEC = pd.DataFrame(sec, columns=["gene", "AUC_nonresp", "P"]); SEC["q_BH7"] = bh(SEC.P.fillna(1).values)
print("\n次要（7 基因，BH）:"); print(SEC.round(4).to_string(index=False))
SEC.to_csv(OUT / "secondary_genes.csv", index=False)

def logit_or(cols):
    X = sm.add_constant(S[cols]); m = sm.Logit(y.astype(float), X).fit(disp=0)
    ci = m.conf_int().loc["PTGS2"]; return np.exp(m.params["PTGS2"]), np.exp(ci[0]), np.exp(ci[1]), m.pvalues["PTGS2"]
print("\nP3 PTGS2 OR/SD（普通逻辑回归）:")
rows3 = []
for name, cols in (("PTGS2", ["PTGS2"]), ("+TIT", ["PTGS2", "TIT"]), ("+CD3", ["PTGS2", "CD3"]), ("+TIT+CD3", ["PTGS2", "TIT", "CD3"])):
    o, l, u, p = logit_or(cols); rows3.append((name, o, l, u, p)); print(f"  {name}: OR={o:.2f} ({l:.2f}-{u:.2f}) P={p:.4g}")
pd.DataFrame(rows3, columns=["model", "OR", "lo", "hi", "P"]).to_csv(OUT / "p3_or_with_TIT.csv", index=False)

# ================= 第 2 部分：Firth
def firth(X, yv, maxit=200, tol=1e-8):
    X = np.asarray(X, float); yv = np.asarray(yv, float); b = np.zeros(X.shape[1])
    for _ in range(maxit):
        eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta)); W = p * (1 - p)
        I = X.T @ (X * W[:, None]); Iinv = np.linalg.inv(I)
        hdiag = W * np.einsum("ij,jk,ik->i", X, Iinv, X)
        U = X.T @ (yv - p + hdiag * (0.5 - p)); step = Iinv @ U
        mx = np.abs(step).max()
        if mx > 5: step *= 5 / mx
        b_new = b + step
        if np.abs(b_new - b).max() < tol: b = b_new; break
        b = b_new
    eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta)); W = p * (1 - p)
    cov = np.linalg.inv(X.T @ (X * W[:, None])); return b, np.sqrt(np.diag(cov))
def firth_or(cols, idx=None):
    D = S if idx is None else S.loc[idx]; yy = y if idx is None else y.loc[idx]
    X = np.column_stack([np.ones(len(D))] + [D[c].values for c in cols]); b, se = firth(X, yy.values)
    k = 1 + cols.index("PTGS2"); return np.exp(b[k]), np.exp(b[k] - 1.96 * se[k]), np.exp(b[k] + 1.96 * se[k]), 2 * norm.sf(abs(b[k] / se[k]))
print("\n=== 第 2 部分：Firth ===")
M = {"M1 PTGS2": ["PTGS2"], "M2 PTGS2+三谱系": ["PTGS2", "髓系", "上皮/肿瘤", "成纤维"]}
fres = {}
for k, c in M.items():
    o, l, u, p = firth_or(c); fres[k] = (o, l, u, p); print(f"  Firth {k}: OR/SD={o:.2f} ({l:.2f}-{u:.2f}) P={p:.4g}")
    # 与普通回归对照
    o2, l2, u2, p2_ = logit_or(c); print(f"     普通: OR/SD={o2:.2f} ({l2:.2f}-{u2:.2f}) P={p2_:.4g}")
# 分层 bootstrap（M2）
idx1 = y.index[y == 1].values; idx0 = y.index[y == 0].values; bs = []
for _ in range(2000):
    ii = np.concatenate([rng.choice(idx1, len(idx1)), rng.choice(idx0, len(idx0))])
    D = S.loc[ii]; yy = y.loc[ii].values
    X = np.column_stack([np.ones(len(D))] + [D[c].values for c in M["M2 PTGS2+三谱系"]])
    try: b, _ = firth(X, yy); bs.append(np.exp(b[1]))
    except Exception: pass
bs = np.array(bs); print(f"  M2 分层 bootstrap（Firth，n={len(bs)}）: OR 中位 {np.median(bs):.2f}，95% 百分位 CI {np.percentile(bs,2.5):.2f}-{np.percentile(bs,97.5):.2f}；OR<1 的比例 {np.mean(bs<1):.4f}")
# VIF
X = sm.add_constant(S[M["M2 PTGS2+三谱系"]]); print("  VIF:", {c: round(variance_inflation_factor(X.values, i), 2) for i, c in enumerate(X.columns) if c != "const"})
# 留一
loo = []
for s in S.index:
    ii = [t for t in S.index if t != s]; loo.append((s, firth_or(M["M2 PTGS2+三谱系"], ii)[0]))
L = pd.DataFrame(loo, columns=["dropped", "OR"]).sort_values("OR"); L.to_csv(OUT / "loo_M2.csv", index=False)
print(f"  留一（M2，Firth）OR 范围 {L.OR.min():.2f}-{L.OR.max():.2f}；最小值对应样本 {L.iloc[0].dropped}")
o, l, u, p = fres["M2 PTGS2+三谱系"]
verdict = "稳健" if (l > 1 and L.OR.min() > 1.5) else "不稳健"
print(f"  预设判据：Firth M2 下界 {l:.2f}>1 且 留一最小 OR {L.OR.min():.2f}>1.5 → {verdict}")
S.assign(y=y).to_csv(OUT / "scores_bulk.csv")
