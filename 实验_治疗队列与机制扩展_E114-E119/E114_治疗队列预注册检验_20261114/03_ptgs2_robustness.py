#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E114/03 · PTGS2 在 GSE132112 的稳健性（非预注册主检验；描述性补充，AI 初算）：
逐探针集结果、均值/中位数两种汇总、每 1 SD 的 logistic 比值比（非应答为事件）、6 个 hub 的并列结果。"""
import gzip, sys, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, norm
sys.argv = [sys.argv[0]] + sys.argv[1:]
GEO = Path(sys.argv[1]); OUT = Path(sys.argv[2])
sm = gzip.open(GEO / "GSE132112_series_matrix.txt.gz", "rt", errors="ignore").read().split("\n")
src = [l for l in sm if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
lab = {}
for s in src:
    m = re.match(r"(LA-NPC_\d+)-(.+)", s.strip('"')); lab[m.group(1)] = 0 if m.group(2) in ("Complete Response", "Partial Response") else 1  # 1 = 非应答(SD)
# 找 PTGS2 的 PSR 探针集
ids = []
with open(GEO / "GPL17585_data.txt", errors="ignore") as fh:
    for l in fh:
        if l.startswith("ID\t"): cols = l.rstrip("\n").split("\t"); iS = cols.index("gene_symbols"); continue
        if l.startswith(("#", "!", "^")): continue
        f = l.rstrip("\n").split("\t")
        if len(f) > iS and f[0].startswith("PSR") and "PTGS2" in re.split(r"\s*///\s*|;|,", f[iS]): ids.append(f[0])
rows = {}
with gzip.open(GEO / "GSE132112_nor_data_probelevel.txt.gz", "rt", errors="ignore") as fh:
    h = fh.readline().rstrip("\n").split("\t"); S = h[1:]
    for l in fh:
        k = l.split("\t", 1)[0]
        if k in ids: rows[k] = np.array(l.rstrip("\n").split("\t")[1:], float)
y = np.array([lab[s] for s in S])
out = []
for k, v in sorted(rows.items()):
    u, p = mannwhitneyu(v[y == 1], v[y == 0], alternative="two-sided"); out.append(dict(probe=k, AUC_nonresponse=u / ((y == 1).sum() * (y == 0).sum()), P=p, mean_nonresp=v[y == 1].mean(), mean_resp=v[y == 0].mean()))
pr = pd.DataFrame(out); pr.to_csv(OUT / "PTGS2_probe_level.csv", index=False)
print(pr.round(4).to_string())
M = np.vstack(list(rows.values()))
for name, g in [("mean", M.mean(0)), ("median", np.median(M, 0))]:
    u, p = mannwhitneyu(g[y == 1], g[y == 0], alternative="two-sided"); print(name, "AUC(非应答)=%.3f P=%.2g" % (u / ((y == 1).sum() * (y == 0).sum()), p))
g = M.mean(0); z = (g - g.mean()) / g.std(ddof=1)
X = np.column_stack([np.ones_like(z), z]); b = np.zeros(2)
for _ in range(50):
    pr_ = 1 / (1 + np.exp(-X @ b)); W = pr_ * (1 - pr_); H = X.T @ (X * W[:, None]); b = b + np.linalg.solve(H, X.T @ (y - pr_))
se = np.sqrt(np.diag(np.linalg.inv(H)))
print("每 +1 SD PTGS2：非应答 OR = %.2f (95%%CI %.2f–%.2f), Wald P=%.2g" % (np.exp(b[1]), np.exp(b[1] - 1.96 * se[1]), np.exp(b[1] + 1.96 * se[1]), 2 * (1 - norm.cdf(abs(b[1] / se[1])))))
print("n probes", len(rows), "n nonresp", y.sum(), "n resp", (1 - y).sum())
