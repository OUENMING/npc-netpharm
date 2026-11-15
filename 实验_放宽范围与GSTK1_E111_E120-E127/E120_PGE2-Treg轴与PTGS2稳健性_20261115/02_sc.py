"""E120/02 · 队列 B（GSE150430 单细胞）：复现 PNAS 的 PTGS2 分组 → Treg 的 TI-Treg 特征比较。预注册见 预注册_先于看数据.md 第 3 部分。
运行：<python> -I 02_sc.py <GSE150430_hg19.txt.gz> <输出目录>
"""
import gzip, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, rankdata, t as tdist

SRC, OUT = sys.argv[1], Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
TIT = ["TNFRSF9", "TNFRSF18", "TNFRSF4", "CTLA4", "CCR8"]; WANT = set(TIT + ["PTGS2", "FOXP3"])
f = gzip.open(SRC, "rt"); cells = f.readline().rstrip("\n").split("\t")
ct = [c.strip('"') for c in f.readline().rstrip("\n").split("\t")[1:]]; assert len(ct) == len(cells)
pt = np.array([c.split("_")[0] for c in cells]); ct = np.array(ct)
rows = {}
for l in f:
    g, rest = l.split("\t", 1)
    if g in WANT: rows[g] = np.array(rest.rstrip("\n").split("\t"), dtype=float)
print("cells", len(cells), "found", sorted(rows), "missing", sorted(WANT - set(rows)))
tum = np.array([p.startswith("P") for p in pt])
MYE = {"Macrophage", "Monocyte", "pDC", "DC1 (CLEC9A+)", "DC2 (CD1C+)", "DC3 (CCR7+)"}
TREG = {"CD4+ Tregs-1", "CD4+ Tregs-2"}
recs = []
for p in sorted(set(pt[tum])):
    m = (pt == p); tr = m & np.isin(ct, list(TREG))
    if tr.sum() < 20: continue
    rec = dict(patient=p, n_cells=int(m.sum()), n_treg=int(tr.sum()), PTGS2_all=rows["PTGS2"][m].mean(),
               myeloid_frac=np.isin(ct[m], list(MYE)).mean(),
               TIT_sc=np.mean([rows[g][tr].mean() for g in TIT]), FOXP3_treg=rows["FOXP3"][tr].mean())
    for g in TIT: rec[g + "_treg"] = rows[g][tr].mean()
    recs.append(rec)
D = pd.DataFrame(recs); D.to_csv(OUT / "sc_per_patient.csv", index=False)
print("\n纳入患者（≥20 个 Treg 细胞）:", len(D)); print(D[["patient", "n_cells", "n_treg", "PTGS2_all", "myeloid_frac", "TIT_sc"]].round(3).to_string(index=False))
n = len(D); k = int(round(n / 3)); o = D.sort_values("PTGS2_all"); lo, hi = o.head(k), o.tail(k)
u = mannwhitneyu(hi.TIT_sc, lo.TIT_sc, alternative="two-sided")
print(f"\nPTGS2 最高 {k} 例 vs 最低 {k} 例 的 Treg TIT_sc：中位 {hi.TIT_sc.median():.3f} vs {lo.TIT_sc.median():.3f}；Mann-Whitney P={u.pvalue:.4g}（U={u.statistic:.0f}）")
r, p = spearmanr(D.PTGS2_all, D.TIT_sc); print(f"全体 Spearman(PTGS2, TIT_sc): rho={r:.3f} P={p:.4g} (n={n})")
def pcor(x, yv, cov):
    R = rankdata; C = np.column_stack([np.ones(len(x)), R(cov)])
    res = lambda v: R(v) - C @ np.linalg.lstsq(C, R(v), rcond=None)[0]
    rx, ry = res(x), res(yv); rr = np.corrcoef(rx, ry)[0, 1]; df = len(x) - 3; tt = rr * np.sqrt(df / (1 - rr * rr)); return rr, 2 * tdist.sf(abs(tt), df)
rr, pp = pcor(D.PTGS2_all.values, D.TIT_sc.values, D.myeloid_frac.values); print(f"偏相关（控制髓系细胞占比）: rho={rr:.3f} P={pp:.4g}")
print("参考：PTGS2 与髓系占比 Spearman rho=%.3f" % spearmanr(D.PTGS2_all, D.myeloid_frac)[0])
print("各 TIT 基因（Treg 内均值）与 PTGS2 的 Spearman:", {g: round(spearmanr(D.PTGS2_all, D[g + "_treg"])[0], 2) for g in TIT})
