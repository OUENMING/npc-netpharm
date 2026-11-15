"""E127/01。运行：python -I 01_sc_tcell.py <GSE150430_hg19.txt.gz> <输出目录>。预注册见 预注册_先于看数据.md"""
import gzip, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata, wilcoxon, spearmanr
SRC, OUT = sys.argv[1], Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
KEY = ["GSTK1", "PDE4D", "IFNG", "GZMB"]
rng = np.random.default_rng(20261116)
f = gzip.open(SRC, "rt"); cells = f.readline().rstrip("\n").split("\t")
ct = np.array([c.strip('"') for c in f.readline().rstrip("\n").split("\t")[1:]]); pt = np.array([c.split("_")[0] for c in cells])
if len(pt) != len(ct): pt = pt[-len(ct):]
istum = np.array([p.startswith("P") for p in pt]); isT = np.array([c.startswith("CD4+") or c.startswith("CD8+") for c in ct]) & istum
idx = np.where(isT)[0]; ngene = np.zeros(len(ct), dtype=np.int32); rows = {}; pool = {}
for l in f:
    g, rest = l.split("\t", 1); v = np.array(rest.rstrip("\n").split("\t"), dtype=np.float32); ngene += (v > 0)
    if g in KEY: rows[g] = v[idx]
    elif rng.random() < 300 / 22000 and g not in KEY: pool[g] = v[idx]
print("T 细胞", len(idx), "患者", len(set(pt[idx])), "关键基因", sorted(rows), "随机池", len(pool))
pT = pt[idx]; ng = ngene[idx].astype(float)
def resid(y, x):
    rx = rankdata(x); ry = rankdata(y); b = np.polyfit(rx, ry, 1); return ry - np.polyval(b, rx)
def pat_rho(a, b, partial):
    out = {}
    for p in sorted(set(pT)):
        m = pT == p
        if m.sum() < 50: continue
        if partial:
            r = spearmanr(resid(a[m], ng[m]), resid(b[m], ng[m])).statistic
        else: r = spearmanr(a[m], b[m]).statistic
        out[p] = r
    return pd.Series(out)
res = []
for g in ["PDE4D", "IFNG", "GZMB"]:
    s_part = pat_rho(rows["GSTK1"], rows[g], True).dropna(); s_raw = pat_rho(rows["GSTK1"], rows[g], False).dropna()
    P = wilcoxon(s_part).pvalue
    res.append(dict(gene=g, n_pat=len(s_part), median_rho_partial=s_part.median(), n_pos_partial=int((s_part > 0).sum()), P_wilcoxon_partial=P,
                    median_rho_raw=s_raw.median(), n_pos_raw=int((s_raw > 0).sum()), P_wilcoxon_raw=wilcoxon(s_raw).pvalue))
R = pd.DataFrame(res)
prim = R[R.gene.isin(["PDE4D", "IFNG"])].P_wilcoxon_partial.values; o = np.argsort(prim); q = prim[o] * 2 / np.arange(1, 3); q = np.minimum.accumulate(q[::-1])[::-1]
qq = np.empty(2); qq[o] = np.minimum(q, 1); R.loc[R.gene.isin(["PDE4D", "IFNG"]), "q_BH_primary"] = qq
pd.set_option("display.width", 200); print(R.round(4).to_string(index=False)); R.to_csv(OUT / "main_results.csv", index=False)
# 随机基因对照
mean = {g: rows[g].mean() for g in ["PDE4D", "IFNG"]}; null = {}
for tgt in ["PDE4D", "IFNG"]:
    m0 = mean[tgt]; cand = [g for g, v in pool.items() if abs(v.mean() - m0) <= 0.25 * m0 and (v > 0).mean() > 0.02]
    vals = [pat_rho(rows["GSTK1"], pool[g], True).median() for g in cand]; obs = R.set_index("gene").loc[tgt, "median_rho_partial"]
    null[tgt] = (len(cand), obs, float(np.mean(np.array(vals) < obs)) if vals else np.nan, np.nanmedian(vals) if vals else np.nan)
    print(f"随机对照 {tgt}: 匹配基因 {len(cand)}，观察 median 偏 rho {obs:.4f}，随机基因中位 {null[tgt][3]:.4f}，分位 {null[tgt][2]:.3f}")
pd.DataFrame(null, index=["n_matched", "obs_rho", "percentile", "null_median"]).T.to_csv(OUT / "random_control.csv")
# 患者间平均水平
pm = pd.DataFrame({g: pd.Series({p: rows[g][pT == p].mean() for p in set(pT) if (pT == p).sum() >= 50}) for g in KEY})
for g in ["PDE4D", "IFNG", "GZMB"]:
    r = spearmanr(pm.GSTK1, pm[g]); print(f"患者平均水平 GSTK1–{g}: rho {r.statistic:.3f} P {r.pvalue:.3g} n={len(pm)}")
