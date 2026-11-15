"""E125/01 · GSTK1 进展信号：免疫浸润校正。预注册见 预注册_先于看数据.md。
运行：python -I 01_adjust.py <项目根目录> <输出目录>"""
import sys, re, gzip
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import spearmanr, rankdata
R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); X = R/"experiments"; OUT.mkdir(exist_ok=True, parents=True); rng = np.random.default_rng(20261117)
with gzip.open(X/"E4_核实思路_20261002/geo_raw/GSE102349_series_matrix.txt.gz", "rt") as f: rows = [l.rstrip("\n").split("\t") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
pid = [re.sub(r"\D", "", x.strip('"')).lstrip("0") for x in rows[0][1:]]; ch = {}
for r in rows[1:]: ch[r[1].strip('"').split(":")[0]] = [x.strip('"').split(": ", 1)[-1] for x in r[1:]]
with gzip.open(X/"E21_患者锚定代谢物_20261003/geo_raw/GSE102349_NPC_mRNA_processed.txt.gz", "rt") as f:
    hdr = next(f).rstrip("\n").split("\t"); rd = {}
    for l in f: p = l.rstrip("\n").split("\t"); rd.setdefault(p[0], np.array([float(v) for v in p[1:]]))
cols = [re.search(r"T0*(\d+)$", c).group(1) for c in hdr[1:]]; pos = {c: i for i, c in enumerate(cols)}
ok = [(i, pos[p]) for i, p in enumerate(pid) if p in pos]; ci = np.array([i for i, _ in ok]); ei = np.array([j for _, j in ok]); colnames = [hdr[1:][j] for _, j in ok]
col = lambda k: np.array([ch[k][i] for i in ci]); tte, ev, sub = col("time to event"), col("event"), col("expression-based subtype")
m = (tte != "N/A") & np.isin(ev, ["Disease progression", "Last follow-up"]); T = tte[m].astype(float); Ev = (ev[m] == "Disease progression").astype(int); S = sub[m]; n = len(T)
lg = lambda g: np.log2(rd[g][ei] + 1); zs = lambda v: (v - v.mean())/(v.std() + 1e-9)
# 基因集成员检查
def gmt(p): return {l.split("\t")[0]: set(l.rstrip("\n").split("\t")[2:]) for l in open(p)}
ch_sets = gmt(X/"E82_补充分析_投稿前_20261101/raw/Charoentong2017.gmt"); est_sets = gmt(X/"E82_补充分析_投稿前_20261101/raw/ESTIMATE_SI_geneset_from_tidyestimate.gmt")
print("GSTK1 在 Charoentong 基因集:", [k for k, v in ch_sets.items() if "GSTK1" in v], "；在 ESTIMATE 基因集:", [k for k, v in est_sets.items() if "GSTK1" in v])
est = pd.read_csv(X/"E82_补充分析_投稿前_20261101/out/03_estimate_GSE102349.csv", encoding="utf-8-sig", index_col=0); ss = pd.read_csv(X/"E82_补充分析_投稿前_20261101/out/03_ssgsea_GSE102349.csv", encoding="utf-8-sig", index_col=0)
est = est.loc[colnames]; ss = ss.T.loc[colnames] if ss.shape[0] < ss.shape[1] else ss.loc[colnames]
tcell_cols = [c for c in ss.columns if "T cell" in c or c in ("Regulatory T cell", "Type 1 T helper cell", "Type 17 T helper cell", "Type 2 T helper cell", "T follicular helper cell")]
print("Tcell ssGSEA 列:", tcell_cols)
V = {"GSTK1": zs(lg("GSTK1")), "CD3": zs(np.mean([zs(lg(g)) for g in ["CD3D", "CD3E", "CD2"]], 0)), "Immune": zs(est["ImmuneSignature"].values.astype(float)), "Stromal": zs(est["StromalSignature"].values.astype(float)), "Tcell": zs(ss[tcell_cols].mean(axis=1).values.astype(float))}
D2 = np.column_stack([(sub == s).astype(float) for s in ["II", "III"]])
def num(k): v = col(k); return np.array([float(x) if re.fullmatch(r"[\d.]+", x) else np.nan for x in v])
V["TIL_intra"], V["TIL_stroma"] = num("% intratumoral tils"), num("% stromal tils")
def cox(cols, mask=None, adjsub=True):
    mask = np.ones(len(sub), bool) if mask is None else mask; mm = mask & m_full
    Xm = np.column_stack([V[c][mm] if c not in () else None for c in cols] + ([D2[mm]] if adjsub else []))
    r = sm.PHReg(T_full[mm], Xm, status=E_full[mm], ties="breslow").fit(disp=0); return r, mm
m_full = np.zeros(len(sub), bool); idxm = np.where((tte != "N/A") & np.isin(ev, ["Disease progression", "Last follow-up"]))[0]; m_full[idxm] = True
T_full = np.where(m_full, np.nan_to_num(np.array([float(x) if re.fullmatch(r"[\d.]+", x) else np.nan for x in tte])), 0); E_full = (ev == "Disease progression").astype(int)
print("n", int(m_full.sum()), "事件", int(E_full[m_full].sum()))
print("\n[Q1] 免疫分数本身与进展（单变量 / 加亚型）")
for k in ["CD3", "Immune", "Stromal", "Tcell", "TIL_intra", "TIL_stroma", "GSTK1"]:
    mk = ~np.isnan(V[k]); a, _ = cox([k], mk, adjsub=False); b, _ = cox([k], mk, adjsub=True); print(f"  {k:10s} 单变量 logHR={a.params[0]:+.3f} P={a.pvalues[0]:.4f} | +亚型 logHR={b.params[0]:+.3f} P={b.pvalues[0]:.4f}")
print("\n[Q3] GSTK1 与浸润分数 Spearman（全 113 例 / 控制亚型）")
xa = zs(lg("GSTK1")); Va = {"CD3": zs(np.mean([zs(lg(g)) for g in ["CD3D", "CD3E", "CD2"]], 0)), "Immune": zs(est["ImmuneSignature"].values.astype(float)), "Stromal": zs(est["StromalSignature"].values.astype(float)), "Tcell": zs(ss[tcell_cols].mean(axis=1).values.astype(float))}
Dall = np.column_stack([np.ones(len(sub))] + [(sub == s).astype(float) for s in ["II", "III"]]); res = lambda v: rankdata(v) - Dall @ np.linalg.lstsq(Dall, rankdata(v), rcond=None)[0]
for k, v in Va.items():
    r1, p1 = spearmanr(xa, v); r2 = np.corrcoef(res(xa), res(v))[0, 1]; print(f"  GSTK1 vs {k:7s}: rho={r1:+.3f} P={p1:.3g} | 控制亚型偏相关 r={r2:+.3f}")
print("\n[Q2] 校正后 GSTK1")
base, _ = cox(["GSTK1"], None, adjsub=False); b0 = base.params[0]; print(f"  单变量 logHR={b0:+.3f} P={base.pvalues[0]:.4f}")
models = {"M0 GSTK1+亚型": ["GSTK1"], "M1 +CD3": ["GSTK1", "CD3"], "M1b +Tcell(ssGSEA)": ["GSTK1", "Tcell"], "M2 +Immune": ["GSTK1", "Immune"], "M3 +Stromal+Immune": ["GSTK1", "Stromal", "Immune"], "M4 +CD3+Stromal+Immune+TIL×2": ["GSTK1", "CD3", "Stromal", "Immune", "TIL_intra", "TIL_stroma"]}
out = []
for name, cs in models.items():
    mk = np.ones(len(sub), bool)
    for c in cs: mk &= ~np.isnan(V[c])
    r, mm = cox(cs, mk); b = r.params[0]; se = r.bse[0]; out.append((name, int(mm.sum()), int(E_full[mm].sum()), b, np.exp(b), np.exp(b-1.96*se), np.exp(b+1.96*se), r.pvalues[0], (b0-b)/abs(b0)*100))
    print(f"  {name:32s} n={mm.sum()} 事件={E_full[mm].sum()} logHR={b:+.3f} HR={np.exp(b):.2f}({np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f}) P={r.pvalues[0]:.4f} 相对单变量缩小 {(b0-b)/abs(b0)*100:.0f}%")
pd.DataFrame(out, columns=["model", "n", "events", "logHR", "HR", "lo", "hi", "P", "shrink_pct_vs_uni"]).to_csv(OUT/"gstk1_adjusted_models.csv", index=False, encoding="utf-8-sig")
print("\n[Q4] GSTK1 残差（浸润+亚型解释不了的部分）对进展")
Z = np.column_stack([np.ones(len(sub)), Va["CD3"], Va["Immune"], Va["Stromal"], (sub == "II").astype(float), (sub == "III").astype(float)]); resid = xa - Z @ np.linalg.lstsq(Z, xa, rcond=None)[0]; print(f"  R²(GSTK1 被浸润+亚型解释) = {1 - resid.var()/xa.var():.3f}"); V["GSTK1res"] = zs(resid)
r, _ = cox(["GSTK1res"], None, adjsub=True); print(f"  残差 +亚型：logHR={r.params[0]:+.3f} HR={np.exp(r.params[0]):.2f} P={r.pvalues[0]:.4f}")
print("\n[Q5] 校准：随机基因经同样 M3 校正后的保留率")
genes = [g for g, v in rd.items() if isinstance(g, str) and (v > 0).mean() >= 0.2]; mean_ex = np.array([np.log2(rd[g][ei] + 1).mean() for g in genes]); dec = pd.qcut(mean_ex, 10, labels=False, duplicates="drop"); gdec = dec[genes.index("GSTK1")]
cand = [g for g, d in zip(genes, dec) if abs(d - gdec) <= 1 and g != "GSTK1"]; rng.shuffle(cand); keep = []; tried = 0
for g in cand:
    V["g"] = zs(lg(g)); a, _ = cox(["g"], None, adjsub=False)
    if a.pvalues[0] < 0.05:
        b, _ = cox(["g", "Stromal", "Immune"], None, adjsub=True); keep.append((g, a.params[0], a.pvalues[0], b.params[0], b.pvalues[0]))
    tried += 1
    if len(keep) >= 300 or tried >= 3000: break
K = pd.DataFrame(keep, columns=["gene", "b_uni", "P_uni", "b_adj", "P_adj"]); K["shrink"] = (K.b_uni.abs() - K.b_adj.abs())/K.b_uni.abs()*100
print(f"  共 {len(K)} 个表达相近且单变量 P<0.05 的随机基因（试了 {tried} 个）：M3 校正后仍 P<0.05 的比例 {(K.P_adj<0.05).mean():.2f}；logHR 缩小比例中位 {K.shrink.median():.0f}%（IQR {K.shrink.quantile(.25):.0f}–{K.shrink.quantile(.75):.0f}%）")
r3, _ = cox(["GSTK1", "Stromal", "Immune"], None, adjsub=True); sh = (abs(b0) - abs(r3.params[0]))/abs(b0)*100; print(f"  GSTK1 的缩小比例 {sh:.0f}%，位于随机基因缩小比例分布的第 {(K.shrink < sh).mean()*100:.0f} 百分位；其校正后 P={r3.pvalues[0]:.4f}")
K.to_csv(OUT/"random_genes_M3.csv", index=False, encoding="utf-8-sig")
