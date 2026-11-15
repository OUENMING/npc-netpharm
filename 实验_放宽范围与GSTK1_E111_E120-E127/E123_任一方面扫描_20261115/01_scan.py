"""E123/01 · 134 核心靶点「任一方面」扫描。预注册见 预注册_先于看数据.md。
运行：python -I 01_scan.py <项目根目录> <输出目录>"""
import sys, re, gzip
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, kruskal, spearmanr
import statsmodels.api as sm
R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True); X = R/"experiments"
core = sorted(open(X/"E81_方案A全套重出_20261101/out/core_A134.txt").read().split())
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o]*n/(np.arange(n)+1); q = np.minimum.accumulate(q[::-1])[::-1]; r = np.empty(n); r[o] = np.minimum(q, 1); return r
fam = {}   # 族 -> DataFrame(gene, stat, P, dir)
# A、B：E114 现成
e114 = X/"E114_治疗队列预注册检验_20261114/out"
a = pd.read_csv(next(e114.glob("GSE132112_IC应答*core134_tests.csv")), encoding="utf-8-sig"); fam["A_诱导化疗应答"] = pd.DataFrame({"gene": a.gene, "stat": 1-a.AUC_pos, "P": a.P, "note": "非应答AUC"})
b = pd.read_csv(next(e114.glob("GSE103611_*core134_tests.csv")), encoding="utf-8-sig"); fam["B_远处转移"] = pd.DataFrame({"gene": b.gene, "stat": 1-b.AUC_pos, "P": b.P, "note": "转移AUC"})
# H：E82
h = pd.read_csv(X/"E82_补充分析_投稿前_20261101/out/02_meta_134.csv", encoding="utf-8-sig"); fam["H_肿瘤vs正常"] = pd.DataFrame({"gene": h.gene, "stat": h.meta_g, "P": h.p, "note": "Hedges g(T-N)"})
# GSE102349
GEO = X/"E4_核实思路_20261002/geo_raw/GSE102349_series_matrix.txt.gz"
with gzip.open(GEO, "rt") as f: rows = [l.rstrip("\n").split("\t") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
pid = [re.sub(r"\D", "", x.strip('"')).lstrip("0") for x in rows[0][1:]]; ch = {}
for r in rows[1:]: ch[r[1].strip('"').split(":")[0]] = [x.strip('"').split(": ", 1)[-1] for x in r[1:]]
print("可用临床字段:", list(ch))
with gzip.open(X/"E21_患者锚定代谢物_20261003/geo_raw/GSE102349_NPC_mRNA_processed.txt.gz", "rt") as f:
    hdr = next(f).rstrip("\n").split("\t"); rowsd = {}
    for l in f:
        p = l.rstrip("\n").split("\t"); rowsd.setdefault(p[0], [float(v) for v in p[1:]])
cols = [re.search(r"T0*(\d+)$", c).group(1) for c in hdr[1:]]; pos = {c: i for i, c in enumerate(cols)}
ok = [(i, pos[p]) for i, p in enumerate(pid) if p in pos]; ci = np.array([i for i, _ in ok]); ei = np.array([j for _, j in ok])
print("匹配患者", len(ok), "/", len(pid))
genes = [g for g in core if g in rowsd and (np.array(rowsd[g]) > 0).mean() >= 0.2]; print("GSE102349 可用核心靶点", len(genes), "/ 134")
Xe = {g: np.log2(np.array(rowsd[g])[ei] + 1) for g in genes}
def col(k): return np.array([ch[k][i] for i in ci])
tte, ev = col("time to event"), col("event"); m = (tte != "N/A") & np.isin(ev, ["Disease progression", "Last follow-up"])
T = tte[m].astype(float); Ev = (ev[m] == "Disease progression").astype(int); print("Cox 样本", m.sum(), "事件", Ev.sum())
rc = []
for g in genes:
    x = Xe[g][m]; x = (x - x.mean())/(x.std()+1e-9)
    try:
        r = sm.PHReg(T, x[:, None], status=Ev, ties="breslow").fit(disp=0); rc.append((g, float(r.params[0]), float(r.pvalues[0])))
    except Exception: pass
fam["C_疾病进展(Cox)"] = pd.DataFrame(rc, columns=["gene", "stat", "P"]).assign(note="log HR/SD")
st = col("clinical stage"); ms = np.isin(st, ["I", "II", "III", "IV"]); adv = np.isin(st, ["III", "IV"])
print("分期", pd.Series(st).value_counts().to_dict())
rd = []
for g in genes:
    x = Xe[g]; u = mannwhitneyu(x[ms & adv], x[ms & ~adv], alternative="two-sided"); rd.append((g, u.statistic/((ms & adv).sum()*(ms & ~adv).sum()), u.pvalue))
fam["D_临床分期(III-IV vs I-II)"] = pd.DataFrame(rd, columns=["gene", "stat", "P"]).assign(note="AUC(晚期)")
sub = col("expression-based subtype"); print("亚型", pd.Series(sub).value_counts().to_dict()); ks = [s for s in sorted(set(sub)) if s not in ("N/A", "")]
re_ = []
for g in genes:
    x = Xe[g]; re_.append((g, *kruskal(*[x[sub == s] for s in ks])))
fam["E_表达亚型"] = pd.DataFrame(re_, columns=["gene", "stat", "P"]).assign(note="KW H")
def num(k): 
    v = col(k); return np.array([float(x) if re.fullmatch(r"[\d.]+", x) else np.nan for x in v])
til = {"stromal": num("% stromal tils"), "intra": num("% intratumoral tils")}; mb = num("non-synonymous mutation burden")
print("TIL 非缺失", {k: int((~np.isnan(v)).sum()) for k, v in til.items()}, "；突变负荷非缺失", int((~np.isnan(mb)).sum()))
rf, rg = [], []
for g in genes:
    x = Xe[g]; best = None
    for k, v in til.items():
        mm = ~np.isnan(v); r, p = spearmanr(x[mm], v[mm]); 
        if best is None or p < best[2]: best = (k, r, p)
    rf.append((g, best[1], min(1, best[2]*2), best[0]))
    mm = ~np.isnan(mb); r, p = spearmanr(x[mm], mb[mm]); rg.append((g, r, p))
fam["F_免疫浸润(TIL%)"] = pd.DataFrame(rf, columns=["gene", "stat", "P", "note"]); fam["G_突变负荷"] = pd.DataFrame(rg, columns=["gene", "stat", "P"]).assign(note="Spearman rho")
allr = []
for k, d in fam.items():
    d = d.dropna(subset=["P"]).copy(); d["family"] = k; d["q_family"] = bh(d.P.values); allr.append(d)
D = pd.concat(allr, ignore_index=True); D["q_global"] = bh(D.P.values); D.to_csv(OUT/"scan_long.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250, "display.max_rows", 200)
S = D.groupby("family").agg(n_genes=("gene", "size"), nominal_P05=("P", lambda s: int((s < 0.05).sum())), frac_P05=("P", lambda s: round((s < 0.05).mean(), 3)), q_family_lt05=("q_family", lambda s: int((s < 0.05).sum())), q_global_lt05=("q_global", lambda s: int((s < 0.05).sum()))); print(S.to_string()); S.to_csv(OUT/"scan_summary_by_family.csv", encoding="utf-8-sig")
hits = D[D.q_family < 0.05].sort_values(["family", "P"]); print("\n宽松命中（族内 q<0.05）:"); print(hits[["family", "gene", "stat", "note", "P", "q_family", "q_global"]].round(5).to_string(index=False))
hits.to_csv(OUT/"hits_family_q05.csv", index=False, encoding="utf-8-sig")
strict = D[D.q_global < 0.05].sort_values("P"); print("\n严格命中（全局 q<0.05）:"); print(strict[["family", "gene", "stat", "P", "q_global"]].round(6).to_string(index=False))
g = D[D.q_family < 0.05].groupby("gene").family.apply(lambda s: "; ".join(sorted(set(s)))); print("\n每个基因命中的族:"); print(g.to_string()); g.to_csv(OUT/"gene_hit_families.csv", encoding="utf-8-sig")
outcome = D[D.family.str[0].isin(list("ABC"))].sort_values("P").head(12); print("\n结局类（A–C）最小 12 个 P:"); print(outcome[["family", "gene", "stat", "P", "q_family", "q_global"]].round(5).to_string(index=False))
