"""E124/01 · GSTK1：全基因组名次、置换校准、多变量 Cox、通路一致性。预注册见 预注册_先于看数据.md。
运行：python -I 01_bulk.py <项目根目录> <输出目录>"""
import sys, re, gzip
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import norm, spearmanr, kruskal
R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); X = R/"experiments"; OUT.mkdir(exist_ok=True, parents=True); rng = np.random.default_rng(20261115)
with gzip.open(X/"E4_核实思路_20261002/geo_raw/GSE102349_series_matrix.txt.gz", "rt") as f: rows = [l.rstrip("\n").split("\t") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
pid = [re.sub(r"\D", "", x.strip('"')).lstrip("0") for x in rows[0][1:]]; ch = {}
for r in rows[1:]: ch[r[1].strip('"').split(":")[0]] = [x.strip('"').split(": ", 1)[-1] for x in r[1:]]
with gzip.open(X/"E21_患者锚定代谢物_20261003/geo_raw/GSE102349_NPC_mRNA_processed.txt.gz", "rt") as f:
    hdr = next(f).rstrip("\n").split("\t"); rd = {}
    for l in f: p = l.rstrip("\n").split("\t"); rd.setdefault(p[0], np.array([float(v) for v in p[1:]]))
cols = [re.search(r"T0*(\d+)$", c).group(1) for c in hdr[1:]]; pos = {c: i for i, c in enumerate(cols)}
ok = [(i, pos[p]) for i, p in enumerate(pid) if p in pos]; ci = np.array([i for i, _ in ok]); ei = np.array([j for _, j in ok])
col = lambda k: np.array([ch[k][i] for i in ci]); tte, ev, sub = col("time to event"), col("event"), col("expression-based subtype")
m = (tte != "N/A") & np.isin(ev, ["Disease progression", "Last follow-up"]); T = tte[m].astype(float); Ev = (ev[m] == "Disease progression").astype(int); S = sub[m]; n = len(T)
print("n", n, "事件", Ev.sum(), "亚型", dict(zip(*np.unique(S, return_counts=True))))
genes = [g for g, v in rd.items() if isinstance(g, str) and (v > 0).mean() >= 0.2]; M = np.array([np.log2(rd[g][ei][m] + 1) for g in genes]); M = (M - M.mean(1, keepdims=True))/(M.std(1, keepdims=True) + 1e-9); print("全基因组可用", len(genes))
def score_z(T, Ev, M):
    """Cox 得分检验 z（β=0，Breslow）：对每个基因向量化。risk sets：T_j >= T_i。"""
    num = np.zeros(M.shape[0]); var = np.zeros(M.shape[0])
    for i in np.where(Ev == 1)[0]:
        r = T >= T[i]; xr = M[:, r]; mu = xr.mean(1); num += M[:, i] - mu; var += ((xr - mu[:, None])**2).mean(1)
    return num/np.sqrt(var + 1e-12)
z0 = score_z(T, Ev, M); P0 = 2*norm.sf(np.abs(z0)); gi = genes.index("GSTK1"); order = np.argsort(P0); rank = int(np.where(order == gi)[0][0]) + 1
print(f"[P1] GSTK1 得分检验 z={z0[gi]:.3f} P={P0[gi]:.2e}；全基因组名次 {rank}/{len(genes)}；P≤{P0[gi]:.2e} 的基因数 {int((P0 <= P0[gi]).sum())}；P<0.05 的基因 {int((P0 < 0.05).sum())}（期望 ≈{0.05*len(genes):.0f}）")
print("全基因组前 10:", [(genes[i], round(float(z0[i]), 2), f"{P0[i]:.1e}") for i in order[:10]])
pd.DataFrame({"gene": genes, "z": z0, "P": P0}).sort_values("P").to_csv(OUT/"genomewide_cox_score.csv", index=False, encoding="utf-8-sig")
# P2 置换
B = 2000; minP = np.zeros(B); minPs = np.zeros(B); gst_unres = np.zeros(B); gst_strat = np.zeros(B)
idx_by = {s: np.where(S == s)[0] for s in np.unique(S)}
for b in range(B):
    p = rng.permutation(n); z = score_z(T[p], Ev[p], M); pv = 2*norm.sf(np.abs(z)); minP[b] = pv.min(); gst_unres[b] = pv[gi]
    q = np.arange(n)
    for s, ix in idx_by.items(): q[ix] = rng.permutation(ix)
    z2 = score_z(T[q], Ev[q], M); pv2 = 2*norm.sf(np.abs(z2)); minPs[b] = pv2.min(); gst_strat[b] = pv2[gi]
fw = (np.sum(minP <= P0[gi]) + 1)/(B + 1); fws = (np.sum(minPs <= P0[gi]) + 1)/(B + 1)
print(f"[P2] 全基因组最小-P 置换 FWER P={fw:.4f}（置换最小 P 的中位 {np.median(minP):.1e}，5% 分位 {np.percentile(minP,5):.1e}）")
print(f"[P2] 亚型内置换（保留亚型）最小-P FWER P={fws:.4f}（中位 {np.median(minPs):.1e}，5% 分位 {np.percentile(minPs,5):.1e}）")
print(f"[P2] 参考：GSTK1 自身在亚型内置换下 P 的分布：P≤{P0[gi]:.1e} 的比例 {np.mean(gst_strat <= P0[gi]):.4f}")
# P3 多变量
def num(k): v = col(k)[m]; return np.array([float(x) if re.fullmatch(r"[\d.]+", x) else np.nan for x in v])
st, it = num("% stromal tils"), num("% intratumoral tils"); D = np.column_stack([(S == s).astype(float) for s in ["II", "III"]])
x = M[gi]; ok3 = ~np.isnan(st) & ~np.isnan(it); zz = lambda v: (v - v[ok3].mean())/v[ok3].std()
for name, Xm, mask in (("GSTK1+亚型", np.column_stack([x, D]), np.ones(n, bool)), ("GSTK1+亚型+stromalTIL+intraTIL", np.column_stack([x, D, zz(st), zz(it)]), ok3), ("GSTK1 单变量(同样本)", x[:, None], ok3)):
    r = sm.PHReg(T[mask], Xm[mask], status=Ev[mask], ties="breslow").fit(disp=0); b0, se, p0 = r.params[0], r.bse[0], r.pvalues[0]
    print(f"[P3] {name}: n={mask.sum()} 事件={Ev[mask].sum()} logHR/SD={b0:.3f}（HR {np.exp(b0):.2f}，95%CI {np.exp(b0-1.96*se):.2f}-{np.exp(b0+1.96*se):.2f}）P={p0:.4f}")
# P4 通路家族
fam = ["GSTK1", "GSTP1", "GSTM1", "GSTT1", "GSTA1", "GSTA2", "MGST1", "MGST2", "MGST3", "GPX4", "GCLC", "GCLM", "GSR", "SLC7A11", "NFE2L2"]; rr = []
for g in fam:
    if g in genes:
        i = genes.index(g); r = sm.PHReg(T, M[i][:, None], status=Ev, ties="breslow").fit(disp=0); rr.append((g, r.params[0], r.pvalues[0], P0[i], (v := rd[g] > 0).mean()))
    else: rr.append((g, np.nan, np.nan, np.nan, np.nan))
F = pd.DataFrame(rr, columns=["gene", "logHR_SD", "P_wald", "P_score", "frac_expr"]); mm = F.P_wald.notna()
p = F.loc[mm, "P_wald"].values; o = np.argsort(p); q = np.minimum.accumulate((p[o]*len(p)/(np.arange(len(p))+1))[::-1])[::-1]; qq = np.empty(len(p)); qq[o] = np.minimum(q, 1); F.loc[mm, "q_BH"] = qq
print("[P4] 通路家族：\n", F.round(4).to_string(index=False)); F.to_csv(OUT/"pathway_family_cox.csv", index=False, encoding="utf-8-sig")
neg = int((F.logHR_SD < 0).sum()); print(f"[P4] 方向为保护（logHR<0）的基因 {neg}/{int(mm.sum())}；GSTK1 之外 P<0.05: {F[(F.gene != 'GSTK1') & (F.P_wald < 0.05)].gene.tolist()}")
# P5 描述：GSTK1 与 TIL、亚型（全部 113 例）
xa = np.log2(rd["GSTK1"][ei] + 1); suba = col("expression-based subtype"); print("[P5] GSTK1 表达（log2）亚型均值:", {s: round(float(xa[suba == s].mean()), 2) for s in ["I", "II", "III"]}, "KW P=%.2g" % kruskal(*[xa[suba == s] for s in ["I", "II", "III"]]).pvalue)
for k in ("% stromal tils", "% intratumoral tils", "non-synonymous mutation burden"):
    v = np.array([float(a) if re.fullmatch(r"[\d.]+", a) else np.nan for a in col(k)]); ok_ = ~np.isnan(v); r_, p_ = spearmanr(xa[ok_], v[ok_]); print(f"[P5] GSTK1 vs {k}: rho={r_:.3f} P={p_:.3g} n={ok_.sum()}")
e82 = pd.read_csv(X/"E82_补充分析_投稿前_20261101/out/02_meta_134.csv", encoding="utf-8-sig").set_index("gene"); print("[P5] E82 肿瘤 vs 正常 GSTK1:", e82.loc["GSTK1", ["g_GSE12452", "g_GSE64634", "g_GSE53819", "g_GSE61218", "g_GSE68799", "meta_g", "ci_lo", "ci_hi", "I2"]].round(2).to_dict())
