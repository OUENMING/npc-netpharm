"""E124/04 · 事后描述：文献里 NPC 的 GST 家族成员 GSTM3（放射敏感性，Br J Cancer 2024）在 GSE102349 里的进展 Cox 与 GSTK1 的相关。预注册外，已标注。
运行：python -I 04_gstm3_check.py <项目根目录>"""
import sys, re, gzip
from pathlib import Path
import numpy as np, statsmodels.api as sm
from scipy.stats import spearmanr
X = Path(sys.argv[1])/"experiments"
with gzip.open(X/"E4_核实思路_20261002/geo_raw/GSE102349_series_matrix.txt.gz", "rt") as f: rows = [l.rstrip("\n").split("\t") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
pid = [re.sub(r"\D", "", x.strip('"')).lstrip("0") for x in rows[0][1:]]; ch = {}
for r in rows[1:]: ch[r[1].strip('"').split(":")[0]] = [x.strip('"').split(": ", 1)[-1] for x in r[1:]]
with gzip.open(X/"E21_患者锚定代谢物_20261003/geo_raw/GSE102349_NPC_mRNA_processed.txt.gz", "rt") as f:
    hdr = next(f).rstrip("\n").split("\t"); rd = {}
    for l in f: p = l.rstrip("\n").split("\t"); rd.setdefault(p[0], np.array([float(v) for v in p[1:]]))
cols = [re.search(r"T0*(\d+)$", c).group(1) for c in hdr[1:]]; pos = {c: i for i, c in enumerate(cols)}
ok = [(i, pos[p]) for i, p in enumerate(pid) if p in pos]; ci = np.array([i for i, _ in ok]); ei = np.array([j for _, j in ok])
col = lambda k: np.array([ch[k][i] for i in ci]); tte, ev, sub = col("time to event"), col("event"), col("expression-based subtype")
m = (tte != "N/A") & np.isin(ev, ["Disease progression", "Last follow-up"]); T = tte[m].astype(float); Ev = (ev[m] == "Disease progression").astype(int); S = sub[m]
D = np.column_stack([(S == s).astype(float) for s in ["II", "III"]])
def z(g): v = np.log2(rd[g][ei][m] + 1); return (v - v.mean())/v.std()
for g in ["GSTK1", "GSTM3", "GSTP1", "GSTM1", "GSTA1", "MGST1", "PDE4D", "PPARG", "CD274", "SIRPA", "CD47", "CALR"]:
    if g not in rd: print(g, "未测"); continue
    a = sm.PHReg(T, z(g)[:, None], status=Ev, ties="breslow").fit(disp=0); b = sm.PHReg(T, np.column_stack([z(g), D]), status=Ev, ties="breslow").fit(disp=0)
    r, p = spearmanr(np.log2(rd["GSTK1"][ei]+1), np.log2(rd[g][ei]+1)) if g != "GSTK1" else (1, 0)
    print(f"{g:7s} 单变量 logHR={a.params[0]:+.3f} P={a.pvalues[0]:.4f} | +亚型 logHR={b.params[0]:+.3f} P={b.pvalues[0]:.4f} | 与GSTK1 rho={r:+.2f} | 表达>0 比例 {(rd[g]>0).mean():.2f}")
