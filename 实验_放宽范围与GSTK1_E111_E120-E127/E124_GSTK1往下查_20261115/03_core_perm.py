"""E124/03 · 事后敏感性（预注册外，已标注）：把置换校准的「挑选范围」限定在 134 核心靶点（E123 实际挑选 GSTK1 的范围）。
运行：python -I 03_core_perm.py <项目根目录>"""
import sys, re, gzip
from pathlib import Path
import numpy as np
from scipy.stats import norm
R = Path(sys.argv[1]); X = R/"experiments"; rng = np.random.default_rng(20261116)
core = open(X/"E81_方案A全套重出_20261101/out/core_A134.txt").read().split()
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
genes = [g for g in core if g in rd and (rd[g] > 0).mean() >= 0.2]; M = np.array([np.log2(rd[g][ei][m] + 1) for g in genes]); M = (M - M.mean(1, keepdims=True))/(M.std(1, keepdims=True) + 1e-9)
def sz(T, Ev, M):
    num = np.zeros(M.shape[0]); var = np.zeros(M.shape[0])
    for i in np.where(Ev == 1)[0]:
        r = T >= T[i]; xr = M[:, r]; mu = xr.mean(1); num += M[:, i] - mu; var += ((xr - mu[:, None])**2).mean(1)
    return num/np.sqrt(var + 1e-12)
P0 = 2*norm.sf(np.abs(sz(T, Ev, M))); gi = genes.index("GSTK1"); print("核心可用", len(genes), "；GSTK1 P=%.2e；核心内名次 %d" % (P0[gi], int((P0 <= P0[gi]).sum())))
idx = {s: np.where(S == s)[0] for s in np.unique(S)}; B = 5000; mu = np.zeros(B); ms = np.zeros(B)
for b in range(B):
    p = rng.permutation(n); mu[b] = (2*norm.sf(np.abs(sz(T[p], Ev[p], M)))).min()
    q = np.arange(n)
    for s, ix in idx.items(): q[ix] = rng.permutation(ix)
    ms[b] = (2*norm.sf(np.abs(sz(T[q], Ev[q], M)))).min()
print("核心 132 基因内最小-P 置换：全置换 FWER P=%.4f；亚型内置换 FWER P=%.4f（亚型内置换最小 P 中位 %.1e，5%%分位 %.1e）" % ((np.sum(mu <= P0[gi])+1)/(B+1), (np.sum(ms <= P0[gi])+1)/(B+1), np.median(ms), np.percentile(ms, 5)))
