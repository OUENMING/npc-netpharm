"""E123/02 · 事后敏感性（预注册外，已标注）：进展 Cox 加表达亚型校正；GSTK1 的稳健性；与 PubMed/连边。
运行：python -I 02_cox_adjust.py <项目根目录> <输出目录>"""
import sys, re, gzip, json, time, urllib.request, urllib.parse
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm
R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); X = R/"experiments"
core = sorted(open(X/"E81_方案A全套重出_20261101/out/core_A134.txt").read().split())
with gzip.open(X/"E4_核实思路_20261002/geo_raw/GSE102349_series_matrix.txt.gz", "rt") as f: rows = [l.rstrip("\n").split("\t") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
pid = [re.sub(r"\D", "", x.strip('"')).lstrip("0") for x in rows[0][1:]]; ch = {}
for r in rows[1:]: ch[r[1].strip('"').split(":")[0]] = [x.strip('"').split(": ", 1)[-1] for x in r[1:]]
with gzip.open(X/"E21_患者锚定代谢物_20261003/geo_raw/GSE102349_NPC_mRNA_processed.txt.gz", "rt") as f:
    hdr = next(f).rstrip("\n").split("\t"); rd = {}
    for l in f: p = l.rstrip("\n").split("\t"); rd.setdefault(p[0], [float(v) for v in p[1:]])
cols = [re.search(r"T0*(\d+)$", c).group(1) for c in hdr[1:]]; pos = {c: i for i, c in enumerate(cols)}
ok = [(i, pos[p]) for i, p in enumerate(pid) if p in pos]; ci = np.array([i for i, _ in ok]); ei = np.array([j for _, j in ok])
col = lambda k: np.array([ch[k][i] for i in ci]); tte, ev, sub = col("time to event"), col("event"), col("expression-based subtype")
m = (tte != "N/A") & np.isin(ev, ["Disease progression", "Last follow-up"]); T = tte[m].astype(float); Ev = (ev[m] == "Disease progression").astype(int); S = sub[m]
print("n", m.sum(), "事件", Ev.sum(), "亚型×事件:", pd.crosstab(S, Ev).to_dict())
D = np.column_stack([(S == s).astype(float) for s in ["II", "III"]])
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o]*n/(np.arange(n)+1); q = np.minimum.accumulate(q[::-1])[::-1]; r = np.empty(n); r[o] = np.minimum(q, 1); return r
rows_ = []
for g in core:
    if g not in rd or (np.array(rd[g]) > 0).mean() < 0.2: continue
    x = np.log2(np.array(rd[g])[ei][m] + 1); x = (x - x.mean())/(x.std()+1e-9)
    try:
        a = sm.PHReg(T, x[:, None], status=Ev, ties="breslow").fit(disp=0); b = sm.PHReg(T, np.column_stack([x, D]), status=Ev, ties="breslow").fit(disp=0)
        rows_.append((g, a.params[0], a.pvalues[0], b.params[0], b.pvalues[0]))
    except Exception as e: pass
C = pd.DataFrame(rows_, columns=["gene", "logHR_uni", "P_uni", "logHR_adj_subtype", "P_adj"]).dropna(); C["q_uni"] = bh(C.P_uni); C["q_adj"] = bh(C.P_adj)
C = C.sort_values("P_uni"); C.to_csv(OUT/"cox_subtype_adjusted.csv", index=False, encoding="utf-8-sig")
print(f"Cox 可算 {len(C)} 个；单变量 P<0.05 {int((C.P_uni<0.05).sum())}；亚型校正后 P<0.05 {int((C.P_adj<0.05).sum())}；q_adj<0.05 {int((C.q_adj<0.05).sum())}")
print(C.head(12).round(5).to_string(index=False))
# GSTK1 稳健性：留一事件、bootstrap
x = np.log2(np.array(rd["GSTK1"])[ei][m] + 1); x = (x - x.mean())/x.std()
rng = np.random.default_rng(20261115); bs = []
for _ in range(1000):
    i = rng.integers(0, len(T), len(T))
    if Ev[i].sum() < 5: continue
    try: bs.append(sm.PHReg(T[i], x[i][:, None], status=Ev[i], ties="breslow").fit(disp=0).params[0])
    except Exception: pass
print("GSTK1 bootstrap logHR 2.5/50/97.5 百分位:", np.round(np.percentile(bs, [2.5, 50, 97.5]), 3), "；HR/SD 中位", round(float(np.exp(np.median(bs))), 3), "；bootstrap 同号比例", round(float((np.array(bs) < 0).mean()), 3))
ev_idx = np.where(Ev == 1)[0]; loo = []
for k in ev_idx:
    keep = np.ones(len(T), bool); keep[k] = False
    try: r = sm.PHReg(T[keep], x[keep][:, None], status=Ev[keep], ties="breslow").fit(disp=0); loo.append((r.params[0], r.pvalues[0]))
    except Exception: pass
loo = np.array(loo); print("GSTK1 留一事件：logHR 范围", loo[:, 0].min().round(3), loo[:, 0].max().round(3), "；最大 P", loo[:, 1].max().round(4))
ed = pd.read_csv(X/"E81_方案A全套重出_20261101/out/gtpm_edges_A134.csv", encoding="utf-8-sig"); print("GSTK1 连到的代谢物:", sorted(set(ed[(ed.target == "GSTK1") & (ed.interaction == "代谢物预测结合靶点")].source)))
def cnt(q):
    u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": q, "retmode": "json", "retmax": 0}); time.sleep(0.4)
    return int(json.load(urllib.request.urlopen(u, timeout=30))["esearchresult"]["count"])
NPC = '("nasopharyngeal carcinoma"[tiab] OR "nasopharyngeal cancer"[tiab] OR "nasopharyngeal"[tiab])'
al = {"GSTK1": '("GSTK1"[tiab] OR "glutathione S-transferase kappa"[tiab])', "IGF1R": '("IGF1R"[tiab] OR "IGF-1R"[tiab])', "ABCC1": '("ABCC1"[tiab] OR "MRP1"[tiab])', "SLC2A1": '("SLC2A1"[tiab] OR "GLUT1"[tiab])', "CD38": '("CD38"[tiab])', "TOP2A": '("TOP2A"[tiab])', "PIK3CA": '("PIK3CA"[tiab])', "PTGS1": '("PTGS1"[tiab] OR "COX-1"[tiab])'}
for g, a in al.items(): print("PubMed", g, "×鼻咽:", cnt(f"{a} AND {NPC}"), "；2026-11-15")
