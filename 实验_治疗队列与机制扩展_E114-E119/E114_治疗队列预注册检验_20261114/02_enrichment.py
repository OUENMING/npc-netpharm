#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E114/02 · 134 核心靶点（及 14 稳健差异基因）对治疗结局关联的整体富集检验（AI 初算；预注册修订 2）。
运行：<python> 02_enrichment.py <geo目录> <输出目录> <核心靶点txt> <Fig6_source_data.csv>
"""
import gzip, sys, re, csv
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

GEO, OUT, CORE, FIG6 = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
rng = np.random.default_rng(20261114)
core = [l.strip() for l in open(CORE, encoding="utf-8") if l.strip()]; assert len(core) == 134
rob = [r["gene"] for r in csv.DictReader(open(FIG6, encoding="utf-8-sig")) if r["robust"] == "True"]; assert len(rob) == 14

def probe2genes(path):
    m = {}
    with open(path, errors="ignore") as fh:
        cols = None
        for l in fh:
            if l.startswith("ID\t"):
                cols = l.rstrip("\n").split("\t"); mode = "A" if "gene_assignment" in cols else "S"
                if mode == "A": iA, iT = cols.index("gene_assignment"), cols.index("probeset_type")
                else: iS = cols.index("gene_symbols")
                continue
            if cols is None or l.startswith(("#", "!", "^")): continue
            f = l.rstrip("\n").split("\t")
            if mode == "A":
                if len(f) <= max(iA, iT) or f[iT] != "main" or f[iA] == "---": continue
                s = {e.split(" // ")[1].strip() for e in f[iA].split(" /// ") if len(e.split(" // ")) > 1}
            else:
                if len(f) <= iS or not f[0].startswith("PSR") or f[iS] == "---": continue
                s = {x.strip() for x in re.split(r"\s*///\s*|;|,", f[iS]) if x.strip()}
            if s: m[f[0]] = list(s)
    return m

def gene_matrix(lines_iter, pmap, nsamp):
    acc, cnt = {}, {}
    for line in lines_iter:
        k, rest = line.split("\t", 1)
        k = k.strip('"')
        gs = pmap.get(k)
        if not gs: continue
        v = np.fromstring(rest, sep="\t")
        if v.size != nsamp: continue
        for g in gs:
            if g in acc: acc[g] += v; cnt[g] += 1
            else: acc[g] = v.copy(); cnt[g] = 1
    return pd.DataFrame({g: acc[g] / cnt[g] for g in acc}).T, cnt

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); r = np.empty(n)
    q = p[o] * n / (np.arange(n) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; r[o] = np.minimum(q, 1); return r

def analyse(expr, pos_idx, neg_idx, name):
    X = expr.values; P = X[:, pos_idx]; N = X[:, neg_idx]
    keep = (P.std(axis=1) + N.std(axis=1)) > 0
    expr = expr[keep]; X = expr.values; P = X[:, pos_idx]; N = X[:, neg_idx]
    U, p = mannwhitneyu(P, N, axis=1, method="asymptotic")
    auc = U / (P.shape[1] * N.shape[1])
    df = pd.DataFrame({"gene": expr.index, "AUC_pos": auc, "P": p, "mean_expr": X.mean(axis=1)}).set_index("gene")
    df["decile"] = pd.qcut(df.mean_expr.rank(method="first"), 10, labels=False)
    inset = [g for g in core if g in df.index]
    df["BH_q_core"] = np.nan
    df.loc[inset, "BH_q_core"] = bh(df.loc[inset, "P"].values)
    obs = int((df.loc[inset, "P"] < 0.05).sum())
    pool = df[~df.index.isin(core)]
    by_dec = {d: pool.index[pool.decile == d].values for d in range(10)}
    decs = df.loc[inset, "decile"].values
    null = np.empty(10000)
    pv = pool["P"]
    for i in range(10000):
        pick = [rng.choice(by_dec[d]) for d in decs]
        null[i] = (pv[pick].values < 0.05).sum()
    emp = (1 + (null >= obs).sum()) / (1 + 10000)
    robset = [g for g in rob if g in df.index]
    summ = dict(cohort=name, n_pos=len(pos_idx), n_neg=len(neg_idx), universe=len(df), core_measured=len(inset), core_P_lt_0_05=obs,
                core_frac=obs / len(inset), null_mean_frac=null.mean() / len(inset), null_p95_count=float(np.percentile(null, 95)), empirical_P=emp,
                core_BH_lt_0_05=int((df.loc[inset, "BH_q_core"] < 0.05).sum()), universe_frac_P_lt_0_05=float((df.P < 0.05).mean()),
                robust14_measured=len(robset), robust14_P_lt_0_05=int((df.loc[robset, "P"] < 0.05).sum()))
    df.loc[inset].sort_values("P").to_csv(OUT / f"{name}_core134_tests.csv", encoding="utf-8-sig")
    return summ, df

summaries = []
# GSE132112
sm = gzip.open(GEO / "GSE132112_series_matrix.txt.gz", "rt", errors="ignore").read().split("\n")
src = [l for l in sm if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
lab = {}
for s in src:
    m = re.match(r"(LA-NPC_\d+)-(.+)", s.strip('"')); lab[m.group(1)] = ("R" if m.group(2) in ("Complete Response", "Partial Response") else "N")
pmap = probe2genes(GEO / "GPL17585_data.txt")
with gzip.open(GEO / "GSE132112_nor_data_probelevel.txt.gz", "rt", errors="ignore") as fh:
    h = fh.readline().rstrip("\n").split("\t"); samples = h[1:]
    expr, cnt = gene_matrix(fh, pmap, len(samples))
expr.columns = samples
pos = [i for i, s in enumerate(samples) if lab[s] == "R"]; neg = [i for i, s in enumerate(samples) if lab[s] == "N"]
s1, d1 = analyse(expr, pos, neg, "GSE132112_IC应答(正类=应答)"); summaries.append(s1); print(s1)
# GSE103611
pmap = probe2genes(GEO / "GPL19251_data.txt")
L = gzip.open(GEO / "GSE103611_series_matrix.txt.gz", "rt", errors="ignore")
meta = []; 
for l in L:
    if l.startswith("!Sample_characteristics_ch1") and "metastasis" in l: meta = [x.strip('"').split(": ")[1] for x in l.rstrip("\n").split("\t")[1:]]
    if l.startswith("!series_matrix_table_begin"): break
hdr = next(L).rstrip("\n").split("\t"); samples = [x.strip('"') for x in hdr[1:]]
expr, cnt = gene_matrix((l for l in L if not l.startswith("!")), pmap, len(samples)); expr.columns = samples
pos = [i for i, m in enumerate(meta) if m == "no"]; neg = [i for i, m in enumerate(meta) if m == "yes"]
s2, d2 = analyse(expr, pos, neg, "GSE103611_远处转移(正类=无转移)"); summaries.append(s2); print(s2)
res = pd.DataFrame(summaries)
res["BH_q_over_2cohorts"] = bh(res.empirical_P.values)
res.to_csv(OUT / "enrichment_summary.csv", index=False, encoding="utf-8-sig")
print(res.T)
