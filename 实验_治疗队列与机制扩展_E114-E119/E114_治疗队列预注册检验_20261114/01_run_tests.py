#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E114/01 · 三个带治疗信息的公开 NPC 队列的预注册检验（AI 初算）。
预注册见 预注册_先于看数据.md。表达量矩阵/注释放在会话 scratchpad（公开数据，可由 GEO 重新下载，URL 见下）：
  GSE224450 series matrix / GSE103611 series matrix + GPL19251 注释 / GSE132112_nor_data_probelevel.txt.gz + GPL17585 注释
运行：<带 numpy/scipy/pandas 的 python> 01_run_tests.py <geo目录> <输出目录>
"""
import gzip, sys, re, csv, math
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, norm

GEO = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
RNG = np.random.default_rng(20261114)   # 自助法固定种子

TIER = {1: ["ADRA2A", "ADRB2"], 2: ["PTGS2", "CA9"], 3: ["ADRA2B", "MMP1", "AKR1C3", "MAOB", "HDAC2"]}
TIER3_EXTRA_224450 = ["AKT1", "CTNNB1", "EGFR"]     # 预注册：仅 GSE224450 面板内与 hub 重合者
ALL_GENES = sorted({g for v in TIER.values() for g in v} | set(TIER3_EXTRA_224450))

def auc_ci(pos, neg, B=2000):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    def auc(p, n):
        r = pd.Series(np.concatenate([p, n])).rank().values
        return (r[:len(p)].sum() - len(p) * (len(p) + 1) / 2) / (len(p) * len(n))
    a = auc(pos, neg)
    bs = [auc(RNG.choice(pos, len(pos)), RNG.choice(neg, len(neg))) for _ in range(B)]
    return a, np.percentile(bs, 2.5), np.percentile(bs, 97.5)

def mde_auc(n1, n2, alpha=0.05, power=0.8):
    """Hanley-McNeil 近似：80% 功效、双侧 α 下可检出的最小 AUC（二分法）。"""
    za, zb = norm.ppf(1 - alpha / 2), norm.ppf(power)
    def se(a):
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        return math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2))
    lo, hi = 0.5001, 0.999
    for _ in range(60):
        mid = (lo + hi) / 2
        if (mid - 0.5) / ((se(0.5) + se(mid)) / 2 if False else se(0.5)) >= za + zb * se(mid) / se(0.5) : hi = mid
        else: lo = mid
    return hi

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); r = np.empty(n)
    q = p[o] * n / (np.arange(n) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; r[o] = np.minimum(q, 1)
    return r

def run(df_expr, labels, pos_label, neg_label, cohort, gene_list_by_tier, note=""):
    """df_expr: gene×sample (log2)；labels: Series sample→label。"""
    rows = []
    s_pos = labels.index[labels == pos_label]; s_neg = labels.index[labels == neg_label]
    n1, n2 = len(s_pos), len(s_neg)
    mde = mde_auc(n1, n2)
    for tier, genes in gene_list_by_tier.items():
        res = []
        for g in genes:
            if g not in df_expr.index:
                rows.append(dict(cohort=cohort, tier=tier, gene=g, status="面板/注释未覆盖", n_pos=n1, n_neg=n2, MDE_AUC=round(mde, 3))); continue
            p = df_expr.loc[g, s_pos].astype(float).values; n = df_expr.loc[g, s_neg].astype(float).values
            u, pv = mannwhitneyu(p, n, alternative="two-sided")
            a, lo, hi = auc_ci(p, n)
            res.append(dict(cohort=cohort, tier=tier, gene=g, status="已检验", n_pos=n1, n_neg=n2, median_pos=np.median(p), median_neg=np.median(n),
                            AUC=a, AUC_lo=lo, AUC_hi=hi, P=pv, MDE_AUC=round(mde, 3)))
        if res:
            q = bh([r["P"] for r in res])
            for r, qq in zip(res, q):
                r["BH_q_in_tier"] = qq
                r["判读"] = "显著关联" if qq < 0.05 else ("提示性线索" if r["P"] < 0.05 else "未检出")
            rows += res
    out = pd.DataFrame(rows); out["note"] = note
    return out

def parse_series(path):
    L = gzip.open(path, "rt", errors="ignore").read().split("\n")
    meta = {}
    for l in L:
        if l.startswith("!Sample_"):
            k = l.split("\t")[0]; meta.setdefault(k, []).append([x.strip('"') for x in l.split("\t")[1:]])
    i = [n for n, l in enumerate(L) if l.startswith("!series_matrix_table_begin")][0]
    hdr = [x.strip('"') for x in L[i + 1].split("\t")]
    return meta, hdr, L[i + 2:]

def gene_probes(annot_path, genes):
    """gene -> set(probe id)。
    GPL19251(HuGene2.0 ST 外显子级)：probeset_type == main，gene_assignment 任一条目符号等于目标基因。
    GPL17585(HTA2.0)：仅取外显子探针集 PSR*.hg.1（排除 JUC 剪接探针），gene_symbols 拆分后等于目标基因。"""
    gp = {g: set() for g in genes}
    with open(annot_path, errors="ignore") as fh:
        cols = None
        for l in fh:
            if l.startswith("ID\t"):
                cols = l.rstrip("\n").split("\t")
                mode = "A" if "gene_assignment" in cols else "S"
                if mode == "A": iA = cols.index("gene_assignment"); iT = cols.index("probeset_type")
                else: iS = cols.index("gene_symbols")
                continue
            if cols is None or l.startswith(("#", "!", "^")): continue
            f = l.rstrip("\n").split("\t")
            if mode == "A":
                if len(f) <= max(iA, iT) or f[iT] != "main" or f[iA] == "---": continue
                syms = {e.split(" // ")[1].strip() for e in f[iA].split(" /// ") if len(e.split(" // ")) > 1}
            else:
                if len(f) <= iS or not f[0].startswith("PSR") or f[iS] == "---": continue
                syms = {x.strip() for x in re.split(r"\s*///\s*|;|,", f[iS]) if x.strip()}
            for g in genes:
                if g in syms: gp[g].add(f[0])
    return gp

results = []

# ---------- GSE224450 ----------
meta, hdr, body = parse_series(GEO / "GSE224450_series_matrix.txt.gz")
titles = meta["!Sample_title"][0]
bor = [c for c in meta["!Sample_characteristics_ch1"] if c[0].startswith("best overall response")][0]
bor = [x.split(": ")[1] for x in bor]
mat = {r.split("\t")[0].strip('"'): [float(x) for x in r.split("\t")[1:]] for r in body if r and not r.startswith("!")}
expr = pd.DataFrame(mat, index=hdr[1:]).T
expr = np.log2(expr + 1)
lab = pd.Series(bor, index=hdr[1:]); pre = [s for s, t in zip(hdr[1:], titles) if t.endswith("_Pre")]
lab_pre = lab[pre].map(lambda x: "PR" if x == "PR" else "nonPR")
print("GSE224450 pre-treatment:", lab_pre.value_counts().to_dict())
tiers_224 = {1: TIER[1], 2: TIER[2], 3: TIER[3] + TIER3_EXTRA_224450}
results.append(run(expr, lab_pre, "PR", "nonPR", "GSE224450_pre(PR vs 非PR)", tiers_224, "NanoString IO360；log2(x+1)"))

# ---------- GSE103611 ----------
if (GEO / "GPL19251_data.txt").exists():
    meta, hdr, body = parse_series(GEO / "GSE103611_series_matrix.txt.gz")
    ms = [c for c in meta["!Sample_characteristics_ch1"] if c[0].startswith("metastasis")][0]
    labs = pd.Series([x.split(": ")[1] for x in ms], index=hdr[1:])
    gp = gene_probes(GEO / "GPL19251_data.txt", ALL_GENES)
    need = {p: g for g, ps in gp.items() for p in ps}
    rows = {}
    for r in body:
        if not r or r.startswith("!"): continue
        k = r.split("\t", 1)[0].strip('"')
        if k in need: rows[k] = [float(x) for x in r.split("\t")[1:]]
    probe_df = pd.DataFrame(rows, index=hdr[1:]).T
    genes_expr = {g: probe_df.loc[[p for p in gp[g] if p in probe_df.index]].mean(axis=0) for g in ALL_GENES if any(p in probe_df.index for p in gp[g])}
    expr = pd.DataFrame(genes_expr).T
    print("GSE103611 genes:", {g: len([p for p in gp[g] if p in probe_df.index]) for g in ALL_GENES}, labs.value_counts().to_dict())
    pd.DataFrame({g: [len([p for p in gp[g] if p in probe_df.index])] for g in ALL_GENES}).to_csv(OUT / "GSE103611_probe_counts.csv", index=False)
    results.append(run(expr, labs, "no", "yes", "GSE103611(无转移 vs 转移)", {k: v for k, v in TIER.items()}, "HuGene2.0 ST 外显子级；基因=main 探针集均值；matrix 已为 log2"))

# ---------- GSE132112 ----------
if (GEO / "GPL17585_data.txt").exists() and (GEO / "GSE132112_nor_data_probelevel.txt.gz").exists():
    meta, hdr, body = parse_series(GEO / "GSE132112_series_matrix.txt.gz") if False else (None, None, None)
    sm = gzip.open(GEO / "GSE132112_series_matrix.txt.gz", "rt", errors="ignore").read().split("\n")
    src = [l for l in sm if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
    labmap = {}
    for s in src:
        s = s.strip('"'); m = re.match(r"(LA-NPC_\d+)-(.+)", s); labmap[m.group(1)] = m.group(2)
    gp = gene_probes(GEO / "GPL17585_data.txt", ALL_GENES)
    need = {p: g for g, ps in gp.items() for p in ps}
    rows = {}
    with gzip.open(GEO / "GSE132112_nor_data_probelevel.txt.gz", "rt", errors="ignore") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for l in fh:
            k = l.split("\t", 1)[0]
            if k in need: rows[k] = [float(x) for x in l.rstrip("\n").split("\t")[1:]]
    probe_df = pd.DataFrame(rows, index=h[1:]).T
    genes_expr = {g: probe_df.loc[[p for p in gp[g] if p in probe_df.index]].mean(axis=0) for g in ALL_GENES if any(p in probe_df.index for p in gp[g])}
    expr = pd.DataFrame(genes_expr).T
    lab = pd.Series({k: ("response" if v in ("Complete Response", "Partial Response") else "nonresponse") for k, v in labmap.items()})
    print("GSE132112:", lab.value_counts().to_dict(), {g: len([p for p in gp[g] if p in probe_df.index]) for g in ALL_GENES})
    pd.DataFrame({g: [len([p for p in gp[g] if p in probe_df.index])] for g in ALL_GENES}).to_csv(OUT / "GSE132112_probe_counts.csv", index=False)
    results.append(run(expr, lab, "response", "nonresponse", "GSE132112(CR+PR vs SD)", {k: v for k, v in TIER.items()}, "HTA2.0 外显子级；基因=main 探针集均值；样本内发现集，不称独立"))

res = pd.concat(results, ignore_index=True)
res.to_csv(OUT / "results_all.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(res[["cohort", "tier", "gene", "status", "n_pos", "n_neg", "AUC", "AUC_lo", "AUC_hi", "P", "BH_q_in_tier", "判读", "MDE_AUC"]].round(4).to_string())
