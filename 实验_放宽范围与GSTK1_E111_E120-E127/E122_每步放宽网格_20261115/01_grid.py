"""E122/01 · 每一步放宽的 36 格网格。预注册见 预注册_先于看数据.md。
运行：python -I 01_grid.py <项目根目录> <geo目录> <输出目录>"""
import sys, re, gzip, itertools
from pathlib import Path
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import mannwhitneyu
R = Path(sys.argv[1]); GEO = Path(sys.argv[2]); OUT = Path(sys.argv[3]); X = R/"experiments"; OUT.mkdir(exist_ok=True, parents=True)
rng = np.random.default_rng(20261115)
pool250 = open(X/"E81_方案A全套重出_20261101/out/sea_stp_A250.txt").read().split(); core134 = set(open(X/"E81_方案A全套重出_20261101/out/core_A134.txt").read().split())
assert len(pool250) == 250 and core134 <= set(pool250)
E = pd.read_csv(OUT/"string_pool250_ge0.15.tsv", sep="\t"); E = E.drop_duplicates(["preferredName_A", "preferredName_B"])
old = pd.read_csv(X/"E81_方案A全套重出_20261101/out/string_edges_A134.tsv", sep="\t"); oldset = {frozenset(x) for x in zip(old.preferredName_A, old.preferredName_B)}
new134 = E[(E.score >= 0.4) & E.preferredName_A.isin(core134) & E.preferredName_B.isin(core134)]; newset = {frozenset(x) for x in zip(new134.preferredName_A, new134.preferredName_B)}
print(f"[核对] 134 核心内 ≥0.4：旧(v12) {len(oldset)} 边，新(v12.5) {len(newset)} 边，重合 {len(oldset & newset)}，仅旧 {len(oldset - newset)}，仅新 {len(newset - oldset)}")
def hubpos(genes, thr):
    sub = E[(E.score >= thr) & E.preferredName_A.isin(genes) & E.preferredName_B.isin(genes) & (E.preferredName_A != E.preferredName_B)]
    G = nx.Graph(); G.add_edges_from(zip(sub.preferredName_A, sub.preferredName_B))
    deg = dict(G.degree()); mnc = {x: max((len(c) for c in nx.connected_components(G.subgraph(G[x]))), default=0) for x in G}
    btw = nx.betweenness_centrality(G); clo = nx.closeness_centrality(G)
    D = pd.DataFrame({"Degree": deg, "MNC": mnc, "Betweenness": btw, "Closeness": clo})
    rk = D.rank(ascending=False, method="min"); mr = rk.mean(axis=1); pos = mr.rank(method="min").astype(int)
    return pd.DataFrame({"mean_rank": mr, "pos": pos}), G.number_of_nodes(), G.number_of_edges()
meta = pd.read_csv(X/"E7_等待期间的论证_20261002/meta_deg_table.csv").set_index("gene")
e82 = pd.read_csv(X/"E82_补充分析_投稿前_20261101/out/02_meta_allgenes.csv", encoding="utf-8-sig").set_index("gene")
DEG = {"D0": set(meta.index[meta.is_DEG == 1]), "D1": set(meta.index[(meta.meta_FDR_BH < 0.05) & (meta.direction_consistent_3of3 == 1)]), "D2": set(e82.index[e82.consistent == True])}
print({k: len(v) for k, v in DEG.items()}, " ∩250:", {k: len(v & set(pool250)) for k, v in DEG.items()}, " ∩134:", {k: len(v & core134) for k, v in DEG.items()})
pools = {"core134": sorted(core134), "pool250": sorted(pool250)}; cells = []; hubs = {}
for pn, genes in pools.items():
    for thr in (0.4, 0.15):
        H, nn, ne = hubpos(set(genes), thr); hubs[(pn, thr)] = H; print(f"网络 {pn} STRING≥{thr}: {nn} 节点 / {ne} 边")
for pn, thr, K, dn in itertools.product(pools, (0.4, 0.15), (10, 20, 30), DEG):
    H = hubs[(pn, thr)]; hk = set(H.index[H.pos <= K]); inter = sorted(hk & DEG[dn] & set(pools[pn]))
    cells.append(dict(pool=pn, string=thr, K=K, deg=dn, n_hub=len(hk), n_inter=len(inter), inter=" ".join(inter)))
C = pd.DataFrame(cells); C.to_csv(OUT/"grid_cells.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250, "display.max_colwidth", 80, "display.max_rows", 100); print(C.drop(columns=[]).to_string(index=False))
from collections import Counter
cnt = Counter(g for s in C.inter for g in s.split()); stab = pd.DataFrame({"gene": list(cnt), "n_cells": list(cnt.values())}).sort_values("n_cells", ascending=False); stab["frac"] = stab.n_cells/36
print("\n稳定性（出现在多少格 /36）:"); print(stab.to_string(index=False))
# ---- GSE132112
sm = gzip.open(GEO/"GSE132112_series_matrix.txt.gz", "rt", errors="ignore").read().split("\n")
src = [l for l in sm if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]; lab = {}
for s in src:
    m = re.match(r"(LA-NPC_\d+)-(.+)", s.strip('"')); lab[m.group(1)] = m.group(2)
allg = set(pool250); p2g = {}
with open(GEO/"GPL17585_data.txt", errors="ignore") as fh:
    cols = None
    for l in fh:
        if l.startswith("ID\t"): cols = l.rstrip("\n").split("\t"); iS = cols.index("gene_symbols"); continue
        if cols is None or l.startswith(("#", "!", "^")): continue
        f = l.rstrip("\n").split("\t")
        if len(f) <= iS or not f[0].startswith("PSR") or f[iS] == "---": continue
        s = {x.strip() for x in re.split(r"\s*///\s*|;|,", f[iS]) if x.strip()} & allg
        if s: p2g[f[0]] = s
sums, ct = {}, {}
with gzip.open(GEO/"GSE132112_nor_data_probelevel.txt.gz", "rt", errors="ignore") as fh:
    h = fh.readline().rstrip("\n").split("\t")[1:]
    for l in fh:
        k, rest = l.split("\t", 1)
        if k in p2g:
            v = np.array(rest.rstrip("\n").split("\t"), dtype=float)
            for g in p2g[k]: sums[g] = sums.get(g, 0) + v; ct[g] = ct.get(g, 0) + 1
expr = pd.DataFrame({g: sums[g]/ct[g] for g in sums}, index=h).T
y = np.array([0 if lab[s] in ("Complete Response", "Partial Response") else 1 for s in h])
print("\nGSE132112 非应答/应答", y.sum(), (1-y).sum(), "；池内已测基因", len(expr), "/ 250；未测:", sorted(allg - set(expr.index)))
def tests(expr, y):
    rows = []
    for g in expr.index:
        v = expr.loc[g].values; a = v[y == 1]; b = v[y == 0]; u = mannwhitneyu(a, b, alternative="two-sided"); rows.append((g, u.statistic/(len(a)*len(b)), u.pvalue))
    return pd.DataFrame(rows, columns=["gene", "AUC_nonresp", "P"]).set_index("gene")
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o]*n/(np.arange(n)+1); q = np.minimum.accumulate(q[::-1])[::-1]; r = np.empty(n); r[o] = np.minimum(q, 1); return r
T = tests(expr, y); T.to_csv(OUT/"GSE132112_pool250_tests.csv", encoding="utf-8-sig")
U = stab.merge(T.reset_index(), on="gene", how="left"); m = U.P.notna(); U.loc[m, "q_union"] = bh(U.loc[m, "P"].values)
U["stable50"] = U.n_cells >= 18; U.to_csv(OUT/"union_stability_outcome.csv", index=False, encoding="utf-8-sig")
print("\n并集族（{} 个基因，{} 个已测）GSE132112:".format(len(U), int(m.sum()))); print(U.round(4).to_string(index=False))
nom = int((U.P < 0.05).sum()); print(f"\n并集内名义 P<0.05: {nom}/{int(m.sum())}（期望 5%≈{0.05*m.sum():.1f}）；池 250 全体名义 P<0.05: {int((T.P<0.05).sum())}/{len(T)} = {100*(T.P<0.05).mean():.1f}%")
# 随机对照：从 250 池已测基因随机抽同样个数，名义 P<0.05 个数 ≥ 实际的频率；以及 q<0.05 个数
pv = T.P.values; n_u = int(m.sum()); sims = np.array([(rng.choice(pv, n_u, replace=False) < 0.05).sum() for _ in range(2000)])
print(f"随机对照：抽 {n_u} 个，名义 P<0.05 个数 ≥ {nom} 的频率 = {(sims >= nom).mean():.3f}（均值 {sims.mean():.2f}）")
print("\n新阳性判据：稳定(≥18格)且非 PTGS2 且 q_union<0.05 ->", U[(U.stable50) & (U.gene != "PTGS2") & (U.q_union < 0.05)].gene.tolist())
print("稳定候选:", U[U.stable50].gene.tolist())
