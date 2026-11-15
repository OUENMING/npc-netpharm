"""E111/01 · hub 扩到 20 ∩ 患者组织差异。预注册见 预注册_先于看数据.md。
运行：python -I 01_hub20.py <项目根目录> <输出目录>"""
import sys, pandas as pd, numpy as np
from pathlib import Path
from scipy.stats import hypergeom
R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
X = R / "experiments"
hs = pd.read_csv(X/"E81_方案A全套重出_20261101/out/hub_scores_A134.csv", encoding="utf-8-sig")
rk = ["rank_MCC", "rank_MNC", "rank_Stress", "rank_Degree", "rank_Radiality"]
hs["mean_rank"] = hs[rk].mean(axis=1); hs = hs.sort_values("mean_rank").reset_index(drop=True); hs["pos"] = hs.mean_rank.rank(method="min").astype(int)
meta = pd.read_csv(X/"E7_等待期间的论证_20261002/meta_deg_table.csv"); deg = meta.set_index("gene")
e82 = pd.read_csv(X/"E82_补充分析_投稿前_20261101/out/02_meta_134.csv", encoding="utf-8-sig").set_index("gene")
ed = pd.read_csv(X/"E81_方案A全套重出_20261101/out/gtpm_edges_A134.csv", encoding="utf-8-sig"); ed = ed[ed.interaction == "代谢物预测结合靶点"]
tiers = pd.read_csv(X/"E59_审计后重算_20261024/gene_tiers.csv", encoding="utf-8-sig").set_index("gene")
dock = pd.read_csv(R/"figures_方案A_20261101/data/Fig5b_docking_matrix.csv", encoding="utf-8-sig", index_col=0); dock_t = {c for c in dock.columns if dock[c].notna().any()}
e114 = pd.read_csv(next((X/"E114_治疗队列预注册检验_20261114/out").glob("GSE132112_IC应答*core134_tests.csv")), encoding="utf-8-sig").set_index("gene")
core = set(hs.gene)
deg_core = {g for g in core if g in deg.index and deg.loc[g, "is_DEG"] == 1}
print("核心", len(core), "；核心内 is_DEG", len(deg_core), sorted(deg_core))
rows = []
for _, r in hs.iterrows():
    g = r.gene; mets = sorted(set(ed[ed.target == g].source))
    rows.append(dict(gene=g, pos=r.pos, mean_rank=round(r.mean_rank, 2), **{k: int(r[k]) for k in rk}, degree=int(r.degree),
        is_DEG=g in deg_core, log2FC=round(deg.loc[g, "mean_log2FC"], 3) if g in deg.index else np.nan, meta_FDR=deg.loc[g, "meta_FDR_BH"] if g in deg.index else np.nan,
        E82_consistent=bool(e82.loc[g, "consistent"]) if g in e82.index else np.nan, E82_I2=round(e82.loc[g, "I2"], 0) if g in e82.index else np.nan,
        n_metabolites=len(mets), metabolites_top5="; ".join(mets[:5]), gutMGene_tier=tiers.loc[g, "category"] if g in tiers.index else "—",
        has_docking=g in dock_t, GSE132112_AUC_nonresp=round(1 - e114.loc[g, "AUC_pos"], 3) if g in e114.index else np.nan,
        GSE132112_P=e114.loc[g, "P"] if g in e114.index else np.nan))
T = pd.DataFrame(rows); T.to_csv(OUT/"hub_ranked_134.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250, "display.max_columns", 30)
for K in (10, 15, 20):
    H = T[T.pos <= K]; a = int(H.is_DEG.sum()); M = len(core); n = len(deg_core)
    p = hypergeom.sf(a-1, M, n, len(H))
    print(f"\nH_{K}：{len(H)} 个靶点（并列并入）；∩ is_DEG = {a}（背景 {n}/{M}={n/M:.1%}）；超几何 P={p:.3g}")
    print(H[["gene", "pos", "mean_rank", "is_DEG", "log2FC", "E82_consistent", "n_metabolites", "gutMGene_tier", "has_docking"]].to_string(index=False))
H20 = T[(T.pos <= 20)]; C = H20[H20.is_DEG].copy(); print("\n候选（H_20∩DEG）:"); print(C.drop(columns=rk + ["metabolites_top5"]).to_string(index=False))
C.to_csv(OUT/"candidates_H20_DEG.csv", index=False, encoding="utf-8-sig")
# 次要：E82 consistent
C2 = H20[H20.E82_consistent == True]; print("\n次要（H_20 ∩ E82 consistent 5/5）:", C2.gene.tolist())
# 对接靶点 vs hub
print("\n对接靶点:", sorted(dock_t), " ∈ H_20:", sorted(dock_t & set(H20.gene)), " ∈ 核心:", sorted(dock_t & core))
for g in sorted(dock_t): 
    r = T[T.gene == g]; print(" ", g, "排名位次", int(r.pos.iloc[0]) if len(r) else "不在核心134")
# 方案：族内 GSE132112 BH（规则 5）
if len(C):
    p = C.GSE132112_P.values; o = np.argsort(p); n = len(p); q = np.minimum.accumulate((p[o]*n/(np.arange(n)+1))[::-1])[::-1]; qq = np.empty(n); qq[o] = np.minimum(q, 1)
    C["q_in_family"] = qq; print("\n族内 GSE132112 BH:"); print(C[["gene", "GSE132112_AUC_nonresp", "GSE132112_P", "q_in_family"]].to_string(index=False))
