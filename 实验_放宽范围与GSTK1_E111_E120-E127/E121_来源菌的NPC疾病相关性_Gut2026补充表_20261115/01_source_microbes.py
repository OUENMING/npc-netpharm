"""E121/01 · 候选代谢物来源菌 × Gut 2026 S4。预注册见 预注册_先于看数据.md。
运行：<python> -I 01_source_microbes.py <S4所在xlsx> <gutMGene csv> <228清单csv> <docking矩阵csv> <输出目录>
"""
import csv, re, sys
from pathlib import Path
import numpy as np, pandas as pd, openpyxl
from scipy.stats import fisher_exact

XLS, GM, L228, DOCK, OUT = sys.argv[1:6]; OUT = Path(OUT); OUT.mkdir(exist_ok=True, parents=True)
norm = lambda s: re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()

wb = openpyxl.load_workbook(XLS, read_only=True, data_only=True); rows = list(wb["Table S4"].iter_rows(values_only=True))
S4 = pd.DataFrame(rows[2:], columns=rows[1]).dropna(subset=["Species"])
S4["sp"] = S4["Species"].map(norm); S4["gn"] = S4["Genus"].map(norm)
S4["status"] = np.where((S4.FDR < 0.05) & (S4.Beta < 0), "耗竭", np.where((S4.FDR < 0.05) & (S4.Beta > 0), "富集", "其他"))
print("S4 物种", len(S4), S4.status.value_counts().to_dict())
sp_status = dict(zip(S4.sp, S4.status))
# 属：S4 中显著物种方向一致才定向
gstat = {}
for g, d in S4.groupby("gn"):
    sig = d[d.status != "其他"].status
    gstat[g] = "耗竭" if len(sig) and (sig == "耗竭").all() else ("富集" if len(sig) and (sig == "富集").all() else "其他")

m228 = pd.read_csv(L228, encoding="utf-8-sig"); cids = set(m228.cid.astype(str)); name_by_cid = dict(zip(m228.cid.astype(str), m228.Name))
G = pd.DataFrame(list(csv.DictReader(open(GM, encoding="utf-8-sig"))))
G["cid"] = G["Metabolite PubChem CID"].astype(str).str.strip(); G["sp"] = G["Gut Microbiota"].map(norm)
G = G[G.cid.isin(cids)].copy()
print("gutMGene 内 228 代谢物的记录", len(G), "；代谢物数", G.cid.nunique())

def analysis(Gsub, level, label):
    """level: species / genus。返回 T1 结果与 T2 表。"""
    if level == "species":
        D = Gsub[Gsub.Rank == "species"].copy(); D["key"] = D.sp; status = sp_status; universe = list(sp_status)
    else:
        D = Gsub[Gsub.Rank == "genus"].copy(); D["key"] = D.sp; status = gstat; universe = list(gstat)
    D = D[D.key.isin(status)]
    prod = set(D.key); U = pd.Series({k: status[k] for k in universe})
    isprod = U.index.isin(prod); dep = (U == "耗竭").values; enr = (U == "富集").values
    a = (isprod & dep).sum(); b = (isprod & ~dep).sum(); c = (~isprod & dep).sum(); d_ = (~isprod & ~dep).sum()
    OR, P = fisher_exact([[a, b], [c, d_]])
    print(f"\n[{label}|{level}] T1 产生物种 {isprod.sum()} / 已检 {len(U)}；耗竭: 产生 {a}/{isprod.sum()} vs 非产生 {c}/{(~isprod).sum()}；OR={OR:.2f} P={P:.4g}")
    rowsT = []
    for cid, d in D.groupby("cid"):
        ks = set(d.key)
        if len(ks) < 3: continue
        mask = U.index.isin(ks); res = {"cid": cid, "metabolite": name_by_cid.get(cid, cid), "n_producers_tested": int(mask.sum()),
                                        "n_depleted": int((mask & dep).sum()), "n_enriched": int((mask & enr).sum())}
        for tag, v in (("depl", dep), ("enr", enr)):
            t = [[(mask & v).sum(), (mask & ~v).sum()], [(~mask & v).sum(), (~mask & ~v).sum()]]
            o, p = fisher_exact(t); res["OR_" + tag] = o; res["P_" + tag] = p
        rowsT.append(res)
    T = pd.DataFrame(rowsT)
    if len(T):
        for tag in ("depl", "enr"):
            p = T["P_" + tag].values; o = np.argsort(p); n = len(p); q = p[o] * n / (np.arange(n) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
            r = np.empty(n); r[o] = np.minimum(q, 1); T["q_" + tag] = r
        T = T.sort_values("P_depl")
    T.to_csv(OUT / f"T2_{label}_{level}.csv", index=False)
    print(f"[{label}|{level}] T2 家族大小 {len(T)}；q_depl<0.05: {(T.q_depl<0.05).sum() if len(T) else 0}；q_enr<0.05: {(T.q_enr<0.05).sum() if len(T) else 0}")
    if len(T): print(T.head(12).round(4).to_string(index=False))
    return D, T

Gh = G[G["human/mouse"] == "human"]
print("=== 主分析：human 标签 ===")
Dh, Th = analysis(Gh, "species", "human")
print("\n=== 次要 (a) 敏感性：全部记录 ===")
Da, Ta = analysis(G, "species", "all")
print("\n=== 次要 (b) 属级（human）===")
Dg, Tg = analysis(Gh, "genus", "human")

print("\n=== 次要 (c) 有对接结果的化合物的产生物种在 S4 的状态 ===")
dock = pd.read_csv(DOCK, encoding="utf-8-sig", index_col=0); names = list(dock.index)
nm = {norm(n): c for c, n in name_by_cid.items()}
for n in names:
    cid = nm.get(norm(n))
    if cid is None: print(" ", n, "：不在 228 清单（名称对不上）"); continue
    d = G[(G.cid == cid) & (G.Rank == "species")]
    st = [(k, sp_status.get(k, "S4 未检")) for k in sorted(set(d.sp))]
    print(" ", n, "(CID", cid, ")：gutMGene 物种级记录", len(d), "→", st[:8])
Dh.to_csv(OUT / "human_species_records_in_S4.csv", index=False)
