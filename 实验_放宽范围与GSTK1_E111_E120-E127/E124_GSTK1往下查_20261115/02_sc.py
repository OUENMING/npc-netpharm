"""E124/02 · GSTK1 在 GSE150430 单细胞里的细胞类型分布（描述）。运行：python -I 02_sc.py <GSE150430_hg19.txt.gz> <输出目录>"""
import gzip, sys
from pathlib import Path
import numpy as np, pandas as pd
SRC, OUT = sys.argv[1], Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
WANT = ["GSTK1", "GSTP1", "KRT19", "EPCAM", "PTPRC", "CD3E", "LYZ", "COL1A1"]
f = gzip.open(SRC, "rt"); cells = f.readline().rstrip("\n").split("\t"); ct = np.array([c.strip('"') for c in f.readline().rstrip("\n").split("\t")[1:]]); pt = np.array([c.split("_")[0] for c in cells]); rows = {}
for l in f:
    g, rest = l.split("\t", 1)
    if g in WANT: rows[g] = np.array(rest.rstrip("\n").split("\t"), dtype=float)
print("cells", len(cells), "found", sorted(rows), "missing", sorted(set(WANT) - set(rows))); tum = np.array([p.startswith("P") for p in pt])
rec = []
for c in sorted(set(ct)):
    m = (ct == c) & tum
    if m.sum() < 50: continue
    rec.append(dict(cell_type=c, n=int(m.sum()), GSTK1_mean=rows["GSTK1"][m].mean(), GSTK1_pos=(rows["GSTK1"][m] > 0).mean(), GSTP1_mean=rows["GSTP1"][m].mean()))
D = pd.DataFrame(rec).sort_values("GSTK1_mean", ascending=False); D.to_csv(OUT/"gstk1_by_celltype.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 200, "display.max_rows", 60); print(D.round(3).to_string(index=False))
mal = np.array(["malignant" in c.lower() or "epithel" in c.lower() or "tumor" in c.lower() for c in ct]); print("上皮/肿瘤类标签:", sorted(set(ct[mal])))
