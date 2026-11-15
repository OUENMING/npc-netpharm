import gzip,sys,collections
import numpy as np,pandas as pd
SRC=sys.argv[1]; OUT=sys.argv[2]
G={"A_PGE2":"PTGS2 PTGS1 PTGES PTGES2 PTGES3 PTGER1 PTGER2 PTGER3 PTGER4 HPGD SLCO2A1",
"B_sensor":"FFAR2 FFAR3 FFAR4 HCAR2 HCAR3 AHR GPBAR1 NR1H4 VDR NR1I2 TLR2 TLR4 TLR5 NOD2",
"C_GPCR":"ADRA2A ADRA2B ADRA2C DRD2 DRD3 DRD4 HRH2 HRH3 HRH4 HTR1A HTR6 MTNR1B GPR32 ADRB2",
"D_net":"AKT1 CTNNB1 EGFR GSK3B PTGS2 SRC MMP1 ADRA2A MMP3 MAOB AKR1C3 CDK1 SLC1A3 MMP12 TOP2A NOX4 HSD17B2 SOST CA9"}
grp={}
for k,v in G.items():
    for g in v.split(): grp.setdefault(g,[]).append(k)
want=set(grp)
f=gzip.open(SRC,"rt")
cells=f.readline().rstrip("\n").split("\t")
ct=f.readline().rstrip("\n").split("\t")[1:]
assert len(ct)==len(cells),(len(ct),len(cells))
pt=[c.split("_")[0] for c in cells]
rows={}
for l in f:
    g,rest=l.split("\t",1)
    if g in want: rows[g]=np.array(rest.rstrip("\n").split("\t"),dtype=float)
print("cells",len(cells),"found genes",len(rows),"missing",sorted(want-set(rows)))
meta=pd.DataFrame({"pt":pt,"ct":ct})
print("细胞类型（全体）:",meta.ct.value_counts().to_dict())
tum=meta.pt.str.startswith("P")
out=[]
for g,v in rows.items():
    d=meta.copy(); d["v"]=v; d["pos"]=(v>0)
    for scope,mask in (("NPC",tum),("normal_N01",~tum)):
        x=d[mask]
        for (p,c),s in x.groupby(["pt","ct"]):
            if len(s)>=20: out.append((g,scope,p,c,len(s),s.v.mean(),s.pos.mean()))
R=pd.DataFrame(out,columns=["gene","scope","patient","celltype","n_cells","mean","frac_pos"])
R.to_csv(OUT+"/per_patient_celltype.csv",index=False)
S=R[R.scope=="NPC"].groupby(["gene","celltype"]).agg(n_pat=("patient","nunique"),median_mean=("mean","median"),median_frac=("frac_pos","median"),min_mean=("mean","min"),max_mean=("mean","max")).reset_index()
S.to_csv(OUT+"/summary_celltype.csv",index=False)
print("\n=== 主要表达细胞类型（NPC，患者中位平均表达最高，纳入≥5 患者）===")
rows2=[]
for g in sorted(rows,key=lambda x:(grp[x][0],x)):
    s=S[(S.gene==g)&(S.n_pat>=5)].sort_values("median_mean",ascending=False)
    if s.empty: continue
    top=s.iloc[0]; second=s.iloc[1] if len(s)>1 else None
    rows2.append((",".join(grp[g]),g,top.celltype,round(top.median_mean,3),round(top.median_frac,3),int(top.n_pat),(second.celltype+":"+str(round(second.median_mean,3))) if second is not None else ""))
T=pd.DataFrame(rows2,columns=["group","gene","top_celltype","median_mean","median_frac_pos","n_patients","second"])
T.to_csv(OUT+"/top_celltype.csv",index=False); pd.set_option("display.width",220); print(T.to_string(index=False))
print("\n=== PTGS2 各细胞类型（NPC）===")
print(S[S.gene=="PTGS2"].sort_values("median_mean",ascending=False).round(3).to_string(index=False))
