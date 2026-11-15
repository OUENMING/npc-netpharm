"""E117：GSE132112 机制模块探针（预注册见同目录）。种子 20261115。"""
import gzip, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
GEO=Path(sys.argv[1]); OUT=Path(sys.argv[2]); rng=np.random.default_rng(20261115)
MOD={
"1 铁死亡防御":"SLC7A11 GPX4 AIFM2 GCLC GCLM GSS NFE2L2 FTH1 HMOX1",
"2 铁死亡易感":"ACSL4 LPCAT3 TFRC NCOA4",
"3 PGE2/COX轴":"PTGS2 PTGES PTGER2 PTGER4 HPGD SLCO2A1",
"4 浆细胞/IgA/B":"JCHAIN IGHA1 IGHA2 MZB1 CD79A CD19 MS4A1 CXCL13 CD38 SDC1 TNFRSF17",
"5 CD8/干扰素/耗竭":"CD8A GZMB PRF1 IFNG PDCD1 LAG3 HAVCR2 TIGIT CXCL9 CXCL10 STAT1",
"6 色氨酸/AhR":"IDO1 IDO2 TDO2 KMO KYNU AHR AHRR CYP1A1 CYP1B1",
"7 多胺":"ODC1 AZIN1 SRM SMS AMD1 SAT1 SMOX PAOX",
"8 菌源信号受体":"NR1H4 GPBAR1 VDR NR1I2 FFAR2 FFAR3 HCAR2 SLC5A8",
"9 缺氧/糖酵解":"CA9 HIF1A EPAS1 VEGFA SLC2A1 LDHA PGK1 BNIP3",
"10 HDAC/表观":"HDAC1 HDAC2 HDAC3 HDAC6 HDAC8 SIRT1"}
MOD={k:v.split() for k,v in MOD.items()}
ALL=sorted({g for v in MOD.values() for g in v})
p2g={}
with open(GEO/"GPL17585_data.txt",errors="ignore") as fh:
    cols=None
    for l in fh:
        if l.startswith("ID\t"): cols=l.rstrip("\n").split("\t"); iS=cols.index("gene_symbols"); continue
        if cols is None or l.startswith(("#","!","^")): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<=iS or not f[0].startswith("PSR") or f[iS]=="---": continue
        s={x.strip() for x in re.split(r"\s*///\s*|;|,",f[iS]) if x.strip()}&set(ALL)
        if s: p2g[f[0]]=s
sm=gzip.open(GEO/"GSE132112_series_matrix.txt.gz","rt",errors="ignore").read().split("\n")
src=[l for l in sm if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
lab={}
for s in src:
    m=re.match(r"(LA-NPC_\d+)-(.+)",s.strip('"')); lab[m.group(1)]=m.group(2)
sums={};cnt={}
with gzip.open(GEO/"GSE132112_nor_data_probelevel.txt.gz","rt",errors="ignore") as fh:
    h=fh.readline().rstrip("\n").split("\t")[1:]
    for l in fh:
        k,rest=l.split("\t",1)
        if k in p2g:
            v=np.array(rest.rstrip("\n").split("\t"),dtype=float)
            for g in p2g[k]: sums[g]=sums.get(g,0)+v; cnt[g]=cnt.get(g,0)+1
expr=pd.DataFrame({g:sums[g]/cnt[g] for g in sums},index=h).T
y=pd.Series([ (0 if lab[s] in("Complete Response","Partial Response") else 1) for s in h],index=h)  # 1=非应答
print("n nonresp/resp:",int(y.sum()),int((1-y).sum()),"genes measured:",expr.shape[0],"of",len(ALL))
print("未测到:",[g for g in ALL if g not in expr.index])
z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1),axis=0)
def auc(x):
    a=x[y==1];b=x[y==0]; u=mannwhitneyu(a,b,alternative="two-sided"); return u.statistic/(len(a)*len(b)),u.pvalue
def ci(x,B=2000):
    a=x[y==1].values;b=x[y==0].values; r=[]
    for _ in range(B):
        aa=rng.choice(a,len(a));bb=rng.choice(b,len(b)); r.append((aa[:,None]>bb[None,:]).mean()+0.5*(aa[:,None]==bb[None,:]).mean())
    return np.percentile(r,[2.5,97.5])
def bh(p):
    p=np.array(p);o=np.argsort(p);q=np.empty(len(p));m=len(p);c=1
    for i in range(m-1,-1,-1):
        c=min(c,p[o[i]]*m/(i+1)); q[o[i]]=c
    return q
rows=[]
for k,gs in MOD.items():
    gm=[g for g in gs if g in z.index]
    sc=z.loc[gm].mean(axis=0); A,p=auc(sc); lo,hi=ci(sc)
    rows.append((k,len(gm),len(gs),A,lo,hi,p))
M=pd.DataFrame(rows,columns=["module","n_measured","n_listed","AUC_nonresponse","ci_lo","ci_hi","P"]); M["q_BH10"]=bh(M.P)
M.to_csv(OUT/"module_results.csv",index=False)
pd.set_option("display.width",200); print(M.round(4).to_string(index=False))
G=[]
for k,gs in MOD.items():
    for g in gs:
        if g in z.index:
            A,p=auc(expr.loc[g]); G.append((k,g,A,p,float(expr.loc[g][y==1].mean()-expr.loc[g][y==0].mean()),cnt[g]))
G=pd.DataFrame(G,columns=["module","gene","AUC_nonresponse","P","meandiff_nonresp_minus_resp_log2","n_probes"]); G["q_BH_all"]=bh(G.P)
G.to_csv(OUT/"gene_results.csv",index=False)
print("\n基因（P<0.05 的，探索性）:\n",G[G.P<0.05].sort_values("P").round(4).to_string(index=False))
print("总基因数",len(G),"P<0.05 个数",int((G.P<0.05).sum()),"期望约",round(0.05*len(G),1))
