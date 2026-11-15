import gzip,re,sys
import numpy as np,pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu
R=Path(sys.argv[1]); rng=np.random.default_rng(20261115)
MOD={"铁死亡防御":"SLC7A11 GPX4 GCLC GCLM NFE2L2 HMOX1".split(),"铁死亡易感":"ACSL4 NCOA4".split()}
OTH="HPGD SLCO2A1 CD79A MS4A1 CD38 SDC1 CD8A HAVCR2 CXCL9 CYP1A1 CYP1B1 ODC1 AZIN1 SRM SMS NR1H4 HIF1A EPAS1 SLC2A1 PGK1 HDAC2 HDAC8 SIRT1".split()
ALL=sum(MOD.values(),[])+OTH
spots={}
for l in open(R/"gpl4381_table.tsv",encoding="utf8",errors="ignore"):
    f=l.rstrip("\n").split("\t")
    if len(f)<8: continue
    for m in re.findall(r"\(([A-Z0-9\-]+)\)",f[7]):
        if m in ALL: spots.setdefault(m,set()).add(f[0])
need={s:g for g,ss in spots.items() for s in ss}
sm=gzip.open(R/"GSE32389_sm.txt.gz","rt",errors="ignore").read().split("\n")
ids=None;rows={};desc=None
for l in sm:
    if l.startswith("!Sample_description") and desc is None: desc=[x.strip('"') for x in l.split("\t")[1:]]
    if l.startswith('"ID_REF"'): ids=[x.strip('"') for x in l.split("\t")[1:]]; continue
    if ids and not l.startswith("!") and l.strip():
        f=l.split("\t"); k=f[0].strip('"')
        if k in need:
            v=[float(x) if x not in("","null","NA") else np.nan for x in f[1:]]
            rows[k]=v
y=pd.Series([1 if d.startswith("RR") else 0 for d in desc],index=ids)  # 1=放疗抵抗
print("抵抗/敏感:",int(y.sum()),int((1-y).sum()),"spots found:",len(rows))
df=pd.DataFrame(rows,index=ids).T
expr=pd.DataFrame({g:df.loc[[s for s in ss if s in df.index]].mean(axis=0) for g,ss in spots.items()}).T
print("基因:",expr.shape[0], "缺失值比例",float(expr.isna().mean().mean().round(3)))
z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1),axis=0)
def auc(x):
    x=x.dropna();a=x[y[x.index]==1];b=x[y[x.index]==0];u=mannwhitneyu(a,b,alternative="two-sided");return u.statistic/(len(a)*len(b)),u.pvalue
def ci(x,B=2000):
    x=x.dropna();a=x[y[x.index]==1].values;b=x[y[x.index]==0].values;r=[]
    for _ in range(B):
        aa=rng.choice(a,len(a));bb=rng.choice(b,len(b));r.append((aa[:,None]>bb[None,:]).mean()+0.5*(aa[:,None]==bb[None,:]).mean())
    return np.percentile(r,[2.5,97.5])
out=[]
for k,gs in MOD.items():
    gm=[g for g in gs if g in z.index]; sc=z.loc[gm].mean(axis=0); A,p=auc(sc);lo,hi=ci(sc); out.append((k,len(gm),A,lo,hi,p))
M=pd.DataFrame(out,columns=["module","n","AUC_resistant","lo","hi","P"]); M["q_BH2"]=[min(1,p*2/(i+1)) for i,p in zip(range(2),M.sort_values("P").P)] if False else None
ps=M.P.values; o=np.argsort(ps); q=np.empty(2); q[o[1]]=ps[o[1]]; q[o[0]]=min(ps[o[0]]*2,q[o[1]]); M["q_BH2"]=q
print(M.round(4).to_string(index=False))
G=[(g,)+auc(expr.loc[g]) for g in expr.index]; G=pd.DataFrame(G,columns=["gene","AUC_resistant","P"]).sort_values("P")
print("逐基因（探索性）:\n",G.round(4).to_string(index=False)); print("P<0.05:",int((G.P<0.05).sum()),"/",len(G),"期望",round(0.05*len(G),1))
M.to_csv(sys.argv[2]+"/gse32389_modules.csv",index=False); G.to_csv(sys.argv[2]+"/gse32389_genes.csv",index=False)
