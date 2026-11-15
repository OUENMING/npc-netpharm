import gzip,re,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu,spearmanr
import statsmodels.api as sm
GEO=Path(sys.argv[1]); rng=np.random.default_rng(20261115)
L={"髓系":"CD14 CD68 CSF1R LYZ CD163 MSR1 ITGAM FCGR3A VCAN S100A8 S100A9","上皮/肿瘤":"EPCAM KRT5 KRT19 KRT17 TP63","成纤维":"COL1A1 COL1A2 DCN LUM FAP"}
L={k:v.split() for k,v in L.items()}; ALL=sorted({g for v in L.values() for g in v}|{"PTGS2"})
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
sm_=gzip.open(GEO/"GSE132112_series_matrix.txt.gz","rt",errors="ignore").read().split("\n")
src=[l for l in sm_ if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
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
y=pd.Series([0 if lab[s] in("Complete Response","Partial Response") else 1 for s in h],index=h)
print("n 非应答/应答",int(y.sum()),int((1-y).sum()),"未测到",[g for g in ALL if g not in expr.index])
z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1),axis=0)
S=pd.DataFrame({k:z.loc[[g for g in gs if g in z.index]].mean(axis=0) for k,gs in L.items()})
S["PTGS2"]=z.loc["PTGS2"]
def auc(x):
    a=x[y==1];b=x[y==0];u=mannwhitneyu(a,b,alternative="two-sided");return u.statistic/(len(a)*len(b)),u.pvalue
def ci(x,B=2000):
    a=x[y==1].values;b=x[y==0].values;r=[]
    for _ in range(B):
        aa=rng.choice(a,len(a));bb=rng.choice(b,len(b));r.append((aa[:,None]>bb[None,:]).mean()+0.5*(aa[:,None]==bb[None,:]).mean())
    return np.percentile(r,[2.5,97.5])
rows=[]
for k in L:
    A,p=auc(S[k]);lo,hi=ci(S[k]);rows.append((k,A,lo,hi,p))
R=pd.DataFrame(rows,columns=["谱系分数","AUC_非应答","lo","hi","P"])
p=R.P.values;o=np.argsort(p);m=3;q=np.empty(3);c=1
for i in range(m-1,-1,-1): c=min(c,p[o[i]]*m/(i+1)); q[o[i]]=c
R["q_BH3"]=q; print(R.round(4).to_string(index=False))
print("\nPTGS2 与谱系分数的 Spearman:")
for k in L: r,pv=spearmanr(S["PTGS2"],S[k]); print(f"  {k}: rho={r:.3f} P={pv:.3g}")
def logit(cols):
    X=sm.add_constant(S[cols]);m=sm.Logit(y.astype(float),X).fit(disp=0)
    b=m.params["PTGS2"];ci_=m.conf_int().loc["PTGS2"];return np.exp(b),np.exp(ci_[0]),np.exp(ci_[1]),m.pvalues["PTGS2"]
for name,cols in (("PTGS2 单独",["PTGS2"]),("+髓系",["PTGS2","髓系"]),("+上皮",["PTGS2","上皮/肿瘤"]),("+成纤维",["PTGS2","成纤维"]),("+三个谱系",["PTGS2","髓系","上皮/肿瘤","成纤维"])):
    o_,l_,h_,pv=logit(cols); print(f"逻辑回归 {name}: PTGS2 OR/SD={o_:.2f} ({l_:.2f}-{h_:.2f}) P={pv:.4g}")
R.to_csv("out/lineage_auc.csv",index=False)
