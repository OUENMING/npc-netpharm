"""E116：GSE132112 非应答 vs 应答信号 × LINCS L1000 化合物共识签名（本地 Fisher）。种子 20261115。"""
import gzip, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, fisher_exact
GEO = Path(sys.argv[1]); OUT = Path(sys.argv[2]); rng = np.random.default_rng(20261115)
# 注释
p2g = {}
with open(GEO/"GPL17585_data.txt", errors="ignore") as fh:
    cols=None
    for l in fh:
        if l.startswith("ID\t"): cols=l.rstrip("\n").split("\t"); iS=cols.index("gene_symbols"); continue
        if cols is None or l.startswith(("#","!","^")): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<=iS or not f[0].startswith("PSR") or f[iS]=="---": continue
        p2g[f[0]]={x.strip() for x in re.split(r"\s*///\s*|;|,", f[iS]) if x.strip()}
sm = gzip.open(GEO/"GSE132112_series_matrix.txt.gz","rt",errors="ignore").read().split("\n")
src=[l for l in sm if l.startswith("!Sample_source_name_ch1")][0].split("\t")[1:]
labmap={}
for s in src:
    m=re.match(r"(LA-NPC_\d+)-(.+)", s.strip('"')); labmap[m.group(1)]=m.group(2)
lab=pd.Series({k:("R" if v in ("Complete Response","Partial Response") else "N") for k,v in labmap.items()})
sums={}; cnt={}
with gzip.open(GEO/"GSE132112_nor_data_probelevel.txt.gz","rt",errors="ignore") as fh:
    h=fh.readline().rstrip("\n").split("\t")[1:]
    for l in fh:
        k,rest=l.split("\t",1)
        if k not in p2g: continue
        v=np.array(rest.rstrip("\n").split("\t"),dtype=float)
        for g in p2g[k]:
            sums[g]=sums.get(g,0)+v; cnt[g]=cnt.get(g,0)+1
expr=pd.DataFrame({g:sums[g]/cnt[g] for g in sums},index=h).T
lab=lab.reindex(h); N=(lab=="N").values; R=(lab=="R").values
print("n N/R:",N.sum(),R.sum(),"genes:",expr.shape[0])
t,_=ttest_ind(expr.values[:,N],expr.values[:,R],axis=1,equal_var=False)
tt=pd.Series(t,index=expr.index).dropna().sort_values()
down=list(tt.index[:150]); up=list(tt.index[-150:][::-1])
pd.DataFrame({"gene":tt.index,"t_nonresponse_vs_response":tt.values}).to_csv(OUT/"signature_t.csv",index=False)
print("PTGS2 t:",tt.get("PTGS2"), "rank from top:", int((tt>tt.get("PTGS2")).sum())+1)
# 库
lib={}
for l in open(OUT/"LINCS_consensus.txt",encoding="utf8"):
    f=l.rstrip("\n").split("\t"); 
    if len(f)<3: continue
    lib[f[0].strip()]=set(x for x in f[2:] if x)
terms=pd.Series(list(lib)); 
comp={}
for tm in lib:
    m=re.match(r"(.*)\s+(Up|Down)$",tm)
    if m: comp.setdefault(m.group(1).strip(),{})[m.group(2)]=lib[tm]
print("compounds with both Up/Down:",sum(1 for c in comp.values() if len(c)==2),"of",len(comp))
U=set(expr.index)&set().union(*lib.values()); print("universe",len(U))
upS=set(up)&U; dnS=set(down)&U
def lp(a,b):
    a=a&U; b=b&U; k=len(a&b); 
    p=fisher_exact([[k,len(a)-k],[len(b)-k,len(U)-len(a)-len(b)+k]],alternative="greater")[1]
    return -np.log10(max(p,1e-300))
rows=[]
for c,d in comp.items():
    if len(d)<2: continue
    rev=lp(upS,d["Down"])+lp(dnS,d["Up"]); conc=lp(upS,d["Up"])+lp(dnS,d["Down"])
    rows.append((c,rev-conc,rev,conc,len(d["Up"]),len(d["Down"])))
df=pd.DataFrame(rows,columns=["compound","reversal_score","rev","conc","nUp","nDown"]).sort_values("reversal_score",ascending=False).reset_index(drop=True)
df["percentile"]=100*(1-(df.index)/(len(df)-1))
df.to_csv(OUT/"compound_reversal_all.csv",index=False)
print("top20:\n",df.head(20)[["compound","reversal_score","percentile"]].to_string())
A=["celecoxib","aspirin","indomethacin","ibuprofen","diclofenac","NS-398","meloxicam","naproxen","rofecoxib"]
B=["quercetin","apigenin","luteolin","kaempferol","sodium butyrate","urolithin"]
low=df.compound.str.lower()
out=[]
for grp,names in (("A_COX",A),("B_ligand",B)):
    for n in names:
        hit=df[low.isin({n.lower(), n.lower().replace(" ","-")})]
        if hit.empty: out.append((grp,n,None,None,None)); print(grp,n,"不在库内")
        else:
            for _,r in hit.iterrows():
                out.append((grp,n,r.compound,r.reversal_score,r.percentile)); print(grp,n,"->",r.compound,round(r.reversal_score,2),"pct",round(r.percentile,1))
o=pd.DataFrame(out,columns=["group","query","matched","score","percentile"]); o.to_csv(OUT/"prespecified_compounds.csv",index=False)
for grp in ("A_COX","B_ligand"):
    g=o[(o.group==grp)&o.matched.notna()].sort_values("percentile",ascending=False).drop_duplicates("query")
    if len(g)==0: continue
    med=g.percentile.median(); k=len(g)
    perm=np.array([np.median(rng.choice(df.percentile.values,k,replace=False)) for _ in range(10000)])
    print(grp,"n in lib",k,"median pct",round(med,1),"perm P(>=)",(perm>=med).mean())
