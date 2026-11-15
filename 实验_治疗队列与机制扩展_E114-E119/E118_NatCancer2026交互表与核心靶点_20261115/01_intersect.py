import sys,pandas as pd,numpy as np
from scipy.stats import fisher_exact
B="<项目根目录>"
core=[l.strip() for l in open(B+"/60_上传文件库/数据_方案A/基因名单/核心靶点_ChEMBL36_134.txt") if l.strip()]
fig6=pd.read_csv(B+"/10_A2_ChEMBL36单引擎版/figures_方案A_20261101/data/Fig6_source_data.csv")
rc=[c for c in fig6.columns if c.lower()=="robust"][0]; gc=[c for c in fig6.columns if c.lower() in("gene","symbol")][0]
deg=fig6[fig6[rc].astype(str).str.lower().isin(["true","1"])][gc].tolist()
hub=["AKT1","CTNNB1","EGFR","GSK3B","PTGS2","SRC"]
print("core",len(core),"deg",len(deg),deg,"hub",hub)
x=pd.read_excel("out/NatCancer2026_supp_tables_1-25.xlsx",sheet_name=None,header=1)
res=[]
for tab,name in (("Supplementary Table 4","mRNA"),("Supplementary Table 5","protein")):
    d=x[tab].dropna(subset=["Gene"]).copy(); d["Gene"]=d["Gene"].astype(str)
    d["inCore"]=d.Gene.isin(core)
    print("\n==",name,len(d),"features; core overlap:",int(d.inCore.sum()))
    print(d[d.inCore][["Gene","DFS_inter_p","DFS_inter_k","OS_inter_p","OS_inter_k","DFS_interaction_classification","OS_interaction_classification"]].round(4).to_string(index=False))
    print("hub in table:",d[d.Gene.isin(hub)].Gene.tolist(),"| DEG in table:",d[d.Gene.isin(deg)].Gene.tolist())
    for ep,col in (("DFS","DFS_interaction_classification"),("OS","OS_interaction_classification")):
        sig=d[col]!="NS"
        a=int((d.inCore&sig).sum());b=int((d.inCore&~sig).sum());c=int((~d.inCore&sig).sum());e=int((~d.inCore&~sig).sum())
        if d.inCore.sum()>=5:
            OR,p=fisher_exact([[a,b],[c,e]]); res.append((name,ep,a,b,c,e,OR,p))
            print(name,ep,"core sig/nonsig",a,b,"| non-core sig/nonsig",c,e,"OR",round(OR,3),"P",round(p,4))
        else: print(name,ep,"core overlap <5, descriptive:",a,b,"| nonCore",c,e)
if res:
    r=pd.DataFrame(res,columns=["table","endpoint","core_sig","core_ns","non_sig","non_ns","OR","P"])
    p=r.P.values;o=np.argsort(p);m=len(p);q=np.empty(m);c=1
    for i in range(m-1,-1,-1): c=min(c,p[o[i]]*m/(i+1)); q[o[i]]=c
    r["q_BH"]=q; print(r.round(4).to_string(index=False)); r.to_csv("out/fisher_results.csv",index=False)
