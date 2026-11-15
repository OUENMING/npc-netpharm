"""E111/02 · PubMed 计数（E-utilities esearch）。运行：python -I 02_pubmed_counts.py out/pubmed_counts.csv"""
import sys, time, json, urllib.request, urllib.parse, csv, datetime
alias = {"PTGS2": '("PTGS2"[tiab] OR "COX-2"[tiab] OR "COX2"[tiab] OR "cyclooxygenase-2"[tiab])', "CDK1": '("CDK1"[tiab] OR "CDC2"[tiab])',
         "MAOB": '("MAOB"[tiab] OR "MAO-B"[tiab] OR "monoamine oxidase B"[tiab])', "MMP1": '("MMP1"[tiab] OR "MMP-1"[tiab])', "MMP3": '("MMP3"[tiab] OR "MMP-3"[tiab])',
         "TOP2A": '("TOP2A"[tiab] OR "topoisomerase II alpha"[tiab] OR "topoisomerase IIalpha"[tiab])', "AKR1C3": '("AKR1C3"[tiab])', "CA9": '("CA9"[tiab] OR "CAIX"[tiab] OR "carbonic anhydrase IX"[tiab] OR "carbonic anhydrase 9"[tiab])',
         "AKT1": '("AKT1"[tiab])', "EGFR": '("EGFR"[tiab])', "SRC": '("SRC"[tiab] OR "c-Src"[tiab])', "GSK3B": '("GSK3B"[tiab] OR "GSK-3beta"[tiab] OR "GSK3beta"[tiab])', "CTNNB1": '("CTNNB1"[tiab] OR "beta-catenin"[tiab])'}
NPC = '("nasopharyngeal carcinoma"[tiab] OR "nasopharyngeal cancer"[tiab] OR "nasopharyngeal"[tiab])'
CLIN = '(chemotherapy[tiab] OR radiotherapy[tiab] OR chemoradiotherapy[tiab] OR radiation[tiab] OR prognosis[tiab] OR survival[tiab] OR resistance[tiab] OR response[tiab])'
met = {"PTGS2": ["quercetin", "apigenin"], "CDK1": ["quercetin"], "MAOB": ["apigenin", "genistein", "acacetin"], "MMP1": ["caffeic acid"], "TOP2A": ["myricetin"], "MMP3": ["quercetin"], "AKR1C3": ["linoleic acid"], "CA9": ["caffeic acid", "chrysin"]}
def cnt(q):
    u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": q, "retmode": "json", "retmax": 0})
    for _ in range(4):
        try: return int(json.load(urllib.request.urlopen(u, timeout=30))["esearchresult"]["count"])
        except Exception as e: time.sleep(2)
    return None
rows = []; today = datetime.date.today().isoformat()
for g, a in alias.items():
    qs = {"gene_only": a, "gene_NPC": f"{a} AND {NPC}", "gene_NPC_clin": f"{a} AND {NPC} AND {CLIN}"}
    for m in met.get(g, []):
        qs[f"gene_NPC_{m}"] = f"{a} AND {NPC} AND ({m}[tiab])"
        qs[f"gene_{m}"] = f"{a} AND ({m}[tiab])"
    for k, q in qs.items():
        rows.append(dict(gene=g, query_type=k, count=cnt(q), date=today, query=q)); time.sleep(0.4)
w = csv.DictWriter(open(sys.argv[1], "w", newline="", encoding="utf-8-sig"), fieldnames=["gene", "query_type", "count", "date", "query"]); w.writeheader(); w.writerows(rows)
import collections; d = collections.defaultdict(dict)
for r in rows: d[r["gene"]][r["query_type"]] = r["count"]
for g, v in d.items(): print(g, v)
