import sys,json,time,urllib.request,urllib.parse
def get(u):
    for i in range(4):
        try: return json.load(urllib.request.urlopen(u,timeout=40))
        except Exception as e: time.sleep(1.5)
    return None
def run(q,n=8,mindate=None):
    p={'db':'pubmed','term':q,'retmax':n,'retmode':'json','sort':'relevance'}
    if mindate: p.update({'mindate':mindate,'maxdate':'2026','datetype':'pdat'})
    d=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urllib.parse.urlencode(p))
    r=d['esearchresult']; ids=r['idlist']; out=[f"## [{r['count']}] {q}"]
    if ids:
        time.sleep(0.4)
        s=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id='+','.join(ids))['result']
        for i in ids:
            x=s[i]; out.append(f"  {i} {x['pubdate'][:4]} {x['source'][:22]} | {x['title'][:120]}")
    time.sleep(0.4); return '\n'.join(out)
if __name__=='__main__':
    qs=[l.strip() for l in open(sys.argv[1],encoding='utf8') if l.strip() and not l.startswith('#')]
    n=int(sys.argv[2]) if len(sys.argv)>2 else 8
    for q in qs: print(run(q,n),flush=True)
