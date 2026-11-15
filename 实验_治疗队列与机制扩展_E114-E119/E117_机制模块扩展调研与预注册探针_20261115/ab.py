import sys,re,urllib.request,time
ids=sys.argv[1].split(',')
x=urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id='+','.join(ids),timeout=60).read().decode('utf8')
for a in re.findall(r'<PubmedArticle>.*?</PubmedArticle>',x,re.S):
    pm=re.search(r'<PMID[^>]*>(\d+)',a).group(1)
    t=re.sub(r'<[^>]+>','',re.search(r'<ArticleTitle>(.*?)</ArticleTitle>',a,re.S).group(1))
    j=re.search(r'<Title>(.*?)</Title>',a).group(1)
    y=re.search(r'<PubDate>.*?<Year>(\d+)',a,re.S); y=y.group(1) if y else ''
    doi=re.search(r'<ArticleId IdType="doi">([^<]+)',a); doi=doi.group(1) if doi else ''
    ab=re.sub(r'<[^>]+>','',' '.join(re.findall(r'<AbstractText[^>]*>.*?</AbstractText>',a,re.S)))
    ab=re.sub(r'\s+',' ',ab)
    ret='RETRACTED' if 'Retracted' in a or 'Retraction' in a[:3000] else ''
    print(f"### {pm} {y} {j} | {doi} {ret}\n{t}\n{ab[:int(sys.argv[2])]}\n")
