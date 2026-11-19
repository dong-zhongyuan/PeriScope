import urllib.request, urllib.parse, json, sys, re, time, xml.etree.ElementTree as ET
UA={'User-Agent':'research-client/1.0'}; B='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(u):
    for i in range(5):
        try: return urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=45).read()
        except Exception as e: err=e; time.sleep(2+3*i)
    print('ERR',u[:80],err); return None
for line in open(sys.argv[1]):
    line=line.strip()
    if not line: continue
    q,*kw=line.split('||'); kw=kw[0] if kw else ''
    b=get(B+'esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':q,'retmode':'json','retmax':4,'sort':'relevance'})); time.sleep(.4)
    if not b: continue
    r=json.loads(b)['esearchresult']; ids=r['idlist']; print('##',q,'n=',r['count'])
    if not ids: continue
    x=get(B+'efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml'})); time.sleep(.4)
    if not x: continue
    for a in ET.fromstring(x).findall('.//PubmedArticle'):
        pid=a.findtext('.//PMID'); t=''.join(a.find('.//ArticleTitle').itertext())[:110]; y=a.findtext('.//PubDate/Year') or (a.findtext('.//PubDate/MedlineDate') or '')[:4]
        ab=' '.join(''.join(e.itertext()) for e in a.findall('.//AbstractText')); s=''
        if kw:
            m=re.search(kw,ab,re.I)
            if m: s=ab[max(0,m.start()-120):m.end()+140]
        print(' ',pid,y,'|',t,'|>',s)
