"""Resolve cited publications against Crossref; cache publisher-deposited metadata.

No credentials are used. Candidate matching is inspected before bibliography
generation. The cache makes the final build independent of network access.
"""
from pathlib import Path
import re, json, html, time
from difflib import SequenceMatcher
from concurrent.futures import ThreadPoolExecutor
import requests

ROOT=Path(__file__).resolve().parent
OVERRIDES={
 1:("Alzheimer's disease: the amyloid cascade hypothesis",'10.1126/science.1566067'),
 3:("Alzheimer's disease",'10.1016/S0140-6736(06)69113-7'),
 10:('Molecular mechanisms of protein aggregation from global fitting of kinetic models','10.1038/nprot.2016.010'),
 12:("Lecanemab in early Alzheimer's disease",'10.1056/NEJMoa2212948'),
 16:('EGCG remodels mature α-synuclein and amyloid-β fibrils and reduces cellular toxicity','10.1073/pnas.0910723107'),
 24:('A Network Diffusion Model of Disease Progression in Dementia','10.1016/j.neuron.2011.12.040'),
 25:('Propagation and spread of pathogenic protein assemblies in neurodegenerative diseases','10.1038/s41593-018-0238-6'),
 18:('Mathematical modeling for the pathogenesis of Alzheimer\'s disease','10.1371/journal.pone.0015176'),
 19:('Alzheimer\'s disease: a mathematical model for onset and progression','10.1093/imammb/dqw003'),
 20:('A discrete mathematical model for the aggregation of β-Amyloid','10.1371/journal.pone.0196402'),
 30:('A molecular chaperone breaks the catalytic cycle that generates toxic Aβ oligomers','10.1038/nsmb.2971'),
 45:('Prion-like spreading of Alzheimer\'s disease within the brain\'s connectome','10.1098/rsif.2019.0356'),
 46:('Donanemab in Early Symptomatic Alzheimer Disease: The TRAILBLAZER-ALZ 2 Randomized Clinical Trial','10.1001/jama.2023.13239'),
}
def norm(t): return re.sub(r'[^a-z0-9]','',html.unescape(t).lower())
def main():
    text=(ROOT/'output/ad_passivation_paper/stage1/final_draft.md').read_text(encoding='utf-8')
    refs={int(n):line for n,line in re.findall(r'^(\d+)\. (.+)$',text.split('## References')[1].split('## Tables')[0],re.M)}
    cache=ROOT/'manuscript/references_verified.json'
    cache.parent.mkdir(exist_ok=True)
    old=json.loads(cache.read_text(encoding='utf-8')) if cache.exists() else {}
    def get(item):
        n,line=item
        if n in OVERRIDES: title,doi=OVERRIDES[n]
        else: title,doi=re.split(r'\(\d{4}\)\. ',line,maxsplit=1)[1].split('. ',1)[0],None
        if str(n) in old and old[str(n)]['query']==title and (not doi or old[str(n)].get('metadata',{}).get('DOI','').lower()==doi.lower()): return str(n),old[str(n)]
        for attempt in range(3):
            try:
                url='https://api.crossref.org/works'+('/'+doi if doi else '')
                r=requests.get(url,params=None if doi else {'query.title':title,'rows':3},timeout=50)
                r.raise_for_status()
                msg=r.json()['message']
                candidates=[msg] if doi else msg['items']
                m=max(candidates,key=lambda c:SequenceMatcher(None,norm(title),norm(c.get('title',[''])[0])).ratio())
                sim=SequenceMatcher(None,norm(title),norm(m.get('title',[''])[0])).ratio()
                return str(n),dict(query=title,similarity=sim,metadata=m,checked='2026-10-15',source='Crossref publisher deposit')
            except Exception as e:
                if attempt==2: return str(n),dict(query=title,error=type(e).__name__)
                time.sleep(attempt+1)
    # Uncited legacy references 36--44 are omitted, rather than silently
    # retaining a bibliography that is disconnected from the text.
    items=[(n,refs.get(n,'')) for n in list(range(1,36))+[45,46]]
    with ThreadPoolExecutor(max_workers=4) as pool:
        data=dict(pool.map(get,items))
    def supplement(item):
        n,d=item
        if 'pubmed' in d:return n,d
        m=d.get('metadata',{})
        if not m.get('DOI'):return n,d
        r=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',
           params={'query':'DOI:'+m['DOI'],'format':'json','resultType':'core'},timeout=50)
        r.raise_for_status()
        candidates=r.json()['resultList']['result']
        if candidates:
            d['pubmed']=candidates[0]
        return n,d
    with ThreadPoolExecutor(max_workers=4) as pool:
        data=dict(pool.map(supplement,data.items()))
    cache.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    for n,d in data.items():
        m=d.get('metadata',{})
        print(n,round(d.get('similarity',0),3),m.get('DOI',d.get('error')),m.get('title',[''])[0])

if __name__=='__main__': main()
