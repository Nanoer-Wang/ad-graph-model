"""Assemble the reviewed manuscript from local, verified source material."""
from pathlib import Path
import json,re,html,shutil
import pandas as pd

ROOT=Path(__file__).resolve().parent

def plain(s):return html.unescape(re.sub('<[^>]+>','',s)).replace('\u2011','-')
def initials(s):return ''.join(x[0].upper() for x in re.findall(r'[^\W\d_]+',s,re.U))

def main():
    source=ROOT/'manuscript';data=ROOT/'results_v3/data'
    body=(source/'revised_body.md').read_text(encoding='utf-8')
    body=body.replace('k_fY','kf Y').replace('k_onP','k_on P').replace('ā_untreated(t = 200)','ā (untreated, t = 200)')
    body=body.replace('12 covariance-based trajectory draws','12 independent parameter draws using the local marginal standard errors')
    body=body.replace('descriptive spread of 12 covariance-based draws','descriptive spread of 12 independent marginal Gaussian parameter draws')
    # Keep the large equilibrium figure near the opening result rather than
    # stranding a short continuation above an otherwise empty page.
    start=body.index('![fig4]');end=body.index('### 4.4',start)
    fig4=body[start:end];body=body[:start]+body[end:]
    at=body.index('Analytic fast-variable elimination gives')
    body=body[:at]+fig4+body[at:]
    start=body.index('![fig6]');end=body.index('### 4.6',start)
    fig6=body[start:end];body=body[:start]+body[end:]
    at=body.index('{{BLISS_RESULTS}}');body=body[:at]+fig6+body[at:]
    body=body.replace('(A) Langmuir occupancy and binding kinetics. (B) Oligomer trajectories and reference interface activity. (C) Oligomer peak and [0,40] integral with descriptive logistic fits.',
        '(A) Equilibrium Langmuir occupancy. (B) Oligomer trajectories at four infusion rates. (C) Peak and integrated oligomer burdens over [0,40], each normalised to its value at the lowest sampled dose.')
    body=body.replace(' Logistic EC50 and Hill values summarise a sampled model curve and are not pharmacological potency estimates.','')
    syn=json.loads((data/'synergy_matched_summary.json').read_text())
    near=syn['nearest_grid'];exact=syn['exact_reference']
    paragraph=(f"With exactly matched single-agent controls, {100*syn['fraction_positive']:.0f}% of the 100 grid points have positive Bliss excess (tolerance 10⁻⁸), with mean excess {syn['mean_excess']:.3f}. Points at baseline clearance have zero excess up to numerical error. At the sampled coordinate nearest the nominal reference, cY = {near['cY']:.4f} and uP = {near['uP']:.6f}, the single-agent benefits are {near['E_P']:.3f} and {near['E_C']:.3f}; the predicted benefit is {near['bliss']:.3f}, the simulated benefit {near['observed']:.3f} and the excess {near['excess']:.3f}. At the exact nominal coordinate cY = 0.50 and uP = 0.012, with cZ still at 0.114, separate simulations give late damage {exact['a200']:.3f}, predicted benefit {exact['bliss']:.3f}, simulated benefit {exact['observed']:.3f} and excess {exact['excess']:.3f}. The lower damage of 0.093 in Table 2 belongs to a different protocol that also enhances fibril clearance.")
    body=body.replace('{{BLISS_RESULTS}}',paragraph.replace('10⁻⁸','10⁻⁶'))
    old=(ROOT/'output/ad_passivation_paper/stage1/final_draft.md').read_text(encoding='utf-8')
    tables=old.split('## Tables\n',1)[1]
    table1=tables.split('**Table 2.**')[0].strip()
    table1=table1.replace('calibrated to reproduce','selected to illustrate').replace('Conservative network Laplacian','Row-normalised burden mixing').replace('Sets the time scale; Xs = λ0/cX = 0.786','Reference monomer scale Xs = λ0/cX = 0.786').replace('Max production at a = (β−1)/2β = 0.39','Maximum supply at a = 0.389')
    table2='''**Table 2.** Reference interventions from onset t = 8. Peak and integrated oligomer burdens use [0,40]. Damage at t = 200 is a late endpoint rather than a demonstrated asymptotic limit. Parameter schedules are specified in Section 3.3. Values are rounded from the archived simulation tables.

| Strategy | Ypeak | AUCY to 40 | Damage at 40 | Damage at 200 |
|---|---|---|---|---|
| Untreated | 0.346 | 4.735 | 0.955 | 0.912 |
| Passivation | 0.080 | 1.467 | 0.302 | 0.663 |
| Clearance enhancement | 0.028 | 0.662 | 0.044 | 0.725 |
| Primary nucleation inhibition | 0.294 | 2.408 | 0.848 | 0.911 |
| Propagation inhibition | 0.352 | 4.761 | 0.954 | 0.912 |
| Combination | 0.014 | 0.448 | 0.024 | 0.093 |
'''
    table3='**Table 3.**'+tables.split('**Table 3.**')[1].split('**Table 4.**')[0]
    table3=table3.replace('full matrix in Fig. 7B and Supplementary Table S1','full matrix in Fig. 7B and the archived PRCC table; p-values are nominal and unadjusted')
    fit=json.loads((data/'fit_result.json').read_text())
    table4='**Table 4.** Synthetic benchmark recovery. SE is the local Gauss–Newton standard error from output-scaled least squares; it is descriptive and does not establish calibrated coverage.\n\n| Parameter | Generating value | Best fit | Local SE | Recovery error |\n|---|---|---|---|---|\n'
    for n,best in fit['best'].items():
        r=fit['recovery'][n]
        table4+=f"| {n} | {r['true']:.4g} | {best:.4g} | {fit['sd'][n]:.3g} | {100*r['rel_err']:+.1f}% |\n"
    body=body.replace('{{TABLES}}','\n\n'.join([table1,table2,table3.strip(),table4]))
    body=body.replace('**Figure 9.** Supplementary Figure S1.','**Figure S1.**')
    # Vancouver numbering follows first occurrence, including ranges.
    order=[]
    def expand(m):
        result=[]
        for chunk in m.group(1).split(','):
            ns=re.split('[–-]',chunk.strip())
            result.extend(range(int(ns[0]),int(ns[-1])+1))
        return result
    pat=re.compile(r'\[([1-9]\d*(?:[–-]\d+)?(?:,[1-9]\d*(?:[–-]\d+)?)*)\]')
    for m in pat.finditer(body):
        for n in expand(m):
            if n not in order:order.append(n)
    numbering={n:i+1 for i,n in enumerate(order)}
    body=pat.sub(lambda m:'['+','.join(str(numbering[n]) for n in expand(m))+']',body)
    cache=json.loads((source/'references_verified.json').read_text(encoding='utf-8'))
    journal_abbr={
      'Nature Reviews Drug Discovery':'Nat Rev Drug Discov','The Lancet':'Lancet',
      "Alzheimer's & Dementia":'Alzheimers Dement','Annual Review of Biochemistry':'Annu Rev Biochem',
      'Nature Reviews Molecular Cell Biology':'Nat Rev Mol Cell Biol',
      'Physical Chemistry Chemical Physics':'Phys Chem Chem Phys',
      'Proceedings of the National Academy of Sciences':'Proc Natl Acad Sci U S A',
      'Nature Protocols':'Nat Protoc','New England Journal of Medicine':'N Engl J Med',
      "The Journal of Prevention of Alzheimer's Disease":'J Prev Alzheimers Dis',
      'Journal of Molecular Biology':'J Mol Biol','Chemical Biology & Drug Design':'Chem Biol Drug Des',
      'Nature Structural & Molecular Biology':'Nat Struct Mol Biol','PLOS ONE':'PLoS One',
      'Mathematical Medicine and Biology: A Journal of the IMA':'Math Med Biol',
      'Journal of Neuroscience':'J Neurosci','The Journal of Neuroscience':'J Neurosci',
      'Nature Neuroscience':'Nat Neurosci','Nature Communications':'Nat Commun','Nature Medicine':'Nat Med',
      'Acta Neuropathologica':'Acta Neuropathol','Journal of Theoretical Biology':'J Theor Biol',
      'Journal of The Royal Society Interface':'J R Soc Interface'}
    ref_lines=[];bib=[];audit=[]
    for n in order:
        d=cache[str(n)]
        if 'metadata' not in d:raise ValueError(f'Unverified reference {n}')
        m=d['metadata'].copy();pm=d.get('pubmed',{});ji=pm.get('journalInfo',{})
        if pm:
            if not m.get('author'):
                m['author']=[dict(family=a.get('lastName',a.get('collectiveName','')),given=a.get('firstName','')) for a in pm.get('authorList',{}).get('author',[]) if a.get('lastName') or a.get('collectiveName')]
            m['volume']=ji.get('volume',m.get('volume',''))
            m['issue']=ji.get('issue',m.get('issue',''))
            m['page']=pm.get('pageInfo',m.get('page',''))
            m['published-print']={'date-parts':[[int(pm['pubYear'])]]}
        authors=m.get('author',[])
        if not authors:raise ValueError(f'Missing authors for {n}')
        names=[plain(a.get('family',''))+' '+initials(a.get('given',''))+(' '+a['suffix'] if a.get('suffix') else '') for a in authors]
        author_text=', '.join(names[:6])+(', et al' if len(names)>6 else '')
        title=plain(m['title'][0]).rstrip('.')
        if m.get('subtitle'):title+=': '+plain(m['subtitle'][0]).rstrip('.')
        # Use the print volume year where available, not an online-first year.
        year=m.get('published-print',m.get('published',m['issued']))['date-parts'][0][0]
        journal=plain(m.get('container-title',[''])[0]);journal=ji.get('journal',{}).get('medlineAbbreviation',journal_abbr.get(journal,journal))
        vol=m.get('volume','');issue=m.get('issue','');page=m.get('page',m.get('article-number',''))
        if not (vol and page): raise ValueError(f'Incomplete volume/pages for {n}')
        doi=m['DOI']
        reference=f'{author_text}. {title}. {journal}. {year};{vol}'+(f'({issue})' if issue else '')+f':{page}. doi: {doi}.'
        ref_lines.append(f'{numbering[n]}. {reference}')
        key=f'ref{numbering[n]}'
        fields=dict(author=' and '.join(plain(a.get('family',''))+', '+plain(a.get('given','')) for a in authors),
            title=title,journal=plain(m['container-title'][0]),year=str(year),volume=str(vol),pages=str(page).replace('-','--'),doi=doi)
        if issue:fields['number']=issue
        bib.append('@article{'+key+',\n'+',\n'.join('  '+k+' = {'+v+'}' for k,v in fields.items())+'\n}')
        audit.append(dict(number=numbering[n],original_id=n,doi=doi,source='https://api.crossref.org/works/'+doi,pubmed_id=pm.get('id'),
            title=title,first_author=names[0],year=year,volume=vol,pages=page,checked='2026-10-15'))
    body=body.replace('{{REFERENCES}}','\n'.join(ref_lines))
    if '{{' in body:raise ValueError('Unresolved manuscript marker')
    (source/'final_draft.md').write_text(body,encoding='utf-8')
    (source/'references.bib').write_text('\n\n'.join(bib)+'\n',encoding='utf-8')
    (source/'reference_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    images=source/'images';images.mkdir(exist_ok=True)
    figs=['fig1_model_schematic','fig2_baseline_dynamics','fig3_spatial_propagation','fig4_bifurcation_bistability',
          'fig5_passivation_dose_response','fig6_therapy_design','fig7_global_sensitivity','fig8_identifiability','figS1_fold_detail']
    for i,name in enumerate(figs):
        dest='figS1.png' if i==8 else f'fig{i+1}.png'
        shutil.copyfile(ROOT/'results_v3/figures'/f'{name}.png',images/dest)
    print(f'Assembled {len(body.split())} words; {len(order)} verified, cited references; 9 figures')

if __name__=='__main__': main()
