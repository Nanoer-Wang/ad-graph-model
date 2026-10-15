"""Upload an allow-listed revision through GitHub Git Data, without git credentials.

The .env token is read only in memory, never printed, embedded in a URL, or
uploaded. Remote updates are fast-forward only. Existing paths outside the
allow list are preserved. A failed upload never updates the branch reference.
"""
from pathlib import Path
import json,base64,hashlib,argparse,time
import requests

ROOT=Path(__file__).resolve().parent
BASE='https://api.github.com/repos/Nanoer-Wang/ad-graph-model'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--publish',action='store_true');args=ap.parse_args()
    token=None
    for line in (ROOT/'.env').read_text(encoding='utf-8').splitlines():
        if '=' in line:
            k,v=line.split('=',1)
            if k.strip()=='GIT_WANG':token=v.strip().strip('"').strip("'")
    if not token:raise RuntimeError('GIT_WANG is unavailable')
    session=requests.Session();session.headers.update(Authorization='Bearer '+token,Accept='application/vnd.github+json',**{'User-Agent':'ad-model-reviewed-revision'})
    def api(method,path,payload=None):
        r=session.request(method,BASE+path,json=payload,timeout=120)
        if not r.ok:raise RuntimeError(f'GitHub {method} {path}: HTTP {r.status_code}')
        return r.json()
    ref=api('GET','/git/ref/heads/main');parent=ref['object']['sha']
    parentcommit=api('GET','/git/commits/'+parent)
    old=api('GET','/git/trees/'+parentcommit['tree']['sha']+'?recursive=1')
    oldfiles={f['path']:f['sha'] for f in old['tree'] if f['type']=='blob'}
    paths=[]
    for pattern in ['src/ad_graph_model/*.py','configs/*.json','tests/*.py','manuscript/*.md','manuscript/*.json','manuscript/*.bib',
                    'manuscript/images/*.png','results_v3/data/*','results_v3/figures/*.pdf','results_v3/figures/*.svg','paper_v3/*.pdf','paper_v3/*.docx','paper_v3/*.txt']:
        paths.extend(p for p in ROOT.glob(pattern) if p.is_file())
    for name in ['README.md','requirements.txt','requirements-revision-lock.txt','make_figures.py','run_pipeline.py','revision_analysis.py',
                 'revision_bliss.py','verify_references.py','assemble_revision.py','render_equations.py','md_to_html_paper.py',
                 'build_pdf.py','build_manuscript_docx.py','make_revision_notes.py','rerun_profile_ci.py','publish_revision.py','build_docx.py',
                 'analysis_network_stats.py','reproduce_revision.py','.gitignore']:
        p=ROOT/name
        if p.exists():paths.append(p)
    # Original inputs required by the lightweight revised analysis.
    paths.extend(p for p in (ROOT/'results_v2/data').glob('*') if p.is_file())
    entries=[]
    for p in sorted(set(paths)):
        rel=p.relative_to(ROOT).as_posix()
        if p.name=='.env' or '.git/' in rel or '__pycache__' in rel:raise RuntimeError('Forbidden upload path')
        raw=p.read_bytes()
        if token.encode() in raw:raise RuntimeError('Credential bytes detected in selected file')
        sha=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        if sha!=oldfiles.get(rel):entries.append((rel,p,sha))
    plan=dict(repository='https://github.com/Nanoer-Wang/ad-graph-model',parent=parent,
              changed_files=len(entries),changed_bytes=sum(p.stat().st_size for _,p,_ in entries),paths=[rel for rel,_,_ in entries])
    (ROOT/'tmp/revision_qa/publish_plan.json').write_text(json.dumps(plan,indent=2))
    print(json.dumps({k:v for k,v in plan.items() if k!='paths'},indent=2),flush=True)
    if not args.publish:return
    if not entries:print('No changes to publish');return
    tree=[]
    for i,(rel,p,expected) in enumerate(entries,1):
        raw=p.read_bytes()
        blob=api('POST','/git/blobs',dict(content=base64.b64encode(raw).decode(),encoding='base64'))
        if blob['sha']!=expected:raise RuntimeError('Remote blob hash mismatch')
        tree.append(dict(path=rel,mode='100644',type='blob',sha=expected))
        if i%15==0:print(f'Uploaded {i}/{len(entries)} verified files',flush=True)
    made=api('POST','/git/trees',dict(base_tree=parentcommit['tree']['sha'],tree=tree))
    commit=api('POST','/git/commits',dict(message='Revise manuscript figures, equations, references and numerical validation',tree=made['sha'],parents=[parent]))
    # No force: an intervening remote change stops publication.
    api('PATCH','/git/refs/heads/main',dict(sha=commit['sha'],force=False))
    verify=api('GET','/git/trees/'+made['sha']+'?recursive=1')
    remote={e['path']:e['sha'] for e in verify['tree'] if e['type']=='blob'}
    if any(remote.get(rel)!=sha for rel,_,sha in entries):raise RuntimeError('Published tree verification failed')
    result=dict(commit=commit['sha'],url='https://github.com/Nanoer-Wang/ad-graph-model/commit/'+commit['sha'],verified_files=len(entries))
    (ROOT/'paper_v3/github_update.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
