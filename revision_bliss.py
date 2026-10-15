"""Matched single-agent controls for every combination-grid coordinate."""
from pathlib import Path
import sys,json,copy,time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from ad_graph_model.config import load_config
from ad_graph_model.simulate import run_simulation

def main():
    out=ROOT/'results_v3/data';cfg=load_config(ROOT/'configs/v2_baseline.json')
    cfg['simulation'].update(t_end=200,n_eval=1201)
    def end(changes):
        c=copy.deepcopy(cfg);c['interventions']=[dict(time=8.,set=changes)]
        res=run_simulation(c)
        return float(res.a[-1].mean())
    base=json.loads((out/'strategy_metrics_long.json').read_text())['Baseline']['a200']
    cys=np.linspace(.25,.6,10);ups=np.logspace(-3.3,-1,10)
    ctrlp=[];ctrlc=[];t0=time.time()
    for u in ups:
        ctrlp.append(end(dict(uP=float(u))))
        print(f'Passivation control {len(ctrlp)}/10 ({time.time()-t0:.0f}s)',flush=True)
    for cy in cys:
        ctrlc.append(end(dict(cY=float(cy))))
        print(f'Clearance control {len(ctrlc)}/10 ({time.time()-t0:.0f}s)',flush=True)
    G=pd.read_csv(out/'combination_grid_long.csv',index_col=0).to_numpy()
    rows=[]
    for i,u in enumerate(ups):
        for j,cy in enumerate(cys):
            ep=1-ctrlp[i]/base;ec=1-ctrlc[j]/base
            pred=ep+ec-ep*ec;obs=1-G[i,j]/base
            rows.append(dict(uP=u,cY=cy,passivation_a200=ctrlp[i],clearance_a200=ctrlc[j],
                E_P=ep,E_C=ec,bliss=pred,observed=obs,excess=obs-pred))
    syn=pd.DataFrame(rows);syn.to_csv(out/'synergy_matched.csv',index=False)
    idx=np.argmin(abs(ups-.012))*len(cys)+np.argmin(abs(cys-.5))
    exactp=end(dict(uP=.012));exactc=end(dict(cY=.5));exactcombo=end(dict(uP=.012,cY=.5))
    ep=1-exactp/base;ec=1-exactc/base;bliss=ep+ec-ep*ec;obs=1-exactcombo/base
    stats=dict(mean_excess=float(syn.excess.mean()),fraction_positive=float((syn.excess>1e-6).mean()),positive_tolerance=1e-6,
        nearest_grid=syn.iloc[idx].to_dict(),exact_reference=dict(uP=.012,cY=.5,cZ=.114,
            E_P=ep,E_C=ec,bliss=bliss,observed=obs,excess=obs-bliss,a200=exactcombo),
        note='Grid changes cY and uP only; cZ remains 0.114. Table 2/window combination also changes cZ to 0.20.')
    (out/'synergy_matched_summary.json').write_text(json.dumps(stats,indent=2))
    print(json.dumps(stats,indent=2))

if __name__=='__main__':main()
