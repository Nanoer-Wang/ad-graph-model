"""Recompute checks affected by the manuscript revision; retain v2 results."""
from pathlib import Path
import json
import shutil
import sys
import numpy as np
import pandas as pd
from scipy.optimize import root, brentq
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'src'))
from ad_graph_model.bifurcation import (fast_equilibrium_q1, qss_damage_rate,
    node_f4, node_jac4, pseudo_potential)
from ad_graph_model.config import load_config
from ad_graph_model.sensitivity import prcc_matrix


def main():
    out = ROOT/'results_v3'/'data'
    shutil.copytree(ROOT/'results_v2'/'data', out, dirs_exist_ok=True)
    p = load_config(ROOT/'configs'/'v2_baseline.json')['parameters']
    kd = p['koff']/p['kon']
    level = lambda th: kd*th/(1-th)
    def folds(pc, x0):
        def fun(x):
            a, th = x
            f = lambda aa: qss_damage_rate(aa, pc, level(th))
            return [f(a), (f(a+1e-5)-f(a-1e-5))/2e-5]
        r = root(fun, x0, tol=1e-10)
        if not r.success or np.max(np.abs(r.fun)) > 1e-8:
            raise RuntimeError('Fold refinement did not converge')
        a, th = r.x
        xyz = fast_equilibrium_q1(a, pc, level(th))
        j = node_jac4(np.r_[xyz,a], pc, level(th))
        return dict(a=float(a), theta=float(th), residual=float(np.max(np.abs(r.fun))),
                    eigenvalues=[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(j)])
    lower, upper = folds(p,[.045,.928]), folds(p,[.27,.99])
    rows=[]
    for cy in np.linspace(.25,.55,13):
        pc=dict(p,cY=float(cy))
        f=folds(pc,[.045,.93-(cy-.25)*1.3])
        xs=p['lambda0']/p['cX']
        ana=1-(cy+p['kf'])*p['cZ']/(xs*(p['ks']*p['kf']+p['ke']*(cy+p['kf'])))
        rows.append(dict(cY=cy,theta_numeric=f['theta'],theta_analytic=ana,a_fold=f['a']))
    pd.DataFrame(rows).to_csv(out/'critical_curve.csv', index=False)
    thetas=np.unique(np.r_[np.linspace(0,.90,70),np.linspace(.90,.999,210),lower['theta']+1e-6,upper['theta']-1e-6])
    rows=[]
    max_fast_res=0.; max_fast_eig=-np.inf
    for th in thetas:
        P=level(th)
        aa=np.linspace(0,1,601)
        ff=np.array([qss_damage_rate(a,p,P) for a in aa])
        for i in np.flatnonzero(ff[:-1]*ff[1:]<0):
            a=brentq(lambda a:qss_damage_rate(a,p,P),aa[i],aa[i+1],xtol=1e-13)
            xyz=fast_equilibrium_q1(a,p,P)
            y4=np.r_[xyz,a]
            lm=float(np.linalg.eigvals(node_jac4(y4,p,P)).real.max())
            rows.append(dict(theta=th,P=P,X=xyz[0],Y=xyz[1],Z=xyz[2],a=a,lmax=lm,stable=lm<0))
    pd.DataFrame(rows).to_csv(out/'bifurcation_theta.csv',index=False)
    roots={}
    for label,th in [('none',0.),('partial',.65),('strong',.95)]:
        aa=np.linspace(0,1,1001)
        ff=np.array([qss_damage_rate(a,p,level(th)) for a in aa])
        pd.DataFrame(dict(a=aa,f=ff,phi=pseudo_potential(aa,ff),theta=th)).to_csv(out/f'potential_{label}.csv',index=False)
        roots[label]=[]
        for i in np.flatnonzero(ff[:-1]*ff[1:]<0):
            a=brentq(lambda a:qss_damage_rate(a,p,level(th)),aa[i],aa[i+1])
            roots[label].append(float(a))
        for a in np.linspace(0,.999,60):
            xyz=fast_equilibrium_q1(a,p,level(th))
            y4=np.r_[xyz,a]
            max_fast_res=max(max_fast_res,float(np.abs(node_f4(y4,p,level(th))[:3]).max()))
            max_fast_eig=max(max_fast_eig,float(np.linalg.eigvals(node_jac4(y4,p,level(th))[:3,:3]).real.max()))
    for name,y0 in [('healthy',[p['lambda0']/p['cX']*1.05,0,0,0]),('diseased',[1.2,.08,.15,.55])]:
        sol=solve_ivp(lambda t,y:node_f4(y,p,level(.95)),[0,2000],y0,
                      t_eval=np.unique(np.r_[np.linspace(0,120,1201),np.linspace(120,2000,1000)]),
                      method='LSODA',rtol=1e-9,atol=1e-11)
        if not sol.success: raise RuntimeError(sol.message)
        pd.DataFrame(dict(time=sol.t,X=sol.y[0],Y=sol.y[1],Z=sol.y[2],a=sol.y[3])).to_csv(out/f'bistability_{name}.csv',index=False)
    fit=json.loads((out/'fit_result.json').read_text())
    local=[]
    for n,best in fit['best'].items():
        sd=fit['sd'][n]
        local.append(dict(parameter=n,best=best,sd_linear=sd,ci_lo=best*np.exp(-1.96*sd/best),ci_hi=best*np.exp(1.96*sd/best)))
    pd.DataFrame(local).to_csv(out/'local_gaussian_intervals.csv',index=False)
    # Correct only p-values; the rank correlations and simulated trajectories
    # are unchanged. There are d-1 controlled variables, not d.
    tab=pd.read_csv(out/'prcc_table.csv')
    from scipy.stats import t as tdist
    d=tab['parameter'].nunique()
    df=400-2-(d-1)
    for col in [c for c in tab if c.startswith('prcc_')]:
        metric=col[5:]
        rr=tab[col].to_numpy()
        tab['p_'+metric]=2*tdist.sf(np.abs(rr)*np.sqrt(df/np.maximum(1-rr*rr,1e-12)),df)
        tab['sig_'+metric]=np.where(tab['p_'+metric]<.05,'*','')
    tab.to_csv(out/'prcc_table.csv',index=False)
    audit=dict(lower_fold=lower,upper_fold=upper,hysteresis_width=upper['theta']-lower['theta'],
               qss_roots=roots,max_fast_equilibrium_residual=max_fast_res,max_fast_eigenvalue=max_fast_eig,
               excluded_legacy_intervals='profile_likelihood_ci.csv: rerun script fixed nuisance parameters at reference values; excluded from revised manuscript',
               local_intervals='Log-Gaussian delta-method intervals from the existing scaled least-squares covariance, descriptive only; no coverage claim')
    (out/'revision_checks.json').write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit,indent=2))

if __name__=='__main__': main()
