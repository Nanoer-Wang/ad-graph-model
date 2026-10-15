"""Uncluttered, data-derived replacement figures for the revised manuscript."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Patch
from matplotlib.lines import Line2D
from .figstyle import C, STRATEGY_COLORS, save_fig, damage_cmap


def label(ax, letter, title):
    ax.set_title(title, fontsize=11, loc='left', pad=11)
    ax.text(-.16,1.085,letter,transform=ax.transAxes,fontsize=14,fontweight='bold')


def legend_below(ax, handles=None, labels=None, ncol=2, y=-.25):
    if handles is None:
        handles,labels=ax.get_legend_handles_labels()
    return ax.legend(handles,labels,loc='upper left',bbox_to_anchor=(-.02,y),
        fontsize=9.2,ncol=ncol,borderaxespad=0,handlelength=2,
        columnspacing=1.1,labelspacing=.7,frameon=False)


def fig1_schematic(outdir, baseline_csv, passivation_csv):
    fig=plt.figure(figsize=(8,8.1))
    gs=fig.add_gridspec(2,2,height_ratios=[1.25,1],hspace=.48,wspace=.45)
    ax=fig.add_subplot(gs[0,:]); net=fig.add_subplot(gs[1,0]); tr=fig.add_subplot(gs[1,1])
    ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
    ax.text(0,6.2,'A  Local reaction and passivation mechanisms',fontsize=12,fontweight='bold')
    def box(x,y,w,h,text,color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.06,rounding_size=.1',
           ec=color,fc='white',lw=1.5,zorder=3))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=11,color=color,zorder=4)
    def arrow(a,b,color=C['dark'],rad=0,style='-|>',ls='-'):
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=12,
           color=color,lw=1.5,connectionstyle=f'arc3,rad={rad}',linestyle=ls,zorder=2))
    box(.45,3.9,1.8,.9,'Monomer $X$',C['X'])
    box(4.05,3.9,1.8,.9,'Oligomer $Y$',C['Y'])
    box(7.65,3.9,1.8,.9,'Fibril $Z$',C['Z'])
    arrow((2.3,4.35),(4,4.35));ax.text(3.15,4.7,'Primary nucleation\n$k_nX^m$',ha='center',fontsize=9.5)
    arrow((5.9,4.35),(7.6,4.35));ax.text(6.75,4.7,'Conversion\n$k_fY$',ha='center',fontsize=9.5)
    arrow((8.55,3.85),(5.4,3.85),C['Y'],rad=-.45)
    ax.text(6.85,3,'Surface nucleation\n$k_sXZ^q(1-\\theta)$',ha='center',fontsize=9.5,color=C['Y'])
    from matplotlib.path import Path as MplPath
    elongation=MplPath([(1.35,4.85),(1.35,5.5),(8.55,5.5),(8.55,4.85)])
    ax.add_patch(FancyArrowPatch(path=elongation,arrowstyle='-|>',mutation_scale=12,lw=1.5,color=C['Z'],zorder=2))
    ax.text(4.95,5.73,'Elongation $k_eXZ(1-\\theta)$',ha='center',fontsize=9.5,color=C['Z'])
    box(.45,1.55,2.6,.85,'Passivator $P$\nOccupancy $\\theta$',C['P'])
    box(4.05,1.55,2.05,.85,'Damage $a$',C['a'])
    arrow((4.95,3.85),(4.95,2.45),C['a'])
    ax.text(4.72,2.87,'Hill toxicity',ha='right',fontsize=9.5,color=C['a'])
    # Passivation targets the reaction labels, with no crossed node arrows.
    ax.text(1.75,2.72,'Free interface: $1-\\theta$',ha='center',fontsize=9.5,color=C['P'])
    ax.text(1.75,1.07,'Suppresses both surface reactions',ha='center',fontsize=9.2,color=C['P'])
    arrow((6.15,1.98),(8.55,1.98),C['purple'])
    ax.text(7.55,2.35,'Network spread',ha='center',fontsize=9.5,color=C['purple'])
    from matplotlib.path import Path as MplPath
    feedback=MplPath([(4,1.8),(3.4,1.8),(3.4,.65),(.05,.65),(.05,4.35),(.4,4.35)])
    ax.add_patch(FancyArrowPatch(path=feedback,arrowstyle='-|>',mutation_scale=12,
        lw=1.4,color=C['a'],linestyle='--',zorder=2))
    ax.text(.05,.33,'Damage feedback on supply: $\\lambda_0(1-a)(1+\\beta a)$',fontsize=10,color=C['a'])
    ax.text(6.1,.33,'Clearance: $c_XX$, $c_YY$, $c_ZZ$',fontsize=9.5,color=C['gray'])
    # B is explicitly illustrative, not an unreported simulation.
    pos=np.array([[0,.4],[.8,.85],[.85,0],[1.7,.75],[1.7,-.1],[2.5,.4]])
    edges=[(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(3,5),(4,5)]
    for i,j in edges:net.plot(pos[[i,j],0],pos[[i,j],1],color='#BFC5CA',lw=1.5,zorder=1)
    dmg=np.array([.9,.65,.45,.25,.1,.02])
    net.scatter(pos[:,0],pos[:,1],s=180+250*dmg,c=dmg,cmap=damage_cmap(),vmin=0,vmax=1,ec='white',lw=1.4,zorder=3)
    net.annotate('Seed',xy=pos[0],xytext=(-.15,1.2),fontsize=10,
       arrowprops=dict(arrowstyle='->',color=C['dark']))
    net.text(1.2,-.7,'Illustrative network\nDarker nodes indicate greater damage',ha='center',fontsize=9.5)
    net.set(xlim=(-.4,2.9),ylim=(-.9,1.4));net.axis('off')
    label(net,'B','Coupling across regions')
    base=pd.read_csv(baseline_csv);pas=pd.read_csv(passivation_csv)
    tr.plot(base.time,base.mean_a,color=C['a'],lw=2,label='Untreated')
    tr.plot(pas.time,pas.mean_a,color=C['P'],lw=2,label='Passivation')
    tr.axvline(8,color=C['gray'],ls=':',lw=1.4)
    tr.set(xlim=(0,40),ylim=(0,1),xlabel='Time (model units)',ylabel='Mean damage')
    tr.text(.03,.93,'Treatment starts at $t=8$',transform=tr.transAxes,fontsize=9.5)
    label(tr,'C','Simulated intervention')
    legend_below(tr,ncol=2,y=-.27)
    fig.subplots_adjust(left=.10,right=.97,bottom=.12,top=.94)
    return save_fig(fig,Path(outdir)/'fig1_model_schematic',formats=('png','pdf','svg'))


def fig4_bifurcation(data_dir,outdir):
    data=Path(data_dir);bif=pd.read_csv(data/'bifurcation_theta.csv');crit=pd.read_csv(data/'critical_curve.csv')
    checks=json.loads((data/'revision_checks.json').read_text())
    low=checks['lower_fold']['theta']; high=checks['upper_fold']['theta'];ana=crit.theta_analytic.iloc[0]
    fig,axs=plt.subplots(2,2,figsize=(8,8.7));a,b,c,d=axs.ravel()
    # Split all branches explicitly; sorting stable roots together draws
    # spurious vertical lines between different attractors.
    for mask,col,ls in [(bif.stable & (bif.a>.27),C['a'],'-'),
                        (bif.stable & (bif.a<.27),C['X'],'-'),(~bif.stable,C['gray'],'--')]:
        sub=bif[mask].sort_values('theta');a.plot(sub.theta,sub.a,color=col,ls=ls,lw=2)
    a.axvline(ana,color=C['P'],ls=':',lw=1.5)
    a.axvspan(low,high,color='#D9EAF2',alpha=.8,zorder=0)
    # Dedicated zoom is an independent supplementary figure, never an inset.
    a.set(xlim=(0,1),ylim=(0,1),xlabel='Interface occupancy $\\theta$',ylabel='Equilibrium damage $a^*$')
    label(a,'A','Equilibrium branches')
    hs=[Line2D([],[],color=C['a'],lw=2),Line2D([],[],color=C['X'],lw=2),
        Line2D([],[],color=C['gray'],lw=2,ls='--'),Line2D([],[],color=C['P'],ls=':'),Patch(fc='#D9EAF2')]
    legend_below(a,hs,['High-damage stable','Low-damage stable','Unstable','Frozen-pool threshold','Bistability'],ncol=1,y=-.29)
    b.plot(crit.cY,crit.theta_analytic,color=C['X'],lw=2,label='Frozen-pool threshold')
    b.plot(crit.cY,crit.theta_numeric,color=C['a'],lw=2,ls='--',label='Low-damage fold')
    b.fill_between(crit.cY,crit.theta_numeric,1,color='#E8F0E7')
    b.fill_between(crit.cY,crit.theta_analytic,crit.theta_numeric,color='#FCEDD9')
    b.set(xlim=(.25,.55),ylim=(0,1),xlabel='Oligomer clearance $c_Y$',ylabel='Occupancy threshold')
    label(b,'B','Clearance and occupancy')
    hh,ll=b.get_legend_handles_labels()
    hh+=[Patch(fc='#E8F0E7'),Patch(fc='#FCEDD9')];ll+=['Low-damage state exists','Boundary separation']
    legend_below(b,hh,ll,ncol=1,y=-.29)
    for name,col in [('healthy',C['X']),('diseased',C['a'])]:
        tr=pd.read_csv(data/f'bistability_{name}.csv')
        c.plot(tr.time,tr.a,color=col,lw=2,label=f'{name.capitalize()} initial state')
    c.set(xlim=(0,2000),ylim=(0,.7),xlabel='Time (model units)',ylabel='Damage $a$')
    label(c,'C','Two attractors at $\\theta=0.95$')
    legend_below(c,ncol=1,y=-.29)
    for name,col in [('none',C['a']),('partial',C['P']),('strong',C['X'])]:
        po=pd.read_csv(data/f'potential_{name}.csv');th=po.theta.iloc[0]
        d.plot(po.a,po.phi,color=col,lw=2,label=f'$\\theta={th:.2f}$')
        for r in checks['qss_roots'][name]:
            slope=np.interp(min(r+1e-4,1),po.a,po.f)-np.interp(max(r-1e-4,0),po.a,po.f)
            d.scatter(r,np.interp(r,po.a,po.phi),s=35,fc=col if slope<0 else 'white',ec=col,zorder=5)
    d.set(xlim=(0,1),xlabel='Damage coordinate $a$',ylabel='Reduced potential $\\Phi(a)$')
    label(d,'D','Continuous quasi-steady reduction')
    legend_below(d,ncol=3,y=-.29)
    d.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.14,wspace=.43,hspace=1.35)
    paths=save_fig(fig,Path(outdir)/'fig4_bifurcation_bistability',formats=('png','pdf','svg'))
    fig,ax=plt.subplots(figsize=(4.8,3.8))
    for mask,col,ls in [(bif.stable & (bif.a>.27),C['a'],'-'),(bif.stable & (bif.a<.27),C['X'],'-'),(~bif.stable,C['gray'],'--')]:
        sub=bif[mask].sort_values('theta');ax.plot(sub.theta,sub.a,color=col,ls=ls,lw=2)
    ax.axvline(low,color=C['gray'],lw=1,ls=':');ax.axvline(high,color=C['gray'],lw=1,ls=':')
    ax.set(xlim=(.92,1),ylim=(0,.55),xlabel='Interface occupancy $\\theta$',ylabel='Equilibrium damage $a^*$')
    ax.set_title(f'Folds: {low:.4f} and {high:.4f}',fontsize=11)
    fig.tight_layout();paths+=save_fig(fig,Path(outdir)/'figS1_fold_detail',formats=('png','pdf','svg'))
    return paths


def fig5_dose(data_dir,outdir):
    data=Path(data_dir);dose=pd.read_csv(data/'dose_response.csv')
    kin=pd.read_csv(data/'dose_kinetics.csv');late=pd.read_csv(data/'dose_response_long.csv')
    checks=json.loads((data/'revision_checks.json').read_text())
    fig,axs=plt.subplots(2,2,figsize=(8,8.5));a,b,c,d=axs.ravel()
    P=np.logspace(-3,1.2,400);a.plot(P,P/(.15+P),color=C['X'],lw=2)
    a.axvline(.15,color=C['gray'],ls=':',lw=1.3);a.scatter(.15,.5,color=C['dark'],s=30,zorder=4)
    a.set(xscale='log',xlim=(.001,16),ylim=(0,1),xlabel='Free passivator $P$',ylabel='Equilibrium occupancy')
    label(a,'A','Langmuir occupancy')
    legend_below(a,[Line2D([],[],color=C['gray'],ls=':')],['$K_D=0.15$'],ncol=1,y=-.27)
    cols=['#D1A95D','#E69F00','#D55E00','#89362A']
    for u,col in zip([.003,.006,.012,.025],cols):
        sub=kin[np.isclose(kin.uP,u)]
        b.plot(sub.time,sub.mean_Y,color=col,lw=2,label=f'$u_P={u:g}$')
    b.set(xlim=(0,40),xlabel='Time (model units)',ylabel='Mean oligomer burden')
    label(b,'B','Oligomer trajectories')
    legend_below(b,ncol=2,y=-.27)
    for metric,col,lab,mk in [('Y_peak',C['Y'],'Peak oligomer','o'),('AUC_Y',C['Z'],'Integrated oligomer','s')]:
        sub=dose.sort_values('uP');relative=sub[metric]/sub[metric].iloc[0]
        c.plot(sub.uP,relative,color=col,lw=1.8,marker=mk,ms=3,label=lab)
    c.set(xscale='log',xlabel='Passivation infusion $u_P$',ylabel='Fraction of lowest-dose value',ylim=(0,1.05))
    label(c,'C','Normalised burden response')
    legend_below(c,ncol=1,y=-.27)
    d.plot(late.theta_ss,late.mean_a_end40,color=C['gray'],lw=2,ls='--',label='$t=40$')
    d.plot(late.theta_ss,late.mean_a_end,color=C['a'],lw=2,label='$t=200$')
    d.axvline(checks['lower_fold']['theta'],color=C['purple'],ls=':',lw=1.4,label='Constant-occupancy fold')
    d.set(xlim=(.1,1),ylim=(0,1),xlabel='Sustained occupancy $\\theta^*$',ylabel='Mean damage at endpoint')
    label(d,'D','Assessment horizon matters')
    legend_below(d,ncol=1,y=-.27)
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.16,wspace=.45,hspace=.83)
    return save_fig(fig,Path(outdir)/'fig5_passivation_dose_response',formats=('png','pdf','svg'))


def fig6_therapy(data_dir,outdir):
    data=Path(data_dir);win=pd.read_csv(data/'treatment_window_long.csv')
    grid=pd.read_csv(data/'combination_grid_long.csv',index_col=0)
    clr=pd.read_csv(data/'clearance_sweep_long.csv');dose=pd.read_csv(data/'dose_response_long.csv')
    base=json.loads((data/'strategy_metrics_long.json').read_text())['Baseline']['a200']
    cys=np.linspace(.25,.6,10);ups=np.logspace(-3.3,-1,10);G=grid.to_numpy()
    fig=plt.figure(figsize=(8,9.0))
    gs=fig.add_gridspec(2,2,hspace=.98,wspace=.48)
    a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,0]);d=fig.add_subplot(gs[1,1])
    cols=['#BCA43E','#E69F00','#D55E00','#A51D35']
    for onset,col in zip([0,8,16,24],cols):
        tr=pd.read_csv(data/f'window_trajectory_long_t{onset}.csv')
        a.plot(tr.time,tr.mean_a,color=col,lw=1.9,label=f'Passivation: $t_0={onset}$')
    tr=pd.read_csv(data/'window_trajectory_long_comb_t0.csv')
    a.plot(tr.time,tr.mean_a,color=C['dark'],lw=2.2,label='Combination: $t_0=0$')
    a.axvline(40,color=C['gray'],ls=':',lw=1.2)
    a.set(xlim=(0,200),ylim=(0,1),xlabel='Time (model units)',ylabel='Mean damage')
    label(a,'A','Transient and late outcomes')
    legend_below(a,ncol=1,y=-.28)
    for s in ['Passivation','Clearance','Combination']:
        sub=win[win.strategy==s].sort_values('onset');col=STRATEGY_COLORS[s]
        b.plot(sub.onset,sub.mean_a_end,color=col,lw=2,label=s)
        b.plot(sub.onset,sub.mean_a_end40,color=col,lw=1.5,ls='--')
    b.axhline(base,color=C['gray'],ls=':',lw=1)
    comb=win[win.strategy=='Combination'];tlast=comb.loc[comb.mean_a_end<=.14,'onset'].max()
    b.axvspan(0,tlast,color=C['Z'],alpha=.09,zorder=0)
    b.set(xlim=(0,32),ylim=(0,1),xlabel='Intervention onset $t_0$',ylabel='Mean damage at endpoint')
    label(b,'B','Onset and assessment horizon')
    h,l=b.get_legend_handles_labels();h += [Line2D([],[],color=C['dark'],lw=2),Line2D([],[],color=C['dark'],ls='--')]
    l += ['$t=200$ (solid)','$t=40$ (dashed)'];legend_below(b,h,l,ncol=1,y=-.28)
    im=c.pcolormesh(cys,ups,G,cmap='YlOrRd',vmin=0,vmax=1,shading='auto')
    c.contour(cys,ups,G,levels=[.15],colors=['#1A1A1A'],linewidths=1.5)
    ir=int(np.argmin(abs(ups-.012)));jr=int(np.argmin(abs(cys-.5)))
    c.scatter(cys[jr],ups[ir],marker='*',s=150,color=C['dark'],ec='white',lw=.8,zorder=4)
    c.set(yscale='log',xlabel='Oligomer clearance $c_Y$',ylabel='Passivation infusion $u_P$')
    label(c,'C','Combination dose landscape')
    # A horizontal colourbar has its own row; it cannot collide with the
    # occupancy axis or panel D. Occupancy mapping is given in the caption.
    cax=c.inset_axes([0,-.40,1,.045]);cb=fig.colorbar(im,cax=cax,orientation='horizontal')
    cb.set_label('Mean damage at $t=200$',fontsize=9.5);cb.ax.tick_params(labelsize=9)
    rows=[]
    for i,u in enumerate(ups):
        jpass=int(np.argmin(abs(dose.uP-u)));ep=1-dose.mean_a_end.iloc[jpass]/base
        for j,cy in enumerate(cys):
            jc=int(np.argmin(abs(clr.cY-cy)));ec=1-clr.mean_a_end.iloc[jc]/base
            bliss=ep+ec-ep*ec;obs=1-G[i,j]/base
            rows.append(dict(uP=u,cY=cy,bliss=bliss,observed=obs,excess=obs-bliss))
    syn=pd.read_csv(data/'synergy_matched.csv')
    syn.to_csv(Path(outdir)/'fig6_synergy_long.csv',index=False)
    d.scatter(syn.bliss,syn.observed,s=21,color=C['X'],alpha=.7,ec='white',lw=.3,zorder=3)
    d.plot([0,1],[0,1],ls='--',color=C['gray'],lw=1.3)
    ir=int(np.argmin(abs(ups-.012)));jr=int(np.argmin(abs(cys-.5)));ref=syn.iloc[ir*len(cys)+jr]
    d.scatter(ref.bliss,ref.observed,marker='*',s=150,color=C['dark'],ec='white',lw=.8,zorder=4)
    d.set(xlim=(-.03,1.03),ylim=(-.03,1.03),xlabel='Bliss-predicted benefit',ylabel='Simulated combination benefit')
    label(d,'D','Interaction relative to Bliss')
    legend_below(d,[Line2D([],[],color=C['gray'],ls='--'),Line2D([],[],marker='*',ls='',color=C['dark'],ms=9)],['Bliss independence','Nearest grid point'],ncol=1,y=-.28)
    # No floating annotation is placed over trajectories or scatter points.
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.22)
    return save_fig(fig,Path(outdir)/'fig6_therapy_design',formats=('png','pdf','svg'))
