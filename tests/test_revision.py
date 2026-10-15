"""Regression checks for numerical and manuscript-integrity corrections."""
import unittest
import sys,json,re
from pathlib import Path
import numpy as np
from scipy.stats import rankdata,t
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from ad_graph_model.model import DEFAULT_PARAMETERS,diffusion
from ad_graph_model.bifurcation import fast_equilibrium_q1,qss_damage_rate,node_f4,node_jac4
from ad_graph_model.sensitivity import prcc_matrix

class RevisionTests(unittest.TestCase):
    def test_fast_equilibria_have_no_spurious_gaps(self):
        p=DEFAULT_PARAMETERS;kd=p['koff']/p['kon']
        for th in [0,.55,.65,.95,.99]:
            P=kd*th/(1-th)
            for a in np.linspace(0,1,41):
                xyz=fast_equilibrium_q1(a,p,P)
                self.assertTrue(np.isfinite(xyz).all());self.assertTrue((xyz>=0).all())
                self.assertLess(abs(node_f4(np.r_[xyz,a],p,P)[:3]).max(),1e-8)
                self.assertTrue(np.isfinite(qss_damage_rate(a,p,P)))

    def test_fold_and_bistable_equilibria(self):
        p=DEFAULT_PARAMETERS;d=json.loads((ROOT/'results_v3/data/revision_checks.json').read_text())
        self.assertAlmostEqual(d['lower_fold']['theta'],.9282103682,places=8)
        self.assertAlmostEqual(d['upper_fold']['theta'],.9901983187,places=8)
        P=p['koff']/p['kon']*.95/.05
        kinds=[]
        for a in d['qss_roots']['strong']:
            xyz=fast_equilibrium_q1(a,p,P);y=np.r_[xyz,a]
            self.assertLess(abs(node_f4(y,p,P)).max(),1e-8)
            kinds.append(np.linalg.eigvals(node_jac4(y,p,P)).real.max()<0)
        self.assertEqual(kinds,[True,False,True])

    def test_transport_conserves_weighted_not_arithmetic_sum(self):
        A=np.array([[0,1,4],[1,0,2],[4,2,0]],float);s=A.sum(1);W=A/s[:,None]
        v=diffusion(np.array([0.,2.,7.]),W)
        self.assertAlmostEqual(s@v,0,places=12)
        self.assertGreater(abs(v.sum()),.1)
        np.testing.assert_allclose(diffusion(np.ones(3),W),0,atol=1e-14)

    def test_prcc_uses_number_of_controlled_covariates(self):
        rng=np.random.default_rng(5);X=rng.normal(size=(80,3));Y=(X[:,0]+rng.normal(size=80))[:,None]
        r,p=prcc_matrix(X,Y);df=80-2-2
        expected=2*t.sf(abs(r[0,0])*np.sqrt(df/(1-r[0,0]**2)),df)
        self.assertAlmostEqual(p[0,0],expected,places=13)

    def test_references_and_manuscript_assets(self):
        text=(ROOT/'manuscript/final_draft.md').read_text(encoding='utf-8')
        refs=json.loads((ROOT/'manuscript/reference_audit.json').read_text(encoding='utf-8'))
        self.assertEqual(len(refs),37);self.assertEqual(len({x['doi'].lower() for x in refs}),37)
        for ref in refs:
            self.assertTrue(ref['first_author']);self.assertTrue(ref['volume']);self.assertTrue(ref['pages'])
        for path in re.findall(r'!\[.*?\]\((.*?)\)',text):self.assertTrue((ROOT/'manuscript'/path).is_file(),path)
        self.assertNotIn('{{',text)

    def test_bliss_controls_are_matched(self):
        import pandas as pd
        s=pd.read_csv(ROOT/'results_v3/data/synergy_matched.csv')
        np.testing.assert_allclose(s.bliss,s.E_P+s.E_C-s.E_P*s.E_C,atol=1e-12)
        np.testing.assert_allclose(s.excess,s.observed-s.bliss,atol=1e-12)
        np.testing.assert_allclose(s.loc[np.isclose(s.cY,.25),'excess'],0,atol=1e-6)

if __name__=='__main__':unittest.main()
