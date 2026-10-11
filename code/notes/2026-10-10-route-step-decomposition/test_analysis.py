"""Intent checks: complete accounting and an identifiable statistical question."""
import unittest

import numpy as np
import pandas as pd

from analyze import (combine_direction_rows, direction_summary, ortholog_background,
                     symmetric_difference)


class AccountingIntent(unittest.TestCase):
    def test_subthreshold_and_boundary_genes_cannot_disappear(self):
        values = np.array([-.02,-.01,.001,.019,.02,.07])
        q = direction_summary(values,.02)
        self.assertEqual([r["count"] for r in q],[1,3,2])
        self.assertNotEqual(q[1]["signed_contribution"],0)
        self.assertAlmostEqual(sum(r["signed_contribution"] for r in q),values.mean(),places=14)

    def test_quantity_and_amplitude_have_separate_meanings(self):
        reference = direction_summary([-.1,-.1,.001,.001],.02)[0]
        more = direction_summary([-.1,-.1,-.1,.001],.02)[0]
        stronger = direction_summary([-.2,-.2,.001,.001],.02)[0]
        n = symmetric_difference(reference,more)
        a = symmetric_difference(reference,stronger)
        self.assertAlmostEqual(n["amplitude_contribution"],0)
        self.assertAlmostEqual(n["number_contribution"],n["difference_signed_contribution"])
        self.assertAlmostEqual(a["number_contribution"],0)
        self.assertAlmostEqual(a["amplitude_contribution"],a["difference_signed_contribution"])

    def test_symmetric_allocation_closes_without_reference_choice(self):
        ref = direction_summary([-.12,-.04,.019,.1,.02],.02)
        comp = direction_summary([-.3,-.08,-.02,.001,.07],.02)
        differences = [symmetric_difference(a,b) for a,b in zip(ref,comp)]
        reversed_differences = [symmetric_difference(b,a) for a,b in zip(ref,comp)]
        for forward, backward in zip(differences,reversed_differences):
            for field in ["number_contribution","amplitude_contribution","difference_signed_contribution"]:
                self.assertAlmostEqual(forward[field],-backward[field])
        self.assertAlmostEqual(sum(d["number_contribution"]+d["amplitude_contribution"] for d in differences),
                               np.mean([-.3,-.08,-.02,.001,.07])-np.mean([-.12,-.04,.019,.1,.02]))

    def test_absent_class_has_no_invented_conditional_amplitude(self):
        ref = direction_summary([.001,.01],.02)[0]
        comp = direction_summary([-.1,.01],.02)[0]
        d = symmetric_difference(ref,comp)
        self.assertFalse(d["symmetric_allocation_identified"])
        self.assertTrue(np.isnan(d["number_contribution"]))
        self.assertTrue(np.isnan(d["amplitude_contribution"]))
        self.assertAlmostEqual(d["difference_signed_contribution"],-.05)

    def test_equal_route_mixture_differs_from_thresholding_route_average(self):
        a = direction_summary([-.1,.1],.02)
        b = direction_summary([.1,.1],.02)
        mixture = combine_direction_rows([a,b])
        self.assertEqual(mixture[0]["count"],.5)
        self.assertAlmostEqual(mixture[0]["fraction"],.25)
        self.assertAlmostEqual(mixture[0]["conditional_signed_mean"],-.1)
        self.assertAlmostEqual(sum(r["signed_contribution"] for r in mixture),.05)
        averaged = direction_summary([0,.1],.02)
        self.assertEqual(averaged[0]["count"],0)

    def test_background_uses_measurement_and_identity_only(self):
        orth = pd.DataFrame(dict(human_gene_id=["h1","h2"],mouse_gene_id=["m1","m2"],h_PS=[False,True]))
        rows=[]
        for route,genes in [("GSE118258",["h1","h2"]),("CIP",["m1","m2"]),("WOT",["m1"])]:
            for gene in genes:
                row=dict(route=route,gene_id=gene,gene_symbol=gene,study_id=route,PAS_number=2,
                         n_triplets=1,n_target_libraries=1,n_source_libraries=1,n_target_dates=1,
                         P_DWUI=.6,M_DWUI=.7,S_DWUI=.5,DWUI_PS=.1,DWUI_PM=-.1,DWUI_MS=.2)
                rows.append(row)
        wide=ortholog_background(pd.DataFrame(rows),orth,("GSE118258","CIP","WOT"))
        self.assertEqual(wide.human_gene_id.tolist(),["h1"])
        duplicate = pd.concat([pd.DataFrame(rows),pd.DataFrame(rows[:1])],ignore_index=True)
        with self.assertRaises(pd.errors.MergeError):
            ortholog_background(duplicate,orth,("GSE118258","CIP","WOT"))

    def test_identity_randomization_cannot_test_a_complete_mean(self):
        h = np.array([[.1,.01,.09],[-.1,.07,-.17],[.2,-.02,.22]])
        m = np.array([[.2,.02,.18],[-.3,.09,-.39],[.1,.05,.05]])
        def step_contrast(x,y):
            delta=y.mean(0)-x.mean(0)
            return delta[2]-delta[1]
        observed=step_contrast(h,m)
        shuffled=h[[2,0,1]]
        self.assertAlmostEqual(step_contrast(shuffled,m),observed)
        # Pairing can test an association, but its null is a different question.
        self.assertNotAlmostEqual(np.mean(h[:,2]*m[:,2]),np.mean(shuffled[:,2]*m[:,2]))


if __name__ == "__main__":
    unittest.main()
