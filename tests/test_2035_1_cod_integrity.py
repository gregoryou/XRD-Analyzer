"""Safe tests requiring only Python stdlib; GUI is tested on Windows."""
from __future__ import annotations
import sys
import unittest
from dataclasses import dataclass, field
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from xrd_kalpha import (predicted_kalpha2_position,find_resolved_kalpha2,
    dataset_kalpha2_candidates,traceability_rows)

@dataclass
class P:
    two_theta:float
    height:float
    fwhm_deg:float=.18
    d_A:float=1.0
    rel_intensity:float=50.0
    phase:str='Si cubic'
    hkl:str='111'
    phase_source:str='cod'

@dataclass
class D:
    peaks:list=field(default_factory=list)
    wavelength_A:float=1.540598
    metadata:dict=field(default_factory=lambda:{
        'wavelength_kalpha1_A':1.540598,'wavelength_kalpha2_A':1.544426,
        'manual_cod_peak_assignments':[{'peak_two_theta':94.9259,'cod_id':2104737,
          'space_group':'F d -3 m','hkl':'511','ref_two_theta':94.9521,'delta_2theta':0.0263}]})

class TestKalpha(unittest.TestCase):
 def test_three_resolved_si_pairs(self):
    p=[]
    for alpha1,alpha2,i1,i2 in ((94.9259,95.2601,433.99,293.25),
        (106.6906,107.0917,340.62,219.06),(114.077,114.5115,473.36,313.72)):
        p.extend([P(alpha1,i1),P(alpha2,i2)])
    d=D(peaks=p)
    result=dataset_kalpha2_candidates(d)
    self.assertEqual(set(result),{1,3,5})
    self.assertEqual([result[i]['parent_index'] for i in (1,3,5)],[0,2,4])
    hdr,rows=traceability_rows(d)
    self.assertEqual(len(rows),6)
    self.assertEqual(rows[0][6],2104737)
    self.assertEqual(rows[1][6],'')
    self.assertIn('Kalpha2', rows[1][11])
 def test_missing_second_line_never_assumed(self):
    d=D(peaks=[P(94.9259,440),P(95.2601,280)])
    del d.metadata['wavelength_kalpha2_A']
    self.assertFalse(dataset_kalpha2_candidates(d))
 def test_changed_primary_wavelength_disables_pair_flagging(self):
    d=D(peaks=[P(94.9259,440),P(95.2601,280)],wavelength_A=1.788965)
    self.assertFalse(dataset_kalpha2_candidates(d))
 def test_false_positive_not_flagged(self):
    d=D(peaks=[P(94.9259,440),P(95.1,800)])
    self.assertFalse(dataset_kalpha2_candidates(d))
 def test_prediction_physics(self):
    v=predicted_kalpha2_position(94.9259,1.540598,1.544426)
    self.assertAlmostEqual(v,95.2366863,places=5)
 def test_no_peak_deletions(self):
    d=D(peaks=[P(94.9259,440),P(95.2601,280)])
    original=[p.two_theta for p in d.peaks]
    traceability_rows(d)
    self.assertEqual(original,[p.two_theta for p in d.peaks])

if __name__=='__main__': unittest.main()
