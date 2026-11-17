import copy
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

os.environ['QT_QPA_PLATFORM'] = 'offscreen'
from test_plot_creator_persistence import xrd, APP, isolate_settings
from xrd_analysis_dialog import ScientificAnalysisDialog
from xrd_profile_analysis import coupled_profile, InstrumentBroadeningModel
import numpy as np


class ScientificWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='xrd-science-')
        self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name); isolate_settings(self.base/'settings')
        home=patch.object(Path,'home',return_value=self.base); home.start();self.addCleanup(home.stop)
        self.messages=[]
        for method in ('information','warning','critical'):
            mock=patch.object(xrd.QMessageBox,method);self.messages.append(mock.start());self.addCleanup(mock.stop)
        self.window=xrd.MainWindow()
        wavelength=1.540598
        self.hkls=['111','200','220','311']
        self.centers=np.rad2deg(2*np.arcsin(wavelength*np.sqrt([3,4,8,11])/(2*3.6)))+.03
        x=np.linspace(35,95,12001)
        y=np.full_like(x,12.)
        rng=np.random.default_rng(721)
        for center in self.centers:
            y+=coupled_profile(x,center,180.,.18,.09,background_reference_deg=center)
        y+=rng.normal(0,.12,len(x))
        self.ds=xrd.Dataset(self.base/'raw.dat','Synthetic FCC','',x,y,wavelength_A=wavelength)
        self.ds.peaks=[xrd.PeakResult(int(np.argmin(abs(x-c))),c,xrd.bragg_d(c,wavelength),100,100,.24,
                                    phase='FCC test',hkl=h) for c,h in zip(self.centers,self.hkls)]
        self.window.datasets=[self.ds];self.window.sample_combo.addItem(self.ds.name)
        self.window.project_dirty=False
        self.dialog=ScientificAnalysisDialog(self.ds,self.window)

    def tearDown(self):
        self.dialog.deleteLater();self.window.project_path=None;self.window.project_dirty=False
        self.window.db.conn.close();self.window.deleteLater();APP.processEvents()

    def fit_with_instrument(self):
        self.dialog.fit_profiles()
        self.assertEqual(len(self.dialog.state['fits']),4)
        model=InstrumentBroadeningModel([30.,100.],[.08,.08],[.02,.02],
            u_gaussian_fwhm_deg=[.001,.001],u_lorentzian_fwhm_deg=[.001,.001],standard='synthetic',optics_id='test')
        self.dialog.state['instrument']=model.to_dict();self.dialog._recalculate_corrections()
        self.dialog._save_state()

    def test_fit_correct_wh_lattice_persistence_and_raw_invariance(self):
        raw=self.ds.y.copy();self.fit_with_instrument()
        first=self.dialog.state['fits']['0']
        self.assertAlmostEqual(first['fit']['two_theta_kalpha1_deg'],self.centers[0],delta=.0002)
        self.assertGreater(first['scherrer']['D_nm'],0)
        self.dialog.calculate_wh();self.dialog.calculate_lattice()
        self.assertIn('williamson_hall',self.dialog.state['results'])
        lattice=self.dialog.state['results']['cubic_lattice']
        self.assertAlmostEqual(lattice['a_A'],3.6,delta=.0001)
        self.assertAlmostEqual(lattice['zero_shift_deg'],.03,delta=.001)
        np.testing.assert_array_equal(raw,self.ds.y)
        expected=copy.deepcopy(self.ds.metadata['scientific_analysis'])
        project=self.base/'science.xrdproj';self.window._write_project(project,silent=True)
        self.dialog.accept()
        with patch.object(xrd.QFileDialog,'getOpenFileName',return_value=(str(project),'')):
            self.window.open_project()
        self.assertFalse(self.messages[2].called,str(self.messages[2].call_args))
        self.assertEqual(self.window.datasets[0].metadata['scientific_analysis'],expected)
        restored=ScientificAnalysisDialog(self.window.datasets[0],self.window)
        self.assertEqual(restored.state['fits'],expected['fits']);restored.deleteLater()
        json.dumps(expected,allow_nan=False)

    def test_raw_profile_used_even_with_processed_data(self):
        self.ds.corrected_y=np.zeros_like(self.ds.y);self.ds.baseline=np.full_like(self.ds.y,1e6)
        self.dialog.fit_profiles()
        self.assertAlmostEqual(self.dialog.state['fits']['0']['fit']['two_theta_kalpha1_deg'],self.centers[0],delta=.0002)
        self.assertNotIn('scherrer',self.dialog.state['fits']['0'])

    def test_radiation_edit_invalidates_fits(self):
        self.dialog.fit_profiles();self.dialog.wave1.setValue(1.6)
        self.assertEqual(self.dialog.state['fits'],{})

    def test_wavelength_main_invalidates_saved_science(self):
        self.ds.metadata['scientific_analysis']={'results':{'old':True}}
        xrd.apply_instrument_profile(self.ds,{'name':'Co','wavelength_A':1.79},False)
        self.assertNotIn('scientific_analysis',self.ds.metadata)

    def test_sin2psi_controls_and_input_guards(self):
        fields=self.dialog.stress_fields
        fields['d0_A'].setText('1.2');fields['xec_half_s2_per_MPa'].setText('0.000006')
        fields['xec_hkl'].setText('311');fields['xec_phase'].setText('FCC test');fields['xec_source'].setText('synthetic XEC')
        for psi in [-60,-40,-20,0,20,40,60]:
            d=1.2*(1+6e-6*250*np.sin(np.deg2rad(psi))**2)
            tt=float(np.rad2deg(2*np.arcsin(1.540598/(2*d))))
            self.dialog._append_stress(dict(psi_deg=psi,phi_deg=0,two_theta_deg=tt,u_two_theta_deg=.0005,hkl='311',phase='FCC test'))
        self.dialog.calculate_stress();result=self.dialog.state['results']['sin2psi']
        self.assertAlmostEqual(result['stress_MPa'],250,places=7)
        fields['xec_half_s2_per_MPa'].clear();self.dialog.calculate_stress()
        self.assertEqual(self.dialog.state['results']['sin2psi']['status'],'INSUFFICIENT_DATA')

    def test_failed_read_restores_existing_project_and_settings(self):
        self.window.project_dirty=False
        self.window.title_edit.setText('Keep me')
        self.window.project_dirty=False
        bad=self.base/'invalid.json'
        bad.write_text(json.dumps({'format':'XRDAnalyzerProject','datasets':[],
                                  'analysis_settings':{'min_snr':'invalid'}}),encoding='utf-8')
        before=list(self.window.datasets)
        with patch.object(xrd.QFileDialog,'getOpenFileName',return_value=(str(bad),'')):
            self.window.open_project()
        self.assertEqual(self.window.datasets,before)
        self.assertEqual(self.window.title_edit.text(),'Keep me')
        self.assertTrue(self.messages[2].called)

    def test_failed_save_prevents_project_replacement(self):
        self.window.project_dirty=True;self.window.project_path=self.base/'keep.xrdproj'
        with patch.object(xrd.QMessageBox,'question',return_value=xrd.QMessageBox.StandardButton.Save), \
             patch.object(self.window,'_write_project',side_effect=OSError('test disk failure')):
            self.assertFalse(self.window._can_replace_project())
        self.assertTrue(self.window.project_dirty);self.assertIs(self.window.datasets[0],self.ds)


if __name__=='__main__':unittest.main()
