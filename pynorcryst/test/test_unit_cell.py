from pynorcrysttk.unit_cell import UnitCell
from pynorcrysttk.unit_cell import PrincipleAxis

import numpy as np
from collections import namedtuple

import unittest
from nose.tools import assert_equal, assert_almost_equal
from numpy.testing import assert_array_almost_equal

class TestUnitCell(unittest.TestCase):

	def setUp(self):
		self.uc = UnitCell()

	def teardownclass(self):
		self.uc = None

	def test_metric_tensor(self):
		#TODO Make more similar to Java test. Add Hexagonal & rhombohedral too.

		#a = 3; angles 90
		self.uc.update_cell(3)
		fake_tensor = np.array([
			[9, 0, 0],
			[0, 9, 0],
			[0, 0, 9]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Cubic tensor incorrectly calculated')

		#a = b = 2; c = 5 angles 90
		self.uc.update_cell(2, 5)
		fake_tensor = np.array([
			[4, 0, 0],
			[0, 4, 0],
			[0, 0, 25]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Tetragonal tensor incorrectly calculated')

		#a = 5; b = 3; c = 2; angles 90
		self.uc.update_cell(5, 3, 2)
		fake_tensor = np.array([
			[25, 0, 0],
			[0, 9, 0],
			[0, 0, 4]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Orthorhombic tensor incorrectly calculated')

		#a = 2; b = 3; c = 5; be = 30 al = ga = 90
		self.uc.update_cell(2, 3, 5, 30)
		fake_tensor = np.array([
			[4, 0, 5*np.sqrt(3)],
			[0, 9, 0],
			[5*np.sqrt(3), 0, 25]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Monoclinic tensor incorrectly calculated')

		#a = 3; b = 5; c = 2; al = 30deg; be = 45 deg; ga = 60 deg
		self.uc.update_cell(3, 5 ,2, 30, 45, 60)
		fake_tensor = np.array([
			[9, 7.5, 6/np.sqrt(2)],
			[7.5, 25, 5*np.sqrt(3)],
			[6/np.sqrt(2), 5*np.sqrt(3), 4]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Triclinic tensor incorrectly calculated')

	def test_specifying_principle_axis(self):
		self.uc.update_cell(2, 3, 5, 30, principle_axis=PrincipleAxis.A)
		fake_tensor = np.array([
			[4, 0, 0],
			[0, 9, 7.5*np.sqrt(3)],
			[0, 7.5*np.sqrt(3), 25]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Monoclinic (principle alpha) tensor incorrectly calculated')
		
		self.uc.update_cell(2, 3, 5, 30, principle_axis=PrincipleAxis.C)
		fake_tensor = np.array([
			[4, 3*np.sqrt(3), 0],
			[3*np.sqrt(3), 9, 0],
			[0, 0, 25]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Monoclinic (principle gamma) tensor incorrectly calculated')

	def test_return_cell(self):
		self.uc.update_cell(3, 5 ,2, 30, 45, 60)

		cell = self.uc.cell()
		assert_equal(3.0, cell.a)
		assert_equal(5.0, cell.b)
		assert_equal(2.0, cell.c)
		assert_almost_equal(30.0, cell.al, places=8)
		assert_almost_equal(45.0, cell.be, places=8)
		assert_almost_equal(60.0, cell.ga, places=8)

	def test_return_cell_volume(self):
		self.uc.update_cell(2,3,5)
		vol = self.uc.get_volume()
		assert_almost_equal(30.0, vol)

	def test_reciprocal_metric_tensor(self):
		#Try orthorhombic case first
		a, b, c, al, be, ga = [2,3,5,90,90,90]
		self.uc.update_cell(a, b, c, al, be, ga)

		fake_tensor = np.array([
			[a**2, 0, 0],
			[0, b**2, 0],
			[0, 0, c**2]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Orthorhombic tensor incorrectly calculated')

		fake_recip_tensor = np.array([
			[b**2*c**2/self.uc.get_volume()**2, 0, 0],
			[0, a**2*c**2/self.uc.get_volume()**2, 0],
			[0, 0, a**2*b**2/self.uc.get_volume()**2]])
		assert_array_almost_equal(fake_recip_tensor, self.uc.reciprocal_metric_tensor, err_msg='Orthorhombic reciprocal tensor incorrectly calculated')

		#Try more complicated triclinic case: This is Anorthoclase
		a, b, c, al, be, ga = [8.28,12.97,7.15,91.05,116.26,90.15]
		al_r, be_r, ga_r = map(np.radians, [al, be, ga])
		self.uc.update_cell(a, b, c, al, be, ga)

		fake_tensor = np.array([
			[a**2, a*b*np.cos(ga_r), a*c*np.cos(be_r)],
			[a*b*np.cos(ga_r), b**2, b*c*np.cos(al_r)],
			[a*c*np.cos(be_r), b*c*np.cos(al_r), c**2]])
		assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Anorthoclase tensor incorrectly calculated')

		#This is the complete expression for the reciprocal metric tensor,
		#(The Reciprocal Lattice, A. Authier, IUCr Pamphlet 4 (1981))
		g12 = a*b*c**2*(np.cos(al_r)*np.cos(be_r)-np.cos(ga_r))/self.uc.get_volume()**2
		g13 = a*b**2*c*(np.cos(al_r)*np.cos(ga_r)-np.cos(be_r))/self.uc.get_volume()**2
		g23 = a**2*b*c*(np.cos(be_r)*np.cos(ga_r)-np.cos(al_r))/self.uc.get_volume()**2
		fake_recip_tensor = np.array([
			[b**2*c**2*np.sin(al_r)**2/self.uc.get_volume()**2, g12, g13],
			[g12, a**2*c**2*np.sin(be_r)**2/self.uc.get_volume()**2, g23],
			[g13, g23, a**2*b**2*np.sin(ga_r)**2/self.uc.get_volume()**2]])
		assert_array_almost_equal(fake_recip_tensor, self.uc.reciprocal_metric_tensor, err_msg='Anorthoclase reciprocal tensor incorrectly calculated')

	def test_vector_magnitude_calc(self):
		a, b ,c, al, be, ga = [2, 3, 5, 90, 90, 90]
		self.uc.update_cell(a, b, c, al, be, ga)

		assert_almost_equal(2, self.uc.find_vector_magnitude([1,0,0])) 
		assert_almost_equal(3, self.uc.find_vector_magnitude([0,1,0])) 
		assert_almost_equal(np.sqrt(4+9), self.uc.find_vector_magnitude([1,1,0]))
		assert_almost_equal(np.sqrt(4+9+25), self.uc.find_vector_magnitude([1,1,1]))

		a, b, c, al, be, ga = [8.28,12.97,7.15,91.05,116.26,90.15]
		self.uc.update_cell(a, b, c, al, be, ga)
		dist = self.uc.find_vector_magnitude([1,1,1])
		assert_almost_equal(15.21688, dist, places=5)

	def test_lattice_plane_dspacing(self):
		#a, b, c, al, be, ga = [8.28,12.97,7.15,91.05,116.26,90.15]
		a, b ,c, al, be, ga = [2, 3, 5, 90, 90, 90]
		self.uc.update_cell(a, b, c, al, be, ga)

		assert_almost_equal(2, self.uc.find_plane_dspacing([1,0,0]))

		a, b, c, al, be, ga = [8.28,12.97,7.15,91.05,116.26,90.15]
		self.uc.update_cell(a, b, c, al, be, ga)
		
		#d-spacing of planes (calculated with PowderCell)
		#(100) 7.42494; (010) 12.9669; (001) 6.41057; (110) 6.41039; (111) 3.84084
		assert_almost_equal(7.42494, self.uc.find_plane_dspacing([1,0,0]), places=5)
		assert_almost_equal(12.96689, self.uc.find_plane_dspacing([0,1,0]), places=5)
		assert_almost_equal(6.41057, self.uc.find_plane_dspacing([0,0,1]), places=5)
		assert_almost_equal(6.41039, self.uc.find_plane_dspacing([1,1,0]), places=5)
		assert_almost_equal(3.84084, self.uc.find_plane_dspacing([1,1,1]), places=5)