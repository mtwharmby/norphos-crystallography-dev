from pynorcrysttk.unit_cell import UnitCell
from pynorcrysttk.unit_cell import PrincipleAxis
from pynorcrysttk.unit_cell import Lattice

import numpy as np
from collections import namedtuple

import unittest
from nose.tools import assert_equal, assert_almost_equal
from numpy.testing import assert_allclose, assert_array_almost_equal

class TestUnitCell(unittest.TestCase):

    def setUp(self):
        #self.uc = UnitCell()
        pass

    def tearDown(self):
        self.uc = None

    def test_metric_tensor(self):
        '''
        metric tensor = [[a**2, a*b*cos(ga), a*c*cos(be)],
                         [a*b*cos(ga), b**2, b*c*cos(al)],
                         [a*c*cos(be), b*c*cos(al), c**2]]
        '''
        #Cubic a = 3
        self.uc = UnitCell(Lattice(3))
        fake_tensor = np.array([
            [9, 0, 0],
            [0, 9, 0],
            [0, 0, 9]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Cubic tensor incorrectly calculated')

        #Rhombohedral a = 3; al = 60
        self.uc.update_cell(Lattice(3, al=60))
        fake_tensor = np.array([
            [9, 4.5, 4.5],
            [4.5, 9, 4.5],
            [4.5, 4.5, 9]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Rhombohedral tensor incorrectly calculated')

        #Hexagonal a = 5; c = 2; ga = 120
        self.uc.update_cell(Lattice(5,c=2,ga=120))
        fake_tensor = np.array([
            [25, -12.5, 0],
            [-12.5, 25, 0],
            [0, 0, 4]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Hexagonal tensor incorrectly calculated')

        #Tetragonal a = b = 2; c = 5
        self.uc.update_cell(Lattice(2, 5))
        fake_tensor = np.array([
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 25]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Tetragonal tensor incorrectly calculated')

        #a = 5; b = 3; c = 2; angles 90
        self.uc.update_cell(Lattice(5, 3, 2))
        fake_tensor = np.array([
            [25, 0, 0],
            [0, 9, 0],
            [0, 0, 4]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Orthorhombic tensor incorrectly calculated')

        #a = 2; b = 3; c = 5; be = 30 al = ga = 90
        self.uc.update_cell(Lattice(2, 3, 5, be=30))
        fake_tensor = np.array([
            [4, 0, 5*np.sqrt(3)],
            [0, 9, 0],
            [5*np.sqrt(3), 0, 25]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Monoclinic tensor incorrectly calculated')

        #a = 3; b = 5; c = 2; al = 30deg; be = 45 deg; ga = 60 deg
        self.uc.update_cell(Lattice(3, 5 ,2, 30, 45, 60))
        fake_tensor = np.array([
            [9, 7.5, 6/np.sqrt(2)],
            [7.5, 25, 5*np.sqrt(3)],
            [6/np.sqrt(2), 5*np.sqrt(3), 4]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Triclinic tensor incorrectly calculated')

    def test_specifying_principle_axis(self):
        self.uc = UnitCell(Lattice(2, 3, 5, 30, principle_axis=PrincipleAxis.A))
        fake_tensor = np.array([
            [4, 0, 0],
            [0, 9, 7.5*np.sqrt(3)],
            [0, 7.5*np.sqrt(3), 25]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Monoclinic (principle alpha) tensor incorrectly calculated')
        
        self.uc.update_cell(Lattice(2, 3, 5, 30, principle_axis=PrincipleAxis.C))
        fake_tensor = np.array([
            [4, 3*np.sqrt(3), 0],
            [3*np.sqrt(3), 9, 0],
            [0, 0, 25]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Monoclinic (principle gamma) tensor incorrectly calculated')

    def test_return_cell(self):
        self.uc = UnitCell(Lattice(3, 5 ,2, 30, 45, 60))

        assert_equal(3.0, self.uc.lattice.a)
        assert_equal(5.0, self.uc.lattice.b)
        assert_equal(2.0, self.uc.lattice.c)
        assert_almost_equal(30.0, self.uc.lattice.al, places=8)
        assert_almost_equal(45.0, self.uc.lattice.be, places=8)
        assert_almost_equal(60.0, self.uc.lattice.ga, places=8)

    def test_return_cell_volume(self):
        self.uc = UnitCell(Lattice(2,3,5))
        assert_almost_equal(30.0, self.uc.volume)

    def test_reciprocal_metric_tensor(self):
        #Try orthorhombic case first
        ortho_lat = Lattice(2,3,5,90,90,90)
        self.uc = UnitCell(ortho_lat)

        fake_tensor = np.array([
            [ortho_lat.a**2, 0, 0],
            [0, ortho_lat.b**2, 0],
            [0, 0, ortho_lat.c**2]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Orthorhombic tensor incorrectly calculated')

        fake_recip_tensor = np.array([
            [ortho_lat.b**2*ortho_lat.c**2/self.uc.volume**2, 0, 0],
            [0, ortho_lat.a**2*ortho_lat.c**2/self.uc.volume**2, 0],
            [0, 0, ortho_lat.a**2*ortho_lat.b**2/self.uc.volume**2]])
        assert_array_almost_equal(fake_recip_tensor, self.uc.reciprocal_metric_tensor, err_msg='Orthorhombic reciprocal tensor incorrectly calculated')

        #Try more complicated triclinic case: This is Anorthoclase
        anorthoclase_lat = Lattice(8.28,12.97,7.15,91.05,116.26,90.15)
        al_r, be_r, ga_r = map(np.radians, [anorthoclase_lat.al, anorthoclase_lat.be, anorthoclase_lat.ga])
        self.uc.update_cell(anorthoclase_lat)

        fake_tensor = np.array([
            [anorthoclase_lat.a**2, anorthoclase_lat.a*anorthoclase_lat.b*np.cos(ga_r), anorthoclase_lat.a*anorthoclase_lat.c*np.cos(be_r)],
            [anorthoclase_lat.a*anorthoclase_lat.b*np.cos(ga_r), anorthoclase_lat.b**2, anorthoclase_lat.b*anorthoclase_lat.c*np.cos(al_r)],
            [anorthoclase_lat.a*anorthoclase_lat.c*np.cos(be_r), anorthoclase_lat.b*anorthoclase_lat.c*np.cos(al_r), anorthoclase_lat.c**2]])
        assert_array_almost_equal(fake_tensor, self.uc.metric_tensor, err_msg='Anorthoclase tensor incorrectly calculated')

        #This is the complete expression for the reciprocal metric tensor,
        #(The Reciprocal Lattice, A. Authier, IUCr Pamphlet 4 (1981))
        g12 = anorthoclase_lat.a*anorthoclase_lat.b*anorthoclase_lat.c**2*(np.cos(al_r)*np.cos(be_r)-np.cos(ga_r))/self.uc.volume**2
        g13 = anorthoclase_lat.a*anorthoclase_lat.b**2*anorthoclase_lat.c*(np.cos(al_r)*np.cos(ga_r)-np.cos(be_r))/self.uc.volume**2
        g23 = anorthoclase_lat.a**2*anorthoclase_lat.b*anorthoclase_lat.c*(np.cos(be_r)*np.cos(ga_r)-np.cos(al_r))/self.uc.volume**2
        fake_recip_tensor = np.array([
            [anorthoclase_lat.b**2*anorthoclase_lat.c**2*np.sin(al_r)**2/self.uc.volume**2, g12, g13],
            [g12, anorthoclase_lat.a**2*anorthoclase_lat.c**2*np.sin(be_r)**2/self.uc.volume**2, g23],
            [g13, g23, anorthoclase_lat.a**2*anorthoclase_lat.b**2*np.sin(ga_r)**2/self.uc.volume**2]])
        assert_array_almost_equal(fake_recip_tensor, self.uc.reciprocal_metric_tensor, err_msg='Anorthoclase reciprocal tensor incorrectly calculated')

    def test_vector_magnitude_calc(self):
        ortho_lat = Lattice(2,3,5,90,90,90)
        self.uc = UnitCell(ortho_lat)
        assert_almost_equal(2, self.uc.find_vector_magnitude([1,0,0])) 
        assert_almost_equal(3, self.uc.find_vector_magnitude([0,1,0])) 
        assert_almost_equal(np.sqrt(4+9), self.uc.find_vector_magnitude([1,1,0]))
        assert_almost_equal(np.sqrt(4+9+25), self.uc.find_vector_magnitude([1,1,1]))

        anorthoclase_lat = Lattice(8.28,12.97,7.15,91.05,116.26,90.15)
        self.uc.update_cell(anorthoclase_lat)
        dist = self.uc.find_vector_magnitude([1,1,1])
        assert_almost_equal(15.21688, dist, places=5)

    def test_lattice_plane_dspacing(self):
        #a, b, c, al, be, ga = [8.28,12.97,7.15,91.05,116.26,90.15]
        ortho_lat = Lattice(2,3,5,90,90,90)
        self.uc = UnitCell(ortho_lat)
        assert_almost_equal(2, self.uc.find_plane_dspacing([1,0,0]))

        anorthoclase_lat = Lattice(8.28,12.97,7.15,91.05,116.26,90.15)
        self.uc.update_cell(anorthoclase_lat)
        #d-spacing of planes (calculated with PowderCell)
        #(100) 7.42494; (010) 12.9669; (001) 6.41057; (110) 6.41039; (111) 3.84084
        assert_almost_equal(7.42494, self.uc.find_plane_dspacing([1,0,0]), places=5)
        assert_almost_equal(12.96689, self.uc.find_plane_dspacing([0,1,0]), places=5)
        assert_almost_equal(6.41057, self.uc.find_plane_dspacing([0,0,1]), places=5)
        assert_almost_equal(6.41039, self.uc.find_plane_dspacing([1,1,0]), places=5)
        assert_almost_equal(3.84084, self.uc.find_plane_dspacing([1,1,1]), places=5)

    def test_reciprocal_lattice_determination(self):
        def calc_reciprocal_length(a, b, ga, vol):
            return a*b*np.sin(np.radians(ga))/vol

        def calc_reciprocal_angle(al, be, ga):
            al_r, be_r, ga_r = map(np.radians, [al, be, ga])
            return np.degrees(np.arccos((np.cos(al_r) * np.cos(be_r) - np.cos(ga_r))/abs(np.sin(al_r) * np.sin(be_r))))

        def generate_reciprocal_lat(real_lat):
            a_star = calc_reciprocal_length(real_lat.b, real_lat.c, real_lat.al, self.uc.volume)
            b_star = calc_reciprocal_length(real_lat.a, real_lat.c, real_lat.be, self.uc.volume)
            c_star = calc_reciprocal_length(real_lat.a, real_lat.b, real_lat.ga, self.uc.volume)
    
            al_star = calc_reciprocal_angle(real_lat.be, real_lat.ga, real_lat.al)
            be_star = calc_reciprocal_angle(real_lat.al, real_lat.ga, real_lat.be)
            ga_star = calc_reciprocal_angle(real_lat.al, real_lat.be, real_lat.ga)

            return Lattice(a_star, b_star, c_star, al_star, be_star, ga_star)

        ortho_lat = Lattice(2,3,5,90,90,90)
        self.uc = UnitCell(ortho_lat)
        
        fake_recip_lat = generate_reciprocal_lat(ortho_lat)
        assert_almost_equal(fake_recip_lat.a, self.uc.reciprocal_lattice.a, places=10)
        assert_almost_equal(fake_recip_lat.b, self.uc.reciprocal_lattice.b, places=10)
        assert_almost_equal(fake_recip_lat.c, self.uc.reciprocal_lattice.c, places=10)
        assert_almost_equal(fake_recip_lat.al, self.uc.reciprocal_lattice.al, places=10)
        assert_almost_equal(fake_recip_lat.be, self.uc.reciprocal_lattice.be, places=10)
        assert_almost_equal(fake_recip_lat.ga, self.uc.reciprocal_lattice.ga, places=10)

        anorthoclase_lat = Lattice(8.28,12.97,7.15,91.05,116.26,90.15)
        self.uc.update_cell(anorthoclase_lat)

        fake_recip_lat = generate_reciprocal_lat(anorthoclase_lat)
        assert_almost_equal(fake_recip_lat.a, self.uc.reciprocal_lattice.a, places=10)
        assert_almost_equal(fake_recip_lat.b, self.uc.reciprocal_lattice.b, places=10)
        assert_almost_equal(fake_recip_lat.c, self.uc.reciprocal_lattice.c, places=10)
        assert_almost_equal(fake_recip_lat.al, self.uc.reciprocal_lattice.al, places=10)
        assert_almost_equal(fake_recip_lat.be, self.uc.reciprocal_lattice.be, places=10)
        assert_almost_equal(fake_recip_lat.ga, self.uc.reciprocal_lattice.ga, places=10)

    def test_direction_cosines_calc(self):
        #Rows in direction cosines matrix should be such that sum of squares of 
        #each element in a row totals 1; orthonormalisation matrix is direction
        #cosines matrix multiplies by magnitude of cell lengths.
        lengths = [2,3,5]
        ortho_lat = Lattice(lengths[0],lengths[1],lengths[2],90,90,90)
        self.uc = UnitCell(ortho_lat)

        sum_sq_orthon_mat = np.sum(np.power(self.uc.orthonormalisation_matrix, 2), axis=1)
        for i in range(3):
            assert_equal(np.matrix(lengths[i]**2), sum_sq_orthon_mat[i])

        #Try a triclinic case
        lengths = [8.28,12.97,7.15]
        anorthoclase_lat = Lattice(lengths[0],lengths[1],lengths[2],91.05,116.26,90.15)
        self.uc.update_cell(anorthoclase_lat)

        sum_sq_orthon_mat = np.sum(np.power(self.uc.orthonormalisation_matrix, 2), axis=1)
        for i in range(3):
            assert_equal(np.matrix(lengths[i]**2), sum_sq_orthon_mat[i])

    def test_cartesian_distance_calc(self):
        ortho_lat = Lattice(2,3,5,90,90,90)
        self.uc = UnitCell(ortho_lat)
        coords = self.uc.find_cartesian_coordinates([2,1,3])
        assert_allclose(np.matrix('4.;3.;15.'), coords)

        #Finding distance across body diagonal of triclinic cell
        anorthoclase_lat = Lattice(8.28,12.97,7.15,91.05,116.26,90.15)
        self.uc.update_cell(anorthoclase_lat)
        coords = self.uc.find_vector_magnitude([1,1,1])
        dist = np.sqrt(np.sum(np.power(coords, 2)))
        assert_almost_equal(15.21688, dist, places=5)
