import numpy as np
from enum import Enum
from collections import namedtuple

from crystal_exceptions import *


class UnitCell(object):

	def __init__(self, lattice=None):
		self.volume = None
		#Angles in radians for convenience
		self.al_r, self.be_r, self.ga_r = None, None, None
		#Tensors
		self.real, self.reciprocal = None, None
	
		#Set the lattices up values
		self.update_cell(lattice)

	def update_cell(self, lattice):
		self.real = self.__evaluate_lattice(lattice)
		self.al_r, self.be_r, self.ga_r = map(np.radians, [self.real.al, 
														   self.real.be,
														   self.real.ga])
		self.metric_tensor = self.__determine_metric_tensor()
		self.volume = self.__calculate_cell_volume()
		self.reciprocal_metric_tensor = np.linalg.inv(self.metric_tensor)

	def __evaluate_lattice(self, lattice):
		a, b, c = lattice.a, lattice.b, lattice.c
		angles = [lattice.al, lattice.be, lattice.ga]
		p_axis = PrincipleAxis.NONE

		#Set any angles with value None to 90 & find out which angles are identical
		for i in range(3):
			if angles[i] is None:
				angles[i] = 90
		angles_compared = [angles[0]==angles[1],
						   angles[0]==angles[2],
						   angles[1]==angles[2]]

		#Analyse the user input
		if False not in angles_compared:
			if (a == b and a == c) or (b is None and c is None):
				b = c = a
				if angles[0] != 90:
					#Rhombohedral
					pass
				else:
					#Cubic
					pass
			elif (a == b and a != c) or (b is None and c is not None):
				#Tetragonal
				b = a
				p_axis = PrincipleAxis.C
			elif (a != b and c is None):
				#Tetragonal (b -> c)
				c = b
				b = a
				p_axis = PrincipleAxis.C
			else:
				#Orthorhombic (a != b != c)
				if (lattice.b is None or lattice.c is None):
					raise LatticeException("Orthorhombic requires all three lengths to be given")
		elif True not in angles_compared:
			#Triclinic
			if (b is None or c is None) and None not in angles:
				raise LatticeException("Triclinic requires three lengths three angles to be given")
		else:
			#Monoclinic
			if lattice.principle_axis != PrincipleAxis.NONE:
				p_axis = lattice.principle_axis
				old_angles = angles
				angles = [90]*3

				if angles_compared[0]:
					#ga different
					angles[p_axis.value] = old_angles[2]
				elif angles_compared[1]:
					#be different
					angles[p_axis.value] = old_angles[1]
				else:
					#al different
					angles[p_axis.value] = old_angles[0]
			else:
				if (lattice.b is None or lattice.c is None):
					raise LatticeException("Monoclinic requires all three lengths to be given")
				if angles_compared[0]:
					#ga different
					p_axis = PrincipleAxis.C
				elif angles_compared[1]:
					#be different
					p_axis = PrincipleAxis.B
				else:
					#al different
					p_axis = PrincipleAxis.A

		return Lattice(a, b, c, angles[0], angles[1], angles[2], p_axis)

	def __determine_metric_tensor(self):
		def __calc_offaxis(a, b, angle):
			result = a * b * np.cos(angle)
			if abs(result) < 1e-10:
				return 0.0
			return result

		#Generate metric tensor and store the real space lattice
		p00 = self.real.a*self.real.a
		p11 = self.real.b*self.real.b
		p22 = self.real.c*self.real.c
		#Off axis
		p01 = __calc_offaxis(self.real.a, self.real.b, self.ga_r)
		p02 = __calc_offaxis(self.real.a, self.real.c, self.be_r)
		p12 = __calc_offaxis(self.real.b, self.real.c, self.al_r)

		return np.array([
			[p00, p01, p02],
			[p01, p11, p12],
			[p02, p12, p22]])

	def __calculate_cell_volume(self):
		return np.sqrt(np.linalg.det(self.metric_tensor))

	def cell(self):
		return self.real

	def find_vector_magnitude(self, vector, cosines_matrix=None):
		#Default to use the metric_tensor
		if (cosines_matrix is None):
			cosines_matrix = self.metric_tensor

		#Check we've been given a real vector
		vector = np.array(vector)
		if vector.shape != (3L,):
			raise MatrixException("Given vector is not a 3x1 matrix")

		#Use a given metric tensor (cosines_matrix) to calculate vector length: |r| = sqrt(r^t . G . r)
		fix_zeros = np.vectorize(lambda x: 0 if (abs(x) <= 1e-10) else x)
		fixed_product = fix_zeros(np.dot(np.dot(np.transpose(vector), cosines_matrix),vector))
		return np.linalg.norm(np.sqrt(fixed_product))

	def find_plane_dspacing(self, plane_indices):
		one_on_d = self.find_vector_magnitude(plane_indices, self.reciprocal_metric_tensor)
		return 1/one_on_d

class PrincipleAxis(Enum):
	__order__ = "A B C" #Needed for python 2.7
	A = 0
	B = 1 #Standard setting (be != 90)
	C = 2
	NONE = -1

Lattice = namedtuple('Lattice', 'a, b, c, al, be, ga, principle_axis')
Lattice.__new__.__defaults__ = (None, None, None, None, None, PrincipleAxis.NONE)