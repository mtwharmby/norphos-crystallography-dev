import numpy as np
from enum import Enum
from collections import namedtuple

Lattice = namedtuple('Lattice', 'a, b, c, al, be, ga, V')


class UnitCell(object):

	def __init__(self):
		self.principle_axis = PrincipleAxis.B

	def update_cell(self, a, b=None, c=None, al=0, be=0, ga=0, **kwargs):
		#Reset principle_axis
		self.principle_axis = PrincipleAxis.B

		#Set correct lattice parameters
		if (b is None):
			#Cubic case, only one parameter given
			b = a
		if (c is None):
			#Tetragonal case, two parameters given, second is c
			c = b
			b = a
			self.principle_axis = PrincipleAxis.C

		#Set correct angles and convert to radians
		angles = [al, be, ga]
		set_angles = map(lambda angle: angle != 0, angles)
		if (sum(set_angles) != 3): #For non-triclinic cases
			given_angles = angles
			angles = [90]*3

			if (sum(set_angles) == 1): #Monoclinic case
				if 'principle_axis' in kwargs:
					self.principle_axis = kwargs['principle_axis']

				angles[self.principle_axis.value] = given_angles[set_angles.index(True)]
		al_r, be_r, ga_r = map(np.radians, angles)

		#Generate metric tensor and store the real space lattice
		self.metric_tensor = np.array([
			[a**2, a*b*np.cos(ga_r), a*c*np.cos(be_r)],
			[a*b*np.cos(ga_r), b**2, b*c*np.cos(al_r)],
			[a*c*np.cos(be_r), b*c*np.cos(al_r), c**2]])
		self.real_space_lattice = Lattice(a, b, c, angles[0], angles[1], angles[2], self.get_volume(True))

		#Calculate reciprocal tensor and reciprocal lattice
		self.reciprocal_metric_tensor = np.linalg.inv(self.metric_tensor)

	def cell(self):
		return self.real_space_lattice

	def get_volume(self, calc=False):
		if calc:
			return np.sqrt(np.linalg.det(self.metric_tensor))
		else:
			return self.real_space_lattice.V

class PrincipleAxis(Enum):
	__order__ = "A B C" #Needed for python 2.7
	A = 0
	B = 1 #Standard setting (be != 90)
	C = 2