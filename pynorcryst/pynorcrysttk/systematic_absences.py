import numpy as np

from crystal_exceptions import *

class SystematicAbsence(object):

	def __init__(self, pattern, condition, divisor):
		self.pattern = np.matrix(pattern)
		self.condition = np.matrix(condition)
		self.divisor = divisor

		if self.condition.size**2 != self.pattern.size:
			#We expect a 3x3 pattern and a 3x1 condition or 4x4 and 4x1 respectively
			raise MatrixException("Systematic absence condition and pattern are not compatible ")

	def is_reflection_type(self, hkl):
		hkl = np.matrix(hkl)
		#If we are dealing with Miller-Bravais 4-index reflections
		if (self.pattern.size == 16) & (hkl.size == 3):
			i_index = hkl[0,0] + hkl[0,1]
			hkl = np.insert(hkl, 2, i_index)

		product = hkl.dot(self.pattern)

		if np.array_equiv(product, hkl):
			return True
		else:
			return False

	def is_reflection_absent(self, hkl):
		if self.is_reflection_type(hkl):
			hkl = np.matrix(hkl).T
			product = self.condition.dot(hkl)
			print "product: "+str(product)
			if product % self.divisor == 0:
				return True
			else:
				return False
		else:
			return False
