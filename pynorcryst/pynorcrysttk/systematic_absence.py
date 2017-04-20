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

		if (self.pattern.size == 9) & (self.pattern.shape != (3,3)):
			self.pattern = self.pattern.reshape((3,3))
		elif (self.pattern.size == 16) & (self.pattern.shape != (4,4)):
			self.pattern = self.pattern.reshape((4,4))
		else:
			raise MatrixException("Expecting square matrix of size 3x3 or 4x4. Got: "+str(self.pattern.size))

	def is_reflection_type(self, hkl):
		"""
		Determine whether the given reflection indices (hkl or hkil) will be 
		affected by this systematic absence. Comparison is made with the 
		configured pattern.

		Parameters
		----------
		hkl : array_like
			  A list or array containing indices hkl or hkil.

		Returns
		-------
		out : bool
			  True if given hkl matches the pattern.
		"""
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
		"""
		Determine whether the reflection with given indices (hkl or hkil) will 
		be systematically absent. Checks first whether reflection is the right 
		type to be affected by this systematic absence.

		Parameters
		----------
		hkl : array_like
			  A list or array containing indices hkl or hkil.

		Returns
		-------
		out : bool
			  True if the reflection is systematically absent.

		"""
		if self.is_reflection_type(hkl):
			hkl = np.matrix(hkl).T
			product = self.condition.dot(hkl)
			if product % self.divisor == 0:
				return True
			else:
				return False
		else:
			return False

	def __str__(self):
		def convert_to_hkl(row, report_none=False):
			if row.size == 3:
				labels = ["h","k","l"]
			else:
				labels = ["h","k","i","l"]
			
			hkl_str = ""
			for i in range(row.size):
				if row[0,i] != 0:
					if row[0,i] < 0:
						hkl_str = hkl_str.rstrip("+")+"-"
					elif abs(row[0,i]) > 1:
						hkl_str = str(row[0,i])
					hkl_str += labels[i]+"+"

			if report_none & (hkl_str == ""):
				return str(0)
			else:
				return hkl_str.rstrip("+")

		absence_str = ""
		for i in range(self.pattern[0].size):
			absence_str += convert_to_hkl(self.pattern.T[i], True)
		absence_str += ": "+convert_to_hkl(self.condition)+"="+str(self.divisor)+"n"
		return absence_str
