from pynorcrysttk.systematic_absences import SystematicAbsence

import unittest
from nose.tools import assert_true, assert_false

class TestSystematicAbsence(unittest.TestCase):

	def setUp(self):
		pass

	def teardown(self):
		self.sys_abs = None

	def test_pattern_match(self):
		#Pattern to match all hkl
		self.sys_abs = SystematicAbsence([[1,0,0],[0,1,0],[0,0,1]], [1,1,1], 2)
		assert_true(self.sys_abs.is_reflection_type([1,0,0]), "Any reflections should be matched")
		assert_true(self.sys_abs.is_reflection_type([4,5,10]), "Any reflections should be matched")
		assert_true(self.sys_abs.is_reflection_type([8,-12,0]), "Any reflections should be matched")

		#Pattern to match hhl type reflections
		self.sys_abs = SystematicAbsence([[1,1,0], [0,0,0],[0,0,1]], [1,1,1], 2)
		assert_true(self.sys_abs.is_reflection_type([2,2,5]), "Should match hhl type reflection")
		assert_true(self.sys_abs.is_reflection_type([-2,-2,5]), "Should match hhl type reflection")
		assert_false(self.sys_abs.is_reflection_type([1,0,0]), "Should not match hkl type reflection")
		assert_false(self.sys_abs.is_reflection_type([1,7,-9]), "Should not match hkl type reflection")

		#Pattern to match h0l type reflections
		self.sys_abs = SystematicAbsence([[1,0,0],[0,0,0],[0,0,1]], [1,1,1], 2)
		assert_true(self.sys_abs.is_reflection_type([1,0,1]), "Should match h0l type reflection")
		assert_false(self.sys_abs.is_reflection_type([1,1,-1]), "Should not match hhl type reflection")

		#Pattern to match h-hl type reflections
		self.sys_abs = SystematicAbsence([[1,-1,0],[0,0,0],[0,0,1]], [1,1,1], 2)
		assert_true(self.sys_abs.is_reflection_type([1,-1,2]), "Should match h-hl type reflections")
		assert_false(self.sys_abs.is_reflection_type([1,1,2]), "Should not match hhl type reflections")

		#Pattern to match h-h0l type reflections (with hkil and hkl input)
		self.sys_abs = SystematicAbsence([[1,-1,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]], [1,1,1,1], 2)
		assert_true(self.sys_abs.is_reflection_type([1,-1,0,4]), "Should match h-h0l type reflections")
		assert_true(self.sys_abs.is_reflection_type([1,-1,4]), "Should match h-h0l type reflections quoted in hkl format")
		assert_false(self.sys_abs.is_reflection_type([1,1,-2,4]), "Should not match hhil type reflections")
		assert_false(self.sys_abs.is_reflection_type([1,1,4]), "Should not match hhil type reflections quoted in hkl format")

	def test_condition_match(self):
		#Body centred lattice absences: (hkl) h+k+l=2n
		self.sys_abs = SystematicAbsence([[1,0,0],[0,1,0],[0,0,1]], [1,1,1], 2)
		assert_true(self.sys_abs.is_reflection_absent([1,1,2]), "112 should be absent in I-centred")
		assert_false(self.sys_abs.is_reflection_absent([1,1,1]), "111 should be absent in I-centred")

		#Rhombohedral (obverse) centering
		self.sys_abs = SystematicAbsence([[1,0,0],[0,1,0],[0,0,1]], [-1,1,1],3)
		assert_true(self.sys_abs.is_reflection_absent([-2,2,2]), "-222 should be absent in R-centred")
		assert_false(self.sys_abs.is_reflection_absent([3,-2,-2]), "3-1-2 should not be absent in R-centred")

		#First test reflections of wrong type ignored, then whether h+l=2n works for h0l
		self.sys_abs = SystematicAbsence([[1,0,0],[0,0,0],[0,0,1]], [1,0,1], 2)
		assert_false(self.sys_abs.is_reflection_absent([1,1,-1]), "Reflection not h0l type")
		assert_true(self.sys_abs.is_reflection_absent([1,0,1]), "101 should be absent with an n-glide perpendicular to b-direction")
		assert_false(self.sys_abs.is_reflection_absent([1,1,1]), "111 should not be absent with an n-glide perpendicular to b-direction")