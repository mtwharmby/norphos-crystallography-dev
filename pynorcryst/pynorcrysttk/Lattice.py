import numpy as np

class Lattice(object):

	def __init__(a, b, c, al, be, ga):
		self.a = a
		self.b = b
		self.c = c
		self.al = np.radians(al)
		self.be = np.radians(be)
		self.ga = np.radians(ga)


		self.metric_tensor = np.array([
			[self.a, self.a*self.b*np.cos(self.ga), self.a*self.c*np.cos(self.be)],
			[]
			[]])
