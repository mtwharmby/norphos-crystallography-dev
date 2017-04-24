from unit_cell import UnitCell

class Crystal(object):

    def __init__(self, lattice, spacegroup=None):
        self.unit_cell = UnitCell(lattice)

