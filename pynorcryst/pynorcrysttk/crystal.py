from unit_cell import UnitCell

class Crystal(object):

    def __init__(self, lattice, spacegroup=None):
        self.unit_cell = UnitCell(lattice)

        #TODO This is not correct, this is a temporary fudge
        self.crystal_system = self.unit_cell.crystal_family
