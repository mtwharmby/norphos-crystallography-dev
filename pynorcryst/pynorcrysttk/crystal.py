from unit_cell import UnitCell

class Crystal(object):

    def __init__(self, lattice, spacegroup=None):
        self.unit_cell = UnitCell(lattice)

class CrystalSystem(Enum):
    __order__ = "triclinic monoclinic orthorhombic tetragonal trigonal hexagonal cubic"
    triclinic = 0
    monoclinic = 1
    orthorhombic = 2
    tetragonal = 3
    trigonal = 4
    hexagonal = 5
    cubic = 6