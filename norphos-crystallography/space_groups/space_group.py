from enum import Enum
from collections import namedtuple

class SpaceGroup(object):

    def __init__(self, number, full_symbol, crystal_system, wyckoff_positions):
        self.number = number
        self.full_symbol = full_symbol
        self.crystal_system = crystal_system
        self.wyckoff_positions = wyckoff_positions

        self.short_symbol = self.__make_short_symbol()

    def __make_short_symbol(self):
        symbol = self.full_symbol.centering
        sg_has_symm = False
        for axis in [self.full_symbol.pri, self.full_symbol.sec, self.full_symbol.tert]:
            if axis is '1':
                sg_has_symm = sg_has_symm | False
                continue
            else:
                sg_has_symm = sg_has_symm | True
            symbol = symbol+axis

        if not sg_has_symm: #Triclinic
            symbol = symbol+'1'

        return symbol


class WyckoffPosition(object):

    def __init__(self, multiplicity, letter, positions, systematic_absences):
        self.multiplicity = multiplicity
        self.letter = letter
        self.positions = positions
        self.systematic_absences = systematic_absences

        self.symbol = self.__make_symbol()

    def __make_symbol(self):
        return str(self.multiplicity)+self.letter


FullSGSymbol = namedtuple('FullSgSymbol', 'centering, pri, sec, tert')


class CrystalSystem(Enum):
    __order__ = "TRICLINIC MONOCLINIC ORTHORHOMBIC TETRAGONAL TRIGONAL RHOMBOHEDRAL HEXAGONAL CUBIC"
    TRICLINIC = 0
    MONOCLINIC = 1
    ORTHORHOMBIC = 2
    TETRAGONAL = 3
    TRIGONAL = 4
    RHOMBOHEDRAL = 5
    HEXAGONAL = 6
    CUBIC = 7