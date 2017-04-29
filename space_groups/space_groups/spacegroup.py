

class SpaceGroup(object):

    def __init__(self, number, full_symbol, crystal_system, wyckoff_positions):
        self.number = number
        self.full_symbol = full_symbol
        self.crystal_system = crystal_system
        self.wyckoff_positions = wyckoff_positions

    def short_symbol(self):
        pass


class WyckoffPosition(object):

    def __init__(self, multiplicity, wyckoff_letter, systematic_absences):
        self.multiplicity = multiplicity
        self.wyckoff_letter = wyckoff_letter
        self.systematic_absences = systematic_absences