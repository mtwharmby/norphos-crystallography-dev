from space_groups.space_group import SpaceGroup
from space_groups.space_group import FullSGSymbol
from space_groups.space_group import WyckoffPosition

import unittest
from nose.tools import assert_true, assert_false, assert_equal

class TestSpaceGroup(unittest.TestCase):

    def test_short_symbol(self):
        self.sg = SpaceGroup(14, FullSGSymbol('P','1','21/c','1'), None, None)
        assert_equal('P21/c',self.sg.short_symbol)

        self.sg = SpaceGroup(2, FullSGSymbol('P','-1','1','1'), None, None)
        assert_equal('P-1', self.sg.short_symbol)

        self.sg = SpaceGroup(1, FullSGSymbol('P','1','1','1'), None, None)
        assert_equal('P1', self.sg.short_symbol)

        self.sg = SpaceGroup(214, FullSGSymbol('I','41','3','2'),None, None)
        assert_equal('I4132', self.sg.short_symbol)



class TestWyckoffPosition(unittest.TestCase):

    def test_wyckoff_symbol(self):
        self.wyk = WyckoffPosition(4, 'e', None, None)
        assert_equal('4e', self.wyk.symbol)