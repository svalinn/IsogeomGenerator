from IsogeomGenerator import vol
import unittest


class TestIsogeom(unittest.TestCase):

    def test_assign_levels(self):

        v = vol.IsoVolume()
        levels = [0.1, 0.2, 0.05]
        exp = sorted(levels)

        v.assign_levels(levels)
        assert(v.levels == exp)
