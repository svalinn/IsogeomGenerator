from IsogeomGenerator import vol
import unittest


class TestIsogeom(unittest.TestCase):

    def test_assign_levels(self):
        """assign predetermined levels (out of order)"""
        v = vol.IsoVolume()
        levels = [0.1, 0.2, 0.05]
        exp = sorted(levels)
        v.assign_levels(levels)
        assert(v.levels == exp)
        assert(v.minN == min(levels))
        assert(v.maxN == max(levels))
        assert(v.N == len(levels))


    def test_generate_levels_linear(self):
        """generate linearly spaced levels"""
        v = vol.IsoVolume()
        N = 6
        minN = 5
        maxN = 15
        exp = [5., 7., 9., 11., 13., 15.]
        v.generate_levels(N, minN, maxN, log=False)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))


    def test_generate_levels_log(self):
        """generate lograthmically spaced levels"""
        v = vol.IsoVolume()
        N = 6
        minN = 1
        maxN = 1e5
        exp = [1., 10., 1.e2, 1.e3, 1.e4, 1.e5]
        v.generate_levels(N, minN, maxN, log=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))

    def test_generate_levels_ratio_1(self):
        """generate levels by ratio, max included"""
        v = vol.IsoVolume()
        N = 5
        minN = 1
        maxN = 625
        exp = [1., 5., 25., 125., 625.]

        v.generate_levels(N, minN, maxN, ratio=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))


    def test_generate_levels_ratio_2(self):
        """generate levels by ratio, max not included"""
        v = vol.IsoVolume()
        N = 5
        minN = 1
        maxN = 700
        exp = [1., 5., 25., 125., 625.]

        v.generate_levels(N, minN, maxN, ratio=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))
