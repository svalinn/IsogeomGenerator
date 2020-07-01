"""tests for the driver script tool.py"""

from IsogeomGenerator import driver
import pytest


# Generate Levels parametrized tests:
# linear: (6, 5, 15, 'lin', [5., 7., 9., 11., 13., 15.]
# log: (6, 1, 1e5, 'log', [1, 10, 1e2, 1e3, 1e4, 1e5])
# ratio, max included: (5, 1, 625, 'ratio', [1., 5., 25., 125., 625.])
# ratio, max not included: (5, 1, 700, 'ratio', [1., 5., 25., 125., 625.])
@pytest.mark.parametrize("N,minN,maxN,mode,exp",
                         [(6, 5, 15, 'lin', [5., 7., 9., 11., 13., 15.]),
                          (6, 1, 1e5, 'log', [1, 10, 1e2, 1e3, 1e4, 1e5]),
                          (5, 1, 625, 'ratio', [1., 5., 25., 125., 625.]),
                          (5, 1, 700, 'ratio', [1., 5., 25., 125., 625.])])
def test_generate_levels(N, minN, maxN, mode, exp):
    """generate levels with different modes"""
    obs = driver.generate_levels(N, minN, maxN, mode=mode)
    assert(obs == exp)
