# test_SloppyCell.py
#
# Bryan Daniels
# 7/24/2013
#
# Tests whether SloppyCell is functional
#

import unittest

from SirIsaac import fittingProblem
import numpy as np
from SloppyCell.ReactionNetworks import *

class TestSloppyCell(unittest.TestCase):

    def test_c_compiling(self):
        m = fittingProblem.CTSNFittingModel(5,['a'],['b'])
        n = m.net
        n.compile()
        
        # check that SloppyCell thinks the compiling worked
        self.assertTrue(n.compiled)
        
        # ensure compiled c code exists
        ccode = n.get_c_code()
        self.assertTrue(ccode in n._c_module_cache)
        self.assertTrue(n._c_module_cache[ccode] is not None)

    def test_CTSN(self):
        m = fittingProblem.CTSNFittingModel(5,['a'],['b'])
        r = m.evaluateVec(np.linspace(0,1,10),'b',[0])
        eps = 1e-5
        return sum(r**2) < eps


if __name__ == '__main__':
    unittest.main()
