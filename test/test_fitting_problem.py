# test_fitting_problem.py
#
# Bryan Daniels
#
# 1.16.2019
#
# Tests for fittingProblem.py
#

import unittest

from SirIsaac.fittingProblem import *
import SloppyCell.ReactionNetworks as scrn

def mock_net():
    """
    Create a simple mock SloppyCell Network object
    """
    return scrn.Network('mocknet')

def mock_data():
    """
    Create a simple mock dataset with a single datapoint
    """
    return [ {'x': {1: (0.,1.) }}, ]

class TestFittingProblem(unittest.TestCase):

    def test_sloppycell_fitting_model_init(self):
        """
        Test that a basic SloppyCell fitting model can be initialized correctly
        """
        net = mock_net()
        mtest = SloppyCellFittingModel(net)
        
        self.assertEqual([], mtest.getParameters())

    def test_fitting_problem_init(self):
        """
        Test that a basic fitting problem can be initialized correctly
        """
        net = mock_net()
        m = SloppyCellFittingModel(net)
        datatest = mock_data()
        indepParamsList = [[]]
        modellisttest = [m]
        ftest = FittingProblem(datatest,modellisttest,
                               indepParamsList=indepParamsList)

        self.assertEqual(None, ftest.getBestModel())