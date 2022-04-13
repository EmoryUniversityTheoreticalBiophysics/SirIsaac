# test_fitting_problem.py
#
# Bryan Daniels
#
# 1.16.2019
#
# Tests for fittingProblem.py
#

import unittest

import SirIsaac.fittingProblem as fp
import SloppyCell.ReactionNetworks as scrn
import numpy as np

def mock_net():
    """
    Create a simple mock SloppyCell Network object
    """
    return scrn.Network('mocknet')

def mock_data():
    """
    Create a mock dataset with a single datapoint
    """
    return [ {'x': {1: (0.,1.) }}, ]
    
def simple_linear_data(N):
    """
    Create a simple example dataset depending linearly on input and time,
    with N datapoints.
    """
    indepParamsList = [ [i,] for i in range(1,N+1) ]
    data = [ {'y': {i[0]: (2.*i[0], 1.)}} for i in indepParamsList ]
    return indepParamsList, data

class TestFittingProblem(unittest.TestCase):

    def test_sloppycell_fitting_model_init(self):
        """
        Test that a basic SloppyCell fitting model can be initialized correctly
        """
        net = mock_net()
        mtest = fp.SloppyCellFittingModel(net)
        
        self.assertEqual([], mtest.getParameters())

    def test_fitting_problem_init(self):
        """
        Test that a basic fitting problem can be initialized correctly
        """
        net = mock_net()
        m = fp.SloppyCellFittingModel(net)
        datatest = mock_data()
        indepParamsList = [[]]
        modellisttest = [m]
        ftest = fp.FittingProblem(datatest,modellisttest,
                                  indepParamsList=indepParamsList)

        self.assertEqual(None, ftest.getBestModel())

    def test_fitAll(self):
        """
        Test that fitAll produces reasonable results on an easy test problem.
        """
        numDatapoints = 3
        complexityList = [0,2]
        
        # ensGen controls the generation of the initial ensemble of
        # parameter starting points.
        totalSteps = 5
        keepSteps = 5
        seeds = (1,1) # use a fixed random seed
        ensTemperature = 100.
        ensGen = fp.EnsembleGenerator( totalSteps, keepSteps,
            temperature=ensTemperature, seeds=seeds )
        
        # Parameters that control when local fitting stops.
        avegtol = 1e-2
        maxiter = 100
        
        # set up simple linear data
        indepParamNames = ['x',]
        indepParamsList,data = simple_linear_data(numDatapoints)
        outputNames = data[0].keys()
        
        p = fp.PowerLawFittingProblem(
                complexityList,
                data,
                outputNames=outputNames,
                indepParamsList=indepParamsList,
                indepParamNames=indepParamNames,
                ensGen=ensGen,
                avegtol=avegtol,
                maxiter=maxiter)
            
        p.fitAll()
        
        # check that the correct number of models have been initialized
        self.assertEqual(len(p.fittingModelNames),len(complexityList))
        self.assertEqual(len(p.numParametersDict),len(complexityList))
        
        # check that some models have indeed been fit
        self.assertTrue(len(p.logLikelihoodDict) > 0)
        
        # loop over models that have been fit
        for modelName in p.logLikelihoodDict:
            # check that all fit models have numerical log-likelihood
            self.assertFalse(np.isnan(p.logLikelihoodDict[modelName]))
        
            # check that costs (distance from data) are roughly what we expect
            self.assertTrue(p.costDict[modelName] < 0.05)
            
            # check we have the correct number of singular values
            self.assertEqual(len(p.singValsDict[modelName]),p.numParametersDict[modelName])
        
