# test_fitting_problem.py
#
# Bryan Daniels
#
# 1.16.2019
#
# Tests for fittingProblem.py
#
'''
import unittest
from SirIsaac.fittingProblem import *
import SloppyCell.ReactionNetworks as scrn

def mock_net():
    """
    from SirIsaac import fittingProblem
import SirIsaac


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
'''
import scipy, pylab
from SirIsaac import fittingProblem


def test_x(self):
    data = scipy.loadtxt('simpleExample_data.txt')
    indepParamsList = [ [ expt[0] ] for expt in data ]
    indepParamsList[:3]
    sirIsaacData = []
    for expt in data:
        sirIsaacData.append( { 'x': { expt[1]: ( expt[2], expt[3] ) } } )
    sirIsaacData[:3]
    outputNames = ['x']
    indepParamNames = ['x_init']
    complexityStepsize = 2 # increase complexity with steps of size 2
    complexityMax = 25 # don't try models with complexity > 25
    complexityList = range(0,complexityMax,complexityStepsize) 

    # ensGen controls the generation of the initial ensemble of 
    # parameter starting points.
    totalSteps = 1e3
    keepSteps = 10
    seeds = (1,1) # use a fixed random seed
    ensTemperature = 100.
    ensGen = fittingProblem.EnsembleGenerator( totalSteps, keepSteps,
        temperature=ensTemperature, seeds=seeds )

    # Parameters that control when local fitting stops.
    avegtol = 1e-2
    maxiter = 100

    # priorSigma controls the width of priors on all parameters
    priorSigma = 3.

    # If you have pypar installed, you can run on multiple processors
    numprocs = 6

    # We'll only use a subset of our data to make the example run faster
    N = 20

    p = fittingProblem.PowerLawFittingProblem( complexityList, 
        sirIsaacData[:N], indepParamsList=indepParamsList[:N], 
        outputNames=outputNames, indepParamNames=indepParamNames, 
        ensGen=ensGen, avegtol=avegtol, maxiter=maxiter,
        priorSigma=priorSigma, numprocs=numprocs, verbose=True )
    p.fitAll()
    #
    fittingProblem.save(p,'simpleExample_savedFittingProblem.data')
