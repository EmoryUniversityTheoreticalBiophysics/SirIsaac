# 3.26.2020
#
# Bryan Daniels
#
# Tests for functions that use parallel processing through MPI.
#

import unittest

import SirIsaac.fittingProblem
import SirIsaac.fakeData

NUMPROCS = 2

def test_data_and_model():
    
    # create simple SloppyCell fitting model with single independent parameter
    complexity = 1
    indepParamNames = ['a']
    outputNames = ['b']
    m = SirIsaac.fittingProblem.CTSNFittingModel(complexity,
        indepParamNames,outputNames)
    
    # create a single set of fake data (with small added noise and constant seed)
    numPoints = 3
    timeInterval = [0,10]
    noiseFracSize = 0.1 # 0.01 hangs for local_fit_parallel??
    data = [SirIsaac.fakeData.noisyFakeData(m.net,
                                            numPoints,
                                            timeInterval,
                                            noiseFracSize=noiseFracSize,
                                            seed=0)]
                                       
    indepParamsList = [[1.] for d in data]
    dataModel = m._SloppyCellDataModel(data,indepParamsList)
                                       
    return m,data,indepParamsList,dataModel

class TestParallel(unittest.TestCase):

    def test_generate_ensemble_parallel(self):
        
        m,data,indepParamsList,dataModel = test_data_and_model()
        
        # set up ensemble generation
        p = m.getParameters()
        totalSteps = 25
        keepSteps = 5
        seeds = (1,1)
        ensGen = SirIsaac.fittingProblem.EnsembleGenerator(totalSteps,keepSteps,
                                                           seeds=seeds)
        
        # first run serially
        ensSerial,_ = ensGen.generateEnsemble(dataModel,p)
        # then run in parallel
        ensParallel,_ = ensGen.generateEnsemble_pypar(NUMPROCS,dataModel,p)
        
        # check that we get the same answer when running in parallel
        for paramIndex in range(len(ensSerial[0])):
            self.assertAlmostEqual(ensSerial[0][paramIndex],
                                   ensParallel[0][paramIndex])

    def test_local_fit_parallel(self):
        
        m,data,indepParamsList,dataModel = test_data_and_model()
        
        # set up fit
        p = m.getParameters()
        
        # first run serially
        fitParamsSerial,_ = m.localFitToData(data,dataModel,startParams=p)
        # then run in parallel
        ens = [p,p,p] # start from multiple (equivalent) points
        outputDictParallel = \
            m.localFitToData_parallel(NUMPROCS,data,dataModel,ens,indepParamsList)
        
        # check that we get the same answer from each parallel instance
        testVar = fitParamsSerial.keys()[0] # just check the first parameter
        varSerial = fitParamsSerial.getByKey(testVar)
        varParallelList = [ result[0].getByKey(testVar) \
                            for result in outputDictParallel.values() ]
        for varParallel in varParallelList:
            self.assertAlmostEqual(varSerial,varParallel)
        
        
