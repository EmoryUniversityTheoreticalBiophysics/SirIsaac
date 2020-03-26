# 3.26.2020
#
# Bryan Daniels
#
# Tests for functions that use parallel processing through MPI.
#

import unittest

import SirIsaac

NUMPROCS = 2

class TestParallel(unittest.TestCase):
    def test_local_fit_parallel(self):
        
        # create simple SloppyCell fitting model
        complexity = 1
        indepParamNames = ['a']
        outputNames = ['b']
        m = SirIsaac.fittingProblem.CTSNFittingModel(complexity,
            indepParamNames,outputNames)
        
        # create some fake data (with small added noise and constant seed)
        numPoints = 3
        timeInterval = [0,10]
        noiseFracSize = 0.1
        data = [SirIsaac.fakeData.noisyFakeData(m.net,
                                                numPoints,
                                                timeInterval,
                                                noiseFracSize=noiseFracSize,
                                                seed=0)]
        
        # set up fit
        indepParamsList = [[1.]]
        dataModel = m._SloppyCellDataModel(data,indepParamsList)
        p = m.getParameters()
        
        # first run serially
        fitParamsSerial,_ = m.localFitToData(data,dataModel,startParams=p)
        # then run in parallel
        ens = [p,p,p] # start from multiple (equivalent) points
        outputDictParallel = \
            m.localFitToData_pypar(NUMPROCS,data,dataModel,ens,indepParamsList)
        
        # check that we get the same answer from each parallel instance
        testVar = fitParamsSerial.keys()[0] # just check the first parameter
        varSerial = fitParamsSerial.getByKey(testVar)
        varParallelList = [ result[0].getByKey(testVar) \
                            for result in outputDictParallel.values() ]
        for varParallel in varParallelList:
            self.assertAlmostEqual(varSerial,varParallel)
        
        
