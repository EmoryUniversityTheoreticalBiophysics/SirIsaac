# 3.26.2020
#
# Bryan Daniels
#
# Tests for functions that use parallel processing through MPI.
#

import unittest

import SirIsaac.fittingProblem
import SirIsaac.fakeData
from SirIsaac.simplePickle import load,save
import subprocess
import os

NUMPROCS = 2
SIRISAACDIR = SirIsaac.fittingProblem.SIRISAACDIR

def example_data_and_model():
    
    # create simple SloppyCell fitting model with single independent parameter
    complexity = 1
    indepParamNames = ['input']
    outputNames = ['output']
    m = SirIsaac.fittingProblem.CTSNFittingModel(complexity,
        indepParamNames,outputNames)
    
    # create a single set of fake data (with small added noise and constant seed)
    numPoints = 3
    timeInterval = [0,10]
    noiseFracSize = 0.1 # 0.01 hangs for local_fit_parallel??
    lenData = 1 # (2 is probably a better test for parallel stuff in generateEnsemble)
    data = [SirIsaac.fakeData.noisyFakeData(m.net,
                                            numPoints,
                                            timeInterval,
                                            noiseFracSize=noiseFracSize,
                                            seed=seed) for seed in range(lenData)]
                                       
    indepParamsList = [[1.+i] for i in range(lenData)]
    dataModel = m._SloppyCellDataModel(data,indepParamsList)
                                       
    return m,data,indepParamsList,dataModel


class TestParallel(unittest.TestCase):

    def test_basic_mpi_functionality(self):
        """
        Test basic MPI functionality
        """
        temp_input_filename = os.path.join(SIRISAACDIR,
                                           "temporary_input_{}.dat".format(os.getpid()))
        temp_output_filename = os.path.join(SIRISAACDIR,
                                            "temporary_output_{}.dat".format(os.getpid()))
        temp_stdout_filename = os.path.join(SIRISAACDIR,
                                            "temporary_stdout_{}.txt".format(os.getpid()))
        
        input_data = {'test':123,'output_filename':temp_output_filename}
        save(input_data,temp_input_filename)
        stdoutFile = open(temp_stdout_filename,'w')
        subprocess.call([ "mpirun","-np",str(NUMPROCS),"python",
                          os.path.join(SIRISAACDIR, "mpi_basic.py"),
                          temp_input_filename ],
                          stderr=stdoutFile,stdout=stdoutFile,
                          env=os.environ)
        output_data = load(temp_output_filename)
        os.remove(temp_input_filename)
        os.remove(temp_output_filename)
        os.remove(temp_stdout_filename)

        self.assertTrue(input_data['test'] == output_data['test_output'])

    def test_local_fit_parallel(self):
        """
        Test localFitToData_parallel
        """
        
        m,data,indepParamsList,dataModel = example_data_and_model()
        
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

    def test_generate_ensemble_parallel(self):
        """
        Test generateEnsemble_parallel
        """

        m,data,indepParamsList,dataModel = example_data_and_model()

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
        ensParallel,_ = ensGen.generateEnsemble_parallel(NUMPROCS,dataModel,p)

        # check that we get the same answer when running in parallel
        for paramIndex in range(len(ensSerial[0])):
            self.assertAlmostEqual(ensSerial[0][paramIndex],
                                   ensParallel[0][paramIndex])

    
        
        
