# generateEnsembleParallel.py
#
# Bryan Daniels
# 7.2.2012
#
# Used in FittingProblem.EnsembleGenerator.generateEnsemble_pypar
#

#!/usr/bin/env python
from numpy import *
import pypar
import time

#from generateFightData import *
from FittingProblem import *
import sys

# for parallel computation supported by SloppyCell
import SloppyCell.ReactionNetworks.RunInParallel as Par


# Constants
MASTER_PROCESS = 0
WORK_TAG = 1
DIE_TAG = 2

MPI_myID = pypar.rank() #Par.my_rank 
num_processors = pypar.size() #Par.num_procs

# read in arguments from command line file name
if len(sys.argv) < 2 or len(sys.argv) > 2:
    print "Usage: python generateEnsembleParallel.py inputDictFile.data"
    exit()
inputDictFile = sys.argv[1]
inputDict = load(inputDictFile)
ensGen = inputDict['ensGen']
dataModel = inputDict['dataModel']
initialParameters = inputDict['initialParameters']
returnCosts = inputDict['returnCosts'] 
scaleByDOF = inputDict['scaleByDOF']
outputFilename = inputDict['outputFilename']
inputDict['outputDict'] = {}
fitFunction = lambda startParamsIndex:                          \
    fittingProblem.localFitToData(fittingData,dataModel,        \
    retall=True,startParams=startParamsList[startParamsIndex])

allOutputsDict = {}

output = ensGen.generateEnsemble(dataModel,initialParameters,   \
    returnCosts=returnCosts,scaleByDOF=scaleByDOF)

# Write data to file
save(output,outputFilename)

