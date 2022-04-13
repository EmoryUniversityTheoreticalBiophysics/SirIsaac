# generateEnsembleParallel.py
#
# Bryan Daniels
# 7.2.2012
#
# Used in FittingProblem.EnsembleGenerator.generateEnsemble_pypar
#

#!/usr/bin/env python
from numpy import *
import time

#from generateFightData import *
from .fittingProblem import *
import sys

# use parallel computation supported by SloppyCell
# (in cost and sensitivity calculations)
import SloppyCell.ReactionNetworks.RunInParallel as Par

# read in arguments from command line file name
if len(sys.argv) < 2:
    print("Usage: python generateEnsembleParallel.py inputDictFile.data [SloppyCell options]")
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

allOutputsDict = {}

output = ensGen.generateEnsemble(dataModel,initialParameters,
    returnCosts=returnCosts,scaleByDOF=scaleByDOF)

# Write data to file
save(output,outputFilename)

