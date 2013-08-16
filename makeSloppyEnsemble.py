# makeSloppyEnsemble.py
#
# Bryan Daniels
# 8.24.2009
#
# 

import sys
import FittingProblem
reload(FittingProblem)
from SloppyCell.ReactionNetworks import *

fpdPhosSigPerfect = Utility.load('0079_fitProb_varying_randomSeed_'             \
    +'PhosphorylationNet_PowerLaw_withInput3_withEnsembleT100000_'              \
    +'steps1000_10_numPoints8_useBest_perfectOnly.dat')
FittingProblem.UpdateOldFitProbDict(fpdPhosSigPerfect)

outputFilename = '0079_PhosphorylationNet_perfectModel_ensembles_'              \
    +'priorSigma100.dat'

# set up ensemble generator
numSteps = 50000 # 10000 # takes 15-30 minutes?
numStepsKept = 5000
sing_val_cutoff = 1e-4
seeds = (1,1)
logParams = True # Was I using log parameters before?  Why are all parameters
                 # positive in my fits?
sloppyEnsGen = FittingProblem.EnsembleGenerator(numSteps,numStepsKept,          \
    sing_val_cutoff=sing_val_cutoff,seeds=seeds,logParams=logParams)

ensembleDict = {}

for seed in fpdPhosSigPerfect.keys():

    # set up the sloppyCell model
    p = fpdPhosSigPerfect[seed]
    p.perfectModel.priorSigma = 100.
    dataModel =                                                                 \
        p.perfectModel._SloppyCellDataModel(p.fittingData,p.indepParamsList)
    initialParameters = p.perfectFitParams

    # generate the ensemble
    ens,costs = sloppyEnsGen.generateEnsemble(dataModel,initialParameters,      \
        returnCosts=True)
        
    # save the result
    ensembleDict[seed] = ens,costs
    Utility.save(ensembleDict,outputFilename)
    