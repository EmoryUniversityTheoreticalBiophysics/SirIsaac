# makeSloppyEnsemble.py
#
# Bryan Daniels
# 8.24.2009
# 9.25.2013
#

import fittingProblem
from SloppyCell.ReactionNetworks import *
import time # for sleep
import os

print "This computer's name is",os.uname()[1]
if (os.uname()[1][:4] == 'node'): # 4.4.2012 emory machines
    print "The current directory is",os.getcwd()
    if os.getcwd().startswith('/star'):
        os.chdir('/star/physics/nemenman/daniels/SirIsaac')
    elif os.getcwd().startswith('/spark'):
        os.chdir('/spark/physics/nemenman/daniels/SirIsaac')
    print "Now the current directory is",os.getcwd()

def makeEnsemble(fpdFilename,numDataPoints,modelName=None,numSteps=1000,        \
    numStepsKept=10,sing_val_cutoff=1e-4,seeds=(100,100),numprocs=1,            \
    verbose=True):
    """
    Generates SloppyCell parameter ensemble and saves to file.
    
    Returns ens,cost:
    ens  = list of ensemble members
    cost = list of costs
    
    modelName (None)        : If None, use maxLogLikelihoodName.
                              If 'Perfect', use perfectModel.
    seeds ((100,100))       : Used to seed the sampling algorithm
    """
    
    # load model
    fpd = Utility.load(fpdFilename)
    fp = fpd[numDataPoints]
    if modelName is None:
        modelName = fp.maxLogLikelihoodName()
    if modelName is 'Perfect':
        m = fp.perfectModel
        initialParameters = fp.perfectFitParams
    else:
        m = fp.fittingModelDict[modelName]
        initialParameters = m.getParameters()
    logParams = m.ensGen.logParams

    # set up output
    outputFilename = fpdFilename[:5] + '_' + str(numDataPoints) +               \
        '_' + modelName + '_ensemble' + str(numStepsKept) + '.dat'

    if verbose:
        print "makeSloppyEnsemble: Generating ensemble for "+outputFilename

    # set up model
    dataModel = m._SloppyCellDataModel(fp.fittingData,fp.indepParamsList)

    # generate the ensemble
    sloppyEnsGen = fittingProblem.EnsembleGenerator(numSteps,numStepsKept,      \
        sing_val_cutoff=sing_val_cutoff,seeds=seeds,logParams=logParams)
    if numprocs == 1:
        output = sloppyEnsGen.generateEnsemble(dataModel,                       \
            initialParameters,returnCosts=True,scaleByDOF=False)
    else:
        output = sloppyEnsGen.generateEnsemble_pypar(numprocs,dataModel,        \
            initialParameters,returnCosts=True,scaleByDOF=False)
        
    # save the result
    ens,ratio,costs = output
    Utility.save(output,outputFilename)

    if verbose:
        print "makeSloppyEnsemble: Acceptance ratio =",ratio

    return ens,ratio,costs

if __name__ == '__main__':

    # Specify the model(s) for which to generate ensemble(s)
    numDataPoints = 300 #200 #100 #52
    filenameList,modelNameList = [],[]
    
    if True: # perfect phosphorylation model
        filenameList.append('k0011_fitProb_varying_numInputs_PhosphorylationNet_'   \
                'PerfectPhosphorylation_withEnsembleT1000_steps10000.0_10_'         \
                'maxiter100_avegtol0.01_noiseFracSize0.1_ratePriorSigma10_'         \
                'seeds0_1.dat')
        #filenameList.append('k0032_fitProb_varying_numInputs_PhosphorylationNet_'   \
        #        'PerfectPhosphorylation_withEnsembleT1000_steps10000.0_10_'         \
        #        'maxiter100_avegtol0.01_noiseFracSize0.1_ratePriorSigma1000.0_'     \
        #        'seeds0_1_restart0006.dat')
        modelNameList.append( 'Perfect' )
    if True: # sigmoidal phosphorylation model
        filenameList.append('k0012_fitProb_varying_numInputs_PhosphorylationNet_'   \
                'CTSN_withEnsembleT1000_steps10000.0_10_maxiter100_avegtol0.01_'    \
                'noiseFracSize0.1_ratePriorSigma10_seeds0_1.dat')
        #filenameList.append('k0009_fitProb_varying_numInputs_PhosphorylationNet_'   \
        #        'CTSN_withEnsembleT1000_steps10000.0_10_useBest_numPoints1_'        \
        #        'maxiter100_avegtol0.01_noClamp_newErrorBars0.1_'                   \
        #        'removeLogForPriors_ratePriorSigma1000.0_seeds0_1_restart0013.dat')
        modelNameList.append( None ) # None => use maxLogLikelihoodName
    if True: # simple phosphorylation model
        filenameList.append('s0004_fitProb_varying_numInputs_PhosphorylationNet_'   \
            'SimplePhosphorylation_withEnsembleT1000_steps10000.0_10_'              \
            'useBest_numPoints1_maxiter100_avegtol0.01_noClamp_newErrorBars0.1_'    \
            'removeLogForPriors_ratePriorSigma1000.0_seeds0_1_restart0001.dat')
        modelNameList.append( 'SimplePhosphorylationModel' )

    # wait for perfect fit to be done
    fp = None
    i = 0
    while not hasattr(fp,'perfectFitParams'):
        i += 1
        Utility.save([i,10000],'131008_temp10000.txt')
        print "makeSloppyEnsemble: Sleeping..."
        time.sleep(60)
        fpd = Utility.load(filenameList[0])
        fp = fpd[numDataPoints]
    
    # set up ensemble generator
    numprocs = 20 #10
    numSteps = 10000 # 1000 takes roughly 60 minutes on 8 processors for perfect 200
    numStepsKept = 1000
    sing_val_cutoff = 0.1 #1e-4 #1

    # create each ensemble
    for fpdFilename,modelName in zip(filenameList,modelNameList):
        ens,ratio,costs = makeEnsemble(fpdFilename,numDataPoints,modelName,         \
            numSteps,numStepsKept,sing_val_cutoff,numprocs=numprocs)

    
    
