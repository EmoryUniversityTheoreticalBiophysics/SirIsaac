# runFittingProblem.py
#
# Bryan Daniels
# 6.30.2009
#
# Code that reproduces the fitting done in the 2013 manuscript.
#
    
import scipy
from SloppyCell.ReactionNetworks import *
import FittingProblem
reload(FittingProblem)
import FakeData
reload(FakeData)
import os,sys,copy
from outputTag import nextFileNumString
import SloppyCellTest

print "This computer's name is",os.uname()[1]
if (os.uname()[1][:4] == 'node'): # 4.4.2012 emory machines
    print "The current directory is",os.getcwd()
    if os.getcwd().startswith('/star'):
        os.chdir('/star/physics/nemenman/daniels/SirIsaac')
    elif os.getcwd().startswith('/spark'):
        os.chdir('/spark/physics/nemenman/daniels/SirIsaac')
    print "Now the current directory is",os.getcwd()
    
def paramsDict(fittingProblem):
    d = {}
    for name in fittingProblem.fittingModelNames:
        params = fittingProblem.fittingModelDict[name].getParameters()
        d[name] = params
    return d

outputDirectory = '.'

# set tag to command-line argument unless there isn't a valid one
if len(sys.argv) < 2:
    prefix = nextFileNumString(outputDirectory)
elif (sys.argv[1].find('&') != -1) or (sys.argv[1].find('>') != -1)         \
    or (sys.argv[1].find('|') != -1):
    prefix = nextFileNumString(outputDirectory)
else:
    prefix = sys.argv[1]

# 9.24.2013 make sure SloppyCell C compiling is working
if not SloppyCellTest.testCcompiling():
    raise Exception, "SloppyCell C compiling not working."

noiseInLog = False # 3.6.2013
usePreviousParams = True #False # 3.6.2013
randomX = True
avegtol = 1e-2 #1. #1e-2 
maxiter = 100 # 200 
verbose = False
numprocs = 11
useDerivs = False
includeEndpoints = False
inputNames = None # default to non-init variables
stopFittingN = 3
numPoints = 1 # take one data point per independent parameter set
valOffsets = None # default; changed for PowerLaw fit to Phosphorylation

# parameters for ensemble generation
useEnsemble = True #False 
totalSteps = 1e4 #100
keepSteps = 10 #20 # 5
ensTemperature = 1e3 #1e5 #1e6 #1e5 # 100000. 100. # changed later for fitting perfect model
sing_val_cutoff = 1.
seeds = (1,1)
if useEnsemble:
    ensGen = FittingProblem.EnsembleGenerator(                                      \
        totalSteps,keepSteps,ensTemperature,sing_val_cutoff,seeds)
else:
    ensGen = None

# () choose data source
#originalString = 'PlanetaryNet'
originalString = 'yeastOscillator'
#originalString = 'PhosphorylationNet'

# () choose fitting model class
#fittingType = 'Polynomial'
#fittingType = 'Laguerre'
#fittingType = 'SimplePhosphorylation'
#fittingType = 'PerfectPhosphorylation'
#fittingType = 'PowerLaw'
fittingType = 'CTSN'


# 8.7.2013 use Planetary network to generate perfect data
if originalString is 'PlanetaryNet':
    
    timeAndNoiseSeed = 0 #0
    ICseed = 1 #1
    switchSigmoid = False
    xiNegative = False
    noiseFracSize = 0.05 #0.01 #0.1
    
    maxNumInputs = 1000 # we'll generate this many random inputs for possible use
    
    scipy.random.seed(ICseed)
    
    varsType = 'rOnly'
    #varsType = 'rAndTheta'
    
    if varsType is 'rAndTheta':
        inputVars = ['r_init','theta_init']
        inputNames = ['r_init']
        outputVars = ['r','theta']
    elif varsType is 'rOnly':
        inputVars = ['r_init']
        inputNames = ['r_init']
        outputVars = ['r']

    inputMin,inputMax = 1.,2.5 #1.,3. # units GM/(v0^2) (1->circle,2->parabola)
    inputList = inputMin + (inputMax-inputMin)*scipy.random.random(maxNumInputs)
    # 5.3.2013 set first two inputs to the two extremes
    inputList[0] = inputMin
    inputList[1] = inputMax

    if varsType is 'rAndTheta':
        theta_init = 2.*scipy.pi # nonzero to avoid problems with power law models
        inputListFull = [ [input,theta_init] for input in inputList ]
    elif varsType is 'rOnly':
        inputListFull = [ [input] for input in inputList ]
    
    timeInterval = [0.,100.] # units GM/(v0^3)
    includeDerivs = False
        
    originalFittingModel = FittingProblem.PlanetaryFittingModel(                    \
        indepParamNames=inputVars,verbose=True,avegtol=avegtol,maxiter=maxiter,     \
        ensGen=ensGen)
    
    rateVars = []
    nonrateVars = []

    ratePriorSigma = 10. #1e3
    nonratePriorSigma = 10.
    
    # 9.5.2013 since we're interested in testing whether we can find the
    # perfect representation, we'll remove the stopping criterion
    stopFittingN = scipy.inf
    
    fakeDataAbs = False # Avoid negative data



# 7.23.2009 use Phosphorylation network to generate perfect data
elif originalString is 'PhosphorylationNet':

    timeAndNoiseSeed = 0 #0
    ICseed = 1 #1
    switchSigmoid = False # False
    xiNegative = False
    noiseFracSize = 0.1
    
    maxNumInputs = 1000 # we'll generate this many random inputs for possible use
    
    scipy.random.seed(ICseed)

    # 9.19.2012 try to avoid crossing zero so much
    # 6.4.2013 also need to avoid having totalPhos_init = 0
    if fittingType is "PowerLaw":
        offset = 1.
        valOffsets = [offset]
    else:
        offset = 0.

    inputVars = ['k23p','totalPhos_init'] # 5.2.2013 added 'totalPhos_init'
    inputLogMin,inputLogMax = -3,3 #-3,4
    inputLogList = inputLogMin + (inputLogMax-inputLogMin)*scipy.random.random(maxNumInputs)
    # 5.3.2013 set first two inputs to the two extremes
    inputLogList[0] = inputLogMin
    inputLogList[1] = inputLogMax
    inputListFull = [ [10**inputLog,offset] for inputLog in inputLogList ]

    timeInterval = [0.,10.] #[0.,3.] #[0.,10.]
    includeDerivs = False

    outputVars = ['totalPhos']

    # 9.8.2013
    # for saving original model using BioNetGen
    #makeOriginalModel = True
    # for machines without BioNetGen
    makeOriginalModel = False

    rateVars = ['k']
    nonrateVars = ['Km','totalPhos']

    ratePriorSigma = 10 #1e3 #10 #1e3 #10. #1e3
    nonratePriorSigma = 10.

    originalModelFilename = 'examplePhosphorylationFittingModel.model'
    if makeOriginalModel:
        n = 5
        rules = [(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(4,5),(5,4)] #[(2,3)]
        originalFittingModel = FittingProblem.PhosphorylationFittingModel(n,rules,      \
            indepParamNames=inputVars[:1],verbose=True,avegtol=avegtol,maxiter=maxiter, \
            ensGen=ensGen,totalOffset=offset)
        # 4.29.2013 set random parameters
        randomParamsSeed = 12345
        scipy.random.seed(randomParamsSeed)
        newParams = {}
        for var in originalFittingModel.getParameters().keys():
            if var.startswith('k'):
                newParams.update({var: abs(scipy.random.normal(scale=nonratePriorSigma))})
            if var.startswith('Km'):
                newParams.update({var: abs(scipy.random.normal(scale=nonratePriorSigma))})
        originalFittingModel.initializeParameters(newParams)
        Utility.save(originalFittingModel,originalModelFilename)
        die
    else:
        # try to load from file
        originalFittingModel = Utility.load(originalModelFilename)
        # 9.9.2013 using different offset for different models, so set it here
        originalFittingModel.net.setInitialVariableValue('totalPhos_offset',offset)

    fakeDataAbs = True # Avoid negative data
    


# yeast oscillator (see RuoChrWol03)
elif originalString is 'yeastOscillator':
  
  #to do 3.22.2012
  originalFittingModel = None #yeastOscillatorFittingModel(inputVars) 
  
  # *****************
  upperRangeMultiple = 1. #1.
  nonzeroMin = True # True
  # *****************
    
  includedIndices = range(3) #range(7) # range(3) 
  allNames = scipy.array(['S1','S2','S3','S4','N2','A3','S4ex']) 
  names = allNames[includedIndices]
  outputVars = names
    
  switchSigmoid = False # False
  xiNegative = False
    
  timesSeed = 0 #0
  noiseSeed = 1 #1
  ICseed = 2 #2
  noiseFracSize = 0.1
  
  fakeDataAbs = False 
  
  # (no SloppyCell perfectModel, so no prior vars)
  rateVars = []
  nonrateVars = []
  
  ratePriorSigma = 10.
  nonratePriorSigma = 10.
  
  def yeastDataFunction(numICs,useDerivs,                                   \
    names=names,timesSeed=timesSeed,noiseSeed=noiseSeed,ICseed=ICseed):
    if (os.uname()[1] != 'star'): # can't do if MATLAB isn't installed
        from simulateYeastOscillator import *
            
    timeInterval = [0.,5.] #[0.,10.] #[0.,1e-5] #[0.,10.] # minutes #
            
    multiplicativeErrorBar = noiseFracSize # 0.1
            
    fittingData,fittingDataDerivs,inputVars,inputList =                     \
        yeastData(numPoints,timeInterval,                                   \
        numICs,useDerivs,includedIndices,timesSeed=timesSeed,               \
        noiseSeed=noiseSeed,ICseed=ICseed,                                  \
        multiplicativeErrorBar=multiplicativeErrorBar,randomX=randomX,      \
        upperRangeMultiple=upperRangeMultiple,nonzeroMin=nonzeroMin)
    
    return fittingData,inputVars,inputList,useDerivs,fittingDataDerivs
    
else:
    raise Exception, "Unrecognized originalString"


# for Laguerre polynomials
degreeListLag = scipy.arange(0,8,1) # (0,10,1)
# polynomials describing parameter dependence have degree half
# that of the degree of the polynomial describing time dependence
polynomialDegreeListListLag = [ (degree/2)*scipy.ones(2+degree,dtype=int)   \
    for degree in degreeListLag ]

# for plain polynomials
degreeListPoly = scipy.arange(0,8,1) # (0,10,1)
# polynomials describing parameter dependence have degree half
# that of the degree of the polynomial describing time dependence
polynomialDegreeListListPoly = [ (degree/2)*scipy.ones(degree,dtype=int)  \
    for degree in degreeListPoly ]

# 4.29.2013 set priors
if (fittingType is 'Polynomial') or (fittingType is 'Laguerre'):
    nonrateVars.extend( ['C','sqrt_abs_alpha','g'] )
elif fittingType is 'PowerLaw':
    rateVars.extend( ['log_alpha','log_beta'] )
    nonrateVars.extend( ['g','h','X'] )
elif fittingType is 'CTSN':
    rateVars.extend( ['w','log_tau'] )
    nonrateVars.extend( ['theta','X'] )
elif (fittingType is 'SimplePhosphorylation')                               \
  or (fittingType is 'PerfectPhosphorylation'):
    pass # we don't include priors for these models
else:
    raise Exception, "Unrecognized fittingType"
# make priorSigma list using ratePriorSigma and nonratePriorSigma
priorSigma = []
for v in rateVars: priorSigma.append( (v,ratePriorSigma) )
for v in nonrateVars: priorSigma.append( (v,nonratePriorSigma) )

fitProbDict = {}

# () optionally restart calculations from a loaded fittingProblemDict
restartDictName = None
#restartDictName = '0062_fitProb_varying_numInputs_yeastOscillator_CTSN_withEnsembleT1000_steps10000.0_10_useBest_numPoints1_maxiter100_avegtol0.01_noClamp_newErrorBars0.1_removeLogForPriors_seeds0_1_2_restart0027.dat'
#restartDictName = 'k0030_fitProb_varying_numInputs_yeastOscillator_CTSN_withEnsembleT1000_steps10000.0_10_maxiter100_avegtol0.01_noiseFracSize0.1_ratePriorSigma10.0_seeds3_4_5_restart0041.dat'
#restartDictName = '0042_fitProb_varying_numInputs_yeastOscillator_CTSN_withEnsembleT1000_steps10000.0_10_useBest_numPoints1_maxiter100_avegtol0.01_noClamp_newErrorBars0.1_removeLogForPriors_ratePriorSigma10.0_seeds6_7_8_restart0038.dat'
if restartDictName is not None:
    fitProbDict = Utility.load(restartDictName)
    i = restartDictName.find('_fitProb_')
    restartStr = '_restart'+restartDictName[i-4:i]
    # try to catch inconsistencies
    if restartDictName.find(fittingType) < 0: raise Exception
    if restartDictName.find(originalString) < 0: raise Exception
    if originalString is "yeastOscillator":
        seedsStr = '_seeds'+str(timesSeed)+'_'+str(noiseSeed)+'_'+str(ICseed)
    elif (originalString is "PhosphorylationNet") or (originalString is "PlanetaryNet"):
        seedsStr = '_seeds'+str(timeAndNoiseSeed)+'_'+str(ICseed)
    else:
        raise Exception
    if restartDictName.find(seedsStr) < 0: raise Exception
    # ***
    # 9.5.2013 temporary change to stopFittingN for planetary fits
    #for key in fitProbDict.keys():
    #    fp = fitProbDict[key]
    #    fp.stopFittingN = stopFittingN
    # ***

# () set up filename for output
fileNumString = prefix
print "runFittingProblem: Output files will start with",fileNumString
configString = '_fitProb_varying_numInputs'                                 \
   +'_'+originalString                                                      \
   +'_'+fittingType                                                         \
   +'_withEnsembleT'+str(int(ensTemperature))                               \
   +'_steps'+str(totalSteps)+'_'+str(keepSteps)                             \
   +'_maxiter'+str(maxiter)                                                 \
   +'_avegtol'+str(avegtol)                                                 \
   +'_noiseFracSize'+str(noiseFracSize)        \
   +'_ratePriorSigma'+str(ratePriorSigma)
if originalString is "yeastOscillator":
    configString += '_seeds'+str(timesSeed)+'_'+str(noiseSeed)+'_'+str(ICseed)
elif (originalString is "PhosphorylationNet") or (originalString is "PlanetaryNet"):
    configString += '_seeds'+str(timeAndNoiseSeed)+'_'+str(ICseed)
elif originalString is "powerLawYeastOscillator":
    configString += '_seeds'+str(timeAndNoiseSeed)+'_'+str(indepParamsSeed)
else:
    raise Exception

if restartDictName is not None:
    configString += restartStr

Utility.save(fitProbDict,fileNumString+configString+'.dat')
saveFilename = fileNumString+configString+'.dat'#+'_partial.dat'

# () set up numIndepParamsList, specifying the lengths of datasets to test
smallerBestParamsDict = {}
if originalString is "PhosphorylationNet":
    deltaNumIndepParams = 5 #2
    maxNumIndepParams = 54
    numIndepParamsList = range(deltaNumIndepParams,maxNumIndepParams,deltaNumIndepParams)
    numIndepParamsList.extend([52,100,200,300,400,500])
elif originalString is "yeastOscillator":
    deltaNumIndepParams = 2
    maxNumIndepParams = 52 #25
    numIndepParamsList = range(deltaNumIndepParams,maxNumIndepParams,deltaNumIndepParams)
elif originalString is "PlanetaryNet":
    deltaNumIndepParams = 10
    maxNumIndepParams = 200
    numIndepParamsList = range(deltaNumIndepParams,maxNumIndepParams,deltaNumIndepParams)
else: raise Exception

# () set up complexityList, specifying which models to test in the model class
useFullyConnected = False #True
if useFullyConnected:
    numModels = 10
    complexityList = scipy.linspace(0,1,numModels) # fraction of possible parameters
else:
    if originalString is "PhosphorylationNet":
        complexityStepsize = 2
        complexityMax = 25 #50
    elif originalString is "yeastOscillator":
        complexityStepsize = 5
        complexityMax = 200
    elif originalString is "PlanetaryNet":
        complexityStepsize = 2
        complexityMax = 25
    else: raise Exception
    complexityMin = 0 #19 #7 #0
    complexityList = scipy.arange(complexityMin,complexityMax,complexityStepsize)

fittingDataDerivs = None
previousParams = None # 4.29.2013 used for fitting perfect model



# () loop over an increasing amount of fittingData and perform fitting
for numIndepParams in numIndepParamsList:

    key = numIndepParams

    # () Create fittingData
    if (originalFittingModel is not None) and (not fitProbDict.has_key(key)):
      # the following need to be set above:
      #     outputVars
      #     inputVars,inputListFull
      #
      # 5.3.2013
      # if inputVars contains variables that are not
      #   defined in the original model, then they are
      #   ignored here, but still passed on to be used
      #   as input to the fittingModels
      #   (eg. 'totalPhos_init' for phosphorylation model)
      originalNet = originalFittingModel.net
      originalNet.compile() # only want to do this once
      
      
      runVars,runList = inputVars,inputListFull[:numIndepParams]
      inputList = inputListFull[:numIndepParams]
      
      # in case inputListFull is too short
      if len(runList) != numIndepParams: raise Exception
  
      fakeData = []
      for i,runVals in enumerate(runList):
        newNet = originalNet.copy()
        for runVar,runVal in zip(runVars,runVals):
          if runVar in newNet.parameters.keys():
            newNet.setInitialVariableValue(runVar,runVal)
          else:
            pass # 5.3.2013 runVal will still be passed on to the fittingModels
        fakeDataSingleRun = {}
        for j,var in enumerate(outputVars):
            
            # 9.24.2013 to match with what I was doing before when len(outputVars)=1
            if len(outputVars) > 1: noiseSeed = int(timeAndNoiseSeed*1e5+(i+1)*1e3+j)
            else: noiseSeed = None
            
            # do individually so every var is measured at the same (random) time
            fakeDataSingleRun.update( FakeData.noisyFakeData(newNet,numPoints,      \
                timeInterval,                                                       \
                seed=int(timeAndNoiseSeed*1e5+i),vars=[var],                        \
                noiseFracSize=noiseFracSize,randomX=randomX,                        \
                includeEndpoints=includeEndpoints,takeAbs=fakeDataAbs,              \
                noiseSeed=noiseSeed,typValOffsets=valOffsets ))
        fakeData.append( fakeDataSingleRun )
    
    elif originalString == 'yeastOscillator':
      print "Using yeast oscillator data."
      fakeData,inputVars,inputList,includeDerivs,fittingDataDerivs =           \
          yeastDataFunction(numIndepParams,useDerivs)

    # () Create fittingProblem p
    if fitProbDict.has_key(key): 
        p = fitProbDict[key]
        p.saveFilename = saveFilename # in case it has changed
    else: # we haven't started 
        kwargs = {  'avegtol': avegtol,
                    'maxiter': maxiter,
                    'ensGen': ensGen,
                    'verbose': verbose,
                    'indepParamNames': inputVars,
                    'indepParamsList': inputList,
                    'perfectModel': copy.deepcopy(originalFittingModel),
                    'saveFilename': saveFilename,
                    'includeDerivs': includeDerivs,
                    'numprocs': numprocs,
                    'smallerBestParamsDict': smallerBestParamsDict,
                    'saveKey': key,
                    'stopFittingN': stopFittingN,
                 }
        if fittingType is 'Laguerre':
          p = FittingProblem.LaguerreFittingProblem(degreeListLag,fakeData,
            outputName=outputVars[0],
            polynomialDegreeListList=polynomialDegreeListListLag,**kwargs)
        elif fittingType is 'Polynomial':
          p = FittingProblem.PolynomialFittingProblem(degreeListPoly,fakeData,
            outputName=outputVars[0],
            polynomialDegreeListList=polynomialDegreeListListPoly,**kwargs)
        elif (fittingType is 'SimplePhosphorylation')                               \
          or (fittingType is 'PerfectPhosphorylation'):
          kwargs.pop('stopFittingN') # only 1 model to fit
          p = FittingProblem.SimplePhosphorylationFittingProblem(fakeData,
            offset=offset,**kwargs)
        elif fittingType is 'PowerLaw':
          p = FittingProblem.PowerLawFittingProblem(complexityList,fakeData,
            outputNames=outputVars,priorSigma=priorSigma,
            fittingDataDerivs=fittingDataDerivs,
            useFullyConnected=useFullyConnected,inputNames=inputNames,**kwargs)
        elif fittingType is 'CTSN':
          kwargs['xiNegative'] = xiNegative
          p = FittingProblem.CTSNFittingProblem(complexityList,fakeData,
            outputNames=outputVars,priorSigma=priorSigma,                       
            inputNames=inputNames,**kwargs)
        else:
            raise Exception, 'No valid fittingType specified.'
        
        fitProbDict[key] = p

    Utility.save(fitProbDict,fileNumString+configString+'.dat')
    
    # () Fit the models in the fittingProblem p
    if fittingType is not 'PerfectPhosphorylation':
        try:
            # fit models
            p.fitAll(usePreviousParams=usePreviousParams)
        except KeyboardInterrupt:
            raise
    Utility.save(fitProbDict,fileNumString+configString+'.dat')

    if fittingType is 'PerfectPhosphorylation':
      if not hasattr(p,'perfectCost'):
        # 4.29.2013 fit perfect model for phosphorylation
        # (see runFitPerfectModelPhos.py)
        p.perfectModel.numprocs = numprocs
        p.perfectModel.priorSigma = priorSigma
        p.perfectModel.speciesNames = ['totalPhos']
        p.perfectModel.ensGen.logParams = True
        p.perfectModel.ensGen.temperature = 10. # 5.1.2013
        p.perfectModel.ensGen.totalSteps = 100. # 5.2.2013
        p.fitPerfectModel(otherStartingPoint=previousParams)
        previousParams = p.perfectFitParams
        Utility.save(fitProbDict,fileNumString+configString+'.dat')
      else:
        previousParams = p.perfectFitParams

    # 4.17.2012
    smallerBestParamsDict = paramsDict(p)
        
    print "runFittingProblem: Done with key", key



