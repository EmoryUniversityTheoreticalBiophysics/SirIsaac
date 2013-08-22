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

noiseFracSize = 0.1 #0.01
noiseInLog = False # 3.6.2013
usePreviousParams = True #False # 3.6.2013
randomX = True #False #True
avegtol = 1e-2 #1. #1e-2
maxiter = 100 # 200
verbose = False
numprocs = 11 
useDerivs = False
restartPhos = False # default
includeEndpoints = False
inputNames = None # default to non-init variables

# for ensemble generation
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
useClampedPreminimization = False

# 8.7.2013 use Planetary network to generate perfect data
if True: 
    
    timeAndNoiseSeed = 0 #0
    ICseed = 1 #1
    switchSigmoid = False 
    
    numPoints = 1
    maxNumInputs = 1000 # we'll generate this many random inputs for possible use
    
    scipy.random.seed(ICseed)
    
    inputVars = ['r_init','theta_init']
    inputMin,inputMax = 1.,2.5 #1.,3. # units GM/(v0^2) (1->circle,2->parabola)
    inputList = inputMin + (inputMax-inputMin)*scipy.random.random(maxNumInputs)
    # 5.3.2013 set first two inputs to the two extremes
    inputList[0] = inputMin
    inputList[1] = inputMax
    theta_init = 2.*scipy.pi # nonzero to avoid problems with power law models
    inputListFull = [ [input,theta_init] for input in inputList ]
    
    timeInterval = [0.,100.] # units GM/(v0^3)
    includeDerivs = False
    
    inputNames = ['r_init']
    outputVars = ['r','theta'] #['r']
    
    originalFittingModel = FittingProblem.PlanetaryFittingModel(                    \
        indepParamNames=inputVars,verbose=True,avegtol=avegtol,maxiter=maxiter,     \
        ensGen=ensGen)
    
    originalString = 'PlanetaryNet'

    maxSVDeig = None
    
    rateVars = []
    nonrateVars = []

    ratePriorSigma = 1e3
    nonratePriorSigma = 10.
    
    fakeDataAbs = True # Avoid negative data



# 7.23.2009 use Phosphorylation network to generate perfect data
if False: 

    timeAndNoiseSeed = 6 #0
    ICseed = 7 #1
    switchSigmoid = False # False <<<<
    
    numPoints = 1
    maxNumInputs = 1000 # we'll generate this many random inputs for possible use
    
    scipy.random.seed(ICseed)

    # 9.19.2012 try to avoid crossing zero so much
    # 6.4.2013 also need to avoid having totalPhos_init = 0
    offset = 1. #0 #1

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

    # 5.8.2013
    # for saving data to use on Emory machines
    #restartPhos = False
    #savePhosDict = True
    #makeOriginalModel = True
    # for running on Emory machines
    restartPhos = True
    savePhosDict = False
    makeOriginalModel = False
    # for runs not on Emory machines
    #restartPhos = False
    #savePhosDict = False
    #makeOriginalModel = True # False if on Emory machines with manual restart

    n = 5
    rules = [(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(4,5),(5,4)] #[(2,3)]
    if makeOriginalModel:
        originalFittingModel = FittingProblem.PhosphorylationFittingModel(n,rules,  \
        indepParamNames=inputVars[:1],verbose=True,avegtol=avegtol,maxiter=maxiter, \
        ensGen=ensGen,totalOffset=offset)
    
    originalString = 'PhosphorylationNet'

    maxSVDeig = None

    rateVars = ['k']
    nonrateVars = ['Km','totalPhos']

    ratePriorSigma = 1e3 #10. #1e3
    nonratePriorSigma = 10.

    if makeOriginalModel:
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
    else:
        originalFittingModel = None

    fakeDataAbs = True # Avoid negative data
    


    
if False: # yeast oscillator (see RuoChrWol03) 
  
  #to do 3.22.2012
  originalFittingModel = None #yeastOscillatorFittingModel(inputVars) 
  originalString = 'yeastOscillator'
  
  # *****************
  numPoints = 1 #10 #100
  upperRangeMultiple = 1. #1.
  nonzeroMin = True # True
  # *****************
    
  includedIndices = range(3) #range(7) # range(3) 
  allNames = scipy.array(['S1','S2','S3','S4','N2','A3','S4ex']) 
  names = allNames[includedIndices]
  outputVars = names
    
  timesSeed = 12 #0
  noiseSeed = 13 #1
  ICseed = 14 #2
  
  fakeDataAbs = False # <--- may want to change to True if rerunning everything
  
  maxSVDeig = None
  
  # (no SloppyCell perfectModel, so no prior vars)
  rateVars = []
  nonrateVars = []
  
  ratePriorSigma = 10.
  nonratePriorSigma = 10.
  
  def yeastDataFunction(numICs,useDerivs,                                   \
    names=names,timesSeed=timesSeed,noiseSeed=noiseSeed,ICseed=ICseed):
    if (os.uname()[1] != 'star'): # can't do if MATLAB isn't installed
        from simulateYeastOscillator import *
            
    # NOTE! if changing timeInterval, make sure you run MATLAB stuff again
    timeInterval = [0.,10.] #[0.,1e-5] #[0.,10.] # minutes #
            
    multiplicativeErrorBar = noiseFracSize # 0.1
            
    fittingData,fittingDataDerivs,inputVars,inputList =                     \
        yeastData(numPoints,timeInterval,                                   \
        numICs,useDerivs,includedIndices,timesSeed=timesSeed,               \
        noiseSeed=noiseSeed,ICseed=ICseed,                                  \
        multiplicativeErrorBar=multiplicativeErrorBar,randomX=randomX,      \
        upperRangeMultiple=upperRangeMultiple,nonzeroMin=nonzeroMin)
    
    return fittingData,inputVars,inputList,useDerivs,fittingDataDerivs
    
    


# for Laguerre polynomials
degreeListLag = scipy.arange(0,8,1) # (0,10,1)
polynomialDegreeListListLag = [ (degree/2)*scipy.ones(3+degree,dtype=int)   \
    for degree in degreeListLag ]

# for plain polynomials
degreeListPoly = scipy.arange(0,8,1) # (0,10,1)
polynomialDegreeListListPoly = [ (degree/2)*scipy.ones(degree+1,dtype=int)  \
    for degree in degreeListPoly ]

# pick fitting model class
#fittingType = 'Polynomial'
#fittingType = 'Laguerre'
#fittingType = 'PowerLaw'
fittingType = 'CTSN'

# 4.29.2013 set priors
if (fittingType is 'Polynomial') or (fittingType is 'Laguerre'):
    raise Exception, "4.29.2013 Need to set priors for "+fittingType
elif fittingType is 'PowerLaw':
    rateVars.extend( ['log_alpha','log_beta'] )
    nonrateVars.extend( ['g','h','X'] )
elif fittingType is 'CTSN':
    if switchSigmoid: # 7.24.2013
        rateVars.extend( ['log_tau'] )
        nonrateVars.extend( ['theta','X','w'] )
    else:
        rateVars.extend( ['w','log_tau'] )
        nonrateVars.extend( ['theta','X'] )
else:
    raise Exception
# make priorSigma list using ratePriorSigma and nonratePriorSigma
priorSigma = []
for v in rateVars: priorSigma.append( (v,ratePriorSigma) )
for v in nonrateVars: priorSigma.append( (v,nonratePriorSigma) )

fitProbDict = {}

# () optionally restart calculations from a loaded fittingProblemDict
restartDictName = None
#restartDictName = 'k0001_fitProb_varying_numInputs_PhosphorylationNet_CTSN_withEnsembleT1000_steps10000.0_10_useBest_numPoints1_maxiter100_avegtol0.01_noClamp_newErrorBars0.1_removeLogForPriors_ratePriorSigma1000.0_seeds0_1.dat'
if restartDictName is not None:
    fitProbDict = Utility.load(restartDictName)
    i = restartDictName.find('_fitProb_')
    restartStr = '_restart'+restartDictName[i-4:i]
    # try to catch inconsistencies
    if restartDictName.find(fittingType) < 0: raise Exception
    if restartDictName.find(originalString) < 0: raise Exception
    if originalString is "yeastOscillator":
        seedsStr = '_seeds'+str(timesSeed)+'_'+str(noiseSeed)+'_'+str(ICseed)
    elif originalString is "PhosphorylationNet":
        seedsStr = '_seeds'+str(timeAndNoiseSeed)+'_'+str(ICseed)
    else:
        raise Exception
    if restartDictName.find(seedsStr) < 0: raise Exception

# () set up filename for output
fileNumString = prefix
print "runFittingProblem: Output files will start with",fileNumString
configString = '_fitProb_varying_numInputs'                                 \
   +'_'+originalString                                                      \
   +'_'+fittingType                                                         \
   +'_withEnsembleT'+str(int(ensTemperature))                               \
   +'_steps'+str(totalSteps)+'_'+str(keepSteps)                             \
   +'_useBest'                                                              \
   +'_numPoints'+str(int(numPoints))                                        \
   +'_maxiter'+str(maxiter)                                                 \
   +'_avegtol'+str(avegtol)                                                 \
   +'_noClamp_newErrorBars'+str(noiseFracSize)+'_removeLogForPriors'        \
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

Utility.save({},fileNumString+configString+'.dat')
saveFilename = fileNumString+configString+'.dat'#+'_partial.dat'

# 5.8.2013 for running on Emory machines that don't have BioNetGen
if restartPhos:
    if switchSigmoid:
        restartPhosDictFilename = 'restartPhosDict_switchSigmoid.data'
    else:
        restartPhosDictFilename = 'restartPhosDict.data'
    restartPhosDict = Utility.load(restartPhosDictFilename)
    phosStr = configString
    if phosStr not in restartPhosDict.keys():
        raise Exception, "restartPhosDict does not have key"+str(phosStr)
    fitProbDict = restartPhosDict[phosStr]

# () set up complexityList, specifying which models to test in the model class
smallerBestParamsDict = {}
deltaNumIndepParams = 10 #2 #1000 #50 #5 #20 #5 #2 
if originalString is "PhosphorylationNet": maxNumIndepParams = 54 #42
elif originalString is "yeastOscillator": maxNumIndepParams = 25
elif originalString is "PlanetaryNet": maxNumIndepParams = 200
else: raise Exception
numIndepParamsList = range(deltaNumIndepParams,maxNumIndepParams,deltaNumIndepParams) 
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
        for var in outputVars:
            # do individually so every var is measured at the same (random) time
            fakeDataSingleRun.update( FakeData.noisyFakeData(newNet,numPoints,      \
                timeInterval,                                                       \
                seed=int(timeAndNoiseSeed*1e5+i),vars=[var],                        \
                noiseFracSize=noiseFracSize,randomX=randomX,                        \
                includeEndpoints=includeEndpoints,takeAbs=fakeDataAbs) )
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
        # below: Laguerre networks, with input variables
        if fittingType is 'Laguerre':
          p = FittingProblem.LaguerreFittingProblem(degreeListLag,fakeData,     \
            outputName=outputVars[0],avegtol=avegtol,maxiter=maxiter,           \
            ensGen=ensGen,                                                      \
           verbose=verbose,indepParamNames=inputVars,indepParamsList=inputList, \
           perfectModel=copy.deepcopy(originalFittingModel),                    \
           polynomialDegreeListList=polynomialDegreeListListLag,                \
           saveFilename=saveFilename,includeDerivs=includeDerivs,               \
           useClampedPreminimization=useClampedPreminimization,                 \
           numprocs=numprocs,smallerBestParamsDict=smallerBestParamsDict,       \
           saveKey=key)
        # below: Polynomial networks, with input variables
        elif fittingType is 'Polynomial':
          p = FittingProblem.PolynomialFittingProblem(degreeListPoly,fakeData,  \
            outputName=outputVars[0],avegtol=avegtol,maxiter=maxiter,           \
            ensGen=ensGen,\
           verbose=verbose,indepParamNames=inputVars,indepParamsList=inputList, \
           perfectModel=copy.deepcopy(originalFittingModel),                    \
           polynomialDegreeListList=polynomialDegreeListListPoly,               \
           saveFilename=saveFilename,includeDerivs=includeDerivs,               \
           useClampedPreminimization=useClampedPreminimization,                 \
           numprocs=numprocs,smallerBestParamsDict=smallerBestParamsDict,       \
           saveKey=key)
        # below: PowerLaw Networks, with input variables
        elif fittingType is 'PowerLaw':
          p = FittingProblem.PowerLawFittingProblem(complexityList,fakeData,    \
            outputNames=outputVars,avegtol=avegtol,maxiter=maxiter,             \
            priorSigma=priorSigma,ensGen=ensGen,                                \
            verbose=verbose,indepParamNames=inputVars,                          \
            indepParamsListList=inputList,                                      \
            perfectModel=copy.deepcopy(originalFittingModel),                   \
            saveFilename=saveFilename,includeDerivs=includeDerivs,              \
            useClampedPreminimization=useClampedPreminimization,                \
            numprocs=numprocs,smallerBestParamsDict=smallerBestParamsDict,      \
            saveKey=key,fittingDataDerivs=fittingDataDerivs,                    \
            useFullyConnected=useFullyConnected,maxSVDeig=maxSVDeig,            \
            inputNames=inputNames)
        elif fittingType is 'CTSN':
          p = FittingProblem.CTSNFittingProblem(complexityList,fakeData,        \
            outputNames=outputVars,avegtol=avegtol,maxiter=maxiter,             \
            priorSigma=priorSigma,ensGen=ensGen,                                \
            verbose=verbose,indepParamNames=inputVars,                          \
            indepParamsListList=inputList,                                      \
            perfectModel=copy.deepcopy(originalFittingModel),                   \
            saveFilename=saveFilename,includeDerivs=includeDerivs,              \
            useClampedPreminimization=useClampedPreminimization,                \
            numprocs=numprocs,smallerBestParamsDict=smallerBestParamsDict,      \
            saveKey=key,switchSigmoid=switchSigmoid,inputNames=inputNames)
        else:
            raise Exception, 'No valid fittingType specified.'
        
        fitProbDict[key] = p

    # () Fit the models in the fittingProblem p
    try:
        # fit models
        p.fitAll(usePreviousParams=usePreviousParams)
        pass
    except KeyboardInterrupt:
        raise
    Utility.save(fitProbDict,fileNumString+configString+'.dat')

    if False and (originalString is "PhosphorylationNet"):
      if not hasattr(p,'perfectCost'):
        # 4.29.2013 fit perfect model for phosphorylation
        # (see runFitPerfectModelPhos.py)
        p.perfectModel.numprocs = numprocs
        p.perfectModel.priorSigma = p.priorSigma
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


# 5.8.2013 creating a file for use on Emory machines
if savePhosDict:
    if switchSigmoid:
        restartPhosDictFilename = 'restartPhosDict_switchSigmoid.data'
    else:
        restartPhosDictFilename = 'restartPhosDict.data'
    restartPhosDict = Utility.load(restartPhosDictFilename)
    phosStr = configString
    restartPhosDict[phosStr] = fitProbDict
    Utility.save(restartPhosDict,restartPhosDictFilename)



