# fitAllParallel.py
#
# Bryan Daniels
# 4.4.2016
#
# Data structure for running in parallel fitting over increasing
# amounts of data and multiple conditions.
#

import scipy
import time, copy, os
from .simplePickle import load,save
from .sloppyCellTest import testCcompiling
from .fittingProblemMultipleCondition import *

def directoryPrefix(fileNumString,conditioni,numTimepoints):
    return fileNumString+'_fitProbs/N'+str(numTimepoints)+'/condition'+str(conditioni)+'/'

def directoryPrefixNonly(fileNumString,numTimepoints):
    return fileNumString+'_fitProbs/N'+str(numTimepoints)+'/'

def createDirectoryStructure(fileNumString,numConditions,numTimepointsList):
    os.mkdir(fileNumString+'_fitProbs/')
    for numTimepoints in numTimepointsList:
      os.mkdir(fileNumString+'_fitProbs/N'+str(numTimepoints))
      for i in range(numConditions):
        os.mkdir(directoryPrefix(fileNumString,i,numTimepoints))

def paramsDict(fittingProblem):
    d = {}
    for name in fittingProblem.fittingModelNames:
        params = fittingProblem.fittingModelDict[name].getParameters()
        d[name] = params
    return d

def saveFitProb(fitProb,saveFilename,fileNumString,conditioni,numTimepoints):
    dirPrefix = directoryPrefix(fileNumString,conditioni,numTimepoints)
    fitProbDict = {numTimepoints: fitProb}
    save(fitProbDict,dirPrefix+saveFilename)

def loadFitProb(saveFilename,fileNumString,conditioni,numTimepoints):
    dirPrefix = directoryPrefix(fileNumString,conditioni,numTimepoints)
    return load(dirPrefix+saveFilename)[numTimepoints]

def loadFitProbData(fileNumString):
    try:
        fitProbData = load(fileNumString+'_fitProbData.dat')
    except (IOError, EOFError):
        print("loadFitProbData: WARNING Unable to load fitProbData file."\
              "Returning None.")
        fitProbData = None
    return fitProbData

def saveFitProbData(fitProbData,fileNumString):
    try:
        save(fitProbData,fileNumString+'_fitProbData.dat')
    except IOError:
        print("saveFitProbData: WARNING Unable to save fitProbData file.")

def setLock(fileNumString):
    save(1,fileNumString+'_fileLocked.dat')

def removeLock(fileNumString):
    os.remove(fileNumString+'_fileLocked.dat')

def waitForUnlocked(fileNumString,maxIter=100):
    lockFilename = fileNumString+'_fileLocked.dat'
    i = 0
    while lockFilename in os.listdir('.'):
        print("waitForUnlocked: Waiting for another process to unlock fitProbData file...")
        # wait a bit
        time.sleep(1.+5.*scipy.rand())

        i += 1
        if i > maxIter:
            raise Exception("Waiting too long for lock on fitProbData file.")

def lockAndLoadFitProbData(fileNumString):
    waitForUnlocked(fileNumString)
    setLock(fileNumString)
    return loadFitProbData(fileNumString)

def saveAndUnlockFitProbData(fitProbData,fileNumString):
    saveFitProbData(fitProbData,fileNumString)
    removeLock(fileNumString)

def updateFitProbData(fitProb,fileNumString,conditioni,numTimepoints,modelj):
    
    fitProbData = lockAndLoadFitProbData(fileNumString)
    
    if fitProbData is not None:
        pDataMultiple = fitProbData[numTimepoints]
        pData = pDataMultiple['fitProbDataList'][conditioni]
        
        modelName = pData['fittingModelNames'][modelj]

        # insert new data for single condition
        pData['logLikelihoodDict'][modelName] = fitProb.logLikelihoodDict[modelName]
        pData['fittingStateDict'][modelName] = 'finished'

        # insert new data for the given model over all conditions
        # if they're all done for that N
        if scipy.all( [ p['fittingStateDict'][modelName] == 'finished' \
                        for p in pDataMultiple['fitProbDataList'] ] ):
            llList = [ p['logLikelihoodDict'][modelName] \
                       for p in pDataMultiple['fitProbDataList'] ]
            pDataMultiple['logLikelihoodDict'][modelName] = scipy.sum( llList )
            
            # check if we are also done fitting models for that N
            # [stop after seeing stopFittingN models with worse logLikelihood]
            orderedLs = []
            stopFittingN = pDataMultiple['stopFittingN']
            for n in pDataMultiple['fittingModelNames']:
                if n in pDataMultiple['logLikelihoodDict']:
                        orderedLs.append(pDataMultiple['logLikelihoodDict'][n])
                if (len(orderedLs) > stopFittingN):
                    if max(orderedLs[-stopFittingN:]) < max(orderedLs):
                        pDataMultiple['fitAllDone'] = True
            if len(orderedLs) == len(pDataMultiple['fittingModelNames']):
                # every fittingModel has been fit
                pDataMultiple['fitAllDone'] = True

    saveAndUnlockFitProbData(fitProbData,fileNumString)

# note: getState and setState are somewhat slow due to sorting
def getState(fitProbData,conditioni,numTimepointsi,modelj):
    numTimepoints = scipy.sort(list(fitProbData.keys()))[numTimepointsi]
    pData = fitProbData[numTimepoints]['fitProbDataList'][conditioni]
    modelName = pData['fittingModelNames'][modelj]
    return pData['fittingStateDict'][modelName]

# note: getState and setState are somewhat slow due to sorting
def setState(fitProbData,conditioni,numTimepointsi,modelj,state):
    numTimepoints = scipy.sort(list(fitProbData.keys()))[numTimepointsi]
    pData = fitProbData[numTimepoints]['fitProbDataList'][conditioni]
    modelName = pData['fittingModelNames'][modelj]
    pData['fittingStateDict'][modelName] = state

def assignWork(fileNumString):

    conditioni = None
    while conditioni is None: # wait for work to be available
        # wait a bit
        time.sleep(1.+scipy.rand())
        
        # load current fitProbData
        fitProbData = lockAndLoadFitProbData(fileNumString)
        if fitProbData is None:
            print("assignWork: Error loading fitProbData")
            removeLock(fileNumString)
        else:
            # find unstarted work to be done
            conditioni,numTimepointsi,modelj = findWork(fitProbData)
            if conditioni is None: removeLock(fileNumString)
    
    # mark work as started
    setState(fitProbData,conditioni,numTimepointsi,modelj,'started')
    
    # save updated fitProbData
    saveAndUnlockFitProbData(fitProbData,fileNumString)
    
    return conditioni,numTimepointsi,modelj

def findWork(fitProbData):
    numTimepointsList = scipy.sort( list(fitProbData.keys()) )
    # loop over numTimepoints
    for numTimepointsi,numTimepoints in enumerate(numTimepointsList):
        
        pDataMultiple = fitProbData[numTimepoints]
        
        # if there's work to be done for this N
        if not pDataMultiple['fitAllDone']:
        
            if numTimepointsi > 0:
                smallerPDataMultiple = fitProbData[numTimepointsList[numTimepointsi-1]]
            else:
                smallerPDataMultiple = None
    
            # pick the first model for which there's work to be done
            modelj = 0
            modelName = pDataMultiple['fittingModelNames'][modelj]
            while modelName in pDataMultiple['logLikelihoodDict']:
              modelj += 1
              modelName = pDataMultiple['fittingModelNames'][modelj]
            # loop over conditions
            for conditioni,pData in enumerate(pDataMultiple['fitProbDataList']):
                # if the model fit to less data has already been fit (if applicable)
                if (smallerPDataMultiple is None) or \
                   (smallerPDataMultiple['fitAllDone']) or \
                   (smallerPDataMultiple['fitProbDataList'][conditioni]\
                                        ['fittingStateDict']\
                                        [modelName] == 'finished'):
                  # and the previous j has been fit (if applicable)
                  if (modelj == 0) or \
                     (pData['fittingStateDict'] \
                           [pData['fittingModelNames'][modelj-1]] == 'finished'):
                    # and the model hasn't already been started
                    if pData['fittingStateDict'][modelName] == 'unstarted':
                        # then this is a model that needs to be fit
                        return conditioni,numTimepointsi,modelj
        
    print("findWork: No work found.")
    return None,None,None

def resetFitProbData(fileNumString):
    """
    Set all 'started' work to 'unstarted'.  (Leave 'finished' alone.)
    """
    fitProbData = lockAndLoadFitProbData(fileNumString)
    for pMultiple in list(fitProbData.values()):
        for p in pMultiple['fitProbDataList']:
            for name in list(p['fittingStateDict'].keys()):
                if p['fittingStateDict'][name] == 'started':
                    p['fittingStateDict'][name] = 'unstarted'
    saveAndUnlockFitProbData(fitProbData,fileNumString)

def setStopFittingN(fileNumString,stopFittingN,resetFitAllDone=True):
    """
    Overwrite stopFittingN values with given value.
    
    resetFitAllDone (True)  : If True, set all fitAllDone to False.
                              If False, leave all fitAllDone alone.
    """
    fitProbData = lockAndLoadFitProbData(fileNumString)
    for pMultiple in list(fitProbData.values()):
        pMultiple['stopFittingN'] = stopFittingN
        if resetFitAllDone: pMultiple['fitAllDone'] = False
    saveAndUnlockFitProbData(fitProbData,fileNumString)

def countFitProbData(fileNumString,printReport=True,printAll=False):
    """
    Print the current status of model fitting.
    """
    fitProbData = loadFitProbData(fileNumString)
    totalSubsets = len(list(fitProbData.values()))
    finishedSubsets = []
    finished,started,unstarted = 0,0,0
    for numTimepoints in scipy.sort(list(fitProbData.keys())):
        line = str(numTimepoints) + ' '
        pMultiple = fitProbData[numTimepoints]
        if pMultiple['fitAllDone']:
            finishedSubsets.append(numTimepoints)
            line += 'done '
        else:
            line += '     '
        subsetUnstarted = True
        for p in pMultiple['fitProbDataList']:
            indFinished,indStarted = 0,0
            for name in p['fittingModelNames']:
                s = p['fittingStateDict'][name]
                if s == 'finished':
                    finished += 1
                    indFinished += 1
                    subsetUnstarted = False
                if s == 'started':
                    started += 1
                    indStarted += 1
                    subsetUnstarted = False
                if s == 'unstarted': unstarted += 1
            line += str(indFinished)
            if indStarted > 0:
                line += '('+str(indStarted)+') '
            else:
                line += '    '
        if (not subsetUnstarted) or printAll:
            if printReport: print(line)

    if printReport:
        print("")
        print("Data subsets:")
        print("  ",len(finishedSubsets),"of",totalSubsets,"finished")
        print("Individual models:")
        print("  ",finished,"finished")
        print("  ",started,"running")
        print("  ",unstarted,"unstarted")
    else:
        return finishedSubsets

def makeFpdLean(fpd):
    """
    Modify in place to create a stripped-down version of fpd 
    that doesn't include the models.
    """
    for N in scipy.sort(list(fpd.keys())):
        fp = fpd[N]
        for f in fp.fittingProblemList:
            f.fittingModelDict = {}
            f.fittingModelList = []

# 10.4.2017
def subsetsWithFits(fileNumString,onlyNew=False):
    """
    Find data subsets (N) that have models that have been fit to
    all conditions.
    
    onlyNew (False)         : Optionally include only subsets that have
                              fits that are not included in the current
                              combined fitProbs.
    """
    fpd = loadFitProbData(fileNumString)
    saveFilename = list(fpd.values())[0]['saveFilename']
    
    Nlist = []
    for N in scipy.sort(list(fpd.keys())):
        # find models that have been fit to all conditions
        if len(fpd[N]['fitProbDataList']) == 1:
            fitModels = list(fpd[N]['fitProbDataList'][0]['logLikelihoodDict'].keys())
        else:
            fitModels = scipy.intersect1d([ list(fp['logLikelihoodDict'].keys()) \
                                            for fp in fpd[N]['fittingProblemList'] ])
        if onlyNew:
            Nfilename = directoryPrefixNonly(fileNumString,N)+'/'+saveFilename
            fileExists = os.path.exists(Nfilename)
            if not fileExists: # no combined file exists
                if len(fitModels) > 0:
                    Nlist.append(N)
            else: # check which fit models are currently included in the saved file
                fpMultiple = load(Nfilename)
                fitModelsSaved = list(fpMultiple.logLikelihoodDict.keys())
                if len(scipy.intersect1d(fitModels,fitModelsSaved)) < len(fitModels):
                    Nlist.append(N)
        else:
            if len(fitModels) > 0:
                Nlist.append(N)
    return Nlist


def combineFitProbs(fileNumString,saveCombined=True,combinedLean=True,
                    reset=False):
    """
    Combine fittingProblems from multiple conditions saved in the 
    parallel file structure into a single fittingProblemDict.
    
    Currently only includes data from models that have been fit to
    all conditions.
    
    saveCombined (True) : Overwrites any current top-level fitProbDict file
                          with a combined fitProbDict containing all 
                          numTimepoints.  Set to False to minimize memory use.
    combinedLean (True) : Combined fpd is saved without models to save memory.
    reset (False)       : If True, overwrite or delete any existing combined 
                          fitProbDicts.  This erases any existing 
                          outOfSampleCost information.
    """
    fitProbData = loadFitProbData(fileNumString)
    saveFilename = list(fitProbData.values())[0]['saveFilename']
    #save({},saveFilename)
    
    if saveCombined: fpdMultiple = {}
    
    fitSubsets = subsetsWithFits(fileNumString)
    subsetsToCombine = subsetsWithFits(fileNumString,onlyNew=not reset)
    
    for numTimepoints in fitSubsets:
    
      Nfilename = directoryPrefixNonly(fileNumString,numTimepoints)+'/'+saveFilename
      fileExists = os.path.exists(Nfilename)
      
      if fileExists and reset:
          # then an old combined file exists -- erase it to reset
          os.remove(Nfilename)
          fileExists = False
          print("combineFitProbs: Reset removed file for numTimepoints =",numTimepoints)
      
      if numTimepoints in subsetsToCombine: # combine
          oldOutOfSampleCostDict = {}
          if fileExists: # grab any out-of-sample cost data
              fpMultiple = load(Nfilename)
              if hasattr(fpMultiple,'outOfSampleCostDict'):
                  oldOutOfSampleCostDict = fpMultiple.outOfSampleCostDict
          
          p = fitProbData[numTimepoints]
          
          fpList = []
          for conditioni in range(len(p['fitProbDataList'])):
            fp = loadFitProb(saveFilename,fileNumString,conditioni,numTimepoints)
            fpList.append(fp)
          # make new multiple condition fitting problem by starting
          # with an empty fitting problem and inserting the fittingProblemList
          saveKey = p['saveKey']
          fp.stopFittingN = p['stopFittingN']
          fpMultiple = FittingProblemMultipleCondition([],[],saveFilename=None,
                                                       saveKey=saveKey,fp0=fp)
          fpMultiple.fittingProblemList = fpList
          fpMultiple.outOfSampleCostDict = oldOutOfSampleCostDict

          # Populate the logLikelihoodDict, etc by running fitAll.
          fpMultiple.fitAll(onlyCombine=True)
          
          if saveCombined: fpdMultiple[numTimepoints] = fpMultiple
          
          save(fpMultiple,Nfilename)
          
          print("combineFitProbs: Done with numTimepoints =",numTimepoints)
      
      else: # no new fits to combine; just load from file
      
          if saveCombined: fpdMultiple[numTimepoints] = load(Nfilename)
          print("combineFitProbs: Done with numTimepoints =",numTimepoints)
      


    if saveCombined:
        if combinedLean:
            makeFpdLean(fpdMultiple)
        save(fpdMultiple,saveFilename[:-4]+'_combined.dat')



def dataSubset(fittingData,numDatapoints,seed=345,maxNumIndepParams=None):
    """
    By default, add one timepoint for each independent parameter first,
    then increase the number of timepoints per independent parameter.
    Timepoints are added randomly for each independent parameter.
    Independent parameters are added in the order of indepParamsList.
    """
    scipy.random.seed(seed)
    subset = []
    numIndepParams = len(fittingData)
    if maxNumIndepParams is None: maxNumIndepParams = numIndepParams
    numDatapoints = int(numDatapoints)
    for i in range(min(numDatapoints,maxNumIndepParams)):
        varNames = scipy.sort( list(fittingData[i].keys()) )
        allTimes = scipy.sort( list(fittingData[i][varNames[0]].keys()) )
        
        possibleIndices = list(range(len(allTimes)))
        scipy.random.shuffle(possibleIndices)
        
        N = numDatapoints/maxNumIndepParams
        if i < numDatapoints%maxNumIndepParams: N += 1
        timeIndices = possibleIndices[:N]
        times = scipy.array(allTimes)[timeIndices]

        s = {}
        for var in varNames:
            s[var] = dict([(t,fittingData[i][var][t]) for t in times])
        subset.append(s)

    return subset


def initializeFitAllParallel(fullFittingProblem,fileNumString,
    deltaNumDatapoints=2,maxTimesPerIndepParam=None,timeOrderSeed=123,
    verbose=True,numIndepParams=None):
    """
    Creates data structure on disk for keeping track of fitting over increasing
    amounts of data and multiple conditions.
    
    After initialization, use runFitAllParallelWorker to run fitting.
    Multiple workers can be run at the same time.
    
    By default, add one timepoint for each independent parameter first,
    then increase the number of timepoints per independent parameter.
    Timepoints are added randomly for each independent parameter.
    Independent parameters are added in the order of indepParamsList.
    
    The amount of data is always kept equal across each condition.
    
    If the length of indepParamsList or the number of timepoints per
    independent parameter varies in the original, the total amount of
    data used will be (#conditions)x(minimum number of indepParams per
    condition)x(minimum number of timepoints per indepParam)
    (that is, NOT ALL DATA WILL BE USED).
    
    fullFittingProblem can be an instance of a FittingProblem or a 
    FittingProblemMultipleCondition.
    
    deltaNumDatapoints (2)      : The change in the number of datapoints
                                  (per condition) between successive fits.
    maxTimesPerIndepParam (None): The maximum number of timepoints used
                                  per independent parameter.
    numIndepParams (None)       : Number of independent parameter
                                  combinations (trials) to use in-sample per
                                  condition.  Defaults to the minimum
                                  number available over conditions.
    """
    # (only one fittingProblem if there are not multiple conditions)
    fittingProblemList = getattr(fullFittingProblem,
                                'fittingProblemList',
                                [fullFittingProblem])
    
    if fullFittingProblem.saveFilename is not None:
        configString = fullFittingProblem.saveFilename[4:-4]
    else:
        configString = ''

    # The length of fittingProblemList[0].fittingData is len(indepParamsList).

    # N is the number of datapoints per condition.
    
    # calculate maxN, the total number of datapoints per condition
    numIndepParamsList,numTimepointsList = [],[]
    for fittingProblem in fittingProblemList:
        numIndepParamsList.append(len(fittingProblem.fittingData))
        for d in fittingProblem.fittingData:
            numTimepointsList.append(len(list(d.values())[0]))
    if numIndepParams is None:
        numIndepParams = min(numIndepParamsList)
    elif numIndepParams > min(numIndepParamsList):
        raise Exception
    minNumTimepoints = min(numTimepointsList)
    if maxTimesPerIndepParam is not None:
        minNumTimepoints = min(minNumTimepoints,maxTimesPerIndepParam)
    maxN = numIndepParams*minNumTimepoints

    Nlist = list(range(deltaNumDatapoints,maxN,deltaNumDatapoints))
    Nlist = Nlist + [maxN]

    createDirectoryStructure(fileNumString,len(fittingProblemList),Nlist)

    # () With each increasing amount of data, make a copy of the fullFittingProblem
    #    that includes only that data.
    fitProbData = {}
    for N in Nlist:
        fitProbDataList = []
        for i,fittingProblem in enumerate(fittingProblemList):
            fittingData = fittingProblem.fittingData
            fittingDataSubset = dataSubset(fittingData,N,seed=timeOrderSeed+i,
                                           maxNumIndepParams=numIndepParams)
            indepParamsListSubset = \
                fittingProblem.indepParamsList[:len(fittingDataSubset)]

            newFittingProblem = copy.deepcopy(fittingProblem)
            newFittingProblem.setData(fittingDataSubset,
                                      indepParamsListSubset,
                                      newFittingProblem.indepParamNames)
            newFittingProblem.saveKey = N
            #fittingProblemListNew.append(newFittingProblem)

            # store each full fittingProblem in separate file
            fitProbDict = { N: newFittingProblem }
            dirPrefix = directoryPrefix(fileNumString,i,N)
            save(fitProbDict,dirPrefix+fileNumString+configString+'.dat')

            # in fitProbData, store only info necessary to decide which
            # fittingProblem to work on next
            fitProb = newFittingProblem
            fittingStateDictInitial = \
                dict( [ (name,'unstarted') for name in fitProb.fittingModelNames ])
            pData = {'logLikelihoodDict': fitProb.logLikelihoodDict,
                     'fittingStateDict': fittingStateDictInitial,
                     'fittingModelNames': fitProb.fittingModelNames,
                     'stopFittingN': fitProb.stopFittingN,
                     'saveFilename': fitProb.saveFilename,
                     'saveKey': N,
                     }
            fitProbDataList.append(pData)

        p = fullFittingProblem
        cp = copy.deepcopy
        pDataMultiple = {'logLikelihoodDict': cp(p.logLikelihoodDict),
                     'fitAllDone': cp(p.fitAllDone),
                     'fittingModelNames': cp(p.fittingModelNames),
                     'fitProbDataList': cp(fitProbDataList),
                     'stopFittingN': cp(p.stopFittingN),
                     'saveFilename': cp(p.saveFilename),
                     'saveKey': cp(p.saveKey),
                     }
        fitProbData[N] = pDataMultiple
        save(fitProbData,fileNumString+'_fitProbData.dat')

        if verbose:
            print("initializeFitAllParallel: Done initializing N =", N)





def runFitAllParallelWorker(fileNumString,endTime=None,verbose=True):
    """
    Each worker node runs this function to look for and perform work.
    
    endTime (None)      : Stop work if endTime hours (wall time) 
                          have elapsed when completing a work unit.  
                          If None, continue indefinitely.
    """

    # check that the fitProbData file exists
    if not fileNumString+"_fitProbData.dat" in os.listdir('.'):
        raise Exception("fitProbData database file not found: "+str(fitProbDatFilename))

    # 9.24.2013 make sure SloppyCell C compiling is working
    if not testCcompiling():
        raise Exception("SloppyCell C compiling not working.")

    if endTime is None: endTime = scipy.inf
    startWallTime = time.time()
    elapsedTimeHours = 0

    while elapsedTimeHours < endTime:
      fitProbData = loadFitProbData(fileNumString)
      saveFilename = list(fitProbData.values())[0]['saveFilename']

      numTimepointsList = scipy.sort(list(fitProbData.keys()))

      # () find a (condition,Np,model) triplet to work on
      conditioni,numTimepointsi,modelj = assignWork(fileNumString)
      numTimepoints = numTimepointsList[numTimepointsi]
      fitProb = loadFitProb(saveFilename,fileNumString,conditioni,numTimepoints)
      
      if verbose:
          print("runFitAllParallelWorker: Assigned work: condition",conditioni,\
            ", numTimepoints",numTimepoints,", model index",modelj)
      
      # set up smallerBestSeenParams
      if (numTimepointsi > 0) and \
         (getState(fitProbData,conditioni,numTimepointsi-1,modelj) == 'finished'):
        smallerFitProb = loadFitProb(saveFilename,fileNumString,conditioni,
                                     numTimepointsList[numTimepointsi-1])
        fitProb.smallerBestParamsDict = paramsDict(smallerFitProb)
      
      # fit the single model
      fitProb.fitAll(maxNumFit=modelj+1)
      
      # save the result in the individual fitProbDict file
      saveFitProb(fitProb,saveFilename,fileNumString,conditioni,numTimepoints)

      # save the result in the more general fitProbData file
      updateFitProbData(fitProb,fileNumString,conditioni,numTimepoints,modelj)

      if verbose:
          print("runFitAllParallelWorker: Finished work.")

      elapsedTimeHours = (time.time() - startWallTime)/3600.


