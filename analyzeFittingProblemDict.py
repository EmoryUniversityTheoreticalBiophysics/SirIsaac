# analyzeFittingProblemDict.py
#
# Bryan Daniels
# 5.10.2012
#

from FittingProblem import *
from plotMatrix import plotMatrix

# the general scheme for calculating over all fpds
def calcForAllFpds(fpdList,func,maxIndex=-3,skip=True,
    addYeastPerfectModel=False,verbose=True):
    dataDict = {}
    for i,fpd in enumerate(fpdList):
        if verbose: print "Fpd",i+1,"of",len(fpdList),":"
        for k in scipy.sort(fpd.keys()):
            fp = fpd[k]
            if addYeastPerfectModel and (fp.perfectModel is None):
                fp.perfectModel = yeastOscillatorFittingModel(fp.indepParamNames)
            mName = fp.maxLogLikelihoodName(maxIndex=maxIndex)
            numDataPts = numDataPoints(fp)
            if (mName is not None) or (not skip):
                result = func(mName,fp)
                if result is not None:
                    if verbose: print "    Done with key",k
                    if not dataDict.has_key(numDataPts): # was using fp.saveKey          
                        dataDict[numDataPts] = []
                    dataDict[numDataPts].append(result)
                else:
                    if verbose: print "    No result.  Skipping key",k
            else:
                if verbose: print "    mName = None.  Skipping key",k
    return dataDict

# 5.2.2013
def split(allFpdsDict):
    """
    Splits an allFpdsDict (eg for bestOutOfSampleCorrs with returnErrors=True)
    """
    N = scipy.shape(allFpdsDict.values()[0][0])[-1]
    allFpdsDictList = []
    for i in range(N):
        d = {}
        for key in allFpdsDict.keys():
            itemList = []
            for item in allFpdsDict[key]:
                #itemList.append(item[i])
                data = (scipy.array(item).T).tolist()
                numItems = len(data)
                if numItems == N:
                    itemList.append(data[i])
                elif numItems == 0: # calculation has not yet been performed
                    itemList.append([])
                else:
                    raise Exception, "Incorrect shape of data"
            d[key] = itemList
        allFpdsDictList.append(d)
    return allFpdsDictList

# 5.31.2012
def numDataPoints(fittingProblem,perDataPoint=True):
    if perDataPoint:
        return len(fittingProblem.fittingData)
    else:
        return scipy.sum([len(d.keys()) for d in fittingProblem.fittingData])

# various things you can calculate
def bestOutOfSampleCorrs(fpdList,type,vars=None,maxIndex=-3,outOfSampleMult=2.,indepParamRanges=None,sampleInLog=False,onlyBest=True,seed=100,timeRange=None,testPerfect=False,returnErrors=True,**kwargs):
    """
    nans (encountered when integration fails) are changed to zeros.
    
    For non-derivative fitting, uses time range of [0,10].
    
    type                    : Can be 'yeast','phos','worm'
    onlyBest (True)         : If False, calculate out-of-sample correlations for all
                              fit models, not just the max likelihood one.
    testPerfect (False)     : If True, test fp.perfectModel.  Overrides onlyBest.
    returnErrors (True)     : If True, also return mean squared errors 
                              (only works for certain cases)
    """
    if hasattr(fpdList[0].values()[0],'fittingDataDerivs'):
        useDerivs = (fpdList[0].values()[0].fittingDataDerivs is not None)
    else:
        useDerivs = False
    if type is 'yeast': #(indepParamRanges is None) and hasattr(fpdList[0].values()[0].perfectModel,'typicalIndepParamRanges'):
        if indepParamRanges is None:
            dummyYeastModel = yeastOscillatorFittingModel(fpdList[0].values()[0].indepParamNames)
            indepParamRanges = dummyYeastModel.typicalIndepParamRanges(outOfSampleMult)
        if vars is None: vars = ['S1','S2','S3']
        addYeastPerfectModel = True #False #True
    elif type is 'phos': #if (indepParamRanges is None):
        indepParamRanges = [[10**-3,10**3]]
        sampleInLog = True
        if vars is None: vars = ['totalPhos']
        addYeastPerfectModel = False
    elif type is 'worm': 
        indepParamRanges = []
        sampleInLog = False
        if vars is None: vars = ['wormSpeed']
        addYeastPerfectModel = False
    else:
        raise Exception, "Unrecognized type: "+str(type)

    def modelCorrsFunc(fp,model):
        c = fp.outOfSampleCorrelation(model,[0,10],vars,indepParamRanges,numTests=100,verbose=False,sampleInLog=sampleInLog,seed=seed,returnErrors=returnErrors)
        if returnErrors:
            return scipy.mean(scipy.nan_to_num(c[0])),scipy.mean(c[1])
        else:
            return scipy.mean(scipy.nan_to_num(c))

    if testPerfect: # 4.17.2013 test fitted perfectModels
      if useDerivs: raise Exception, "testPerfect + useDerivs not yet supported."
      corrsFunc = lambda mName,fp: modelCorrsFunc(fp,fp.perfectModel)
    elif onlyBest:
      if not useDerivs: # not fitting derivatives
        corrsFunc = lambda mName,fp: modelCorrsFunc(fp,fp.fittingModelDict[mName])
      else: # fitting derivatives 1.10.2013
        corrsFunc = lambda mName,fp: scipy.nan_to_num( fp.outOfSampleCorrelation_deriv(fp.fittingModelDict[mName],vars,indepParamRanges,numTests=100,sampleInLog=sampleInLog,seed=seed) ) 
      skip = True
    else: # 7.25.2012 calculate for all fit models
      if not useDerivs: # not fitting derivatives
        corrsFunc = lambda mName,fp: \
          [ modelCorrsFunc(fp,fp.fittingModelDict[name]) for name in orderedFitNames(fp) ]
      else: # fitting derivatives 1.10.2013
        if sampleInLog: raise Exception, "sampleInLog not supported"
        indepParamsSeed = seed
        scipy.random.seed(seed)
        timeSeed = scipy.random.randint(1e6)
        corrsFunc = lambda mName,fp: [ scipy.nan_to_num( fp.outOfSampleCorrelation_deriv(fp.fittingModelDict[name],vars,indepParamRanges,numTests=100,timeSeed=timeSeed,indepParamsSeed=indepParamsSeed,timeRange=timeRange) ) for name in orderedFitNames(fp) ]
      skip = False
    return calcForAllFpds(fpdList,corrsFunc,maxIndex=maxIndex,addYeastPerfectModel=addYeastPerfectModel,**kwargs)

def bestModelCostPerMeasurement(fpdList,testPerfect=False,onlyBest=True,**kwargs):
    """
    testPerfect (False)     : If True, look for fp.perfectCost.
    
    Returns 2*cost/numMeasurements.  (2*cost should be chi^2.)
    """
    if testPerfect:
        def costFunc(mName,fp):
            if hasattr(fp,'perfectCost'):
                return 2.*fp.perfectCost/len(fp.fittingData)
            else:
                print "bestModelCostPerMeasurement: No perfectCost.  Returning None."
                return None
        kwargs['skip'] = False # don't skip if we haven't fit the fittingModels
    elif onlyBest:
        costFunc = lambda mName,fp: 2.*fp.costDict[mName]/len(fp.fittingData)
    else: # calculate for all models
        costFunc = lambda mName,fp: \
            [ 2.*fp.costDict[name]/len(fp.fittingData) for name in orderedFitNames(fp) ]
    return calcForAllFpds(fpdList,costFunc,**kwargs)

def bestModelNumParams(fpdList,**kwargs):
    numParamsFunc = lambda mName,fp: fp.numParametersDict[mName]
    return calcForAllFpds(fpdList,numParamsFunc,**kwargs)
    
def bestModelEffectiveNumParams(fpdList,**kwargs):
    effectiveNumParamsFunc = lambda mName,fp: fp.numStiffSingValsDict[mName]
    return calcForAllFpds(fpdList,effectiveNumParamsFunc,**kwargs)
    
def logLikelihoods(fpdList,**kwargs):
    llFunc = lambda mName,fp: [ fp.newLogLikelihoodDict[name] for name in orderedFitNames(fp) ]
    return calcForAllFpds(fpdList,llFunc,**kwargs)
    
def totalNumFunctionCalls(fpdList,**kwargs):
    """
    Sum of all cost calls plus grad calls over ALL models tested.
    """
    totalFuncCallsFunc = lambda mName,fp: scipy.sum( [ scipy.sum(fp.fittingModelDict[name].numCostCallsList) + scipy.sum(fp.fittingModelDict[name].numGradCallsList) for name in fp.newLogLikelihoodDict.keys() ] )
    return calcForAllFpds(fpdList,totalFuncCallsFunc,skip=False,**kwargs)

# updated 10.10.2013
def totalNumEvaluations(fpdList,stopFittingN=scipy.inf,**kwargs):
    """
    Total number of times daeint was called.
    
    Sum of
       (cost calls) * (number of measurements)
       (grad calls + ensemble steps) * (1 + number of measurements)
    over ALL models tested.
    
    (I'm not sure that the '1 + ' is necessary for the gradient
    calculation, but that seems to be the number of times
    SloppyCell calls daeint.)
    
    stopFittingN (inf)          : Only include models up to stopFittingN
                                  past the maxLogLikelihoodName
    """
    def totalFuncCallsFunc(mName,fp):
        if hasattr(fp,'newLogLikelihoodDict'):
            keyList = orderedFitNames(fp,stopFittingN=stopFittingN)
        else:
            keyList = []
        ND = len(fp.fittingData)
        mList = [ fp.fittingModelDict[name] for name in keyList ]
        return scipy.sum( [                                             \
           ND * scipy.sum(m.numCostCallsList)                           \
         + ND*(1+len(m.getParameters())) *                              \
            ( scipy.sum(m.numGradCallsList) + m.ensGen.totalSteps )
         for m in mList ] )
    return calcForAllFpds(fpdList,totalFuncCallsFunc,skip=False,**kwargs)
    
def totalWallTimesHours(fpdList,includedNames=None,**kwargs):
    """
    Sum of all ensemble times plus minimization times over ALL models tested.
    
    (Time taken when running in series.)
    """
    if includedNames is None: 
        includedNamesFunc = lambda name: True
    else:
        includedNamesFunc = lambda name: name in includedNames
    totalWallTimeFunc = lambda mName,fp: scipy.sum( [ includedNamesFunc(name)*( scipy.sum(fp.fittingModelDict[name].ensTimeSecondsList) + scipy.sum(fp.fittingModelDict[name].minimizationTimeSecondsList) ) for name in fp.newLogLikelihoodDict.keys() ] )/3600.
    return calcForAllFpds(fpdList,totalWallTimeFunc,skip=False,**kwargs)

def parallelWallTimesHours(fpdList,includedNames=None,**kwargs):
    """
    Sum of all ensemble times plus max of minimization times over ALL models tested.
    
    (Time taken when running in parallel.)
    """
    if includedNames is None: 
        includedNamesFunc = lambda name: True
    else:
        includedNamesFunc = lambda name: name in includedNames
    totalWallTimeFunc = lambda mName,fp: scipy.sum( [ includedNamesFunc(name)*( scipy.sum(fp.fittingModelDict[name].ensTimeSecondsList) + max(fp.fittingModelDict[name].minimizationTimeSecondsList) ) for name in fp.newLogLikelihoodDict.keys() ] )/3600.
    return calcForAllFpds(fpdList,totalWallTimeFunc,skip=False,**kwargs)


def totalEnsembleTimesHours(fpdList,includedNames=None,**kwargs):
    """
    Sum of all ensemble times over ALL models tested.
    """
    if includedNames is None: 
        includedNamesFunc = lambda name: True
    else:
        includedNamesFunc = lambda name: name in includedNames
    totalEnsTimeFunc = lambda mName,fp: scipy.sum( [ includedNamesFunc(name)*( scipy.sum(fp.fittingModelDict[name].ensTimeSecondsList) ) for name in fp.newLogLikelihoodDict.keys() ] )/3600.
    return calcForAllFpds(fpdList,totalEnsTimeFunc,skip=False,**kwargs)
    
def totalMinimizationTimesHours(fpdList,**kwargs):
    """
    Sum of all minimization times over ALL models tested.
    """
    totalMinTimeFunc = lambda mName,fp: scipy.sum( [ scipy.sum(fp.fittingModelDict[name].minimizationTimeSecondsList) for name in fp.newLogLikelihoodDict.keys() ] )/3600.
    return calcForAllFpds(fpdList,totalMinTimeFunc,skip=False,**kwargs)
    
def totalMinimizationTimesSeconds(fpdList,**kwargs):
    """
    Sum of all minimization times over ALL models tested.
    """
    totalMinTimeFunc = lambda mName,fp: scipy.sum( [ scipy.sum(fp.fittingModelDict[name].minimizationTimeSecondsList) for name in fp.logLikelihoodDict.keys() ] )
    return calcForAllFpds(fpdList,totalMinTimeFunc,skip=False,**kwargs)
    
def bestModelConvFlags(fpdList,maxIndex=-3,**kwargs):
    convFlagFunc = lambda mName,fp: fp.fittingModelDict[mName].convFlag
    return calcForAllFpds(fpdList,convFlagFunc,maxIndex=maxIndex,**kwargs)
    
def bestModelMeanAbsGradient(fpdList,**kwargs):
    """
    Mean magnitude of the gradient with respect to parameters.
    (Should be small if we've converged)
    """
    meanAbsGradFuncByName = lambda mName,fp: meanAbsGradFunc(fp.fittingModelDict[mName],fp)
    return calcForAllFpds(fpdList,meanAbsGradFuncByName,**kwargs)

# 9.26.2012
def bestIndices(fpdList,**kwargs):
    """
    Index of best ensemble member for all tested models.
    
    As of 9.26.2012, the order of ensemble starting points is:
        Index 0 : Fit parameters from model with fewer data points,
                  if using smallerBestParamsDict.
        Index -1: Initial starting point of Monte Carlo ensemble steps.
        Others  : Increasing index means fewer Monte Carlo
                  steps from initial starting point.
    """
    bestIndex = lambda mName,fp: [ fp.fittingModelDict[n].bestIndex for n in orderedFitNames(fp) ]
    return calcForAllFpds(fpdList,bestIndex,**kwargs)

# 9.26.2012
def acceptanceRatios(fpdList,**kwargs):
    a = lambda mName,fp: [ fp.fittingModelDict[n].acceptanceRatio for n in orderedFitNames(fp) ]
    return calcForAllFpds(fpdList,a,**kwargs)

def meanAbsGradFunc(m,fp):
    m2 = m._SloppyCellDataModel(fp.fittingData,fp.indepParamsList)
    res = m2.res(m.getParameters())
    j = scipy.array(m2.jacobian_sens(m.getParameters()))
    grad = scipy.mat(res)*scipy.mat(j)
    return scipy.mean(abs(2.*grad))
    
def plotAllFpdsDict(dataDict,marker='o',ls='',color='b',label=None,     \
    plotMeans=False,errorBars=True,filterNans=False,                    \
    returnData=False,makePlot=True,**kwargs):
    """
    filterNans (False)      : If True, NaNs are removed from the data
                              before taking the mean.
    """
    kwargs['label'] = label
    kList,meanList,stdList = [],[],[]
    for i,k in enumerate(scipy.sort(dataDict.keys())):
        vals = dataDict[k]
        if i==1: # include label for legend only once
            kwargs.pop('label')
        if plotMeans:
            if filterNans: fVals = filter(lambda v: not scipy.isnan(v),vals)
            else: fVals = vals
            meanVal,stdVal = scipy.mean(fVals),scipy.std(fVals,ddof=1)
            if makePlot:
              if errorBars:
                pylab.errorbar([k],[meanVal],yerr=[stdVal],             \
                    marker=marker,ls=ls,color=color,**kwargs)
              else:
                pylab.plot([k],[meanVal],marker=marker,ls=ls,           \
                    color=color,**kwargs)
            kList.append(k)
            meanList.append(meanVal)
            stdList.append(stdVal)
            
        else: # scatter plot of all values
            if makePlot:
              pylab.plot(scipy.repeat(k,len(vals)),vals,marker=marker,  \
                ls=ls,color=color,**kwargs)
            #kList.extend(scipy.repeat(k,len(vals)))
            #meanList.extend(vals)
            #stdList.extend(scipy.repeat(None,len(vals)))
            kList.append(k)
            meanList.append(vals)
            stdList.append(None)

            
    if returnData:
        return kList,meanList,stdList
        
# 7.25.2012
def orderedFitNames(fp,stopFittingN=scipy.inf):
    """
    stopFittingN (inf)      : Only include names of models up to
                              stopFittingN past the maxLogLikelihoodName
    """
    # make list of names
    names = []
    if hasattr(fp,'newLogLikelihoodDict'):
      for name in fp.fittingModelNames:
        if fp.newLogLikelihoodDict.has_key(name):
            names.append(name)
    #print "orderedFitNames: debug: ",names

    # determine how many we should return based on stopFittingN
    if stopFittingN < scipy.inf:
        maxName = fp.maxLogLikelihoodName()
    else:
        maxName = None
    try:
        n = names.index(maxName) + stopFittingN + 1
    except ValueError: # maxName wasn't in the list
        n = len(names)

    return names[:n]

# 7.25.2013
def plotAllFpdsDict2D(dataDict,fpdList=None,newFigure=True,defaultValue=0,\
    index=0,**kwargs):
    """
    matrix plot, Number of measurements vs. model number
    
    index (0)               : Index of fittingProblem to plot
    fpdList (None)          : If given, use to indicate the selected model
    """
    maxNumModels = max( [ max([ len(data) for data in d ]) for d in dataDict.values() ] )
    mat = defaultValue * scipy.ones((maxNumModels,len(dataDict.keys())))
    for i,key in enumerate( scipy.sort(dataDict.keys()) ):
        data = dataDict[key][index]
        for j,val in enumerate(data):
            mat[j,i] = val

    plotMatrix(mat,X=scipy.sort(dataDict.keys()),Y=range(maxNumModels),**kwargs)
                
    pylab.xlabel('Number of measurements $N$')
    pylab.ylabel('Model index')

    # indicate the selected model
    if fpdList is not None:
        fpd = fpdList[index]
        for key in scipy.sort(fpd.keys()):
            fp = fpd[key]
            name = fp.maxLogLikelihoodName()
            if name is not None:
                bestModelIndex = fp.fittingModelNames.index(name)
                pylab.plot([key],[bestModelIndex],'ro')
    
        
            

# 1.24.2013
def plotOutOfSampleCorrelationVsMeasurements(fpdList,varList,seed=100,
    newFigure=True,ls=':'):
    """
    Plots mean over fpds in fpdList.
    """
    cW = Plotting.ColorWheel()
    if newFigure:
        pylab.figure()
    for var in varList:
        d = bestOutOfSampleCorrs(fpdList,vars=[var],seed=seed)
        keyList = list( scipy.sort(d.keys()) )
        corrs = [ scipy.mean(d[k]) for k in keyList ]
        # remove exact zeros (which come from evaluation errors)
        while 0. in corrs:
            i = corrs.index(0.)
            corrs.pop(i)
            keyList.pop(i)
        # keep same colors each time you plot
        cWfmt = cW.next()
        if newFigure: label = var
        else: label = '_nolegend_'
        pylab.plot(keyList,corrs,marker=cWfmt[1],ls=ls,color=cWfmt[0],label=label)
    
    # make pretty
    pylab.xlabel('Number of derivative data points')
    pylab.ylabel('Out-of-sample correlation')
    pylab.legend(loc=4)