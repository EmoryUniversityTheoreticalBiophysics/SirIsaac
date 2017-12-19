# simulateYeastOscillator.py
# 
# Bryan Daniels
# 8.19.2011
#
# Runs yeast oscillator using Abhishek Soni's MATLAB code.

from MATLAB import *

import scipy
import scipy.io
import os
import pylab # for pylab.find
import simplePickle
import copy

#codeDir = "~/Research/SloppySimplification/yeastOscillator"
codeDir = "./yeastOscillator"

def simulateYeastOscillator(times,temperature,returnParams=True,    \
    initialConditions=None):
    """
    Returns times,14-dimensional data for each time (7 concentrations,7 derivatives).
    (S1,S2,S3,S4,N2,A3,S4ex) and their derivatives in the same order
    (glucose,GAP/DAP,BPG,Pyruvate/Acetaldehyde,NADH,ATP,Ext. Pyruvate/Acetaldehyde)
    
    returnParams (True)         : Also return the values of the
                                  parameters, some of which are
                                  temperature-dependent.
                                  (J0,k1,k2,k3,k4,k5,k6,k,kappa,q,K1,psi)
    initialConditions (None)    : Optionally pass a 7-dimensional vector
                                  consisting of initial conditions for 
                                  the 7 concentrations.  If None, use
                                  default values.
    """
    randFileNo = str(os.getpid())
    
    if initialConditions is None:
        # 2.22.2012 took from run_ruoff_model_python_interface
        #% initial conditions from run_Ruoff_model.m
        #% Set initial conditions as per Table 2 in the paper. 
        #y0(1) = 1.187;  % S1; 
        #y0(2) = 0.193;  % S2;
        #y0(3) = 0.050;  % S3;
        #y0(4) = 0.115;  % S4;
        #y0(5) = 0.077;  % N2;
        #y0(6) = 2.475;  % A3;
        #y0(7) = 0.077;  % S4ex;
        initialConditions =                                     \
            [1.187,0.193,0.050,0.115,0.077,2.475,0.077]
    elif len(initialConditions) != 7:
        raise Exception,'initialConditions should have length 7.'
    
    outputFilenameList =                                                    \
      ["run_ruoff_model_python_interface_output"+randFileNo+".txt",         \
       "run_ruoff_model_python_interface_output_params"+randFileNo+".txt"]
    
    # 1.10.2013 for calling with len(times)=1
    desiredTimes = None
    if len(times)==1:
        if times[0] == 0.:
            times = [0.,1.e-5]
            desiredTimes = [0]
        else:
            times = [0.,times[0]]
            desiredTimes = [-1]
    
    timesFilename = "times_temp"+randFileNo+".data"
    scipy.savetxt(timesFilename,times)
    
    initialConditionsCallString =                               \
        ''.join([','+str(ic) for ic in initialConditions])
    
    callStr = "run_ruoff_model_python_interface('"              \
        +randFileNo+"',"+str(temperature)+                      \
        initialConditionsCallString+");"
    
    try:
        allData,params = callMATLAB(callStr,                    \
            outputFilenameList,codeDirList=[codeDir])
        os.remove(timesFilename)
    except:
        os.remove(timesFilename)
        raise
    
    allData = scipy.transpose(allData)
    
    # 1.10.2013
    if desiredTimes is None: desiredTimes = range(len(allData[0]))
    
    if returnParams:
        return allData[0,desiredTimes],allData[1:,desiredTimes],params
    else:
        return allData[0,desiredTimes],allData[1:,desiredTimes]
    

# 1.9.2013 taken from runFittingProblem.py
def yeastData(numPoints,timeInterval,                                   \
              numICs,useDerivs=True,includedIndices=range(7),           \
              timesSeed=3,noiseSeed=4,ICseed=5,                         \
              multiplicativeErrorBar=0.1,upperRangeMultiple=1.,         \
              randomX=True,inputVars=None,nonzeroMin=True,              \
              yeastSavedDataFile=None,includeNoise=True):
    """
    upperRangeMultiple (1.)         : Each range of initial conditions is
                                      expanded by this factor by increasing
                                      the upper limit.
    nonzeroMin (True)               : Use the lower limit for initial
                                      conditions given in SchValJen11;
                                      otherwise use 0. as lower limit.
    """
    usePreloadedData = True #<<<<<<<*******************
    
    allNames = scipy.array(['S1','S2','S3','S4','N2','A3','S4ex']) 
    names = allNames[includedIndices]
    varList = names
    
    defaultTemperature = 286.5 #changed 10.16.2013 (was 288 K before)
    
    # Table 2 of RuoChrWol03
    defaultICs = scipy.array([1.187,0.193,0.050,0.115,0.077,2.475,0.077]) # mM
    
    # choose input variable(s)
    if False: # use temp. as input
        inputVars = ['Temperature']
        inputList = [[278],[288],[293]]
        inputDescriptors = [inputs[0] for inputs in inputList]
    elif False: # only a single temperature value
        inputVars = ['Temperature']
        inputList = [[defaultTemperature]]
        inputDescriptors = [inputs[0] for inputs in inputList]
    elif inputVars is None: # use varying initial conditions on all 7 species
        inputVars = [name+"_init" for name in names]  
        print "inputVars =",inputVars      
        # taken from SchValJen11 Table 2
        ICranges = scipy.array(                                             \
                   [[0.15,1.60],[0.19,2.16],                                \
                    [0.04,0.20],[0.10,0.35],                                \
                    [0.08,0.30],[0.14,2.67],[0.05,0.10]] )[includedIndices] # mM
        
        # as I vary the number of ICs, I want the same sequence of random ICs
        scipy.random.seed(ICseed)
        randICmults = []
        for i in range(numICs):
            randICmults.append(scipy.rand(len(includedIndices)))
        
        #randomICs = scipy.rand(numICs,len(includedIndices))*                \
        #    (ICranges[:,1]-ICranges[:,0]) + ICranges[:,0]
        randomICs = scipy.array(randICmults)*upperRangeMultiple*             \
            (ICranges[:,1]-nonzeroMin*ICranges[:,0]) + nonzeroMin*ICranges[:,0]
        inputList = randomICs
        inputDescriptors = range(len(inputList))
    else:
        raise Exception, "Changing inputVars not yet implemented"
    
    
    # as I vary the number of ICs, I want the same sequence of lists of times
    scipy.random.seed(timesSeed)
    timesRandList = []
    for i in range(numICs):
        timesRandList.append(scipy.rand(numPoints))
    
    # all other random stuff is done now except for the noise
    # as I vary the number of ICs, I want the same sequence of noise
    scipy.random.seed(noiseSeed)
    
    # set up preloaded data dict
    if yeastSavedDataFile is None:
        yeastSavedDataFile = "yeastExampleData.fittingData"
    yeastDict = {}
    if usePreloadedData:
        try:
            #yeastTimes,loadedICs,yeastRawData =                             \
            #    simplePickle.load(yeastDataFilename)
            yeastDict = simplePickle.load(yeastSavedDataFile)
        except IOError:
            print "yeastData: Could not load "+yeastSavedDataFile
            usePreloadedData = False
        

    yeastOscillatorData = {}
    yeastOscillatorDataDerivs = {}
    for inputs,inputDescriptor,timesRand in zip(inputList,inputDescriptors,timesRandList):
        stopAfterRunningMATLAB = False
        
        # deal with independent parameters
        # if it's temperature (only)
        if inputVars is ['Temperature']:
            temperature = inputs[0]
        else:
            temperature = defaultTemperature
        
        initialConditions = copy.deepcopy(defaultICs)
        for inputVar,input in zip(inputVars,inputs):
            # if it's initial conditions
            if inputVar[-5:] == "_init":
                index = pylab.find(allNames==inputVar[:-5])[0]
                initialConditions[index] = input
        print "initialConditions =",initialConditions
        
        eps = (timeInterval[1]-timeInterval[0])/1000. #1e-3 #1e-4 # (minutes) resolution of timepoints
        integratorTimes = scipy.arange(timeInterval[0],timeInterval[1]+eps,eps)
        if not randomX:
            desiredTimes = scipy.linspace(timeInterval[0],timeInterval[1],numPoints)
        else:
            desiredTimes = timesRand * (timeInterval[1]-timeInterval[0])        \
                + timeInterval[0]

        yeastDataKey = copy.deepcopy(                                           \
            (tuple(initialConditions),inputDescriptor,tuple(desiredTimes)))
        
        runMATLAB = not usePreloadedData
        if usePreloadedData:
                if yeastDict.has_key(yeastDataKey):
                    yeastTimes,loadedICs,yeastRawData = yeastDict[yeastDataKey]
                    if not scipy.all(loadedICs == initialConditions):
                        print "loadedICs =",loadedICs
                        print "initialConditions =",initialConditions
                        raise Exception, "loadedICs != initialConditions"
                else:
                    print "yeastData: "+yeastSavedDataFile+" does not "         \
                        "contain the necessary data."
                    runMATLAB = True
        if runMATLAB:
            print "yeastData: Running simulateYeastOscillator."
            yeastTimes,yeastRawData,yeastParams =                               \
                simulateYeastOscillator(integratorTimes,temperature,            \
                                        initialConditions=initialConditions)
            # the integrator returns more timepoints than I want
            desiredIndices = [ pylab.find(abs(yeastTimes-time)<eps/2)[0]        \
                              for time in desiredTimes ]
            yeastRawData = yeastRawData[:,desiredIndices]
            yeastTimes = yeastTimes[:,desiredIndices]
            
            yeastDict[yeastDataKey] = (yeastTimes,initialConditions,yeastRawData)

        if False: # OLD! use temp. as input, calculate v1 as output, mult. error bars
            S1data = yeastRawData[0]
            A3data = yeastRawData[5]
            k1,K1,q = yeastParams[1],yeastParams[10],yeastParams[9]
            v1data = k1*S1data*A3data/(1.+(A3data/K1)**q)
            yeastOscillatorData[inputDescriptor] = {}
            yeastOscillatorData[inputDescriptor]['v1'] =                        \
                dict( zip(times, zip(v1data,multiplicativeErrorBar*v1data)) )
            varList = ['v1']
        
        if useDerivs: # use derivatives as output (7-dimensional), const. error bars
            # 2.22.2012 from Table 2 of SchValJen11
            # 11.29.2012 using these even though they are stddevs of values, not derivs
            stddevs = [0.4872,0.6263,0.0503,0.0814,0.0379,0.7478,0.0159] # mM
            stddevs = scipy.array(stddevs)[includedIndices]
            
            #derivNames = [ 'ddt_'+name for name in names ]
            derivDataList = yeastRawData[7:][includedIndices]
            yeastOscillatorDataDerivs[inputDescriptor] = {}
            for derivName,derivData,stddev in zip(names,derivDataList,stddevs):
                constErrorBar = multiplicativeErrorBar*stddev
                yeastOscillatorDataDerivs[inputDescriptor][derivName] =         \
                    dict( zip(yeastTimes,                                         \
                              zip(derivData,scipy.ones_like(derivData)*constErrorBar)) )
        #varList = names
        
        # 10.3.2011
        if True: # use all species as output (7-dimensional), const. error bars
            # 2.22.2012 from Table 2 of SchValJen11
            stddevs = [0.4872,0.6263,0.0503,0.0814,0.0379,0.7478,0.0159] # mM
            
            #includeNoise = False #True # <<<<****** 12.5.2012 to avoid negative values
            verbose = True
            
            dataList = scipy.array(yeastRawData)[:7][includedIndices]
            stddevs = scipy.array(stddevs)[includedIndices]
            
            # 11.17.2011 reorder stuff 
            #names,dataList = names[::-1],dataList[::-1]
            
            yeastOscillatorData[inputDescriptor] = {}
            for name,data,stddev in zip(names,dataList,stddevs):
                #constErrorBar = multiplicativeErrorBar*scipy.mean(data) #*data[0]
                constErrorBar = multiplicativeErrorBar*stddev
                if includeNoise:
                    noise = scipy.random.normal(scale=constErrorBar,size=len(data))
                    data = data + noise
                yeastOscillatorData[inputDescriptor][name] =                    \
                    dict( zip(yeastTimes,                                         \
                              zip(data,scipy.ones_like(data)*constErrorBar )) )
    
    #if usePreloadedData:
    # save data for future use via usePreloadedData
    simplePickle.save(yeastDict,yeastSavedDataFile)
    
    #Plotting.plot(yeastTimes,v1data,'o-',label=str(temperature))
    
    if (not usePreloadedData) and stopAfterRunningMATLAB: die
    
    fittingData = [ yeastOscillatorData[d] for d in inputDescriptors ]
    if useDerivs:
        fittingDataDerivs = [ yeastOscillatorDataDerivs[d] for d in inputDescriptors ]
    else:
        fittingDataDerivs = None
    return fittingData,fittingDataDerivs,inputVars,inputList


