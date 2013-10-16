# runPowerLawDerivFit.py
#
# Bryan Daniels
# 6.28.2012
#
# Trying for a proof-of-concept of 
# alternating log-linear power-law fitting
# using 7-dimensional yeast oscillator model.
#
# ********************************************
# Outdated as of 10.16.2013, since we're now
# using parameters that match SchValJen11
# instead of using T = 288 K
# (see Ruoff_model_original.m).
# ********************************************
#

import FittingProblem
reload(FittingProblem)
from simulateYeastOscillator import *
from simplePickle import load,save
import pylab

# taken from runFittingProblem.yeastDataFunction
defaultICs = scipy.array([1.187,0.193,0.050,0.115,0.077,2.475,0.077]) # mM
ICranges = scipy.array(                                             \
    [[0.15,1.60],[0.19,2.16],                                       \
     [0.04,0.20],[0.10,0.35],                                       \
     [0.08,0.30],[0.14,2.67],[0.05,0.10]] ) # mM

# load/produce yeast derivative data
loadData = True
#dataFilename = '120719_yeastSpeciesAndDerivsData_100_numSpecies'+str(numSpecies)+'.data'
ICseed = 0
timesSeed = 1

# set up random initial conditions
numICs = 10
scipy.random.seed(ICseed)
randICmults = []
for i in range(numICs):
    randICmults.append(scipy.rand(7))
ICsList = scipy.array(randICmults)*                                 \
    2.*(ICranges[:,1]-ICranges[:,0]) #+ ICranges[:,0]

# set up random measurement times
maxTime = 10. #1e-5 #10.
scipy.random.seed(timesSeed)
randTimes = scipy.rand(numICs)*maxTime
singleTimepoint = True # single time point for each IC

useYeast = True
usePhosphorylation = False
makeFittingData = True # use False if you're getting data from a previous run

if useYeast:
  numVisibleSpecies = 3
  numSpecies = 4
  dataFilename = '120725_yeastSpeciesAndDerivsData_10_numSpecies'+str(numVisibleSpecies)+'.data'
  visibleSpeciesNames = ['S1','S2','S3','S4','N2','A3','S4ex'][:numVisibleSpecies]
  indepParamNames = [ n+'_init' for n in visibleSpeciesNames ]
  if loadData:
    speciesData,derivsData = load(dataFilename)
  elif not loadData:
    # set up list of random initial conditions
    # taken from runFittingProblem.yeastDataFunction

    # run MATLAB simulator
    speciesData = []
    derivsData = []
    i = 0
    for ICs,time in zip(ICsList,randTimes):
        i += 1
        print "runPowerLawDerivFit: Running yeast simulation",i,"of",len(ICsList)
        T = 288 # K
        initialConditions = defaultICs
        initialConditions[:numSpecies] = ICs[:numSpecies]
        if singleTimepoint: # single time point for each IC
            times,yeastRawData,params = simulateYeastOscillator([0,time],T,         \
                initialConditions=initialConditions)
            speciesData.append( yeastRawData[:7,-1] )
            derivsData.append( yeastRawData[-7:,-1] )
        else: # 7.25.2012 keep all timepoints
            times,yeastRawData,params = simulateYeastOscillator([0,maxTime],T,      \
                initialConditions=initialConditions)
            speciesData.extend( scipy.transpose(yeastRawData[:7]) )
            derivsData.extend( scipy.transpose(yeastRawData[-7:]) )
    if True:
        speciesData = (scipy.transpose(speciesData)[:numSpecies]).T
        derivsData = (scipy.transpose(derivsData)[:numSpecies]).T
    else:
        speciesData = speciesData[:numSpecies]
        derivsData = derivsData[:numSpecies]
    save((speciesData,derivsData),dataFilename)
    
    # 7.19.2012 try using all data from single IC
    #speciesData = yeastRawData[:numSpecies]
    #derivsData = yeastRawData[-7:(-7+numSpecies)]
  
  offset = 10 # 8.22.2012 staying away from zero is better?
  indepParamsList = ICsList + offset
  speciesData = speciesData + 10
  
  # (one data point for each IC)
  sigmas = 0.1*scipy.ones_like(speciesData)
  fittingData =       [ dict( zip( visibleSpeciesNames,                         \
      [ dict( [(randTime,(data,sigma))] )                                       \
                for data,sigma in zip(dataRow,sigmaRow)] ) )                    \
            for randTime,dataRow,sigmaRow in zip(randTimes,speciesData,sigmas) ]
  fittingDataDerivs = [ dict( zip( visibleSpeciesNames,                         \
      [ dict( [(randTime,data)] )                                               \
                for data in dataRow ] ) )                                       \
            for randTime,dataRow in zip(randTimes,derivsData) ]

elif not useYeast:
    times = scipy.linspace(0,10,100)
    numSpecies = 3
    indepParamsList = [[]]
    indepParamNames = []
    
    if True:    
        # 7.10.2012 try silly simple case
        visibleSpeciesNames = [ 'Species0' ] #,'Species1'
        speciesData = scipy.array([ 1. - scipy.exp(-times) + scipy.exp(-2.*times) ])#, scipy.exp(-times) ])
        derivsData = scipy.array([ scipy.exp(-times) - 2.*scipy.exp(-2*times) ])#, -scipy.exp(-times) ])
    elif False:
        # 7.19.2012 oscillatory test
        visibleSpeciesNames = [ 'Species0']# ,'Species1' ]
        speciesData = scipy.array( [ 4.+scipy.sin(times), 4.+scipy.cos(times) ] )
        derivsData = scipy.array( [ scipy.cos(times), -scipy.sin(times) ] )
    elif usePhosphorylation:
        # 8.23.2012 phosphorylation
        dir = '/Users/bdaniels/Dropbox/Private/SloppySimplification/'
        fpd = load(dir+'0040_fitProb_varying_numInputs_PhosphorylationNet_CTSN_withEnsembleT1000_steps10000.0_10_useBest_numPoints1_maxiter100_avegtol0.01_noClamp_newErrorBars_removeLogForPriors_seeds0_1.dat')
        fp = fpd[12]
        fittingData = fp.fittingData
        visibleSpeciesNames = ['totalPhos']
        numVisibleSpecies = len(visibleSpeciesNames)
        makeFittingData = False
        indepParamsList = fp.indepParamsList
        indepParamNames = fp.indepParamNames
        # calculate derivatives (numerically for now....)
        fittingDataDerivs = []
        for indepParams,conditionData in zip(fp.indepParamsList,fittingData):
            conditionDerivs = {}
            for var in conditionData.keys():
                varDerivs = {}
                for time in conditionData[var].keys():
                    delta = 1e-5
                    val1,val2 = fp.perfectModel.evaluateVec(                    \
                        [time-delta,time+delta],var,indepParams)
                    deriv = (val2-val1)/(2.*delta)
                    varDerivs[time] = deriv
                conditionDerivs[var] = varDerivs
            fittingDataDerivs.append(conditionDerivs)
        
    if makeFittingData:
        # (many data points for a single IC)
        sigmas = 0.1*scipy.ones_like(speciesData)
        fittingData =       [ dict( zip( visibleSpeciesNames,                       \
            [ dict( zip(times,zip(data,sigma)) )                                    \
                for data,sigma in zip(speciesData,sigmas)] ) ) ]
        fittingDataDerivs = [ dict( zip( visibleSpeciesNames,                       \
            [ dict( zip(times,data) ) for data in derivsData ] ) ) ]

    
#speciesNames = visibleSpeciesNames + [ 'X_'+str(i) for i in range(numSpecies-len(visibleSpeciesNames)) ]

    
# () set up power-law fitting model
# want fully connected numSpecies-dimensional model
#networkList = [ [ 5,dict([(j,2) for j in range(i)]+[(j,2)                           \
#                  for j in range(i+1,numSpecies)]) ] for i in range(numSpecies) ]
p = FittingProblem.PowerLawFittingModel_FullyConnected(numSpecies,                  \
    outputNames=visibleSpeciesNames,indepParamNames=indepParamNames)

# () do the fit
numLinearIter = 1 #10 #20 #4e2
maxiter = 100 #100
paramsList = p.fitToDataDerivs(fittingData,fittingDataDerivs,                  \
                  numLinearIter=int(numLinearIter),                                 \
                  verbose=True,maxiter=maxiter,seed=1,                              \
                  indepParamsList=indepParamsList)
                  


# () compare integrated behavior
colors = ['g','b','r','k','y'][:len(visibleSpeciesNames)] + list( scipy.repeat('gray',numSpecies-len(visibleSpeciesNames)) )
T = 288 #K
pylab.figure()
if useYeast: # for yeast
    # set initial conditions to the last example
    p.net.set_var_ics(dict(zip(p.speciesNames,indepParamsList[-1])))
    times,yeastRawData,params = simulateYeastOscillator([0,10],T,initialConditions=ICsList[-1])
    #pylab.plot(times,speciesData[0])
    for i,color in zip(range(numVisibleSpecies),colors):
        pylab.plot(times,yeastRawData[i]+offset,color=color,ls=':')
elif usePhosphorylation:
    # 8.23.2012 how should I set the initial condition more generally?
    #defaultIC = 1.
    #p.net.set_var_ics(dict(zip(visibleSpeciesNames,defaultIC*scipy.ones(len(visibleSpeciesNames)))))
    #p.net.setOptimizables(dict(zip(fp.indepParamNames,indepParamsList[-1])))
    #perfectData = fp.perfectModel.evaluateVec(times,visibleSpeciesNames,indepParamsList[-1])
    #for i,color in zip(range(numVisibleSpecies),colors):
    #    pylab.plot(times,perfectData[i],color=color,ls=':')
    p.plotResults(fittingData,indepParamsList)
else: # simple examples
    p.net.set_var_ics(dict(zip(p.speciesNames,speciesData[:,0])))
    for i,color in zip(range(len(speciesData)),colors):
        pylab.plot(times,speciesData[i],color=color,ls=':',label=p.speciesNames[i]+" data")
    times = scipy.linspace(0,10,100)
for speciesName,color in zip(p.speciesNames,colors): #range(numSpecies)
    pylab.plot(times,p.evaluateVec(times,speciesName,indepParamsList[-1]),color=color,label=speciesName+" fit") #'X_'+str(i)

pylab.legend(loc=4)
pylab.xlabel("Time")

die

# as a function of maxReplaceValue
maxReplaceValues = scipy.linspace(0,1,21)
costList = [ p.fitToDataDerivs(speciesData,derivsData,numiter=100,      \
             maxReplaceValue=maxReplaceValue,setModelParams=False,      \
             verbose=False)                                             \
            for maxReplaceValue in maxReplaceValues ]




