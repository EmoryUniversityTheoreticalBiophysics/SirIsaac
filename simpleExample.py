
# coding: utf-8

# *Note: This file is provided in two formats: 
# Python (simpleExample.py) and a Jupyter iPython 
# notebook (simpleExample.ipynb).  The 
# iPython notebook opens in a web browser and 
# includes plots in an interactive format.  To 
# open the .ipynb file, run:*
#     
#     jupyter notebook simpleExample.ipynb
# 
# *To run the .py file in iPython at the command line, run:*
# 
#     ipython --pylab
#     %run simpleExample.py
#     show()

# simpleExample.ipynb
# -------------------
# 
# - Bryan Daniels
# 
# - 1.29.2014
# - updated 1.18.2019
# 
# - Uses a simple 1D harmonic oscillator example to demonstrate usage of SirIsaac.

# In[1]:

import scipy, pylab
from SirIsaac import fittingProblem


# Load example data
# -----------------

# In the example data file, we have four columns, each with 100 data points, listing:
# 
# * Initial condition *x_init*
# * Measurement time *t*
# * Measurement value *x*
# * Measurement uncertainty (standard deviation)

# In[2]:

data = scipy.loadtxt('simpleExample_data.txt')


# We now put this in a format compatible with SirIsaac.  First we make a list of input values (in this case initial conditions):
# 
# 

# In[3]:

indepParamsList = [ [ expt[0] ] for expt in data ]


# In[4]:

indepParamsList[:3]


# Next, we have a corresponding list of data taken at each of those input values, in the format below.  In this case, we only have one variable *x*.  (Note: In general, multiple timepoints could be also be measured at each input value; in all of our examples, we measure all variables at a single timepoint per input value.)

# In[5]:

# [ {'var1': { time0: ( value, uncertainty ) },
#    'var2': { time0: ( value, uncertainty ) },
#     ... },
#   {'var1': { time1: ( value, uncertainty ) },
#    'var2': { time1: ( value, uncertainty ) },
#     ... },
#   ... ]


# In[6]:

sirIsaacData = []
for expt in data:
    sirIsaacData.append( { 'x': { expt[1]: ( expt[2], expt[3] ) } } )


# In[7]:

sirIsaacData[:3]


# Finally, SirIsaac will need to know what to call the input and output values.  In this case, the input corresponds to the initial value of *x*.  The way to indicate this to SirIsaac is by using the name 'x_init', where 'x' is the name of the corresponding variable.
# 
# Here we have one input and one output:

# In[8]:

outputNames = ['x']
indepParamNames = ['x_init']


# Create SirIsaac FittingProblem
# ------------------------------

# We'll attempt to fit a model in the power law class.  To do this, we'll create an instance of a PowerLawFittingProblem.  Here we set up its arguments and create it:

# In[9]:

# complexityList lists which models in the model class may be tested.
# (Note that by default SirIsaac will still stop once 3 models have 
#  smaller estimated log-likelihood.)
complexityStepsize = 2 # increase complexity with steps of size 2
complexityMax = 25 # don't try models with complexity > 25
complexityList = range(0,complexityMax,complexityStepsize) 

# ensGen controls the generation of the initial ensemble of 
# parameter starting points.
totalSteps = 1e3
keepSteps = 10
seeds = (1,1) # use a fixed random seed
ensTemperature = 100.
ensGen = fittingProblem.EnsembleGenerator( totalSteps, keepSteps,
    temperature=ensTemperature, seeds=seeds )

# Parameters that control when local fitting stops.
avegtol = 1e-2
maxiter = 100

# priorSigma controls the width of priors on all parameters
priorSigma = 3.

# If you have pypar installed, you can run on multiple processors
numprocs = 10

# We'll only use a subset of our data to make the example run faster
N = 20

p = fittingProblem.PowerLawFittingProblem( complexityList, 
    sirIsaacData[:N], indepParamsList=indepParamsList[:N], 
    outputNames=outputNames, indepParamNames=indepParamNames, 
    ensGen=ensGen, avegtol=avegtol, maxiter=maxiter,
    priorSigma=priorSigma, numprocs=numprocs, verbose=True )


# Run parameter fitting
# ---------------------

# The bulk of computation time is used to fit the parameters of each model to the data.  Uncomment the following lines to run the parameter fitting, which takes a few hours using 10 processors.  Or skip ahead to load a version that has already been fit.

# In[11]:

## Uncomment to run parameter fitting.
#p.fitAll()
#
#fittingProblem.save(p,'simpleExample_savedFittingProblem.data')


# In[10]:

# Load saved version of fittingProblem that has already been fit.
p = fittingProblem.load('simpleExample_savedFittingProblem.data')


# Analyze the selected model
# --------------------------

# Here we plot predicted timecourses from the selected model for the first 10 in-sample initial conditions, using plotBestModelResults:

# In[11]:

pylab.figure(figsize=(20,2))
p.plotBestModelResults(plotInitialConditions=True,indices=range(10));


# And now for out-of-sample data:

# In[12]:

pylab.figure(figsize=(20,2))
m = p.getBestModel()
m.plotResults(sirIsaacData[20:30],indepParamsList[20:30],
              plotInitialConditions=True,plotFittingData=True);


# We can look at the selected model's parameters:

# In[13]:

m = p.getBestModel()
print m.getParameters()


# The following will use SloppyCell to output a latex file with the ODEs describing the selected model:

# In[14]:

m = p.getBestModel()
fittingProblem.IO.eqns_TeX_file(m.net,filename='simpleExample_selectedModel.tex')


# More details
# ------------

# We can examine the dynamics of the hidden nodes as well using plotResults.

# In[15]:

pylab.figure(figsize=(20,6))
m = p.getBestModel()
m.plotResults(p.fittingData[:10],p.indepParamsList[:10],
              plotInitialConditions=True,plotHiddenNodes=True);


# We have access to raw trajectories using evaluateVec.  Here we use this to plot a projection of trajectories in phase space for the first in-sample initial conditions:

# In[16]:

pylab.figure(figsize=(4,4))
times = scipy.linspace(0,1,1000)
xdata = m.evaluateVec(times,'x',p.indepParamsList[0])
X1data = m.evaluateVec(times,'X_1',p.indepParamsList[0])
fittingProblem.Plotting.plot(xdata,X1data)
pylab.xlabel('x')
pylab.ylabel('X_1')


# We can also look at other models that SirIsaac fit in searching for the best one.  In this case, 'Model 6' was selected because it has the largest estimated log-likelihood:

# In[17]:

for name in p.fittingModelNames:
  if name in p.logLikelihoodDict.keys():
    print name, ': #species =',len(p.fittingModelDict[name].speciesNames),                ', #params =',p.numParametersDict[name],                ', L =', p.logLikelihoodDict[name]
print
print 'Selected model:',p.maxLogLikelihoodName()


# A model with more parameters fits in-sample data better but out-of-sample data worse:

# In[18]:

pylab.figure(figsize=(20,2))
m2 = p.fittingModelDict['Model 9']
m2.plotResults(sirIsaacData[:10],indepParamsList[:10],
              plotInitialConditions=True,plotFittingData=True);


# In[19]:

pylab.figure(figsize=(20,2))
m2.plotResults(sirIsaacData[30:40],indepParamsList[30:40],
              plotInitialConditions=True,plotFittingData=True);


# Also potentially useful is the Hessian at the best-fit parameters:

# In[20]:

hess = p.HessianDict[p.maxLogLikelihoodName()]
u,singVals,vt = scipy.linalg.svd( hess )
scipy.sort(singVals)


# Other details about what happened during parameter fitting are stored within each fittingModel:

# In[21]:

m = p.getBestModel()
print "Acceptance ratio for initial parameter ensemble =",m.acceptanceRatio
c = sum(scipy.array(m.currentResiduals(p.fittingData,p.indepParamsList,includePriors=False))**2)
print "Sum of squared residuals at best-fit (without priors) =",c
print "Convergence flags for local fits:",m.convFlagList
print "Number of cost evaluations for local fits:",m.numCostCallsList
print "Number of gradient evaluations for local fits:",m.numGradCallsList


# Finally, since in this case we know the function used to create the data, we can compare:

# In[22]:

pylab.figure(figsize=(20,2))
indicesToPlot = range(5)
axArray = p.plotBestModelResults(plotInitialConditions=True,indices=indicesToPlot)

# compare to model that generated the data
f = lambda x0,t: 1.5 + 0.5*scipy.sin(4.*scipy.pi*t + scipy.arcsin(2.*x0 - 3.))
for i,indepParams in enumerate(scipy.array(indepParamsList)[indicesToPlot]):
    times = scipy.linspace(0,1,100)
    x0 = indepParams[0]
    fittingProblem.Plotting.sca(axArray[0][i])
    fittingProblem.Plotting.plot(times,f(x0,times),'k:')


# In[ ]:



