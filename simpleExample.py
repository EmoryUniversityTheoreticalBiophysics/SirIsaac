# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

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

# <markdowncell>

# simpleExample.ipynb
# -------------------
# 
# - Bryan Daniels
# 
# - 1.29.2014
# 
# - Uses a simple 1D harmonic oscillator example to demonstrate usage of SirIsaac.

# <codecell>

import scipy, pylab
from SirIsaac.fittingProblem import *

# <markdowncell>

# Load example data
# -----------------

# <markdowncell>

# In the example data file, we have four columns, each with 100 data points, listing:
# 
# * Initial condition *x_init*
# * Measurement time *t*
# * Measurement value *x*
# * Measurement uncertainty (standard deviation)

# <codecell>

data = scipy.loadtxt('simpleExample_data.txt')

# <markdowncell>

# We now put this in a format compatible with SirIsaac.  First we make a list of input values (in this case initial conditions):
# 

# <codecell>

indepParamsList = [ [ expt[0] ] for expt in data ]

# <codecell>

indepParamsList[:3]

# <markdowncell>

# Next, we have a corresponding list of data taken at each of those input values, in the format below.  In this case, we only have one variable *x*.  (Note: In general, multiple timepoints could be also be measured at each input value; in all of our examples, we measure all variables at a single timepoint per input value.)

# <codecell>

# [ {'var1': { time0: ( value, uncertainty ) },
#    'var2': { time0: ( value, uncertainty ) },
#     ... },
#   {'var1': { time1: ( value, uncertainty ) },
#    'var2': { time1: ( value, uncertainty ) },
#     ... },
#   ... ]

# <codecell>

sirIsaacData = []
for expt in data:
    sirIsaacData.append( { 'x': { expt[1]: ( expt[2], expt[3] ) } } )

# <codecell>

sirIsaacData[:3]

# <markdowncell>

# Finally, SirIsaac will need to know what to call the input and output values.  In this case, the input corresponds to the initial value of *x*.  The way to indicate this to SirIsaac is by using the name 'x_init', where 'x' is the name of the corresponding variable.
# 
# Here we have one input and one output:

# <codecell>

outputNames = ['x']
indepParamNames = ['x_init']

# <markdowncell>

# Create SirIsaac FittingProblem
# ------------------------------

# <markdowncell>

# We'll attempt to fit a model in the power law class.  To do this, we'll create an instance of a PowerLawFittingProblem.  Here we set up its arguments and create it:

# <codecell>

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
ensGen = EnsembleGenerator( totalSteps, keepSteps,
    temperature=ensTemperature, seeds=seeds )

# Parameters that control when local fitting stops.
avegtol = 1e-2
maxiter = 100

# priorSigma controls the width of priors on all parameters
priorSigma = 3.

# If you have pypar installed, you can run on multiple processors
numprocs = 1

# We'll only use a subset of our data to make the example run faster
N = 20

p = PowerLawFittingProblem( complexityList, 
    sirIsaacData[:N], indepParamsList=indepParamsList[:N], 
    outputNames=outputNames, indepParamNames=indepParamNames, 
    ensGen=ensGen, avegtol=avegtol, maxiter=maxiter,
    priorSigma=priorSigma, numprocs=numprocs, verbose=True )

# <markdowncell>

# Run parameter fitting
# ---------------------

# <markdowncell>

# The bulk of computation time is used to fit the parameters of each model to the data.  Uncomment the following line to run the parameter fitting, which takes a few hours using 10 processors.  Or skip ahead to load a version that has already been fit.

# <codecell>

# p.fitAll()
#
# save(p, 'simpleExample_savedFittingProblem.data')

# <codecell>

p = load('simpleExample_savedFittingProblem.data')

# <markdowncell>

# Analyze the selected model
# --------------------------

# <markdowncell>

# Here we plot predicted timecourses from the selected model for the first 10 in-sample initial conditions, using plotBestModelResults:

# <codecell>

pylab.figure(figsize=(20,2))
p.plotBestModelResults(plotInitialConditions=True,indices=range(10));
pylab.show()

# <markdowncell>

# And now for out-of-sample data:

# <codecell>

pylab.figure(figsize=(20,2))
m = p.getBestModel()
m.plotResults(sirIsaacData[20:30],indepParamsList[20:30],
              plotInitialConditions=True,plotFittingData=True);
pylab.show()

# <markdowncell>

# We can look at the selected model's parameters:

# <codecell>

m = p.getBestModel()
print m.getParameters()

# <markdowncell>

# The following will use SloppyCell to output a latex file with the ODEs describing the selected model:

# <codecell>

m = p.getBestModel()
IO.eqns_TeX_file(m.net,filename='simpleExample_selectedModel.tex')

# <markdowncell>

# More details
# ------------

# <markdowncell>

# We can examine the dynamics of the hidden nodes as well using plotResults.

# <codecell>

pylab.figure(figsize=(20,6))
m = p.getBestModel()
m.plotResults(p.fittingData[:10],p.indepParamsList[:10],
              plotInitialConditions=True,plotHiddenNodes=True);
pylab.show()

# <markdowncell>

# We have access to raw trajectories using evaluateVec.  Here we use this to plot a projection of trajectories in phase space for the first in-sample initial conditions:

# <codecell>

pylab.figure(figsize=(4,4))
times = scipy.linspace(0,1,1000)
xdata = m.evaluateVec(times,'x',p.indepParamsList[0])
X1data = m.evaluateVec(times,'X_1',p.indepParamsList[0])
Plotting.plot(xdata,X1data)
pylab.xlabel('x')
pylab.ylabel('X_1')
pylab.show()

# <markdowncell>

# We can also look at other models that SirIsaac fit in searching for the best one.  In this case, 'Model 7' was selected because it has the largest estimated log-likelihood:

# <codecell>

for name in p.fittingModelNames:
  if name in p.logLikelihoodDict.keys():
    print name, ': #species =',len(p.fittingModelDict[name].speciesNames),\
                ', #params =',p.numParametersDict[name],\
                ', L =', p.logLikelihoodDict[name]
print
print 'Selected model:',p.maxLogLikelihoodName()

# <markdowncell>

# A model with more parameters fits in-sample data better but out-of-sample data worse:

# <codecell>

pylab.figure(figsize=(20,2))
m2 = p.fittingModelDict['Model 10']
m2.plotResults(sirIsaacData[:10],indepParamsList[:10],
              plotInitialConditions=True,plotFittingData=True);
pylab.show()

# <codecell>

pylab.figure(figsize=(20,2))
m2.plotResults(sirIsaacData[30:40],indepParamsList[30:40],
              plotInitialConditions=True,plotFittingData=True);
pylab.show()

# <markdowncell>

# Also potentially useful is the Hessian at the best-fit parameters:

# <codecell>

hess = p.HessianDict[p.maxLogLikelihoodName()]
u,singVals,vt = scipy.linalg.svd( hess )
scipy.sort(singVals)

# <markdowncell>

# Other details about what happened during parameter fitting are stored within each fittingModel:

# <codecell>

m = p.getBestModel()
print "Acceptance ratio for initial parameter ensemble =",m.acceptanceRatio
c = sum(scipy.array(m.currentResiduals(p.fittingData,p.indepParamsList,includePriors=False))**2)
print "Sum of squared residuals at best-fit (without priors) =",c
print "Convergence flags for local fits:",m.convFlagList
print "Number of cost evaluations for local fits:",m.numCostCallsList
print "Number of gradient evaluations for local fits:",m.numGradCallsList

# <markdowncell>

# Finally, since in this case we know the function used to create the data, we can compare:

# <codecell>

pylab.figure(figsize=(20,2))
indicesToPlot = range(5)
axArray = p.plotBestModelResults(plotInitialConditions=True,indices=indicesToPlot)

# compare to underlying known model
f = lambda x0,t: 1.5 + 0.5*scipy.sin(4.*scipy.pi*t + scipy.arcsin(2.*x0 - 3.))
for i,indepParams in enumerate(scipy.array(indepParamsList)[indicesToPlot]):
    times = scipy.linspace(0,1,100)
    x0 = indepParams[0]
    Plotting.sca(axArray[0][i])
    Plotting.plot(times,f(x0,times),'k:')
pylab.show()
