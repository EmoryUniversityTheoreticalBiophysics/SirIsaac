#!/usr/bin/env python
# coding: utf-8

# *Note: This file is provided in two formats: 
# Python (simpleExample.py) and a Jupyter 
# notebook (simpleExample.ipynb).  The 
# Jupyter notebook opens in a web browser and 
# includes plots in an interactive format.  To 
# open the .ipynb file, run:*
#     
#     jupyter notebook simpleExample.ipynb
# 
# *To run the .py file in iPython at the command line, run:*
# 
#     ipython
#     %run simpleExample.py
#     plt.show()

# simpleExample.ipynb
# -------------------
# 
# - Bryan Daniels
# 
# - Last updated 2022/6/13
# 
# - Uses a simple exponential decay example to demonstrate usage of SirIsaac

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
from SirIsaac import fittingProblem


# Load example data
# -----------------

# In the example data file, we have five columns, each with 100 data points, listing:
# 
# * Initial condition *x_init*
# * Input condition *tau*
# * Measurement time *t*
# * Measurement value *x*
# * Measurement uncertainty (standard deviation)

# In[2]:


data = np.loadtxt('simpleExample_data.csv',delimiter=',')


# We now put this in a format compatible with SirIsaac.  First we make a list of input values (in this case initial conditions):
# 
# 

# In[3]:


indepParamsList = [ [ expt[0],expt[1] ] for expt in data ]


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
    sirIsaacData.append( { 'x': { expt[2]: ( expt[3], expt[4] ) } } )


# In[7]:


sirIsaacData[:3]


# Finally, SirIsaac will need to know what to call the input and output values.  In this case, the input corresponds to the initial value of *x*.  The way to indicate this to SirIsaac is by using the name 'x_init', where 'x' is the name of the corresponding variable.
# 
# Here we have two inputs and one output:

# In[8]:


outputNames = ['x']
indepParamNames = ['x_init','tau']


# Create SirIsaac FittingProblem
# ------------------------------

# We'll attempt to fit a model in the "continuous time sigmoidal network" (CTSN) class.  To do this, we'll create an instance of a CTSNFittingProblem.  Here we set up its arguments and create it:

# In[9]:


# complexityList lists which models in the model class may be tested.
# (Note that by default SirIsaac will still stop once 3 models have 
#  smaller estimated log-likelihood.)
complexityStepsize = 2 # increase complexity with steps of size 2
complexityMax = 25 # don't try models with complexity > 25
complexityList = range(0,complexityMax,complexityStepsize) 

# ensGen controls the generation of the initial ensemble of 
# parameter starting points.
# (We use a small number for totalSteps here so that the example
# runs faster; a more typical value for totalSteps is 1e3)
totalSteps = 20 #1e3 
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

# If you have mpi4py installed, you can run on multiple processors
numprocs = 4 #10

# We'll only use a subset of our data to make the example run faster
N = 20

p = fittingProblem.CTSNFittingProblem( complexityList, 
    sirIsaacData[:N], indepParamsList=indepParamsList[:N], 
    outputNames=outputNames, indepParamNames=indepParamNames, 
    ensGen=ensGen, avegtol=avegtol, maxiter=maxiter,
    priorSigma=priorSigma, numprocs=numprocs, verbose=True )


# Run parameter fitting
# ---------------------

# The bulk of computation time is used to fit the parameters of each model to the data.  Uncomment the following lines to run the parameter fitting, which takes an hour or two using 4 processors.  Or skip ahead to load a version that has already been fit.

# In[16]:


## Uncomment to run parameter fitting.
#p.fitAll()
#
#fittingProblem.save(p,'simpleExample_savedFittingProblem.pkl')


# In[10]:


# Load saved version of fittingProblem that has already been fit.
p = fittingProblem.load('simpleExample_savedFittingProblem.pkl')


# Analyze the selected model
# --------------------------

# Here we plot predicted timecourses from the selected model for the first 10 in-sample initial conditions, using plotBestModelResults:

# In[11]:


plt.figure(figsize=(20,2))
p.plotBestModelResults(plotInitialConditions=True,indices=range(10));


# And now for out-of-sample data:

# In[12]:


plt.figure(figsize=(20,2))
m = p.getBestModel()
m.plotResults(sirIsaacData[20:30],indepParamsList[20:30],
              plotInitialConditions=True,plotFittingData=True);


# We can look at the selected model's parameters:

# In[13]:


m = p.getBestModel()
print(m.getParameters())


# The following will use SloppyCell to output a latex file with the ODEs describing the selected model:

# In[14]:


m = p.getBestModel()
fittingProblem.IO.eqns_TeX_file(m.net,filename='simpleExample_selectedModel.tex')


# More details
# ------------

# We can examine the dynamics of the hidden nodes as well using plotResults.

# In[15]:


plt.figure(figsize=(20,4))
m = p.getBestModel()
m.plotResults(p.fittingData[:10],
              p.indepParamsList[:10],
              plotInitialConditions=True,
              plotHiddenNodes=True);


# We have access to raw trajectories using evaluateVec.  Here we use this to plot a projection of trajectories in phase space for the first in-sample initial conditions:

# In[17]:


plt.figure(figsize=(4,4))
times = np.linspace(0,1,1000)
xdata = m.evaluateVec(times,'x',p.indepParamsList[0])
X2data = m.evaluateVec(times,'X_2',p.indepParamsList[0])
fittingProblem.Plotting.plot(xdata,X2data)
plt.xlabel('x')
plt.ylabel('X_2');


# We can also look at other models that SirIsaac fit in searching for the best one.  In this case, 'Model 4' was selected because it has the largest estimated log-likelihood:

# In[18]:


for name in p.fittingModelNames:
    if name in p.logLikelihoodDict.keys():
        print('{}: #species = {}, #params = {}, L = {}'.format(name,
              len(p.fittingModelDict[name].speciesNames),
              p.numParametersDict[name],
              p.logLikelihoodDict[name]))
print()
print('Selected model:',p.maxLogLikelihoodName())


# We can also look at predictions of a model with more parameters:

# In[19]:


plt.figure(figsize=(20,2))
m2 = p.fittingModelDict['Model 7']
m2.plotResults(sirIsaacData[:10],indepParamsList[:10],
              plotInitialConditions=True,plotFittingData=True);


# Also potentially useful is the Hessian at the best-fit parameters:

# In[20]:


hess = p.HessianDict[p.maxLogLikelihoodName()]
u,singVals,vt = np.linalg.svd( hess )
np.sort(singVals)


# Other details about what happened during parameter fitting are stored within each fittingModel:

# In[21]:


m = p.getBestModel()
print("Acceptance ratio for initial parameter ensemble = {}".format(m.acceptanceRatio))
c = sum(np.array(m.currentResiduals(p.fittingData,p.indepParamsList,includePriors=False))**2)
print("Sum of squared residuals at best-fit (without priors) = {}".format(c))
print("Convergence flags for local fits: {}".format(m.convFlagList))
print("Number of cost evaluations for local fits: {}".format(m.numCostCallsList))
print("Number of gradient evaluations for local fits: {}".format(m.numGradCallsList))


# Finally, since in this case we know the function used to create the data, we can compare:

# In[22]:


plt.figure(figsize=(20,2))
indicesToPlot = range(10)
axArray = p.plotBestModelResults(plotInitialConditions=True,indices=indicesToPlot)

# compare to model that generated the data (see simpleExample_makeData.py)
f = lambda x0,tau,t: x0 * np.exp(-times/tau)

for i,indepParams in enumerate(np.array(indepParamsList)[indicesToPlot]):
    times = np.linspace(0,1,100)
    x0,tau = indepParams
    fittingProblem.Plotting.sca(axArray[0][i])
    fittingProblem.Plotting.plot(times,f(x0,tau,times),'k:')


# In[ ]:




