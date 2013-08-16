# SloppyCellTest.py
#
# Bryan Daniels
# 7/24/2013
#
# Tests whether SloppyCell is functional
#

import os
print "This computer's name is",os.uname()[1]
if (os.uname()[1][:4] == 'node'): # 4.4.2012 emory machines
    print "The current directory is",os.getcwd()
    if os.getcwd().startswith('/star'):
        os.chdir('/star/physics/nemenman/daniels/SloppySimplification')
    elif os.getcwd().startswith('/spark'):
        os.chdir('/spark/physics/nemenman/daniels/SloppySimplification')
    print "Now the current directory is",os.getcwd()    
print ""

import FittingProblem
import scipy
from SloppyCell.ReactionNetworks import *

print ""
print "FittingProblem.py (including SloppyCell) imported."
print ""

print "Testing evaluation of CTSNFittingModel..."
m = FittingProblem.CTSNFittingModel(5,['a'],['b'])
r = m.evaluateVec(scipy.linspace(0,1,10),'b',[0])
if sum(r**2)==0.:
    print "CTSNFittingModel successfully evaluated."
else:
    print "Error in evaluating CTSNFittingModel."

print ""
print "Testing evaluation of model from restartPhosDict..."
restartPhosDictFilename = 'restartPhosDict_switchSigmoid.data'
restartPhosDict = Utility.load(restartPhosDictFilename)
fitProbDict = restartPhosDict.values()[0]
fp = fitProbDict.values()[0]
m = fp.fittingModelList[0]
r = m.evaluateVec(scipy.linspace(0,1,10),'totalPhos',\
    fp.indepParamsList[0])
print "Sum(r)   =",sum(r)
print "Expected = 8.1894469483"
print ""

