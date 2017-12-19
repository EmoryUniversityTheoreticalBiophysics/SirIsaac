# SloppyCellTest.py
#
# Bryan Daniels
# 7/24/2013
#
# Tests whether SloppyCell is functional
#

import os,sys
print "This computer's name is",os.uname()[1]
if (os.uname()[1][:4] == 'node'): # 4.4.2012 emory machines
    print "The current directory is",os.getcwd()
    if os.getcwd().startswith('/star'):
        os.chdir('/star/physics/nemenman/daniels/SirIsaac')
    elif os.getcwd().startswith('/spark'):
        os.chdir('/spark/physics/nemenman/daniels/SirIsaac')
    print "Now the current directory is",os.getcwd()    
print ""

import FittingProblem
import scipy
from SloppyCell.ReactionNetworks import *

def testCcompiling(verbose=True):
    m = FittingProblem.CTSNFittingModel(5,['a'],['b'],
        verbose=verbose)
    n = m.net
    n.compile()
    
    if not n.compiled:
        if verbose:
          print "SloppyCellTest: testCcompiling: FAIL: "        \
                "m.net.compiled=False"
        return False
    ccode = n.get_c_code()
    if n._c_module_cache.has_key(ccode):
        if n._c_module_cache[ccode] is not None:
            return True
        else:
            if verbose:
              print "SloppyCellTest: testCcompiling: FAIL: "    \
                    "No compiled c code found. (2)"
            return False
    else:
        if verbose:
          print "SloppyCellTest: testCcompiling: FAIL: "        \
                "No compiled c code found. (1)"
        return False
        

def testCTSN(verbose=True):
    m = FittingProblem.CTSNFittingModel(5,['a'],['b'],
        verbose=verbose)
    r = m.evaluateVec(scipy.linspace(0,1,10),'b',[0])
    eps = 1e-5
    return sum(r**2) < eps



if __name__ == '__main__':

    print ""
    print "FittingProblem.py (including SloppyCell) imported."
    print ""

    print "Testing evaluation of CTSNFittingModel..."
    if testCTSN():
        print "CTSNFittingModel successfully evaluated."
    else:
        print "Error in evaluating CTSNFittingModel."
    print ""

    print "Testing C compiling..."
    if testCcompiling():
        print "SloppyCell C compiling successful."
    else:
        print "Error in SloppyCell C compiling."

    if False:
        print ""
        print "Testing evaluation of "                          \
              "examplePhosphorylationModel..."
        originalModelFilename =                                 \
            'examplePhosphorylationFittingModel.model'
        m = Utility.load(originalModelFilename)
        r = m.evaluateVec(scipy.linspace(0,1,10),'totalPhos',   \
                          [0.001,1.0])
        print "Sum(r)   =",sum(r)
        print "Expected = 10.1951142162"
        print ""

        # 9.24.2013
        print ""
        print "Testing creation of perfect data..."
        net = m.net
        params = net.GetParameters()
        timeInterval = [0,10]
        vars = ['totalPhos']
        randomX = True
        numPoints = 1
        scipy.random.seed(123)
        data = PerfectData.discrete_data(net,params,numPoints,  \
                    timeInterval,vars=vars,random=randomX)
        print "Output   =",data.values()[0].values()[0][0]
        print "Expected = 1.98956858449"

        

