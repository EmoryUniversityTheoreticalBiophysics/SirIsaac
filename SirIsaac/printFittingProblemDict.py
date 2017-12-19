# printFittingProblemDict.py
#
# Bryan Daniels
# 8.5.2012
#


from FittingProblem import *
from analyzeFittingProblemDict import orderedFitNames
import sys

def die():
    print "Use: python printFittingProblemDict.py fittingProblemDictFilename.dat"
    sys.exit()

if len(sys.argv) != 2:
    die()
fpdFilename = sys.argv[1]

fpd = load(fpdFilename)

print ""
print fpdFilename
print ""

for k in scipy.sort(fpd.keys()):
    
    print ""
    print "-------------- Key",k,"--------------"
    fp = fpd[k]
    names = orderedFitNames(fp) 
    
    print ""
    print "Log-likelihoods"
    for name in names:
        print name,fp.logLikelihoodDict[name]
    print ""
    
    print ""
    print "Number of parameters dictionary"
    for name in names:
        print name,fp.numParametersDict[name]
    print ""
    
    print ""
    print "Effective number of parameters dictionary"
    for name in names:
        print name,fp.numStiffSingValsDict[name]
    print ""