# simplePickle.py
#
# Bryan Daniels
# 1.24.2012
#

import io, sys
import cPickle

# include some old name variants for back-compatability
import fittingProblem
sys.modules['FittingProblem'] = fittingProblem

def save(obj,filename):
    fout = io.open(filename,'wb')
    cPickle.dump(obj,fout,2)
    fout.close()

def load(filename):
    fin = io.open(filename,'rb')
    obj = cPickle.load(fin)
    fin.close()
    return obj
