# simplePickle.py
#
# Bryan Daniels
# 1.24.2012
#

import io, sys
import pickle

# include some old name variants for back-compatability
#import fittingProblem
#sys.modules['FittingProblem'] = fittingProblem

def save(obj,filename):
    fout = io.open(filename,'wb')
    pickle.dump(obj,fout,2)
    fout.close()

def load(filename):
    fin = io.open(filename,'rb')
    obj = pickle.load(fin)
    fin.close()
    return obj
