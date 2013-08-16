# simplePickle.py
#
# Bryan Daniels
# 1.24.2012
#

import io
import cPickle

def save(obj,filename):
    fout = io.open(filename,'wb')
    cPickle.dump(obj,fout,2)
    fout.close()

def load(filename):
    fin = io.open(filename,'rb')
    obj = cPickle.load(fin)
    fin.close()
    return obj