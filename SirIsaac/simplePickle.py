# simplePickle.py
#
# Bryan Daniels
# 1.24.2012
#

import io, sys
import pickle

if sys.version_info.major > 2:
    # Python 3
    import pickle
else:
    # Python 2
    import cPickle as pickle

def save(obj,filename):
    fout = io.open(filename,'wb')
    # Note: we currently save using backward-compatible protocol 2
    pickle.dump(obj,fout,2)
    fout.close()

def load(filename):
    fin = io.open(filename,'rb')
    try:
        obj = pickle.load(fin)
    except UnicodeDecodeError:
        # try using backward-compatible encoding in case
        # the file was saved using python 2
        fin.close()
        fin = io.open(filename,'rb')
        obj = pickle.load(fin,encoding='bytes')
    fin.close()
    return obj
