# linalgTools.py
#
# Bryan Daniels
# 9.7.2012
#
# 
#

import pylab
import scipy
import scipy.linalg

# 7.10.2012
def svdInverse(mat,maxEig=1e10,minEig=1e-10): #1e10,1e-10
    u,w,vt = scipy.linalg.svd(mat)
    if any(w==0.):
        raise ZeroDivisionError("Singular matrix.")
    wInv = w ** -1
    largeIndices = pylab.find( abs(wInv) > maxEig )
    if len(largeIndices) > 0: print("svdInverse:",len(largeIndices),"large singular values out of",len(w))
    wInv[largeIndices] = maxEig*scipy.sign(wInv[largeIndices])
    smallIndices = pylab.find( abs(wInv) < minEig )
    if len(smallIndices) > 0: print("svdInverse:",len(smallIndices),"small singular values out of",len(w))
    wInv[smallIndices] = minEig*scipy.sign(wInv[smallIndices])
    return scipy.dot( scipy.dot(vt.T,scipy.diag(wInv)), u.T )