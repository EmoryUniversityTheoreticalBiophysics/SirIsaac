# GaussianPrior.py
#
# Bryan Daniels
# 7.02.2009
#
# Residual class to be used with SloppyCell; implements a simple
# Gaussian prior.

import SloppyCell.Residuals
import numpy as np

class GaussianPrior(SloppyCell.Residuals.Residual):
    """
    Like Residual.PriorInLog, but without the log.
    
    The associated prior is
        exp[ -(1/2) ( (paramValue - bestPVal)/sigmaPVal )^2 ]
    """
    def __init__(self, key, pKey, bestPVal, sigmaPVal):
        SloppyCell.Residuals.Residual.__init__(self, key)
        self.pKey = pKey
        self.bestPVal = bestPVal
        self.sigmaPVal = sigmaPVal

    def GetValue(self, predictions, internalVars, params):
        return ( params.getByKey(self.pKey) - self.bestPVal )/self.sigmaPVal

    def dp(self, predictions, internalVars, params):
        return {self.pKey: 1./self.sigmaPVal}

    def dy(self, predictions, internalVars, params):
        return {}
    
    def dintVars(self, predictions, internalVars, params):  # was dintVar
        return {}
        
      
# 11.23.2011
class GaussianPriorExp(SloppyCell.Residuals.Residual):
    """
    The associated prior is
        exp[ -(1/2) ( (exp(paramValue) - exp(bestPVal))/sigmaPVal )^2 ]
    """
    def __init__(self, key, pKey, bestPVal, sigmaPVal):
        SloppyCell.Residuals.Residual.__init__(self, key)
        self.pKey = pKey
        self.bestPVal = bestPVal
        self.sigmaPVal = sigmaPVal

    def GetValue(self, predictions, internalVars, params):
        return ( np.exp( params.getByKey(self.pKey) )                    \
            - np.exp(self.bestPVal) )/self.sigmaPVal

    def dp(self, predictions, internalVars, params):
        return {self.pKey:                                                  \
            np.exp( params.getByKey(self.pKey) )/self.sigmaPVal}

    def dy(self, predictions, internalVars, params):
        return {}
    
    def dintVars(self, predictions, internalVars, params):  # was dintVar
        return {}
        
        
        
        
        
        
