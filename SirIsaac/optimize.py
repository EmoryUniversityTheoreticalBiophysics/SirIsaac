# optimize.py
#
# Bryan Daniels
# 5.8.2014
#
# Branching modified versions from
# SloppyCell.Optimization.
#

import numpy as np
import SloppyCell.lmopt as lmopt

def fmin_lm_log_params(m, params, *args, **kwargs):
    """
    Minimize the cost of a model using Levenberg-Marquardt 
    in terms of log parameters.
    """
    jac = lambda lp: np.asarray(m.jacobian_log_params_sens(lp))
    sln = lmopt.fmin_lm(f=m.res_log_params, x0=np.log(params),
        fprime=jac, *args, **kwargs)
    return ( np.exp(sln[0]), ) + sln[1:]

def fmin_lm(m, params, *args, **kwargs):
    """
    Minimize the cost of a model using Levenberg-Marquardt.
    """
    jac = lambda p: np.asarray(m.jacobian_sens(p))
    sln = lmopt.fmin_lm(f=m.res, x0=params, fprime=jac,
                        *args, **kwargs)
    return sln
