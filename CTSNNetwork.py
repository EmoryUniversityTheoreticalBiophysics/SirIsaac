# CTSNNetwork.py
#
# Bryan Daniels
# 7.29.2009
#
# A SloppyCell implementation of CTSNs (continuous-time sigmoidal networks).
# 
# (modeled after PowerLawNetwork.py)


from SloppyCell.ReactionNetworks import *
from scipy import random, log
import copy
import GaussianPrior

def CTSN_List(networkList,speciesNames=None,                                    \
    logParams=True,netid='CTSN',switchSigmoid=False,xiNegative=False):
    """
    Defines a CTSN based on a connection list.
    
    A SloppyCell implementation of CTSNs:
    
    d X_i / d t  =  
      1/tau_i * ( -X_i + sum_j=1^n w_i_j xi(y_j+theta_j) )
    
    tau is given the default value of 1, and xinit, theta, and w
    parameters are by default 0.
    
    Right now, inputs come into the sum as w_i_j*input_j.
    (Not sure if this is "correct"...)
    
    networkList  : list of the form
  
    [ [nodeType, { connectFrom: connectType, connectFrom: connectType, ...}], ... ]
        nodeType    : integer between 0 and 4 (the number of optimizable 
                      parameters specifying the node's behavior; 0 for input node)
        connectFrom : integer index of node to get connection from
        connectType : integer, either 1 or 2 (the number of parameters specifying
                      the connection)
                      
    speciesNames    : list of species names (length n+m).  If None,
                      species are named X_i for i in range(n+m).
    logParams       : if True, the time constants tau are written 
                      as log_tau (to facilitate parameter searches)
    switchSigmoid   : If True, use sigmoid(sum) instead of sum(sigmoid) in
                      each node's ODE rate rule.  See notes 7.8.2013.
    
    (The X_js beyond n (up to n+m) are constant inputs.)
    """
    
    n = len(networkList)
    #m = 0
    
    # the order in which to add parameters
    order = dict( zip(['xinit','wself','tau','theta'], range(5)) )
    orderConnect = dict( zip(['w'], range(1)) )
    
    net = Network(netid, name='CTSN')
    net.addCompartment('Comp',name='Compartment')
    
    net.addParameter('n', n, isOptimizable=False)
    #net.addParameter('m', m, isOptimizable=False)
    
    defaultParam = 0.
    defaultLogParam = 1.
    defaultW = 0.
    #defaultExpParam = 0.
    
    if speciesNames is None:
        speciesNames = [ 'X_'+str(i) for i in range(n) ]
    
    # add parameters
    for i in range(n):
      
      nodeType, connectionDict = networkList[i]
      
      if nodeType != 0: # if it's not an input node
      
        notLog = not logParams
        
        net.addParameter('wself_'+str(i), defaultW,                             \
            isOptimizable=order['wself']<nodeType)
        net.addParameter('theta_'+str(i), defaultParam,                         \
            isOptimizable=order['theta']<nodeType)
        net.addParameter('tau_'+str(i), defaultLogParam,                        \
            isOptimizable=(notLog and order['tau']<nodeType))
        
        if logParams:
            net.addParameter('log_tau_'+str(i), log(defaultLogParam),           \
                isOptimizable=order['tau']<nodeType,typicalValue=1.)
            net.addAssignmentRule('tau_'+str(i),'exp(log_tau_'+str(i)+')')
        
        # connect to others
        for j in connectionDict.keys():
            net.addParameter('w_'+str(i)+'_'+str(j), defaultW,                  \
                isOptimizable=orderConnect['w']<connectionDict[j])
        
        net.addParameter(speciesNames[i]+'_init', defaultParam,                 \
            isOptimizable=order['xinit']<nodeType)
            
    # add species
    for i in range(n):
      nodeType, connectionDict = networkList[i]
      if nodeType != 0: # if it's not an input node
        net.addSpecies( speciesNames[i], 'Comp', speciesNames[i]+'_init' )
      else: # it is an input node
        # add as a parameter if it's not already there
        if speciesNames[i] not in net.parameters.keys():
          net.addParameter( speciesNames[i], 0., isOptimizable=False )
    
    # reaction rate rules
    for i in range(n):
      nodeType, connectionDict = networkList[i]
      if (nodeType != 0) and not switchSigmoid: # default
        sum = ''
        # always connect to yourself
        if xiNegative:
          sum += 'wself_'+str(i)                                                  \
              +' / (1. + exp(-'+speciesNames[i]+' - theta_'+str(i)+')) + '
        else: # prior to 12.19.2013
            sum += 'wself_'+str(i)                                                  \
                +' / (1. + exp('+speciesNames[i]+' + theta_'+str(i)+')) + '
        for j in connectionDict.keys():
          if networkList[j][0] != 0: # the connection is not from an input node
            if xiNegative: 
                sum += 'w_'+str(i)+'_'+str(j)                                       \
                  +' / (1. + exp(-'+speciesNames[j]+' - theta_'+str(j)+')) + '
            else: # prior to 12.19.2013
                sum += 'w_'+str(i)+'_'+str(j)                                       \
                    +' / (1. + exp('+speciesNames[j]+' + theta_'+str(j)+')) + '
          else:  # it is an input node.  XXX How should I do this?
            sum += 'w_'+str(i)+'_'+str(j)+' * '+speciesNames[j]+' + '
        # 3.30.2012 trying having tau only divide the decay term
        #net.addRateRule( speciesNames[i],                                       \
        #    '1./tau_'+str(i)+'*( -'+speciesNames[i]+' + '+sum[:-3]+' )')
        net.addRateRule( speciesNames[i],                                       \
            '1./tau_'+str(i)+'*( -'+speciesNames[i]+') + '+sum[:-3] )
      elif (nodeType !=0) and switchSigmoid: # 7.8.2013
        # new version proposed by Ilya
          sum = ''
          # always connect to yourself
          sum += 'wself_'+str(i)+'*('+speciesNames[i]+' + theta_'+str(i)+') + '
          for j in connectionDict.keys():
              if networkList[j][0] != 0: # the connection is not from an input node
                  sum += 'w_'+str(i)+'_'+str(j)                                     \
                      +'*('+speciesNames[j]+'+ theta_'+str(j)+') + '
              else:  # it is an input node.  XXX How should I do this?
                  sum += 'w_'+str(i)+'_'+str(j)+' * '+speciesNames[j]+' + '
          sigmoidSum = '1. / (1. + exp('+sum[:-3]+'))'
          net.addRateRule( speciesNames[i],                                       \
                          '1./tau_'+str(i)+'*( -'+speciesNames[i]+') + '+sigmoidSum )
      else: # it's an input node
        pass

    return net

    
    
def setRandomParameters(net,seed=None,randFunc=random.random):
    """
    Sets parameters to random values given by the function randFunc (by
    default, uniformly distributed on [0,1) ).
    """
    random.seed(seed)
    net.setOptimizables( randFunc(len(net.GetParameters())) )
    return net.GetParameters()



