# PowerLaw_Network.py
#
# Bryan Daniels
# 06.18.2009
#
# A SloppyCell implementation of Savageau's power law representation.
#


from SloppyCell.ReactionNetworks import *
from scipy import random, log
import copy
import GaussianPrior
    
def PowerLaw_Network_List(networkList,speciesNames=None,                        \
    logParams=True,netid='PowerLaw_Network',logParamsInit=False,                \
    includeRegularizer=False,regStrength=0.,useDeltaGamma=False):
    """
    Defines a power-law network based on a connection list.
    
    A SloppyCell implementation of Savageau's power law representation:
    
    d X_i / d t  =  
      alpha_i prod_{j=0}^{n+m-1} X_j^g_i_j - beta_i prod_{j=0}^{n+m-1} X_j^h_i_j
    
    alpha and beta parameters are given the default value of 1, and g and h
    parameters are by default 0.  There are in general a maximum of
    2n(n+m) + 3n + m parameters, including initial conditions.
    
    networkList  : list of the form
  
    [ [nodeType, { connectFrom: connectType, connectFrom: connectType, ...}], ... ]
        nodeType    : integer between 0 and 5 (the number of optimizable 
                      parameters specifying the node's behavior; 0 for input node)
        connectFrom : integer index of node to get connection from
        connectType : integer, either 1 or 2 (the number of parameters specifying
                      the connection)
                      
    speciesNames    : list of species names (length n+m).  If None,
                      species are named X_i for i in range(n+m).
    logParams       : if True, the optimizable multiplicative constants
                      (alpha and beta) and initial values are written 
                      as log_alpha, log_beta, etc. (to facilitate 
                      parameter searches)
    includeRegularizer (True)   : See notes 1.30.2013.  Regularization forces
                                  solutions to behave nicely (not go to zero
                                  or infinity).
    regStrength (0.)            : See notes 1.30.2013.  Sets the strength
                                  of regularization.
    useDeltaGamma (False)       : 2.13.2013 If True, use different parameterization:
                                  delta = alpha, gamma = alpha * beta
                                  (This was usually True prior to 2.13.2013.)
    
    (The X_js beyond n (up to n+m) are constant inputs.)
    """
    
    n = len(networkList)
    #m = 0
    
    # the order in which to add parameters
    if useDeltaGamma:
        order = dict( zip(['xinit','g','gamma','h','delta'], range(5)) )
    else:
        order = dict( zip(['xinit','g','beta','h','alpha'], range(5)) )
    orderConnect = dict( zip(['g','h'], range(2)) )
    
    #fullConnectionList = copy.deepcopy(connectionList)
    # add "self-connections"
    #for i in range(n):
    #  if (i==0) or (i>numIndepParams): # if not an input node
    #    fullConnectionList[i] += [i]
    
    net = Network(netid, name='Power Law Network')
    net.addCompartment('Comp',name='Compartment')
    
    net.addParameter('n', n, isOptimizable=False)
    #net.addParameter('m', m, isOptimizable=False)
    
    if includeRegularizer:
        net.addParameter('regStrength',regStrength,isOptimizable=False)
    
    defaultParam = 1.
    defaultExpParam = 0.
    
    if speciesNames is None:
        speciesNames = [ 'X_'+str(i) for i in range(n) ]
    
    # add parameters
    for i in range(n):
      
      nodeType, connectionDict = networkList[i]
      
      if nodeType != 0: # if it's not an input node
      
        notLog = not logParams
        
        if useDeltaGamma:
            multNames = ['delta','gamma']
        else:
            multNames = ['alpha','beta']
        
        for multName in multNames:
            net.addParameter(multName+'_'+str(i), defaultParam,                 \
                isOptimizable=(notLog and order[multName]<nodeType))

        if logParams:
          for multName in multNames:
            net.addParameter('log_'+multName+'_'+str(i), log(defaultParam),     \
                isOptimizable=order[multName]<nodeType)
            net.addAssignmentRule(multName+'_'+str(i),                          \
                'exp(log_'+multName+'_'+str(i)+')')
        
        # always connect to yourself
        net.addParameter('g_'+str(i)+'_'+str(i), defaultExpParam,               \
            isOptimizable=order['g']<nodeType)
        net.addParameter('h_'+str(i)+'_'+str(i), defaultExpParam,               \
            isOptimizable=order['h']<nodeType)
        
        # connect to others
        for j in connectionDict.keys():
            net.addParameter('g_'+str(i)+'_'+str(j), defaultExpParam,           \
                isOptimizable=orderConnect['g']<connectionDict[j])
            net.addParameter('h_'+str(i)+'_'+str(j), defaultExpParam,           \
                isOptimizable=orderConnect['h']<connectionDict[j])
        
        net.addParameter(speciesNames[i]+'_init', defaultParam,                 \
            isOptimizable=(notLog and order['xinit']<nodeType))
        if logParamsInit:
            net.addParameter('log_'+speciesNames[i]+'_init', log(defaultParam), \
                isOptimizable=order['xinit']<nodeType)
            net.addAssignmentRule(speciesNames[i]+'_init',\
                'exp(log_'+speciesNames[i]+'_init)')
    
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
      if nodeType != 0: 
        product1,product2 = '',''
        # always connect to yourself
        product1 += speciesNames[i]+'**g_'+str(i)+'_'+str(i)+'*'
        product2 += speciesNames[i]+'**h_'+str(i)+'_'+str(i)+'*'
        if includeRegularizer:
            # 1.30.2013
            product1 += 'exp(regStrength*1./'+speciesNames[i]+')*'
            product2 += 'exp(regStrength*'+speciesNames[i]+')*'
        for j in connectionDict.keys():
            product1 += speciesNames[j]+'**g_'+str(i)+'_'+str(j)+'*'
            product2 += speciesNames[j]+'**h_'+str(i)+'_'+str(j)+'*'
            #if includeRegularizer:
            #    # 1.30.2013
            #    product1 += 'exp(1./speciesNames[j])*'
            #    product2 += 'exp(speciesNames[j])*'
        if useDeltaGamma:
            net.addRateRule( speciesNames[i],                                   \
                'delta_'+str(i)+'*( '+product1[:-1]+' - '                       \
               +'gamma_'+str(i)+'*( '+product2[:-1]+' ) )' )
        else:
            net.addRateRule( speciesNames[i],                                   \
                'alpha_'+str(i)+'*( '+product1[:-1]+' )- '                      \
               +'beta_'+str(i)+'*( '+product2[:-1]+' )  ' )
        # 06.22.2009 Ilya playing with parameterizations
        #net.addRateRule( speciesNames[i],                                      \
        #   'alpha_'+str(i)+'*( '+product1[:-1]+' )*( 1. - '                    \
        #   +'beta_'+str(i)+'*( '+product2[:-1]+' ) )' )
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


