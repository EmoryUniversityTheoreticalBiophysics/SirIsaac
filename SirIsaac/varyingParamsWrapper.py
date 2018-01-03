# varyingParamsWrapper.py
#
# Bryan Daniels
# 7.14.2009
#
# A wrapper for SloppyCell networks that makes their parameters
# vary on a given input signal.
#

import copy

def VaryingParamsNet_Polynomial(SloppyCellNetwork,degreeList,inputName='input'):
    """
    degreeList          : should be the length of the number of parameters in
                          the original SloppyCellNetwork, specifying the degree
                          of the polynomial that describes each parameter's
                          dependence on the input
    
    Note!  This messes with the network a fair bit...
    Check to be sure everything is correct.
    """
    net = SloppyCellNetwork.copy()
    paramList = net.GetParameters()
        
    if len(degreeList) != len(paramList):
        raise Exception, "VaryingParamsNet_Polynomial: "                        \
            +"Length of degreeList not equal to length of paramList."
    
    # remove current rate rules, assignment rules, and 
    # non-optimizable variables 
    # (we add them back after the new assignment rules, 
    # in case the old depend on the new ones)
    oldRateRules = copy.deepcopy(net.rateRules)
    oldRules = copy.deepcopy( net.assignmentRules )
    oldVariables = filter( lambda var: not var.is_optimizable,                  \
        copy.deepcopy(net.variables) )
    for name in oldRateRules.keys():
        net.remove_component(name)
    for name in oldRules.keys():
        net.remove_component(name)
    for var in oldVariables:
        net.remove_component(var.id)
    
    # add input parameter
    net.addParameter(inputName,0,isOptimizable=False)
    
    for id,value,degree in zip(paramList.keys(),paramList.values(),degreeList):
        
        # make the variable non-optimizable
        net.set_var_optimizable(id,False)
        
        # make new optimizable variables, and build up the polynomial
        # (the zeroth order term should be set to the original value)
        polynomial = ''
        for i in range(degree+1):
          newName = id+'_'+str(i)
          defaultValue = 0
          if i==0:
            defaultValue = value
          net.addParameter(newName, defaultValue, isOptimizable=True)
          polynomial += newName+'*'+inputName+'**'+str(i)+' + '
        polynomial = polynomial[:-3]  # shave off the last ' + '
        
        # update old assignment rules involving this variable
        #for var,rule in zip(oldRules.keys(),oldRules.values()):
        #  net.addAssignmentRule(var, rule.replace(id,'('+polynomial+')'))
        
        # set up polynomial dependence
        net.addAssignmentRule(id, polynomial)
        
    for var in oldVariables:
        net._add_variable(var)
    for name,rule in zip(oldRules.keys(),oldRules.values()):
        net.addAssignmentRule(name,rule)
    for name,rateRule in zip(oldRateRules.keys(),oldRateRules.values()):
        net.addRateRule(name,rateRule)
        
    return net
        
        
        
        
    
