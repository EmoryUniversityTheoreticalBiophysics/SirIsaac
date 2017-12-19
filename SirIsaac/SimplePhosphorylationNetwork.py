# SimplePhosphorylationNetwork.py
#
# Bryan Daniels
# 9.6.2013
#
# SloppyCell network that does not have any ODEs, but simply
# consists of a simple model for the phosphorylation example.
# (Branched from PolynomialNetwork.py)
#

from SloppyCell.ReactionNetworks import *

def SimplePhosphorylationNetwork(outputName='totalPhos',inputName='k23p',       \
    netid='SimplePhosphorylationNet',offset=1.,offsetName='totalPhos_init'):
    """
    SloppyCell ('network') model that does not have any ODEs, but simply
    consists of a simple model for the phosphorylation example.
    
    Currently produces one output, with 5 parameters (a,b,c,d,t0):
        
    output(t) = 
        offset + [ a + b/2(1+tanh( (log(input)-d)/c )) ] * [ 1 - exp(-t/t0) ]
    
    with default parameters
    a = 2, b = 0.3, c = 1, d = 1, t0 = 0.5
    """
    
    net = Network(netid, name='Simple Phosphorylation network')
    net.addCompartment('Comp',name='Compartment')

    # Do I need this?  Unfortunately, I think I do...
    net.addSpecies('dummy','Comp',0.)
    net.addRateRule('dummy','0.')
    
    # add optimizable parameters (and input parameter 'k')
    params = ['a','b','c','d','t0',inputName,offsetName]
    defaultParams = [2.,0.3,1.,1.,0.5,1.,offset]
    for param,defaultParam in zip(params,defaultParams):
        net.addParameter( param, defaultParam, isOptimizable=True )
        
    net.addSpecies( outputName, 'Comp', 0. )
    
    functionStr = offsetName +                                              \
                '+ ( a + b/2.*(1.+tanh( (log('+inputName+')-d)/c )) ) '     \
                '* ( 1. - exp(-time/t0) )'
    
    net.addAssignmentRule( outputName, functionStr )
        
    return net

