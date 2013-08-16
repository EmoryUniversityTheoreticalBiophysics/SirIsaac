# PolynomialNetwork.py
#
# Bryan Daniels
# 7.20.2009
#
# SloppyCell network that does not have any ODEs, but simply
# consists of an nth degree polynomial.
# (Branched from LaguerreNetwork.py)
#

from SloppyCell.ReactionNetworks import *
import scipy.special

def PolynomialNetwork(degree,outputName='output',netid='PolynomialNet'):
    """
    SloppyCell ('network') model that does not have any ODEs, but simply
    consists of an nth degree polynomial.
    
    Currently produces one output, with number of parameters = degree+1:
        
    output(t) = sum_{i=0}^{degree} g_i x^i
    """
    
    net = Network(netid, name='Polynomial (network) model')
    net.addCompartment('Comp',name='Compartment')
    
    net.addParameter( 'degree', degree, isOptimizable=False )
    
    # Do I need this?  Unfortunately, I think I do...
    net.addSpecies('dummy','Comp',0.)
    net.addRateRule('dummy','0.')
    
    # add optimizable parameters
    for i in range(degree+1):
        net.addParameter( 'g'+str(i), 0., isOptimizable=True )
        
    net.addSpecies( outputName, 'Comp', 0. )
    
    polynomial = ''
    for i in range(degree+1):
        polynomial += 'g'+str(i)+'*time**'+str(i)+' + ' 
    
    net.addAssignmentRule( outputName, polynomial[:-2] )
        
    return net
    
    
	