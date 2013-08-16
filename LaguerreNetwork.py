# LaguerreNetwork.py
#
# Bryan Daniels
# 6.12.2009
#
# SloppyCell ('network') model that does not have any ODEs, but simply
# consists of an nth degree Laguerre polynomial [times e^(-t/2)], which
# should be a set of (nearly) orthogonal basis functions.  (They're
# exactly orthogonal when you have infinite data out to infinity.)
#

from SloppyCell.ReactionNetworks import *
import scipy.special

# Note: there may be a more natural way of implementing this
# in SloppyCell that doeesn't call this model a ReactionNetwork.
# Maybe ask Ryan about this...

def LaguerreNetwork(degree,outputName='output',netid='PolynomialNet'):
    """
    SloppyCell ('network') model that does not have any ODEs, but simply
    consists of an nth degree Laguerre polynomial [times e^(-t/2)], which
    should be a set of (nearly) orthogonal basis functions.  (They're
    exactly orthogonal when you have infinite data out to infinity.)
    
    Currently produces one output, with degree+3 parameters:
        
    output(t) = C + sum_{i=0}^{degree} g_i L_i(2t/alpha) exp(-x/alpha)
    
    (The optimizable time constant parameter is squared to 
     avoid problems with negative alpha.)
    """
    
    net = Network(netid, name='Laguerre (network) model')
    net.addCompartment('Comp',name='Compartment')
    
    net.addParameter( 'degree', degree, isOptimizable=False )
    
    # Do I need this?  Unfortunately, I think I do...
    net.addSpecies('dummy','Comp',0.)
    net.addRateRule('dummy','0.')
    
    # add optimizable parameters
    net.addParameter( 'C', 0., isOptimizable=True )
    net.addParameter( 'alpha', 1., isOptimizable=False )
    net.addParameter( 'sqrt_abs_alpha', 1., isOptimizable=True )
    for i in range(degree+1):
        net.addParameter( 'g'+str(i), 0., isOptimizable=True )
        
    # square the optimizable parameter to avoid problems with negative alpha
    net.addAssignmentRule( 'alpha', 'sqrt_abs_alpha**2' )
        
    net.addSpecies( outputName, 'Comp', 0. )
    
    polynomial = ''
    for i in range(degree+1):
        polynomial += 'g'+str(i)                                                \
          +'*( '+poly2str(scipy.special.laguerre(i),'2.*time/alpha')+' ) + ' 
    
    net.addAssignmentRule( outputName, 'C + exp(-time/alpha)*'                  \
        +'( '+polynomial[:-2]+')' )
        
    return net
    
def poly2str(p, variableString):
    string = ''
    coeffArray = p.c
    degree = len(coeffArray) - 1
    for i in range(degree+1):
        string += str(coeffArray[i])+'*('+variableString+')**'+str(degree-i)+' + '
    return string[:-2]
    
    
    
    
    
	