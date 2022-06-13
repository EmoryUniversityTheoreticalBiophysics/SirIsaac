# simpleSinusoidalNetwork.py
#
# Bryan Daniels
# 5.16.2015
#
# SloppyCell network that does not have any ODEs, but
# consists of a simple sinusoidal model.
# (Branched from polynomialNetwork.py)
#

from SloppyCell.ReactionNetworks import *

def SimpleSinusoidalNetwork(outputNameList=[],
    netid='SimpleSinusoidalNet'):
    """
    SloppyCell ('network') model that does not have any ODEs, but
    implements a simple sinusoidal model.
    
    Produces outputs with names given in outputNameList:
        
    output0(t) = y0_1 + A_1 sin(omega t + phi)
    output1(t) = y0_2 + A_2 sin(omega t + phi + phi_2)
    output2(t) = y0_3 + A_3 sin(omega t + phi + phi_3)
    
    with default parameters
    y0_i = 0, A_i = 1, omega = 8 pi / 5, [outputName]_init = 1., phi_i = pi/2.
    
    and phi calculated to match initial conditions (approximately
    so if the initial conditions cannot be matched exactly).
    """
    
    if len(outputNameList) != 3:
        raise Exception("Only 3-dimensional network currently supported.")
    
    net = Network(netid, name='Simple sinusoidal network')
    net.addCompartment('Comp',name='Compartment')

    # Do I need this?  Unfortunately, I think I do...
    net.addSpecies('dummy','Comp',0.)
    net.addRateRule('dummy','0.')

    net.add_func_def('arctan2',('y','x'),
                     '2.*atan( y / (sqrt(x**2 + y**2) + x) )')

    # Add species and parameters
    for i,outputName in enumerate(outputNameList):
        net.addSpecies( outputName, 'Comp', 0. )
        # add optimizable parameters
        istr = str(i+1)
        params = [outputName+'_init','A_'+istr,'y0_'+istr]
        defaultParams = [1.,1.,0.]
        for param,defaultParam in zip(params,defaultParams):
            net.addParameter( param, defaultParam, isOptimizable=True )
    net.addParameter('phi_2', 3.14159/2., isOptimizable=True )
    net.addParameter('phi_3', 3.14159/2., isOptimizable=True )
    net.addParameter('omega', 8.*3.14159/5., isOptimizable=True )

    # Add functions (first calculate phi from initial conditions)
    net.addSpecies('x1','Comp',0.)
    net.addSpecies('x2','Comp',0.)
    net.addSpecies('x3','Comp',0.)
    x1functionStr = '( '+outputNameList[0]+'_init - y0_1 ) / A_1'
    x2functionStr = '( '+outputNameList[1]+'_init - y0_2 ) / A_2'
    x3functionStr = '( '+outputNameList[2]+'_init - y0_3 ) / A_3'
    net.addAssignmentRule('x1',x1functionStr)
    net.addAssignmentRule('x2',x2functionStr)
    net.addAssignmentRule('x3',x3functionStr)

    net.addSpecies('cproj','Comp',0.)
    cFunctionStr = '(x1 * sin(phi_2 - phi_3) + x2 * sin(phi_3) - x3 * sin(phi_2))'      \
                   ' / (sin(phi_2 - phi_3)**2 + sin(phi_2)**2 + sin(phi_3)**2)'
    net.addAssignmentRule('cproj',cFunctionStr)

    net.addSpecies('x1tilde','Comp',0.)
    net.addSpecies('x2tilde','Comp',0.)
    x1tildeStr = 'x1 + cproj * sin(phi_3 - phi_2)'
    x2tildeStr = 'x2 - cproj * sin(phi_3)'
    net.addAssignmentRule('x1tilde',x1tildeStr)
    net.addAssignmentRule('x2tilde',x2tildeStr)

    net.addSpecies('phi','Comp',0.)
    phiStr = 'arctan2( x1tilde*sin(phi_2), x2tilde - x1tilde*cos(phi_2) )'
    net.addAssignmentRule('phi',phiStr)

    # Add functions
    for i,outputName in enumerate(outputNameList):
        istr = str(i+1)
        if i == 0:
            functionStr = 'y0_'+istr+' + A_'+istr+' * sin(omega*time + phi)'
        else:
            functionStr = 'y0_'+istr+' + A_'+istr+                                      \
                          ' * sin(omega*time + phi + phi_'+istr+')'
        net.addAssignmentRule( outputName, functionStr )

    return net

