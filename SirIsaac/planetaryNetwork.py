# PlanetaryNetwork.py
#
# Bryan Daniels
# 8.6.2013
#
# A SloppyCell implementation of planetary motion.
#

from SloppyCell.ReactionNetworks import *

def Planetary_net(r_init=1,theta_init=0,netid='Planetary'):
    """
    A SloppyCell implementation of planetary dynamics.
    
    Units of distance are rc = GM/(v0^2)
    Units of time are t0 = rc/v0 = GM/(v0^3)
    where G  = gravitational constant
    M  = mass of sun
    v0 = initial speed of object
    (assumed to be moving perpendicular
    to the line connecting it to the sun)
    """
    
    # make new SloppyCell network
    net = Network(netid, name='Planetary')
    net.addCompartment('Comp',name='Compartment')

    net.addParameter('r_init',r_init)
    net.addParameter('theta_init',theta_init)

    net.addSpecies('r','Comp','r_init')
    net.addSpecies('drdt','Comp',0)
    net.addSpecies('theta','Comp','theta_init')

    net.addRateRule('r','drdt')
    net.addRateRule('drdt','1./r**2 * (r_init**2/r - 1)')
    net.addRateRule('theta','r_init/r**2')

    return net



















