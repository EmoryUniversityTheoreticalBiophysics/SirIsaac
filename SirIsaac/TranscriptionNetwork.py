# TranscriptionNetwork.py
#
# Bryan Daniels
# 06.09.2009
# 07.05.2009 added input variables
#
# A SloppyCell implementation of small transcription regulation networks.
#


from SloppyCell.ReactionNetworks import *

def TranscriptionNetworkZiv(netid='TranscriptionNetwork'):
    """
    Creates SloppyCell transcription network with 4 species as in 
    ZivNemWig07.  Each of three transcription factors (X, Y, and Z) 
    is regulated (either activated or inhibited) by exactly one 
    other transcription factor.  GFP is always 
    down-regulated by transcription factor 3.
    
    For now, just implements the first circuit from ZivNemWig07
    (number 1 in Figure 1) until I figure out a better way of
    enumerating the possibilities.
    """
    
    net = Network(netid, name='Transcription Network')
    net.addCompartment('Comp',name='Compartment')
    
    net.addParameter( 'n', 2, isOptimizable=False ) # Hill coefficient
    
    nameList = ['TFX','TFY','TFZ','GFP']
    initialConditions = [0.,0.,0.,1.]
    
    suffixList = ['X','Y','Z','G']
    connectionList = ['TFX','TFX','TFY','TFZ']  # which TF regulates each
    signList = [-1,-1,-1,-1]                    # +1 activation, -1 inhibition
    
    
    # transcription factor species
    for name,IC in zip(nameList,initialConditions):
        net.addSpecies( name, 'Comp', IC )
    
    # decay rates (rG set to definite number in paper)
    defaultr = 1.
    for suffix in suffixList[:-1]:
        net.addParameter( 'r'+suffix, defaultr, isOptimizable=True )
    net.addParameter( 'rG', defaultr, isOptimizable=False) # <---- need to set this
    
    # Michaelis constants
    defaultK = 1.
    for suffix in suffixList:
        net.addParameter( 'K'+suffix, defaultK, isOptimizable=True )
    
    # range parameters
    defaulta = 1.
    for suffix in suffixList:
        net.addParameter( 'a'+suffix, defaulta, isOptimizable=True )
    
    # leak parameter
    defaulta0 = 1.
    net.addParameter( 'a0', defaulta0, isOptimizable=True )
    
    # input parameters
    defaults = 1.
    for suffix in suffixList:
        net.addParameter( 's'+suffix, defaults, isOptimizable=False )
    
    # reaction rate rules
    for name,suffix,connection,sign in                                      \
        zip(nameList,suffixList,connectionList,signList):
      if sign == +1 :   # activated
        net.addRateRule( name,                                              \
            '-r'+suffix+'*'+name+' + a0 + a'+suffix+                        \
            '*('+connection+'/s'+suffix+')**n /'+'(K'+suffix+'**n + ('      \
            +connection+'/s'+suffix+')**n)' )
        #net.addRateRule( name,                                             \
        #   '-r'+suffix+'*'+name+' + a0 + a'+suffix+'*'+connection+'**n /'  \
        #   +'(K'+suffix+'**n + '+connection+'**n)' )
      else:             # inhibited
        net.addRateRule( name,                                              \
            '-r'+suffix+'*'+name+' + a0 + a'+suffix+'*K'+suffix+'**n /'     \
            +'(K'+suffix+'**n + ('+connection+'/s'+suffix+')**n)' )
        #net.addRateRule( name,                                             \
        #   '-r'+suffix+'*'+name+' + a0 + a'+suffix+'*K'+suffix+'**n /'     \
        #   +'(K'+suffix+'**n + '+connection+'**n)' )
    
    return net
	