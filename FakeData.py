# fakeData.py
#
# Bryan Daniels
# 7.20.2009
#
# Make fake data compatible with SloppyCell.

from SloppyCell.ReactionNetworks import *
import scipy

# (originally from runTranscriptionNetwork.py)
def noisyFakeData(net,numPoints,timeInterval,                                   \
        vars=None,noiseFracSize=0.1,seed=None,params=None,randomX=True,         \
        includeEndpoints=True,takeAbs=False,noiseSeed=None,typValOffsets=None):
    """
    Adds Gaussian noise to data: 
        mean 0, stdev noiseFracSize*("typical value" of variable)
        (By default, the "typical value" is the maximum value for that variable;
         see SloppyCell.ReactionNetworks.PerfectData.update_typical_vals)
        
    randomX             : if False, the data points are distributed evenly over
                          the interval.  if True, they are spread randomly and
                          evenly over each variable.
    includeEndpoints    : if True, the initial and final time are included as
                          part of the numPoints points 
                          (not sure if this works right with randomX=False)
    takeAbs (False)     : 5.1.2013 If True, take the absolute value of the 
                          data (to avoid having negative data).
    seed (None)         : Random seed for selection of timepoints.
    noiseSeed (None)    : Random seed for addition of Gaussian noise.
    typValOffsets (None): List (of the same length as vars) consisting of offsets
                          to be subtracted from var_typical_val before
                          calculating noiseSize (useful when using an offset to
                          fit powerLaw models).  Defaults to zeros.
    """
    if seed is not None: scipy.random.seed(seed)
    
    if vars is None:
        vars = net.dynamicVars.keys()
    if params is None:
        params = net.GetParameters()
        
    PerfectData.update_typical_vals([net],[timeInterval])
    
    if includeEndpoints:
        numPoints -= 2
    
    data = PerfectData.discrete_data(net,params,numPoints,timeInterval,         \
        vars=vars,random=randomX)
        
    if includeEndpoints:
        traj = net.integrate(timeInterval)
        for var in vars:
          for time in timeInterval:
            data[var][time] = ( traj.get_var_val(var,time), 0. )

    if noiseSeed is not None: scipy.random.seed(noiseSeed)

    if typValOffsets is None: typValOffsets = scipy.zeros(len(vars))

    for var,offset in zip(data.keys(),typValOffsets):
        noiseSize = noiseFracSize * ( net.get_var_typical_val(var) - offset )
        for key in data[var].keys():
            old = data[var][key]
            if noiseSize > 0:
                new0 = old[0] + scipy.random.normal(0.,noiseSize)
                if takeAbs: new0 = abs(new0)
                new = (new0, noiseSize)
            else:
                new = (old[0], 0.)
            data[var][key] = new
    
    return data

def noisyFakeDataFromData(data,numPoints,varName,noiseFracSize=0.1,seed=None):
    
    scipy.random.seed(seed)
    
    n = len(data)
    typicalSize = scipy.average(data)
    noiseSize = noiseFracSize * typicalSize
    
    fakeDataDict = {}
    
    for i in range(numPoints):
      xVal = scipy.random.randint(0,n)
      yVal = data[xVal]
      fakeDataDict[ xVal ] =                                                    \
        ( yVal + scipy.random.normal(0.,noiseSize), noiseSize )
        
    return {varName: fakeDataDict}