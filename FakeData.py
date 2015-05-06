# fakeData.py
#
# Bryan Daniels
# 7.20.2009
#
# Make fake data compatible with SloppyCell.

from SloppyCell.ReactionNetworks import *
import scipy

# (originally from runTranscriptionNetwork.py)
def noisyFakeData(net,numPoints,timeInterval,
        vars=None,noiseFracSize=0.1,seed=None,params=None,randomX=True,
        includeEndpoints=True,takeAbs=False,noiseSeed=None,typValOffsets=None,
        trueNoiseRange=None):
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
    trueNoiseRange (None): If None, the reported error bars on the data are
                          the same size as the stdev used to add the noise.
                          Otherwise, a constant fractional noise
                          for each variable is chosen from 'trueNoiseRange'
                          and used to add the noise, while 'noiseFracSize' 
                          is reported as the stdev in the data.
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

    if trueNoiseRange is None:
        noiseFracSizeList = [noiseFracSize for var in data.keys()]
    elif len(trueNoiseRange) == 2:
        if noiseSeed is not None: scipy.random.seed(noiseSeed+1)
        a,b = trueNoiseRange
        noiseFracSizeList = scipy.random.uniform(a,b,len(data.keys()))
    else:
        raise Exception, "Unrecognized form of trueNoiseRange"

    if noiseSeed is not None: scipy.random.seed(noiseSeed)

    if typValOffsets is None: typValOffsets = scipy.zeros(len(vars))

    for var,offset,trueNoiseFracSize in zip(data.keys(),typValOffsets,noiseFracSizeList):
        trueNoiseSize = trueNoiseFracSize * ( net.get_var_typical_val(var) - offset )
        reportedNoiseSize = noiseFracSize * ( net.get_var_typical_val(var) - offset )
        for key in data[var].keys():
            old = data[var][key]
            if trueNoiseSize > 0:
                new0 = old[0] + scipy.random.normal(0.,trueNoiseSize)
                if takeAbs: new0 = abs(new0)
                new = (new0, reportedNoiseSize)
            else:
                new = (old[0], reportedNoiseSize)
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