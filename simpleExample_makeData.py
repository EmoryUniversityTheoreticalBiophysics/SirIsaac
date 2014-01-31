# simpleExample_makeData.py
#
# Bryan Daniels
# 1.29.2014
#
# Construct fake data to be used in
# simpleExample.ipynb.
#
# Data taken at times t sampled uniformly between
# 0 and 1 and initial conditions x0 sampled
# uniformly between 1 and 2.
#
# Output is a sinusoid varying between 1 and 2
# with the given initial condition x0.
#
# x = 3/2 + 1/2 sin( 4 pi t + arcsin(2x0 - 3) ) + N(0,0.1)
#

import scipy

numSamples = 100
noiseStd = 0.1

scipy.random.seed(10000)
times = scipy.random.rand(numSamples)
noises = scipy.random.normal(scale=noiseStd,size=numSamples)
x0s = 1. + scipy.random.rand(numSamples)
xs = 1.5 + 0.5 * scipy.sin( 4.*scipy.pi*times \
                            + scipy.arcsin(2.*x0s - 3.) )
xsNoisy = xs + noises

# in the output file,
# column 1: initial condition
# column 2: measurement time
# column 3: measurement value
# column 4: measurement uncertainty (standard deviation)

noiseStds = scipy.repeat(noiseStd,100)

data = zip( x0s, times, xsNoisy, noiseStds )

scipy.savetxt('simpleExample_data.txt',data)



