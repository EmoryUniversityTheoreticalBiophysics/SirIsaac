# simpleExample_makeData.py
#
# Bryan Daniels
# 1.29.2014
# updated 2022/6/13
#
# Construct fake data to be used in
# simpleExample.ipynb.
#
# Data taken at times t sampled uniformly between
# 0 and 1, initial conditions x0 sampled
# uniformly between 1 and 2, and decay rate tau
# sampled uniformly between 0 and 1.
#
# Output is an exponential decay with the given
# initial condition x0 and decay rate tau:
#
# x = x_0 * exp(-t/tau) + N(0,0.1)
#

import numpy as np

numSamples = 100
noiseStd = 0.1

np.random.seed(10000)
times = np.random.rand(numSamples)
noises = np.random.normal(scale=noiseStd,size=numSamples)
x0s = 1. + np.random.rand(numSamples)
taus = np.random.rand(numSamples)
xs = x0s * np.exp(-times/taus)
xsNoisy = xs + noises

# in the output file,
# column 1: initial condition x0
# column 2: decay rate tau
# column 3: measurement time
# column 4: measurement value
# column 5: measurement uncertainty (standard deviation)

noiseStds = np.repeat(noiseStd,100)

data = list(zip( x0s, taus, times, xsNoisy, noiseStds ))

np.savetxt('simpleExample_data.csv',data,delimiter=',')



