# localFitParallel.py
#
# Bryan Daniels
# 1.23.2012
# 3.21.2012 copied from sampleParallel.py
#
# Uses design pattern at 
# http://www.shocksolution.com/2010/04/managing-a-pool-of-mpi-processes-with-python-and-pypar/
# to run local fits in parallel.  
# Used in FittingProblem.SloppyCellFittingModel.localFitToData_pypar
#

#!/usr/bin/env python
from numpy import *
import pypar
import time

#from generateFightData import *
from fittingProblem import *
import sys

# for parallel computation supported by SloppyCell
#import SloppyCell.ReactionNetworks.RunInParallel as Par

# 7.3.2012 disable SloppyCell's parallel stuff
# see SloppyCell's __init__.py
sc = IO.SloppyCell
modules = [ sc,sc.ReactionNetworks,Ensembles,Dynamics,
            Collections,PerfectData ]
import socket
for module in modules:
    module.HAVE_PYPAR = False
    module.num_procs = 1
    module.my_rank = 0
    module.my_host = socket.gethostname()

# Constants
MASTER_PROCESS = 0
WORK_TAG = 1
DIE_TAG = 2

MPI_myID = pypar.rank() #Par.my_rank 
num_processors = pypar.size() #Par.num_procs

if num_processors < 2:
    raise Exception, "Pypar has failed to initialize more than one processor."

# read in arguments from command line file name
if len(sys.argv) < 2 or len(sys.argv) > 2:
    print "Usage: python localFitParallel.py inputDictFile.data"
    exit()
inputDictFile = sys.argv[1]
inputDict = load(inputDictFile)
fittingProblem = inputDict['fittingProblem']
fittingData = inputDict['fittingData']
indepParamsList = inputDict['indepParamsList']
dataModel = inputDict['dataModel'] 
startParamsList = inputDict['startParamsList']
outputFilename = inputDict['outputFilename']
inputDict['outputDict'] = {}
num_work_processors = num_processors - 1 # can we make this num_processors?
fitFunction = lambda startParamsIndex:                          \
    fittingProblem.localFitToData(fittingData,dataModel,        \
    retall=True,startParams=startParamsList[startParamsIndex])
allOutputsDict = {}



### Master Process ###
if MPI_myID == MASTER_PROCESS:
    from simplePickle import save

    num_processors = pypar.size()
    print "Master process found " + str(num_processors) + " worker processors."
    
    # Create a list of dummy arrays to pass to the worker processes
    #work_size = 10
    #work_array = range(0,work_size)
    #for i in range(len(work_array)):
    #    work_array[i] = arange(0.0, 10.0)
    
    # list of startParams indices to pass to the workers
    work_array = range(len(startParamsList))
    work_size = len(work_array)
    
    # Dispatch jobs to worker processes
    work_index = 0
    num_completed = 0
    
    
    
    
    # Start all worker processes
    for i in range(1, min(num_processors, work_size+1)):
        pypar.send(work_index, i, tag=WORK_TAG)
        pypar.send(work_array[work_index], i)
        print "Sent work index " + str(work_index) + " to processor " + str(i)
        work_index += 1

    # Receive results from each worker, and send it new data
    for i in range(num_processors, work_size+1):
        results, status = pypar.receive(source=pypar.any_source, tag=pypar.any_tag, return_status=True)
        # save results
        result_index = results.pop()
        allOutputsDict[result_index] = results
        index = status.tag
        proc = status.source
        num_completed += 1
        # start next
        pypar.send(work_index, proc, tag=WORK_TAG)
        pypar.send(work_array[work_index], proc)
        print "Sent work index " + str(work_index) + " to processor " + str(proc)
        work_index += 1

    # Get results from remaining worker processes
    while num_completed < work_size: #-1
        results, status = pypar.receive(source=pypar.any_source, tag=pypar.any_tag, return_status=True)
        # save results
        result_index = results.pop()
        allOutputsDict[result_index] = results
        num_completed += 1
    
    # Shut down worker processes
    for proc in range(1, num_processors):
        print "Stopping worker process " + str(proc)
        pypar.send(-1, proc, tag=DIE_TAG)
        
    # Write data to file
    save(allOutputsDict,outputFilename)

else:
    ### Worker Processes ###
    continue_working = True
    while continue_working:
        
        work_index, status =  pypar.receive(                        \
            source=MASTER_PROCESS, tag=pypar.any_tag,               \
            return_status=True)
        
        if status.tag == DIE_TAG:
            continue_working = False
        else:
            work_array, status = pypar.receive(                     \
                source=MASTER_PROCESS, tag=pypar.any_tag,           \
                return_status=True)
            work_index = status.tag
            
            # Code below simulates a task running
            #time.sleep(random.random_integers(low=0, high=5))
            #result_array = work_array.copy()
            
            #from simplePickle import save
            #save(work_index,"temporary_debug_"+str(work_array)+".data")
            
            startTime = time.clock()
            
            # do the work
            workerResults = list( fitFunction(work_array) )
            
            minimizationTimeSeconds = time.clock() - startTime
            workerResults.append(minimizationTimeSeconds)
            workerResults.append(work_array)
            
            pypar.send(workerResults, destination=MASTER_PROCESS, tag=work_index)
#### while
#### if worker

pypar.finalize()
print "MPIEND"
