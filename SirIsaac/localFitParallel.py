# localFitParallel.py
#
# Bryan Daniels
# 1.23.2012
# 3.21.2012 copied from sampleParallel.py
#
# Uses design pattern at 
# http://www.shocksolution.com/2010/04/managing-a-pool-of-mpi-processes-with-python-and-pypar/
# to run local fits in parallel.  
# Used in FittingProblem.SloppyCellFittingModel.localFitToData_parallel
#

#!/usr/bin/env python
from numpy import *
from mpi4py import MPI
import time

from .fittingProblem import *
import sys

## 7.3.2012 disable SloppyCell's parallel stuff
## see SloppyCell's __init__.py
sc = IO.SloppyCell
modules = [ sc,sc.ReactionNetworks,Ensembles,Dynamics,
            Collections,PerfectData ]
import socket
for module in modules:
    module.HAVE_MPI = False
    module.num_procs = 1
    module.my_rank = 0
    module.my_host = socket.gethostname()

# Constants
MASTER_PROCESS = 0
DIE_INDEX = -1

comm = MPI.COMM_WORLD
MPI_myID = comm.Get_rank()
num_processors = comm.Get_size()

if num_processors < 2:
    raise Exception("mpi4py has failed to initialize more than one processor.")

# read in arguments from command line file name
if len(sys.argv) < 2:
    print("Usage: python localFitParallel.py inputDictFile.data [SloppyCell options]")
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
fitFunction = lambda startParamsIndex:                          \
    fittingProblem.localFitToData(fittingData,dataModel,        \
    retall=True,startParams=startParamsList[startParamsIndex])
allOutputsDict = {}



### Master Process ###
if MPI_myID == MASTER_PROCESS:
    from .simplePickle import save

    num_processors = comm.Get_size()
    print(("localFitParallel: "\
          "Master process found {} other worker processors".format(num_processors-1)))
    
    # list of startParams indices to pass to the workers
    work_array = list(range(len(startParamsList)))
    work_size = len(work_array)
    
    # Dispatch jobs to worker processes
    work_index = 0
    num_completed = 0
    
    # Start all worker processes
    for i in range(1, min(num_processors, work_size+1)):
        comm.send(work_index, dest=i)
        comm.send(work_array[work_index], dest=i)
        print(("localFitParallel: Sent work index {} to worker {}".format(work_index,i)))
        work_index += 1

    # Receive results from each worker, and send it new data
    # 3.26.2020 For now, wait for results from each worker in order.
    #           It may be possible to be more fancy and get results as they are done.
    for i in range(num_processors, work_size+1):
        proc = 1 + (i-1) % (num_processors - 1)
        results = comm.recv(source=proc)
        # save results
        result_index = num_completed
        allOutputsDict[result_index] = results
        num_completed += 1
        # start next
        comm.send(work_index, dest=proc)
        comm.send(work_array[work_index], dest=proc)
        print(("localFitParallel: Sent work index {} to worker {}".format(work_index,proc)))
        work_index += 1

    # Get results from remaining worker processes
    while num_completed < work_size:
        proc += 1
        if proc > (num_processors - 1): proc = 1
        results = comm.recv(source=proc)
        # save results
        result_index = num_completed
        allOutputsDict[result_index] = results
        num_completed += 1
        
    assert(num_completed == work_size)
    
    # Shut down worker processes
    for proc in range(1, num_processors):
        print(("localFitParallel: Stopping worker process {}".format(proc)))
        comm.send(DIE_INDEX, dest=proc)
        
    # Write data to file
    save(allOutputsDict,outputFilename)

else:
    ### Worker Processes ###
    continue_working = True
    while continue_working:
        
        work_index =  comm.recv(source=MASTER_PROCESS)
        
        if work_index == DIE_INDEX:
            continue_working = False
        else:
            work_array = comm.recv(source=MASTER_PROCESS)
            print(("localFitParallel: "\
                  "worker {} working on work_array {}...".format(MPI_myID,work_array)))
            
            startTime = time.clock()
            
            # do the work
            workerResults = list( fitFunction(work_array) )
            
            minimizationTimeSeconds = time.clock() - startTime
            workerResults.append(minimizationTimeSeconds)
            workerResults.append(work_array)
            
            comm.send(workerResults, dest=MASTER_PROCESS)
            
            print(("localFitParallel: "\
                  "worker {} finished with work_array {}".format(MPI_myID,work_array)))
#### while
#### if worker


