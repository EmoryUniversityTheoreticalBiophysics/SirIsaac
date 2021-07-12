# mpi_basic.py
#
# Bryan Daniels
# 2020/04/20
#
# Basic MPI functionality.  As of 2020/06/03, only used by test_parallel.
#
# Run like so:
#   mpirun -np 3 python mpi_test.py test_input_filename
#
# RUN apt-get update & apt-get install python2-dev build-essential libssl-dev libffi-dev libxml2-dev libxslt1-dev zlib1g-dev

import sys
from simplePickle import load,save

from mpi4py import MPI
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_procs = comm.Get_size()

if __name__ == '__main__':
    
    if len(sys.argv) < 2 or len(sys.argv) > 2:
        print "Usage: mpirun -np [numprocs] python mpi_test.py test_input_filename"
        exit()
    test_input_filename = sys.argv[1]
    
    if num_procs < 2:
        raise Exception("No worker processes detected.")

    if my_rank == 0:
    
        # master process
        
        file_data = load(test_input_filename)
    
        # send work
        for worker in range(1,num_procs):
            comm.send(("file_data['test']",
                       {'comp':{1:2},'file_data':file_data}),
                      dest=worker)
            
        # get results
        for worker in range(1,num_procs):
            msg = comm.recv(source=worker)
            print("mpi_test: Worker {} said {}".format(worker,msg))
            file_data['test_output'] = msg
            
        # stop workers
        for worker in range(1, num_procs):
            comm.send(SystemExit(), dest=worker)
            
        save(file_data,file_data['output_filename'])
            
    while my_rank != 0:
        # worker process
        
        # Wait for a message
        message = comm.recv(source=0)
        
        if isinstance(message, SystemExit):
            sys.exit()
        
        command, msg_locals = message
        locals().update(msg_locals)
        
        result = eval(command)
        comm.send(result, dest=0)
        

