from mpi4py import MPI

import numpy as np
import glob
import csv
from collections import defaultdict
from pathlib import Path
import os
from timeit import default_timer as timer
from scipy.io import savemat
import h5py

def fill_op(xmem, ymem, dt):
    x = np.frombuffer(xmem)
    y = np.frombuffer(ymem)
    y[~np.isfinite(y)] = x[~np.isfinite(y)]
    y[:] = y


FILL = MPI.Op.Create(fill_op, commute=True)

def readCSVDict(filename):
    out = []
    with open(filename, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=' ')
        out=csv_reader.fieldnames.copy()
        out = [x for x in out if x.strip()]

    return out
    
def readCSV(filename):
    out = defaultdict(list)

    with open(filename, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=' ')
        for line in csv_reader:
            for key, value in line.items():
                out[key].append(value)

    return out

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("[MPI environment] rank = %d, size = %d.\n" % (rank, size))
prefix = "results_"
postfix = ".dat"
infiles = glob.glob("postProcessing/"+prefix+"*"+postfix)
if rank==0:
    print("Found "+ str(len(infiles))+" dat files.")
inputArray = readCSVDict(infiles[0])
if rank == 0:
    print(inputArray)
    
N = len(infiles)
num_probes = len(readCSV(infiles[0])[inputArray[0]])

# Initial output fields
output_file = "eta_withfields.h5"
corrupted = 0

nbsteps = np.arange(N)
nsub = len(nbsteps[rank::size])
#print(nbsteps[rank::size])

offset = np.zeros(1, dtype='int64')
comm.Exscan([np.array([nsub], dtype='int64'), MPI.INT64_T], [offset, MPI.INT64_T])
nbsteps = offset + np.arange(nsub)
print(nbsteps)

fields = {
    "t": np.nan + np.zeros((nsub,)),
}
for arr in inputArray:
    fields[arr] = np.nan + np.zeros((nsub, num_probes))

if rank == 0:
    t_start = timer()

with h5py.File(output_file, 'w', driver='mpio', comm=comm) as file:
    # Create a dataset collectively
    dset = file.create_dataset('t', shape=(N, ), dtype='float64')
    
    for arr in inputArray:
        dset = file.create_dataset(arr, shape=(N, num_probes, ), dtype='float64')
    
    
with h5py.File(output_file, 'r+', driver='mpio', comm=comm) as file:
    
    for i, ii in zip(nbsteps, range(nsub)):
        # display some information
        print(str(i) + "/" + str(N))
        
        # assign time
        fields["t"][ii] = float(infiles[i].rstrip(postfix).split(prefix)[-1])

        # Read optional field data
        data = readCSV(infiles[i])

        for arr in inputArray:
            fields[arr][ii, :] = data[arr]

    # assign time
    file['t'][offset[0]:offset[0]+nsub] = fields["t"]
    print("[%d]: writing array t into h5 file." % (rank, ))

    for arr in inputArray:
        file[arr][offset[0]:offset[0]+nsub,:] = fields[arr]
        print("[%d]: writing matrix %s into h5 file." % (rank, arr))

comm.Barrier()

if rank==0:
    t_end = timer()
    t_total = t_end - t_start

    print('  Total time:\t\t\t  %6.2f seconds' % (t_total))

MPI.Finalize()
