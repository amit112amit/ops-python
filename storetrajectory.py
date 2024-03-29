import collections
import itertools
import os
from math import log, sqrt

import numpy
import h5py
import pandas

from ops import Model, Solver

# -------------------------------------------------------------------------------
# Get batch job task id if any.
# -------------------------------------------------------------------------------
try:
    taskid = os.environ['SGE_TASK_ID']
except KeyError:
    taskid = 1

# -------------------------------------------------------------------------------
# Initialization
# -------------------------------------------------------------------------------
vtkfile = 'T7.vtk'
polydatafile = 'VTKFile-' + str(taskid) + '.h5'
outfile = 'DetailedOutput-' + str(taskid) + '.h5'
statefile = 'SimulationState-' + str(taskid) + '.dat'

savefrequency = 1200  # DO NOT INCREASE FURTHER due to h5py chunking size limits

# Create the model
ops = Model()

# New simulation or restore an old one from a saved state.
if os.path.isfile(statefile):
    ops.restoreSavedState(statefile)
elif os.path.isfile(vtkfile):
    ops.initializeFromVTKFile(vtkfile)
else:
    print('FileNotFound :', vtkfile, 'or', statefile)
    quit()

# Read the simulation schedule and identify starting state
schedule = pandas.read_csv('schedule-{}.dat'.format(taskid), sep='\t')
vtkcount = ops.vtkfilecount
step = ops.timestep
skiprows, totalsteps = next(itertools.dropwhile(
    lambda x: x[1] <= step, enumerate(schedule.ViterMax.cumsum())))
rowstep = schedule.ViterMax[skiprows] - (totalsteps - step)

# Create solver object
solver = Solver(ops)

# -------------------------------------------------------------------------------
# Prepare output containers and the HDF5 file
# -------------------------------------------------------------------------------
N = ops.numberofpoints
output = collections.deque(maxlen=savefrequency*N*3)

if not os.path.isfile(outfile):
    with h5py.File(outfile, 'w') as hfile:
        hfile.create_dataset('Positions', dtype='float32', shape=(savefrequency, N*3),
                             maxshape=(schedule.ViterMax.sum(), N*3),
                             chunks=(savefrequency, N*3), compression='gzip')
    outputsize = 0
else:
    with h5py.File(outfile, 'a') as hfile:
        outputsize = len(hfile['Positions'])

    # -------------------------------------------------------------------------------
    # Loop over every row in `schedule` starting from row `rowid`
    # -------------------------------------------------------------------------------
for index, row in itertools.islice(schedule.iterrows(), skiprows, None):
    ops.fvk = row.Gamma
    ops.morseWellWidth = log(2.0)*100.0/row.PercentStrain
    ops.constraint = row.AreaConstr
    if index == 0 and rowstep == 0:
        ops.brownCoeff = 0.0
        ops.viscosity = 0.0
        solver.solve()
        ops.saveInitialPosition()

    ops.viscosity = row.Alpha
    ops.brownCoeff = sqrt(2*row.Alpha / row.Beta)

    rowsteps = int(row.ViterMax)
    printstep = int(row.PrintStep)

    ops.updatePreviousX()

    # ---------------------------------------------------------------------------
    # Inner loop
    # ---------------------------------------------------------------------------
    for i in range(rowstep, rowsteps):
        ops.generateParallelKicks()
        ops.lagrangeCoeff = 10.0
        ops.penaltyCoeff = 1000.0

        # -----------------------------------------------------------------------
        # Augmented Lagrangian loop
        # -----------------------------------------------------------------------
        notdone = True
        maxiter = 10
        while notdone and maxiter > 0:
            solver.solve()
            ops.uzawaUpdate()
            maxiter -= 1
            notdone = not ops.constraintSatisfied()

        ops.applyKabschAlgorithm()
        ops.updateTriangles()

        # -----------------------------------------------------------------------
        # Collect statistics for output
        # -----------------------------------------------------------------------
        output.append(ops.positionvectors.flatten())

        # -----------------------------------------------------------------------
        # Save visualization, simulation state and/or the output statistics
        # -----------------------------------------------------------------------
        if (step % savefrequency) == savefrequency - 1:
            ops.timestep = step
            ops.vtkfilecount = vtkcount
            ops.writeSimulationState(statefile)
            with h5py.File(outfile, 'a') as hfile:
                hfile['Positions'].resize((outputsize + savefrequency, N*3))
                data = numpy.asarray(output).squeeze()
                hfile['Positions'].write_direct(data, numpy.s_[:], numpy.s_[
                                                outputsize:outputsize + savefrequency, :])
            output.clear()
            outputsize += savefrequency

        if i % printstep == 0 and printstep <= rowsteps:
            vtkcount += 1

        ops.updatePreviousX()
        step += 1

    # Reset rowstep for next row
    rowstep = 0

# ---------------------------------------------------------------------------
# Save any leftover output data
# ---------------------------------------------------------------------------
if len(output) > 0:
    with h5py.File(outfile, 'a') as hfile:
        hfile['Positions'].resize((outputsize + len(output), N*3))
        data = numpy.asarray(output).squeeze()
        hfile['Positions'].write_direct(data, numpy.s_[:], numpy.s_[
                                        outputsize:outputsize + len(output), :])
