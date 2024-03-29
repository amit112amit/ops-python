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

savefrequency = 250000  # DO NOT INCREASE FURTHER due to h5py chunking

# Create the model
ops = Model()

# New simulation or restore an old one from a saved state.
if os.path.isfile(statefile):
    ops.restoreSavedState(statefile)
    ops.vtkfilecount += 1
elif os.path.isfile(vtkfile):
    ops.initializeFromVTKFile(vtkfile)
else:
    print('FileNotFound :', vtkfile, 'or', statefile)
    quit()

# Read the simulation schedule and identify starting state
schedule = pandas.read_csv('schedule-{}.dat'.format(taskid), sep='\t')
if schedule.ViterMax.sum() < savefrequency:
    savefrequency = int(schedule.ViterMax.sum())

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
output = {'Volume': collections.deque(maxlen=savefrequency),
          'MSD': collections.deque(maxlen=savefrequency),
          'MSDT': collections.deque(maxlen=savefrequency),
          'RMSAngleDeficit': collections.deque(maxlen=savefrequency)}

if not os.path.isfile(outfile):
    with h5py.File(outfile, 'w') as hfile:
        for key in output:
            hfile.create_dataset(key, dtype='float32', shape=(savefrequency,),
                                 maxshape=(schedule.ViterMax.sum(),),
                                 chunks=(savefrequency,), compression='gzip')
    outputsize = 0
else:
    with h5py.File(outfile, 'a') as hfile:
        outputsize = len(hfile['Volume'])

# -------------------------------------------------------------------------------
# Loop over every row in `schedule` starting from row `rowid`
# -------------------------------------------------------------------------------
with h5py.File(polydatafile, 'a') as vtkhfile:
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
            radialmsd, tangentialmsd = ops.bothMSD()
            output['Volume'].append(ops.volume)
            output['MSD'].append(radialmsd)
            output['MSDT'].append(tangentialmsd)
            output['RMSAngleDeficit'].append(ops.rmsAngleDeficit)

            # -----------------------------------------------------------------------
            # Save visualization, simulation state and/or the output statistics
            # -----------------------------------------------------------------------
            if (step % savefrequency) == savefrequency - 1:
                ops.timestep = step
                ops.vtkfilecount = vtkcount
                ops.writeSimulationState(statefile)
                with h5py.File(outfile, 'a') as hfile:
                    for key, data in output.items():
                        hfile[key].resize((outputsize + savefrequency,))
                        hfile[key].write_direct(numpy.array(data), numpy.s_[:], numpy.s_[outputsize:])
                        data.clear()
                outputsize += savefrequency

            points, normals, cells = ops.polyDataParts()
            pointskey = '/T{0}/Points'.format(vtkcount)
            normalskey = '/T{0}/Normals'.format(vtkcount)
            cellskey = '/T{0}/Polygons'.format(vtkcount)
            vtkhfile.attrs['TimeSteps'] = vtkcount + 1
            vtkhfile.create_dataset(pointskey, data=points, dtype='float32')
            vtkhfile.create_dataset(normalskey, data=normals, dtype='float32')
            vtkhfile.create_dataset(cellskey, data=cells, dtype='uint32')
            vtkcount += 1

            ops.updatePreviousX()
            step += 1

        # Reset rowstep for next row
        rowstep = 0

# ---------------------------------------------------------------------------
# Save any leftover output data
# ---------------------------------------------------------------------------
if len(output['Volume']) > 0:
    with h5py.File(outfile, 'a') as hfile:
        for key, data in output.items():
            hfile[key].resize((outputsize + len(data),))
            hfile[key].write_direct(numpy.array(data), numpy.s_[:], numpy.s_[outputsize:])
