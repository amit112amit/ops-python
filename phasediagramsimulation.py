import argparse
import collections
import itertools
import os
import logging
import sys
import time
from math import log, sqrt

sys.path.append('/usr/lib/python3.7/site-packages')

import numpy
import h5py
import pandas

from ops import Model, Solver

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# -------------------------------------------------------------------------------
# Check if there is an input argument or else look for an environment variable
# -------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Array task id.', required=False)
args = parser.parse_args()

if args.input is not None:
    taskid = int(args.input)
    logging.info('Command line argument provided `taskid=%d`.', taskid)
else:
    try:
        taskid = os.environ['SGE_TASK_ID']
        logging.info('Environment variable provided `taskid=%d`.', taskid)
    except KeyError:
        logging.info('No command line input or SGE_TASK_ID environment variable found. Using `taskid`=1.')
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
    logging.info('Restoring old simulation from %s...', statefile)
    ops.restoreSavedState(statefile)
elif os.path.isfile(vtkfile):
    logging.info('Starting a new simulation using %s.', vtkfile)
    ops.initializeFromVTKFile(vtkfile)
else:
    raise FileNotFoundError('Unable to find {} or {} '
            'in the current directory.'.format(vtkfile, statefile))

# Read the simulation schedule and identify starting state
schedule = pandas.read_csv('schedule-{}.dat'.format(taskid), sep='\t')
if schedule.ViterMax.sum() < savefrequency:
    savefrequency = int(schedule.ViterMax.sum())
    logging.info('Decreasing `savefrequency` to %d.', savefrequency)

vtkcount = ops.vtkfilecount
step = ops.timestep
skiprows, totalsteps = next(itertools.dropwhile(
    lambda x: x[1] <= step, enumerate(schedule.ViterMax.cumsum())))
rowstep = schedule.ViterMax[skiprows] - (totalsteps - step)

logging.info('Simulation will start from `rowstep=%d` and `step=%d`.', rowstep, step)

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
    logging.info('Creating new data file %s.', outfile)
    with h5py.File(outfile, 'w') as hfile:
        for key in output:
            hfile.create_dataset(key, dtype='float32', shape=(savefrequency,),
                                 maxshape=(schedule.ViterMax.sum(),),
                                 chunks=(savefrequency,), compression='gzip')
    outputsize = 0
else:
    logging.info('Appending output data to %s.', outfile)
    with h5py.File(outfile, 'a') as hfile:
        outputsize = len(hfile['Volume'])

# -------------------------------------------------------------------------------
# Loop over every row in `schedule` starting from row `rowid`
# -------------------------------------------------------------------------------
for index, row in itertools.islice(schedule.iterrows(), skiprows, None):
    ops.fvk = row.Gamma
    ops.morseWellWidth = log(2.0)*100.0/row.PercentStrain
    ops.constraint = row.AreaConstr
    logging.info('Setting gamma=%f, Morse-width=%f, area constraint=%f.', row.Gamma,
            ops.morseWellWidth, row.AreaConstr)
    if index == 0 and rowstep == 0:
        logging.info('Relaxing the system at zero temperature...')
        ops.brownCoeff = 0.0
        ops.viscosity = 0.0
        solver.solve()
        ops.saveInitialPosition()

    logging.info('Continuing with prescribed alpha=%f, beta=%f', row.Alpha, row.Beta)
    ops.viscosity = row.Alpha
    ops.brownCoeff = sqrt(2*row.Alpha / row.Beta)

    rowsteps = int(row.ViterMax)
    printstep = int(row.PrintStep)

    #We need to ensure that the `printstep` is compatible with `savefrequency`
    try:
        assert(savefrequency % printstep == 0)
    except AssertionError:
        logging.warning('`printstep` should be perfect divisor of %d.', savefrequency)
        printstep = savefrequency // 125
        logging.warning('Resetting `printstep` to %d', printstep)

    # Prepare containers to hold 3d shell structure
    polypartssize = savefrequency // printstep
    polyparts = collections.deque(maxlen=polypartssize)

    ops.updatePreviousX()

    starttime = time.time()
    logging.info('Starting inner loop. This will take a long time...')
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
            logging.info('Writing simulation state to file. Completed %d steps.', step)
            ops.timestep = step
            ops.vtkfilecount = vtkcount
            ops.writeSimulationState(statefile)
            with h5py.File(outfile, 'a') as hfile:
                for key, data in output.items():
                    hfile[key].resize((outputsize + savefrequency,))
                    hfile[key].write_direct(numpy.array(data), numpy.s_[:], numpy.s_[outputsize:])
                    data.clear()
            outputsize += savefrequency
            # Write the polydata parts to file
            with h5py.File(polydatafile, 'a') as hfile:
                hfile.attrs['TimeSteps'] = vtkcount + 1
                for subdict in polyparts:
                    for key, data in subdict.items():
                        hfile.create_dataset(key, data=data, dtype='float32')
            polyparts.clear()

        if i % printstep == 0 and printstep <= rowsteps:
            points, normals, cells = ops.polyDataParts()
            pointskey = '/T{0}/Points'.format(vtkcount)
            normalskey = '/T{0}/Normals'.format(vtkcount)
            cellskey = '/T{0}/Polygons'.format(vtkcount)
            polyparts.append({pointskey:points, normalskey:normals, cellskey:cells})
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

if len(polyparts) > 0:
    # Write the polydata parts to file
    with h5py.File(polydatafile, 'a') as hfile:
        hfile.attrs['TimeSteps'] = vtkcount + 1
        for subdict in polyparts:
            for key, data in subdict.items():
                hfile.create_dataset(key, data=data, dtype='float32')

endtime = time.time()
logging.info('Simulation completed. Last inner loop took %f seconds.', endtime - starttime)
