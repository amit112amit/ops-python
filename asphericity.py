import os
from math import log, sqrt

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
outfile = 'T7_OPS-' + str(taskid) + '.dat'
percentstrain = 15

# Create the model
ops = Model()

try:
    ops.initializeFromVTKFile(vtkfile)
except FileNotFoundError:
    print('Unable to locate', vtkfile, 'in the working directory!')

# Read the simulation schedule
with open('schedule-{}.dat'.format(taskid)) as inpfile:
    allgamma = inpfile.read().split('\n')[:-1]
    schedule = [float(g) for g in allgamma]

vtkcount = ops.vtkfilecount

# Create solver object
solver = Solver(ops)

# Create output data object
output = {'Gamma': [],
          'Asphericity': [],
          'Radius': [],
          'Volume': [],
          'Area': [],
          'MorseEnergy': [],
          'NormalityEn': [],
          'CircularityEn': [],
          'TotalOPSEnergy': [],
          'RMSAngleDeficit': []
          }

# -------------------------------------------------------------------------------
# Start solving the zero temperature problem
# -------------------------------------------------------------------------------
ops.lagrangeCoeff = 0.0
ops.penaltyCoeff = 0.0
ops.brownCoeff = 0.0
ops.viscosity = 0.0
ops.morseWellWidth = log(2.0)*100.0/percentstrain

ops.updatePreviousX()

for index, gamma in enumerate(schedule):
    ops.fvk = gamma
    solver.solve()
    ops.updateTriangles()

    # -----------------------------------------------------------------------
    # Collect statistics for output
    # -----------------------------------------------------------------------
    output['Gamma'].append(gamma)
    output['Asphericity'].append(ops.asphericity)
    output['Radius'].append(ops.averageRadius)
    output['Volume'].append(ops.volume)
    output['Area'].append(ops.area)
    output['MorseEnergy'].append(ops.morseEnergy)
    output['NormalityEn'].append(ops.normalityEnergy)
    output['CircularityEn'].append(ops.circularityEnergy)
    output['TotalOPSEnergy'].append(ops.totalEnergy)
    output['RMSAngleDeficit'].append(ops.rmsAngleDeficit)

    # -----------------------------------------------------------------------
    # Save visualization
    # -----------------------------------------------------------------------
    ops.vtkfilecount = vtkcount
    points, normals, cells = ops.polyDataParts()
    pointskey = '/T{0}/Points'.format(vtkcount)
    normalskey = '/T{0}/Normals'.format(vtkcount)
    cellskey = '/T{0}/Polygons'.format(vtkcount)
    with h5py.File(polydatafile, 'a') as hfile:
        hfile.attrs['TimeSteps'] = vtkcount + 1
        hfile.create_dataset(pointskey, data=points, dtype='float32')
        hfile.create_dataset(normalskey, data=normals, dtype='float32')
        hfile.create_dataset(cellskey, data=cells, dtype='uint32')
    vtkcount += 1

# -----------------------------------------------------------------------
# Write output data to file
# -----------------------------------------------------------------------
outputdata = pandas.DataFrame(output)
outputdata.to_csv(outfile, sep='\t', index=False)
