import itertools
import os
from math import log, sqrt

import pandas

from opsmodel import OPSModel, Solver

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
statefile = 'SimulationState-' + str(taskid) + '.dat'

ops = OPSModel()

if os.path.isfile(statefile):
    ops.restoreSavedState(statefile)
elif os.path.isfile(vtkfile):
    ops.initializeFromVTKFile(vtkfile)
else:
    print('FileNotFound :', vtkfile, 'or', statefile)
    quit()

solver = Solver(ops)

savefrequency = 500
vtkcount = ops.vtkfilecount
step = ops.timestep

schedule = pandas.read_csv('schedule-{}.dat'.format(taskid), sep=' ')
skiprows, totalsteps = next(itertools.dropwhile(
    lambda x: x[1] <= step, enumerate(schedule.ViterMax.cumsum())))
rowstep = schedule.ViterMax[skiprows] - (totalsteps - step)

outputrows = []

# -------------------------------------------------------------------------------
# Loop over every row in `schedule` starting from row `rowid`
# -------------------------------------------------------------------------------
for index, row in itertools.islice(schedule.iterrows(), skiprows, None):
    ops.fvk = row.Gamma
    ops.morseWellWidth = log(2.0)*100.0/row.PercentStrain
    ops.constraint = row.AreaConstr
    if index is 0 and rowstep is 0:
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
        # Save visualization or simulation state
        # -----------------------------------------------------------------------
        if (step % savefrequency) == savefrequency - 1:
            ops.timestep = step
            ops.vtkfilecount = vtkcount
            ops.writeSimulationState(statefile)

        if i % printstep is 0 and printstep <= rowsteps:
            vtkfile = 'T7-relaxed-' + str(vtkcount) + '.vtk'
            ops.printVTKFile(vtkfile)
            vtkcount += 1

        # -----------------------------------------------------------------------
        # Collect statistics for output
        # -----------------------------------------------------------------------
        radialmsd, tangentialmsd = ops.bothMSD()
        outputrows.append({'Volume': ops.volume, 'MSD': radialmsd, 'MSDT': tangentialmsd,
                           'RMSAngleDeficit': ops.rmsAngleDeficit})

        ops.updatePreviousX()
        step += 1

    # Reset rowstep for next row
    rowstep = 0

pandas.DataFrame(outputrows).to_excel('Output.xlsx')
