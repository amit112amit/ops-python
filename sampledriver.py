from opsmodel import OPSModel, Solver
from math import log

ops = OPSModel()

ops.initializeFromVTKFile('T7.vtk')
namesuffix = ops.getVTKFileSuffix()
step = ops.getTimeStep()
ops.setViscosity(0.0)
ops.setBrownCoeff(0.0)
ops.setMorseDistance(1.0)
ops.setMorseWellWidth(100.0/15 * log(2.0))
ops.setFVK(1000.0)

solver = Solver(ops)
solver.solve()

ops.printVTKFile('Relaxed.vtk')
