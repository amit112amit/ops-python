from numpy import dot, cross
from numpy.linalg import norm
from math import log

import sys
sys.path.append('/home/amit/WorkSpace/UCLA/ops-python/build')
from ops import Model, Solver

opsmodel = Model()

opsmodel.initializeFromVTKFile('T7.vtk')
namesuffix = opsmodel.vtkfilecount
step = opsmodel.timestep
opsmodel.viscosity = 0.0
opsmodel.brownCoeff = 0.0
opsmodel.lagrangeCoeff = 0.0
opsmodel.penaltyCoeff = 0.0
opsmodel.morseDistance = 1.0
opsmodel.morseWellWidth = 100.0/15 * log(2.0)
opsmodel.fvk = 1000.0

solver = Solver(opsmodel)
solver.maxiter = 1000
solver.solve()

opsmodel.updateTriangles()
opsmodel.printVTKFile('Relaxed.vtk')

coords, normals, triangles = opsmodel.polyDataParts()

print('Total energy = ', opsmodel.totalEnergy)
print('Volume = ', opsmodel.volume)