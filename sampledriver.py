from opsmodel import OPSModel, Solver
from math import log

ops = OPSModel()

ops.initializeFromVTKFile('T7.vtk')
namesuffix = ops.vtkfilecount
step = ops.timestep
ops.viscosity = 0.0
ops.brownCoeff = 0.0
ops.morseDistance = 1.0
ops.morseWellWidth = 100.0/15 * log(2.0)
ops.fvk = 1000.0

solver = Solver(ops)
solver.maxiter = 1000
solver.solve()

ops.printVTKFile('Relaxed.vtk')

coords, normals, triangles = ops.polyDataParts()

print('Total energy = ', ops.totalEnergy)
print('Volume = ', ops.volume)
