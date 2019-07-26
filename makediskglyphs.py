from vtkmodules.vtkIOLegacy import vtkPolyDataReader, vtkPolyDataWriter
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.numpy_interface import dataset_adapter as dsa
import numpy

rd = vtkPolyDataReader()
rd.SetFileName('T7.vtk')
rd.Update()

pd = rd.GetOutput()
pdw = dsa.WrapDataObject(pd)

# Get the number of neighbors for each point
pd.BuildLinks()
valence = []
neighbors = vtkIdList()

for i in range(pd.GetNumberOfPoints()):
    pd.GetPointCells(i, neighbors)
    valence.append(neighbors.GetNumberOfIds())

# Add valence and position vectors as pointdata
pdw.PointData.append(numpy.array(valence), 'Valence')
pdw.PointData.append(pdw.Points, 'PositionVectors')

wr = vtkPolyDataWriter()
wr.SetFileName('TryThis.vtk')
wr.SetInputData(pd)
wr.Write()

