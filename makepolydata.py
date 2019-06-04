import argparse

import pandas
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter

parser = argparse.ArgumentParser(prog='makepolydata', description='Construct VTK Polydata from points, normals and cells in a HDF5 file.')
parser.add_argument('-i', '--input', help='The HDF5 file with points, normals and cells data.')
parser.add_argument('-o', '--output', help='Base name for the output VTK files.')
args = parser.parse_args()

with pandas.HDFStore(args.input, 'r') as store:
    for path, subgroups, subkeys in store.walk():
        if len(subkeys) is 3:
            polydata = vtkPolyData()

            points = store['{}/Points'.format(path)].to_numpy()
            vtkpoints = vtkPoints()
            vtkpoints.SetData(numpy_to_vtk(points))
            polydata.SetPoints(vtkpoints)

            normals = store['{}/Normals'.format(path)].to_numpy()
            polydata.GetPointData().SetNormals(numpy_to_vtk(normals))

            cells = store['{}/Polygons'.format(path)].to_numpy()
            triangles = vtkCellArray()
            for i, j, k in cells:
                triangles.InsertNextCell(3)
                triangles.InsertCellPoint(i)
                triangles.InsertCellPoint(j)
                triangles.InsertCellPoint(k)
            polydata.SetPolys(triangles)

            vtkfilename = '{}-{}.vtk'.format(args.output, path[2:])
            wr = vtkPolyDataWriter()
            wr.SetFileName(vtkfilename)
            wr.SetInputData(polydata)
            wr.Write()
