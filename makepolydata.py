import argparse

import h5py
import numpy
import pandas

from ops import cgalconvexhull
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter


def writepolydata(filename, points, normals, cells):
    polydata = vtkPolyData()
    vtkpoints = vtkPoints()
    vtkpoints.SetData(numpy_to_vtk(points))
    polydata.SetPoints(vtkpoints)
    polydata.GetPointData().SetNormals(numpy_to_vtk(normals))
    triangles = vtkCellArray()
    for i, j, k in cells:
        triangles.InsertNextCell(3)
        triangles.InsertCellPoint(i)
        triangles.InsertCellPoint(j)
        triangles.InsertCellPoint(k)
    polydata.SetPolys(triangles)
    wr = vtkPolyDataWriter()
    wr.SetFileName(filename)
    wr.SetInputData(polydata)
    wr.Write()


def redoconvexhull(filename, points, normals):
    polydata = vtkPolyData()
    vtkpoints = vtkPoints()
    vtkpoints.SetData(numpy_to_vtk(points))
    polydata.SetPoints(vtkpoints)
    polydata.GetPointData().SetNormals(numpy_to_vtk(normals))
    triangles = vtkCellArray()
    cells = cgalconvexhull(points)
    for i, j, k in cells:
        triangles.InsertNextCell(3)
        triangles.InsertCellPoint(i)
        triangles.InsertCellPoint(j)
        triangles.InsertCellPoint(k)
    polydata.SetPolys(triangles)
    wr = vtkPolyDataWriter()
    wr.SetFileName(filename)
    wr.SetInputData(polydata)
    wr.Write()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='makepolydata', description='Construct VTK Polydata from points, normals and cells in a HDF5 file.')
    parser.add_argument(
        '-i', '--input', help='The HDF5 file with points, normals and cells data.')
    parser.add_argument(
        '-o', '--output', help='Base name for the output VTK files.')
    args = parser.parse_args()

    with h5py.File(args.input, 'r') as hfile:
        for key, value in hfile.items():
            filename = ''.join([args.output, '-', key[1:], '.vtk'])
            points = numpy.empty(value['Points'].shape, dtype=numpy.float64)
            value['Points'].read_direct(points, numpy.s_[:,:], numpy.s_[:, :])
            normals = numpy.empty(value['Normals'].shape, dtype=numpy.float64)
            value['Normals'].read_direct(normals, numpy.s_[:,:], numpy.s_[:, :])
            cells = value['Polygons'][:]
            writepolydata(filename, points, normals, cells)
            filename = ''.join(['New', args.output, '-', key[1:], '.vtk'])
            redoconvexhull(filename, points, normals)
