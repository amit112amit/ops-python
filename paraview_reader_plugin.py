""" Reference:
https://gitlab.kitware.com/paraview/paraview/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py
"""

from paraview.util.vtkAlgorithm import *
import numpy

#------------------------------------------------------------------------------
# A reader example.
#------------------------------------------------------------------------------
def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

# To add a reader, we can use the following decorators
#   @smproxy.source(name="PythonCSVReader", label="Python-based CSV Reader")
#   @smhint.xml("""<ReaderFactory extensions="csv" file_description="Numpy CSV files" />""")
# or directly use the "@reader" decorator.
@smproxy.reader(name="PythonH5Reader", label="Python-based HDF5 Reader",
                extensions="h5", file_description="HDF5 files")
class PythonH5Reader(VTKPythonAlgorithmBase):
    """A reader that reads a HDF5 file. The data is treated as a temporal dataset."""
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkPolyData')
        self._filename = None
        self._hfile = None
        self._timesteps = None

        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))

    def _get_raw_data(self, requested_time=None):
        if self._hfile is not None:
            if requested_time is not None:
                path = 'T' + str(requested_time)
                return self._hfile[path]
            return self._hfile['T0']

        if self._filename is None:
            # Note, exceptions are totally fine!
            raise RuntimeError("No filename specified")

        self._hfile = File(self._filename, 'r')
        self._timesteps = list(range(0, self._hfile.attrs['TimeSteps']))
        self._arrayselection.AddArray('Normals')
        return self._get_raw_data(requested_time)

    def _get_timesteps(self):
        self._get_raw_data()
        return self._timesteps if self._timesteps is not None else None

    def _get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self._get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

    def _get_array_selection(self):
        return self._arrayselection

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="csv", file_description="Numpy CSV files")
    def SetFileName(self, name):
        """Specify filename for the file to read."""
        if self._filename != name:
            self._filename = name
            self._ndata = None
            self._timesteps = None
            self.Modified()

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._get_timesteps()

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())

        timesteps = self._get_timesteps()
        if timesteps is not None:
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
        from vtkmodules.numpy_interface import dataset_adapter as dsa
        from vtkmodules.util.numpy_support import numpy_to_vtk

        data_time = self._get_update_time(outInfoVec.GetInformationObject(0))
        raw_data = self._get_raw_data(data_time)

        points = numpy.empty(raw_data['Points'].shape, dtype=numpy.float64)
        normals = numpy.empty(raw_data['Normals'].shape, dtype=numpy.float64)
        raw_data['Points'].read_direct(points, numpy.s_[:, :], numpy.s_[:, :])
        raw_data['Normals'].read_direct(normals, numpy.s_[:, :], numpy.s_[:, :])
        cells = raw_data['Polygons'][:]
        triangles = vtkCellArray()
        for i, j, k in cells:
            triangles.InsertNextCell(3)
            triangles.InsertCellPoint(i)
            triangles.InsertCellPoint(j)
            triangles.InsertCellPoint(k)

        output = dsa.WrapDataObject(vtkPolyData.GetData(outInfoVec))
        output.Points = points
        output.VTKObject.SetPolys(triangles)
        output.VTKObject.GetPointData().SetNormals(numpy_to_vtk(normals))

        if data_time is not None:
            output.GetInformation().Set(output.DATA_TIME_STEP(), data_time)
        return 1
