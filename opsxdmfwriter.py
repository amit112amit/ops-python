"""
Create XDMF files for visualization in Paraview.
"""
from ops import cgalconvexhull
import h5py
import numpy


def makeshellandwrite(infile, outfile, timesteps=None):
    """
    Reads particle positions from `infile` for specified time-steps, performs
    meshing and writes the mesh to a XDMF file.

    Parameters:
    -----------
    infile: name of the input HDF5 file containing TxN array where T is number
            of time-steps and N is number of 3 times the number of particles.
    outfile: name of the output XDMF file. The actual data will be stored in a
             HDF5 file with same name.

    timesteps: if `int` then it denotes the step size e.g. write every 1000
               time-steps to file
               if `tuple` then write only the specified time-steps to file
               So to write only the Nth time step we should give (N,)
               if string 
               if `None` write all time-steps to output file.
    """
    hdf5filename = outfile.split('.')[0] + '.h5'

    with h5py.File(infile, 'r') as inpfile:
        # Get shape of the input dataset
        numtimesteps = inpfile['Positions'].shape[0]
        rowlength = inpfile['Positions'].shape[1]
        numparticles = int(rowlength/3)
        row = numpy.empty(rowlength)

        # Decide what timesteps to write
        if isinstance(timesteps, tuple):
            stepstowrite = timesteps
        else:
            step = 1 if timesteps is None else timesteps
            stepstowrite = range(0, numtimesteps, step)

        with open(outfile, 'w') as xdmffile:
            # Write the XDMF header tags
            xdmffile.write(
                    '<?xml version="1.0"?>\n'
                    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
                    '<Xdmf Version="3.0">\n'
                    '  <Domain>\n'
                    '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n'
                    )
            with h5py.File(hdf5filename, 'w') as hdf5file:
                hdf5file.create_group('Geometry')
                hdf5file.create_group('Topology')
                for i in stepstowrite: 
                    # Read coordinates and make a mesh
                    inpfile['Positions'].read_direct(row, numpy.s_[i, :], numpy.s_[:])
                    coordinates = row.reshape(numparticles, 3)
                    triangles = cgalconvexhull(coordinates)
                    # Write the coordinates and mesh to HDF5 file
                    hdf5file.create_dataset('Geometry/{}'.format(i), data=coordinates,
                            compression="gzip", compression_opts=9)
                    hdf5file.create_dataset('Topology/{}'.format(i), data=triangles,
                            compression="gzip", compression_opts=9)
                    # Write a tag to XDMF file
                    xdmffile.write(
                            '      <Grid Name="Mesh" GridType="Uniform">\n'
                            '        <Topology NumberOfElements="{0}" TopologyType="Triangle" NodesPerElement="3">\n'
                            '          <DataItem Dimensions="{0} 3" NumberType="Int" Format="HDF">{1}:/Topology/{2}</DataItem>\n'
                            '        </Topology>\n'
                            '        <Geometry GeometryType="XYZ">\n'
                            '          <DataItem Dimensions="{3} 3" Format="HDF">{1}:/Geometry/{2}</DataItem>\n'
                            '        </Geometry>\n'
                            '        <Time Value="{2}" />\n'
                            '      </Grid>\n'.format(triangles.shape[0], hdf5filename, i, numparticles))
                xdmffile.write(
                        '    </Grid>\n'
                        '  </Domain>\n'
                        '</Xdmf>')


def writexdmf(inpfile, outfile):
    """
    Given a HDF5 file containing geometry, topology and point normals write a
    XDMF file.

    Parameter:
    ----------

    inpfile: a HDF5 file containing a structure as
            "/T{t}/[Points|Polygons|Normals]" where t is time-step
    outfile: string containing name of the output XDMF file.
    """
    with open(outfile, 'w') as xdmffile:
        # Write the XDMF header tags
        xdmffile.write(
                '<?xml version="1.0"?>\n'
                '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
                '<Xdmf Version="3.0">\n'
                '  <Domain>\n'
                '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n'
                )
        with h5py.File(inpfile, 'r') as hdf5file:
            numkeys = len(hdf5file.keys()) - 1    # One key has 'Asphericity'
            for key in range(numkeys):
                # key has the format "T<N>" where N is 0,1,2,3... without angular brackets.
                value = hdf5file['/T{0}'.format(key)]
                # Write a tag to XDMF file
                xdmffile.write(
                        '      <Grid Name="Mesh" GridType="Uniform">\n'
                        '        <Topology NumberOfElements="{0}" TopologyType="Triangle" NodesPerElement="3">\n'
                        '          <DataItem Dimensions="{0} 3" NumberType="Int" Format="HDF">{1}:/T{2}/Polygons</DataItem>\n'
                        '        </Topology>\n'
                        '        <Geometry GeometryType="XYZ">\n'
                        '          <DataItem Dimensions="{3} 3" Format="HDF">{1}:/T{2}/Points</DataItem>\n'
                        '        </Geometry>\n'
                        '        <Attribute Name="PointNormals" AttributeType="Vector" Center="Node">\n'
                        '          <DataItem Dimensions="{3} 3" Format="HDF">{1}:/T{2}/Normals</DataItem>\n'
                        '        </Attribute>\n'
                        '        <Time Value="{2}" />\n'
                        '      </Grid>\n'.format(value['Polygons'].shape[0], inpfile, key, value['Points'].shape[0]))
            # Close the outermost XML tags
            xdmffile.write(
                    '    </Grid>\n'
                    '  </Domain>\n'
                    '</Xdmf>')
