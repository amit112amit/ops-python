#include <string.h>
#include <array>
#include <Eigen/Dense>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkIdFilter.h>
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include "CGAL_ConvexHull_Traits.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
typedef std::pair<K::Point_3, unsigned> Point_with_info;
typedef CGAL::Surface_mesh<Point_with_info> Surface_mesh;

typedef Eigen::Matrix<double_t, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3dR;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3iR;

MatrixX3iR cgalconvexhull(Eigen::Ref<MatrixX3dR> positions);
MatrixX3iR vtkdelaunay(Eigen::Ref<MatrixX3dR> points);
MatrixX3iR cgaldelaunay(Eigen::Ref<MatrixX3dR> positions);