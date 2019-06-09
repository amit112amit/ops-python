#include "ConvexHull.h"

MatrixX3iR cgalconvexhull(Eigen::Ref<MatrixX3dR> positions)
{
    Surface_mesh sm;
    auto N = positions.rows();

    // Copy coordinates
    MatrixX3dR points(N, 3);
    points = positions;

    // Recenter and project points to unit sphere
    points = points.rowwise() - points.colwise().mean();
    points.rowwise().normalize();

    // Insert the projected points in a CGAL vertex_with_info vector
    std::vector<Point_with_info> verts;
    for (auto j = 0; j < N; ++j)
    {
        verts.push_back(std::make_pair(K::Point_3(points(j, 0), points(j, 1), points(j, 2)), j));
    }

    // Find the spherical convex hull
    CGAL::convex_hull_3(verts.begin(), verts.end(), sm, OPS::CH_traits_for_point_with_info());
    auto origids = sm.points();

    MatrixX3iR triangles(sm.number_of_faces(), 3);

    //Iterate over all triangles related to the center vertex
    auto i = 0;
    for (const auto &f : sm.faces())
    {
        auto h = sm.halfedge(f);
        auto f0 = origids[sm.target(h)].second;
        auto f1 = origids[sm.target(sm.next(h))].second;
        auto f2 = origids[sm.target(sm.next(sm.next(h)))].second;
        triangles.row(i) << f0, f1, f2;
        ++i;
    }

    return triangles;
}

MatrixX3iR vtkdelaunay(Eigen::Ref<MatrixX3dR> points)
{
    vtkNew<vtkPoints> pts;
    vtkNew<vtkDoubleArray> unitPData;
    vtkNew<vtkPolyData> unitSphere;
    vtkNew<vtkDataSetSurfaceFilter> dssf;
    vtkNew<vtkIdFilter> idf;
    vtkNew<vtkDelaunay3D> d3D;
    vtkNew<vtkIdList> pointIds;

    size_t N = points.rows();
    MatrixX3dR unitP(N, 3);
    unitP = points;

    // Make the origin as the center and project the points to a unit sphere
    unitP = unitP.rowwise() - unitP.colwise().mean();
    unitP.rowwise().normalize();

    // We need to do this hack to get vtkDelaunay3D to behave well
    unitP *= 100;
    unitP.array().round();
    unitP /= 100;

    void *voidPtr = (void *)unitP.data();
    unitPData->SetVoidArray(voidPtr, 3 * N, 1);
    unitPData->SetNumberOfComponents(3);
    pts->SetData(unitPData);

    unitSphere->SetPoints(pts);

    idf->SetIdsArrayName("PointIds");
    idf->PointIdsOn();
    idf->SetInputData(unitSphere);
    d3D->SetInputConnection(idf->GetOutputPort());
    dssf->SetInputConnection(d3D->GetOutputPort());
    dssf->Update();

    auto finalpd = dssf->GetOutput();
    auto cells = finalpd->GetPolys();
    auto origIds = finalpd->GetPointData()->GetArray("PointIds");

    cells->InitTraversal();

    MatrixX3iR triangles(cells->GetNumberOfCells(), 3);
    auto i = 0;
    while (cells->GetNextCell(pointIds))
    {
        for (int j = 0; j < 3; j++)
        {
            triangles(i, j) = (int)origIds->GetTuple1(pointIds->GetId(j));
        }
        ++i;
    }
    return triangles;
}

MatrixX3iR cgaldelaunay(Eigen::Ref<MatrixX3dR> positions)
{
    auto N = positions.rows();

    // Copy coordinates
    MatrixX3dR points(N, 3);
    points = positions;

    // Reset the center of the sphere to origin by translating
    points = points.rowwise() - points.colwise().mean();

    // Project points to unit sphere
    points.rowwise().normalize();

    // Insert the projected points in a CGAL vertex_with_info vector
    std::vector<std::pair<Point, unsigned>> verts;
    for (auto j = 0; j < N; ++j)
    {
        verts.push_back(std::make_pair(Point(points(j, 0), points(j, 1), points(j, 2)), j));
    }

    Delaunay dt(verts.begin(), verts.end());

    // Get the connectivity of outer faces
    auto inf = dt.infinite_vertex();
    std::vector<Delaunay::Cell_handle> cellhandles;
    dt.incident_cells(inf, std::back_inserter(cellhandles));
    MatrixX3iR triangles(cellhandles.size(), 3);
    auto i = 0;
    for (const auto &ch : cellhandles)
    {
        auto infindex = ch->index(inf);
        for(auto j = 0; j < 3; ++j){
            auto vi = Delaunay::vertex_triple_index(infindex, j);
            triangles(i, j) = ch->vertex(vi)->info();
        }
        ++i;
    }

    return triangles;
}