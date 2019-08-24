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
        verts.push_back(std::make_pair(K::Point_3(points(j, 0), points(j, 1),
                                                  points(j, 2)),
                                       j));
    }

    // Find the spherical convex hull
    CGAL::convex_hull_3(verts.begin(), verts.end(), sm, CHT(Pmap()));
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
    auto pts = vtkSmartPointer<vtkPoints>::New();
    auto unitPData = vtkSmartPointer<vtkDoubleArray>::New();
    auto unitSphere = vtkSmartPointer<vtkPolyData>::New();
    auto dssf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    auto idf = vtkSmartPointer<vtkIdFilter>::New();
    auto d3D = vtkSmartPointer<vtkDelaunay3D>::New();
    auto pointIds = vtkSmartPointer<vtkIdList>::New();

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
        verts.push_back(std::make_pair(Point(points(j, 0), points(j, 1),
                                             points(j, 2)),
                                       j));
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
        for (auto j = 0; j < 3; ++j)
        {
            auto vi = Delaunay::vertex_triple_index(infindex, j);
            triangles(i, j) = ch->vertex(vi)->info();
        }
        ++i;
    }

    return triangles;
}

void shellstats(Eigen::Ref<MatrixdR> data, Eigen::Ref<VectorXi> neighbors,
                Eigen::Ref<VectorXd> asphericity, Eigen::Ref<VectorXd> msd,
                Eigen::Ref<VectorXd> radius, Eigen::Ref<VectorXd> volume)
{
    size_t N = int(data.cols() / 3);
    Eigen::Map<MatrixX3dR> initial(&data(0, 0), N, 3);
    for (size_t i = 0; i < data.rows(); ++i)
    {
        Eigen::Map<MatrixX3dR> current(&data(i, 0), N, 3);
        Eigen::VectorXd R(N);
        R = current.rowwise().norm();
        auto R0 = R.mean();
        radius(i) = R0;
        asphericity(i) = ((R.array() - R0).square()).sum();
        asphericity(i) /= (N * R0 * R0);
        auto triangles = cgalconvexhull(current);
        volume(i) = 0.0;
        for (size_t j = 0; j < triangles.rows(); ++j)
        {
            auto a = triangles(j, 0);
            auto b = triangles(j, 1);
            auto c = triangles(j, 2);
            volume(i) +=
                0.166666667 * (current.row(a).dot(
                                  current.row(b).cross(current.row(c))));
        }

        msd(i) = 0;
        for (size_t j = 0; j < N; ++j)
        {
            Vector3d xi, xj, diff, xi_diff, xj_diff;
            Vector3d xi0, xj0, xi1, xj1, ni0, nj0;

            auto nn = neighbors(j);

            xi0 = initial.row(j);
            xi1 = current.row(j);

            xj0 = initial.row(nn);
            xj1 = current.row(nn);

            xi_diff = (xi1 - xi0);
            xj_diff = (xj1 - xj0);

            diff = xi_diff - xj_diff;
            msd(i) += diff.dot(diff);
        }
        msd(i) /= N;
    }
}