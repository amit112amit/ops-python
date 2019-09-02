#include "ExportConvexHull.h"

//! For calling from Julia/Python as a C library
//! The caller must ensure that there is enough memory allocated
//! for the triangles array. Minimum size is 60 + (N - 12)*6
//! integers
void cgalconvexhullc(size_t N, float *coordinates, int triangles[])
{
    Eigen::Map<Matrix3Xf> positions(coordinates, 3, N);
    Surface_mesh sm;

    // Copy coordinates
    Matrix3Xf points(3, N);
    points = positions;

    // Recenter and project points to unit sphere
    points = points.colwise() - points.rowwise().mean();
    points.colwise().normalize();

    // Insert the projected points in a CGAL vertex_with_info vector
    std::vector<Point_with_info> verts;
    for (auto j = 0; j < N; ++j)
    {
        verts.push_back(std::make_pair(K::Point_3(points(0, j), points(1, j),
                                                  points(2, j)),
                                       j));
    }

    // Find the spherical convex hull
    CGAL::convex_hull_3(verts.begin(), verts.end(), sm, CHT(Pmap()));
    auto origids = sm.points();

    //Iterate over all triangles related to the center vertex
    auto i = 0;
    std::cout<< "Number of triangles = " << sm.number_of_faces() << std::endl;
    for (const auto &f : sm.faces())
    {
        auto h = sm.halfedge(f);
        auto f0 = origids[sm.target(h)].second;
        auto f1 = origids[sm.target(sm.next(h))].second;
        auto f2 = origids[sm.target(sm.next(sm.next(h)))].second;
        //std::cout<< f0 <<", " << f1 << ", " << f2 << std::endl;
        triangles[i++] = int(f0);
        triangles[i++] = int(f1);
        triangles[i++] = int(f2);
    }

    return;
}