#include <string.h>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef std::pair<K::Point_3, unsigned> Point_with_info;
typedef CGAL::Surface_mesh<Point_with_info> Surface_mesh;
typedef CGAL::First_of_pair_property_map<Point_with_info> Pmap;
typedef CGAL::Extreme_points_traits_adapter_3<Pmap, CGAL::Convex_hull_traits_3<K>> CHT;

typedef Eigen::Matrix3Xf Matrix3Xf;

extern "C"{
    void cgalconvexhullc(size_t N, float* coords, int triangles[]);
}