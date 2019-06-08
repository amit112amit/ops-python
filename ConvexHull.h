#if !defined(__CONVEXHULL_H__)
#define __CONVEXHULL_H__

#include <vector>
#include <utility>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_traits_2.h>

namespace OPS{

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef std::pair<K::Point_3, unsigned> Point_with_info;
typedef CGAL::Surface_mesh<Point_with_info> Surface_mesh;

template <class F>
struct Forward_functor
    : public F
{
    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q) const
    {
        return static_cast<const F *>(this)->operator()(p.first, q.first);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
    {
        return static_cast<const F *>(this)->operator()(p.first, q.first, r.first);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s) const
    {
        return static_cast<const F *>(this)->operator()(p.first, q.first, r.first, s.first);
    }
};

template <class F>
struct Forward_functor_xy
    : public F
{
    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q) const
    {
        K::Point_2 a(p.first.x(), p.first.y());
        K::Point_2 b(q.first.x(), q.first.y());
        return static_cast<const F *>(this)->operator()(a, b);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
    {
        K::Point_2 a(p.first.x(), p.first.y());
        K::Point_2 b(q.first.x(), q.first.y());
        K::Point_2 c(r.first.x(), r.first.y());
        return static_cast<const F *>(this)->operator()(a, b, c);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s) const
    {
        K::Point_2 a(p.first.x(), p.first.y());
        K::Point_2 b(q.first.x(), q.first.y());
        K::Point_2 c(r.first.x(), r.first.y());
        K::Point_2 d(s.first.x(), s.first.y());
        return static_cast<const F *>(this)->operator()(a, b, c, d);
    }
};

template <class F>
struct Forward_functor_yz
    : public F
{
    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q) const
    {
        K::Point_2 a(p.first.y(), p.first.z());
        K::Point_2 b(q.first.y(), q.first.z());
        return static_cast<const F *>(this)->operator()(a, b);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
    {
        K::Point_2 a(p.first.y(), p.first.z());
        K::Point_2 b(q.first.y(), q.first.z());
        K::Point_2 c(r.first.y(), r.first.z());
        return static_cast<const F *>(this)->operator()(a, b, c);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s) const
    {
        K::Point_2 a(p.first.y(), p.first.z());
        K::Point_2 b(q.first.y(), q.first.z());
        K::Point_2 c(r.first.y(), r.first.z());
        K::Point_2 d(s.first.y(), s.first.z());
        return static_cast<const F *>(this)->operator()(a, b, c, d);
    }
};

template <class F>
struct Forward_functor_xz
    : public F
{
    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q) const
    {
        K::Point_2 a(p.first.x(), p.first.z());
        K::Point_2 b(q.first.x(), q.first.z());
        return static_cast<const F *>(this)->operator()(a, b);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
    {
        K::Point_2 a(p.first.x(), p.first.z());
        K::Point_2 b(q.first.x(), q.first.z());
        K::Point_2 c(r.first.x(), r.first.z());
        return static_cast<const F *>(this)->operator()(a, b, c);
    }

    template <class Point_3>
    bool operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s) const
    {
        K::Point_2 a(p.first.x(), p.first.z());
        K::Point_2 b(q.first.x(), q.first.z());
        K::Point_2 c(r.first.x(), r.first.z());
        K::Point_2 d(s.first.x(), s.first.z());
        return static_cast<const F *>(this)->operator()(a, b, c, d);
    }
};

class Point_with_info_triple
{
protected:
    typedef K::FT FT;
    typedef K::Vector_3 Vector_3;
    typedef Point_with_info Point_3;

    Point_3 p_, q_, r_;

public:
    Point_with_info_triple() {}

    Point_with_info_triple(const Point_3 &p, const Point_3 &q, const Point_3 &r)
        : p_(p), q_(q), r_(r)
    {
    }

    const K::Point_3 &p() const { return p_.first; }
    const K::Point_3 &q() const { return q_.first; }
    const K::Point_3 &r() const { return r_.first; }
};

class Point_with_info_triple_has_on_positive_side_3
{
public:
    typedef Point_with_info Point_3;
    typedef Point_with_info_triple Plane_3;
    bool
    operator()(const Plane_3 &pl, const Point_3 &p) const
    {
        typename K::Orientation_3 o;
        return (o(pl.p(), pl.q(), pl.r(), p.first) == CGAL::POSITIVE);
    }

    typedef bool result_type;
};

class Point_with_info_triple_construct_orthogonal_vector_3
{
public:
    typedef K::Vector_3 Vector_3;
    typedef Point_with_info_triple Plane_3;

    Vector_3 operator()(const Plane_3 &plane) const
    {
        K::Construct_orthogonal_vector_3 construct_orthogonal_vector_3;
        return construct_orthogonal_vector_3(plane.p(), plane.q(), plane.r());
    }
};

class Point_with_info_triple_oriented_side_3
{
public:
    typedef Point_with_info Point_3;
    typedef Point_with_info_triple Plane_3;
    typedef CGAL::Oriented_side result_type;

    result_type
    operator()(const Plane_3 &pl, const Point_3 &p) const
    {
        K::Orientation_3 o;
        CGAL::Orientation ori = o(pl.p(), pl.q(), pl.r(), p.first);
        if (ori > 0)
            return CGAL::ON_POSITIVE_SIDE;
        if (ori < 0)
            return CGAL::ON_NEGATIVE_SIDE;
        return CGAL::ON_ORIENTED_BOUNDARY;
    }
};

class Point_with_info_triple_less_signed_distance_to_plane_3
{
public:
    typedef Point_with_info Point_3;
    typedef Point_with_info_triple Plane_3;

    typedef bool result_type;

    bool
    operator()(const Plane_3 &h, const Point_3 &p, const Point_3 &q) const
    {
        const K::Point_3 &hp = h.p();
        const K::Point_3 &hq = h.q();
        const K::Point_3 &hr = h.r();
        K::Less_signed_distance_to_plane_3 less_signed_distance_to_plane_3;
        return less_signed_distance_to_plane_3(hp, hq, hr, p.first, q.first);
    }
};

class CH_traits_for_point_with_info
{
public:
    typedef K::Segment_3 Segment_3;
    typedef K::Triangle_3 Triangle_3;
    typedef K::Vector_3 Vector_3;
    typedef CGAL::Convex_hull_traits_3<K, Surface_mesh> Base;
    typedef Point_with_info Point_3;
    typedef Point_with_info_triple Plane_3;
    typedef Surface_mesh Polyhedron_3;

    typedef Base::Construct_segment_3 Construct_segment_3;
    typedef Base::Construct_ray_3 Construct_ray_3;
    typedef Base::Construct_triangle_3 Construct_triangle_3;
    typedef Base::Construct_centroid_3 Construct_centroid_3;

    class Construct_plane_3
    {
    public:
        Plane_3 operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
        {
            return Plane_3(p, q, r);
        }
    };

    typedef Point_with_info_triple_construct_orthogonal_vector_3 Construct_orthogonal_vector_3;

    typedef Forward_functor<K::Equal_3> Equal_3;
    typedef Forward_functor<K::Orientation_3> Orientation_3;
    typedef Forward_functor<K::Collinear_3> Collinear_3;
    typedef Forward_functor<K::Coplanar_3> Coplanar_3;
    typedef Forward_functor<K::Less_distance_to_point_3> Less_distance_to_point_3;

    typedef Point_with_info_triple_has_on_positive_side_3 Has_on_positive_side_3;
    typedef Point_with_info_triple_less_signed_distance_to_plane_3
        Less_signed_distance_to_plane_3;

    // required for degenerate case of all points coplanar
    class Traits_xy_3
    {
    public:
        typedef Point_with_info Point_2;
        typedef CGAL::convex_hull_traits_2<K> Base;
        typedef Forward_functor_xy<Base::Less_xy_2> Less_xy_2;
        typedef Forward_functor_xy<Base::Less_yx_2> Less_yx_2;
        typedef Forward_functor_xy<Base::Left_turn_2> Left_turn_2;
        typedef Forward_functor_xy<Base::Equal_2> Equal_2;
        struct Orientation_2
        {
            CGAL::Orientation
            operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const
            {
                K::Point_2 a(p.first.x(), p.first.y());
                K::Point_2 b(q.first.x(), q.first.y());
                K::Point_2 c(r.first.x(), r.first.y());
                return Base::Orientation_2()(a, b, c);
            }
        };

        Less_xy_2 less_xy_2_object() const
        {
            return Less_xy_2();
        }
        Less_yx_2 less_yx_2_object() const
        {
            return Less_yx_2();
        }
        Left_turn_2 left_turn_2_object() const
        {
            return Left_turn_2();
        }
        Equal_2 equal_2_object() const
        {
            return Equal_2();
        }
        Orientation_2 orientation_2_object() const
        {
            return Orientation_2();
        }
    };

    class Traits_yz_3
    {
    public:
        typedef Point_with_info Point_2;
        typedef CGAL::convex_hull_traits_2<K> Base;
        typedef Forward_functor_yz<Base::Less_xy_2> Less_xy_2;
        typedef Forward_functor_yz<Base::Less_yx_2> Less_yx_2;
        typedef Forward_functor_yz<Base::Left_turn_2> Left_turn_2;
        typedef Forward_functor_yz<Base::Equal_2> Equal_2;
        struct Orientation_2
        {
            CGAL::Orientation
            operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const
            {
                K::Point_2 a(p.first.y(), p.first.z());
                K::Point_2 b(q.first.y(), q.first.z());
                K::Point_2 c(r.first.y(), r.first.z());
                return Base::Orientation_2()(a, b, c);
            }
        };

        Less_xy_2 less_xy_2_object() const
        {
            return Less_xy_2();
        }
        Less_yx_2 less_yx_2_object() const
        {
            return Less_yx_2();
        }
        Left_turn_2 left_turn_2_object() const
        {
            return Left_turn_2();
        }
        Equal_2 equal_2_object() const
        {
            return Equal_2();
        }
        Orientation_2 orientation_2_object() const
        {
            return Orientation_2();
        }
    };

    class Traits_xz_3
    {
    public:
        typedef Point_with_info Point_2;
        typedef CGAL::convex_hull_traits_2<K> Base;
        typedef Forward_functor_xz<Base::Less_xy_2> Less_xy_2;
        typedef Forward_functor_xz<Base::Less_yx_2> Less_yx_2;
        typedef Forward_functor_xz<Base::Left_turn_2> Left_turn_2;
        typedef Forward_functor_xz<Base::Equal_2> Equal_2;
        struct Orientation_2
        {
            CGAL::Orientation
            operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const
            {
                K::Point_2 a(p.first.x(), p.first.z());
                K::Point_2 b(q.first.x(), q.first.z());
                K::Point_2 c(r.first.x(), r.first.z());
                return Base::Orientation_2()(a, b, c);
            }
        };

        Less_xy_2 less_xy_2_object() const
        {
            return Less_xy_2();
        }
        Less_yx_2 less_yx_2_object() const
        {
            return Less_yx_2();
        }
        Left_turn_2 left_turn_2_object() const
        {
            return Left_turn_2();
        }
        Equal_2 equal_2_object() const
        {
            return Equal_2();
        }
        Orientation_2 orientation_2_object() const
        {
            return Orientation_2();
        }
    };

    Traits_xy_3 construct_traits_xy_3_object() const
    {
        return Traits_xy_3();
    }
    Traits_yz_3 construct_traits_yz_3_object() const
    {
        return Traits_yz_3();
    }
    Traits_xz_3 construct_traits_xz_3_object() const
    {
        return Traits_xz_3();
    }

    typedef Base::Construct_vector_3 Construct_vector_3;

    // for postcondition checking
    typedef K::Ray_3 Ray_3;
    typedef Forward_functor<K::Has_on_3> Has_on_3;
    typedef Point_with_info_triple_oriented_side_3 Oriented_side_3;
    typedef Forward_functor<K::Do_intersect_3> Do_intersect_3;

    Construct_segment_3
    construct_segment_3_object() const
    {
        return Construct_segment_3();
    }

    Construct_ray_3
    construct_ray_3_object() const
    {
        return Construct_ray_3();
    }

    Construct_plane_3
    construct_plane_3_object() const
    {
        return Construct_plane_3();
    }

    Construct_triangle_3
    construct_triangle_3_object() const
    {
        return Construct_triangle_3();
    }

    Construct_centroid_3
    construct_centroid_3_object() const
    {
        return Construct_centroid_3();
    }

    Construct_orthogonal_vector_3
    construct_orthogonal_vector_3_object() const
    {
        return Construct_orthogonal_vector_3();
    }

    Collinear_3
    collinear_3_object() const
    {
        return Collinear_3();
    }

    Coplanar_3
    coplanar_3_object() const
    {
        return Coplanar_3();
    }

    Has_on_3
    has_on_3_object() const
    {
        return Has_on_3();
    }

    Less_distance_to_point_3
    less_distance_to_point_3_object() const
    {
        return Less_distance_to_point_3();
    }

    Has_on_positive_side_3
    has_on_positive_side_3_object() const
    {
        return Has_on_positive_side_3();
    }

    Oriented_side_3
    oriented_side_3_object() const
    {
        return Oriented_side_3();
    }

    Equal_3
    equal_3_object() const
    {
        return Equal_3();
    }

    Do_intersect_3
    do_intersect_3_object() const
    {
        return Do_intersect_3();
    }

    Less_signed_distance_to_plane_3
    less_signed_distance_to_plane_3_object() const
    {
        return Less_signed_distance_to_plane_3();
    }

    Orientation_3
    orientation_3_object() const
    {
        return Orientation_3();
    }

    Construct_vector_3
    construct_vector_3_object() const
    {
        return Construct_vector_3();
    }
};

}
#endif //__CONVEXHULL_H__