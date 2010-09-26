///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef CGAL_includes_and_types_H
#define CGAL_includes_and_types_H

#include "../../../../../mepp/Polyhedron/polyhedron.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for specific 2D polygon intersections in the plane
#include <CGAL/Gmpzf.h> // exact floating point type
#include <CGAL/Gmpz.h> // arbitrary precision integer based on the GNU Multiple Precision Arithmetic Library
#include <CGAL/Gmpq.h> // arbitrary precision rational number based on the GNU Multiple Precision Arithmetic Library.
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Filtered_extended_homogeneous.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for exact computation without spending too much time
#include  <CGAL/Cartesian.h>
#include  <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include  <CGAL/Quotient.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <CGAL/Nef_polyhedron_2.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for Delaunay
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/intersections.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for 2D arrangement computation
#include  <CGAL/Arr_segment_traits_2.h> // line segments are used to build 2D arrangements
#include  <CGAL/Arrangement_2.h>
#include <CGAL/Arr_accessor.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Nef_polyhedron_2 methods
///////////////////////////////////////////////////////////////////////////
//typedef double Number_type;
//typedef CGAL::Simple_cartesian<Number_type> Kernel_N;
//////////////////////////////////////////////////////////////////////////
//typedef int FT;
//typedef CGAL::Gmpz FT; // quicker, but not enough precion to give correct results!!
typedef CGAL::Gmpq FT; // for Cartesian, we must use Gmpq
//typedef CGAL::Gmpzf FT; // exact floating point type
//////////////////////////////////////////////////////////////////////////
typedef CGAL::Extended_cartesian< FT > Extended_kernel;
//typedef CGAL::Extended_homogeneous< RT > Extended_kernel;
//typedef CGAL::Filtered_extended_homogeneous< RT > Extended_kernel; // I have some bugs with this one

typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron; // the kernel must be a model of the concept |ExtendedKernelTraits_2| (otherwise it will generate an error!!)
typedef Nef_polyhedron::Point Point;
typedef Nef_polyhedron::Line  Line;

typedef Extended_kernel::Vector_3 Vector3Dex;
typedef Nef_polyhedron::Explorer Explorer; // read-only exploration!!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Delaunay methods
typedef  CGAL::Quotient<CGAL::MP_Float>             Number_type_delaunay; // it works!!
//typedef  CGAL::Gmpq                               Number_type;

typedef CGAL::Lazy_exact_nt<Number_type_delaunay>       Lazy_exact_nt_type_delaunay;
typedef CGAL::Cartesian<Lazy_exact_nt_type_delaunay>    k_delaunay;

//typedef CGAL::Cartesian<long double> number_type_delaunay;
//typedef CGAL::Cartesian<CGAL::Gmpq> number_type_delaunay; // long double sinon avec double ça plante la contruction de la triangulation
typedef CGAL::Point_2< k_delaunay > Point_dt;
typedef CGAL::Segment_2< k_delaunay > Segment_dt;
typedef CGAL::Delaunay_triangulation_2< k_delaunay > Delaunay;
typedef Delaunay::Triangle Triangle_dt;
typedef Delaunay::Vertex_handle vh;
typedef Delaunay::Finite_vertices_iterator v_it;
typedef Delaunay::Vertex_circulator v_cir;
typedef Delaunay::Finite_faces_iterator f_it;

//#define USE_SEGMENT // to use the part dedicated to segments
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ARRANGEMENTS
/*
The Arr_segment_traits_2<Kernel> class is very efficient for maintaining arrangements of a large number
of intersecting line segments, especially if it is instantiated with the appropriate geometric kernel.
Using Cartesian<Gmpq> as the kernel type is generally a good choice; the coordinates of the segment
endpoints are represented as exact rational numbers, and this ensures the robustness and correctness
of any computation. However, if the Gmp library is not installed, it is possible to use the Quotient<MP_Float>
number-type, provided by the support library of Cgal, which is somewhat less efficient.

*/
typedef  CGAL::Quotient<CGAL::MP_Float>         Number_type; // it works!!
//typedef  CGAL::Gmpq                             Number_type;

typedef CGAL::Lazy_exact_nt<Number_type>        Lazy_exact_nt_type;
typedef CGAL::Cartesian<Lazy_exact_nt_type>     Kernel_ar;


typedef  CGAL::Arr_segment_traits_2<Kernel_ar>  Traits_2;
typedef  Traits_2::Point_2                      Point_ar;
typedef  Traits_2::X_monotone_curve_2           Segment_ar;
typedef  Traits_2::Curve_2                      Curve_ar;
typedef  CGAL::Arrangement_2<Traits_2>          Arrangement;
typedef  Arrangement::Vertex                    Vertex_ar;
typedef  Arrangement::Halfedge                  Halfedge_ar;
typedef  Arrangement::Face                      Face_ar;        // represents a 2-dimensional cell in the subdivision
typedef  Arrangement::Vertex_iterator           Vertex_ar_it;
typedef  Arrangement::Edge_iterator             Edge_ar_it;
typedef  Arrangement::Face_iterator             Face_ar_it;
typedef  Arrangement::Face_const_iterator       Face_ar_const_it;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// general
#include <assert.h>

#endif
