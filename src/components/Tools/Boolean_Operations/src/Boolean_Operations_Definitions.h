#ifndef BOOLEAN_OPERATIONS_DEFINITIONS_H
#define BOOLEAN_OPERATIONS_DEFINITIONS_H

#include "../../../../mepp/Polyhedron/polyhedron.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_exact_nt.h>

#include <mepp_config.h>
#define BOOLEAN_OPERATIONS_DEBUG

enum Bool_Op {UNION, INTER, MINUS, XOR};

typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>											num_type;
//typedef double											num_type;

typedef CGAL::Simple_cartesian<num_type>										Exact_Kernel;

typedef CGAL::Vector_3<Exact_Kernel>											Vector_exact;
typedef CGAL::Point_3<Exact_Kernel>												Point3d_exact; 


inline Point3d_exact point_to_exact(Point3d p)
{
	return Point3d_exact(p.x(),p.y(),p.z());
}

inline Point3d point_to_double(Point3d_exact pe)
{
	return Point3d(to_double(pe.x()),to_double(pe.y()),to_double(pe.z()));
}

inline Vector_exact Compute_Normal_direction(Halfedge_handle he)
{
	return CGAL::cross_product(	point_to_exact(he->next()->vertex()->point())			- point_to_exact(he->vertex()->point()),
								point_to_exact(he->next()->next()->vertex()->point())	- point_to_exact(he->vertex()->point()));
}
	
#ifdef BOOLEAN_OPERATIONS_DEBUG
inline double tr(double n)
{
	return floor(n*1000)/1000;
}
#endif // BOOLEAN_OPERATIONS_DEBUG
	
#endif // BOOLEAN_OPERATIONS_DEFINITIONS_H
