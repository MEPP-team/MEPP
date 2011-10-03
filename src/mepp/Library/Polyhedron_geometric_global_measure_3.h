/***************************************************************************
Polyhedron_geometric_global_measure_3.h  -  Geometric global
measures for polyhedra.
----------------------------------------------------------------------------
begin                : March 2009
copyright            : (C) 2008 Hichem Barki, LIRIS M2DisCo
email                : hichem.barki@liris.cnrs.fr
***************************************************************************/
#ifndef _POLYHEDRON_GEOMETRIC_GLOBAL_MEASURE_3_H
#define _POLYHEDRON_GEOMETRIC_GLOBAL_MEASURE_3_H

#include <cmath>

// CGAL Stuff
#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include "Polyhedron_topological_transformation_3.h"

using CGAL::to_double;

namespace Polyhedron_geometric_global_measure_3
{
	template <class Polyhedron_3>
	bool is_convex_without_coplanar_facets(const Polyhedron_3& P)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits								Traits;
		typedef typename Traits::FT									Number_type;
		typedef typename Traits::Triangle_3								Triangle_3;
		typedef typename Polyhedron_3::Facet_const_iterator						Facet_const_it;
		typedef typename Polyhedron_3::Edge_const_iterator						Edge_const_it;
		typedef typename Traits::Vector_3								Vector_3;

		for (Edge_const_it eit = P.edges_begin(); eit != P.edges_end(); ++eit)
		{
			Vector_3 edge(eit->opposite()->vertex()->point(), eit->vertex()->point());
			Vector_3 cur_normal = CGAL::cross_product(edge, Vector_3(eit->opposite()->vertex()->point(), eit->next()->vertex()->point()));
			Vector_3 opp_normal = CGAL::cross_product(-edge, Vector_3(eit->vertex()->point(), eit->opposite()->next()->vertex()->point()));
			Vector_3 cp = CGAL::cross_product(opp_normal, cur_normal);
			if (cp == CGAL::NULL_VECTOR)
			{
				return false;
			}
			else
			{
				if (cp * edge > 0)
				{
					return false;
				}
			}
		}
		return true;		
	}
	template <class Polyhedron_3>
	bool is_convex(const Polyhedron_3& P)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits								Traits;
		typedef typename Traits::FT									Number_type;
		typedef typename Traits::Triangle_3								Triangle_3;
		typedef typename Polyhedron_3::Facet_const_iterator						Facet_const_it;
		typedef typename Polyhedron_3::Edge_const_iterator						Edge_const_it;
                typedef typename Traits::Vector_3                                                               Vector_3; // MT

		for (Edge_const_it eit = P.edges_begin(); eit != P.edges_end(); ++eit)
		{
			Vector_3 edge(eit->opposite()->vertex()->point(), eit->vertex()->point());
			Vector_3 cur_normal = CGAL::cross_product(edge, Vector_3(eit->opposite()->vertex()->point(), eit->next()->vertex()->point()));
			Vector_3 opp_normal = CGAL::cross_product(-edge, Vector_3(eit->vertex()->point(), eit->opposite()->next()->vertex()->point()));
			if (edge * CGAL::cross_product(opp_normal, cur_normal) > 0)
			{
				return false;
			}
		}
		return true;		
	}
	template <class Polyhedron_3>
	typename Polyhedron_3::Traits::FT get_surface(const Polyhedron_3& P)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits									Traits;
		typedef typename Traits::FT												Number_type;
		typedef typename Traits::Triangle_3										Triangle_3;
		typedef typename Polyhedron_3::Facet_const_iterator						Facet_const_it;
		typedef typename Polyhedron_3::Halfedge_const_iterator					Halfedge_const_it;
		Number_type				surface = 0;
		Polyhedron_3 temp = P;
		Polyhedron_topological_transformation_3::convex_triangulate_polyhedron<Polyhedron_3>(temp);
		for (Facet_const_it fit = temp.facets_begin(); fit != temp.facets_end(); ++fit)
		{
			Halfedge_const_it hh = fit->facet_begin();
			Triangle_3 t(hh->vertex()->point(), hh->next()->vertex()->point(), hh->next()->next()->vertex()->point());
			surface += Number_type(sqrt(to_double(t.squared_area())));
		}
		return surface;
	}
	template <class Polyhedron_3>
	typename Polyhedron_3::Traits::FT get_volume(const Polyhedron_3& P)
	{
		CGAL_assertion(&P != NULL);
		typedef typename Polyhedron_3::Traits									Traits;
		typedef typename Traits::FT												Number_type;
		typedef typename Traits::Tetrahedron_3									Tetrahedron_3;
		typedef typename Polyhedron_3::Facet_const_iterator						Facet_const_it;
		typedef typename Polyhedron_3::Halfedge_const_iterator					Halfedge_const_it;
		Number_type volume = 0;
		Polyhedron_3 temp = P;
		Polyhedron_topological_transformation_3::convex_triangulate_polyhedron<Polyhedron_3>(temp);
		for (Facet_const_it fit = temp.facets_begin(); fit != temp.facets_end(); ++fit)
		{
			Halfedge_const_it hh = fit->facet_begin();
			Tetrahedron_3 t(CGAL::ORIGIN, hh->vertex()->point(), hh->next()->vertex()->point(), hh->next()->next()->vertex()->point());
			volume += t.volume();
		}
		return volume;
	}
	

} // end of namespace Polyhedron_geometric_global_measure_3

#endif // end of _POLYHEDRON_GEOMETRIC_GLOBAL_MEASURE_3_H