/***************************************************************************
Polyhedron_topological_transformation_3.h - Transformations affecting only 
the polyhedron's topology
----------------------------------------------------------------------------
begin                : March 2009
copyright            : (C) 2008 Hichem Barki, LIRIS M2DisCo
email                : hichem.barki@liris.cnrs.fr
***************************************************************************/

#ifndef _POLYHEDRON_TOPOLOGICAL_TRANSFORMATION_3_H
#define _POLYHEDRON_TOPOLOGICAL_TRANSFORMATION_3_H

#include <CGAL/basic.h>
#include <CGAL/circulator.h>

namespace Polyhedron_topological_transformation_3
{
	template <typename Polyhedron_3>
	inline void convex_triangulate_polyhedron(Polyhedron_3& P)
	{
		// Triangulate a polyhedron, all facets must be convex,
		// otherwise the result will not be correct
		// Uses Euler modifiers
		
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Facet_handle										Facet_handle;
		typedef typename Polyhedron_3::Halfedge_around_facet_circulator					Halfedge_facet_circulator;
		typedef typename Polyhedron_3::Halfedge_handle									Halfedge_handle;

		Halfedge_facet_circulator hfc;
		for (Facet_handle fit = P.facets_begin(); fit != P.facets_end(); ++fit)
		{

			unsigned long split_times = (unsigned long) CGAL::circulator_size(fit->facet_begin()) - 3;

			if (split_times > 0) // The facet is not triangular
			{
				Halfedge_handle hh = fit->facet_begin();
				for (int i = 0; i != split_times; ++i)
				{
					hh = P.split_facet(hh, hh->next()->next());
				}
			}
		}
	}

} // end of namespace Polyhedron_topological_transformation_3

#endif // end of _POLYHEDRON_TOPOLOGICAL_TRANSFORMATION_3_H