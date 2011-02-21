///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////
#ifndef CGAL_Example_ITEMS_H
#define CGAL_Example_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class P, class Norm, class Plane>
class CGAL_Example_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		CGAL_Example_Facet() {}
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class CGAL_Example_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		CGAL_Example_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>
class CGAL_Example_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		CGAL_Example_Vertex() {}
};


template <class kernel, class items>
class CGAL_Example_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		CGAL_Example_Polyhedron() {}
};

#endif

#endif // CGAL_Example_ITEMS_H
