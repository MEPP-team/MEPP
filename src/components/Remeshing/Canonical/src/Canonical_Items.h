#ifndef Canonical_ITEMS_H
#define Canonical_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_Canonical

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class Norm>
class Canonical_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:
		Canonical_Facet() {}
		int Facet_Flag_S;
		
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Canonical_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Canonical_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>
class Canonical_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:
		Canonical_Vertex() {}
		int Vertex_Flag_S;
		int Vertex_Number_S;
		
		//For valence_driven
		int Vertex_Sign_S;
		
};


template <class kernel, class items>
class Canonical_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Canonical_Polyhedron() {}
};

#endif

#endif // Canonical_ITEMS_H
