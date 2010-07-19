#ifndef Boolean_Operations_ITEMS_H
#define Boolean_Operations_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"
#include "Boolean_Operations_Polyhedron.h"

template <class Refs, class T, class Norm>
class Boolean_Operations_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		Boolean_Operations_Facet() {}

		//subdivision :
		bool Issub;
		//operations booleennes :
		bool IsExt;
		bool IsOK;
		FacetId Label;
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Boolean_Operations_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Boolean_Operations_Halfedge() {}
		
		//subdivision :
		bool Isnew;
		//operations booleennes :
		HalfedgeId Label;
};


template <class Refs, class T, class P, class Norm>
class Boolean_Operations_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		Boolean_Operations_Vertex() {}

		//subdivision :
		bool Isnew;
		//operations booleennes :
		VertexId Label;
};


template <class kernel, class items>
class Boolean_Operations_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Boolean_Operations_Polyhedron() {}
};

#endif // BOOLEAN_OPERATIONS

#endif // Boolean_Operations_ITEMS_H
