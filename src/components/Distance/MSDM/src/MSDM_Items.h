#ifndef MSDM_ITEMS_H
#define MSDM_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class P, class Norm, class Plane>
class MSDM_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		MSDM_Facet() {}
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MSDM_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		MSDM_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>
class MSDM_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		MSDM_Vertex() {}
		float CourbureCoVariance;
		float CourbureVariance;
		float CourbureMoyenne;
		float MSDM_Local;


};


template <class kernel, class items>
class MSDM_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		MSDM_Polyhedron() {}
};

#endif

#endif // MSDM_ITEMS_H
