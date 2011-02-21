#ifndef VSA_ITEMS_H
#define VSA_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class P, class Norm, class Plane>
class VSA_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		VSA_Facet() {}

		//label issu de la segmentation
		int LabelVSA;		
		
		};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class VSA_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		VSA_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>
class VSA_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		VSA_Vertex() {}

        
};


template <class kernel, class items>
class VSA_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		VSA_Polyhedron() {}

		int NbFaceLabel;//nombre de régions issues de la segmentation
};

#endif

#endif // VSA_ITEMS_H
