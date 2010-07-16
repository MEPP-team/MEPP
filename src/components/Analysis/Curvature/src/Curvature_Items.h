#ifndef Curvature_ITEMS_H
#define Curvature_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class Norm>
class Curvature_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		Curvature_Facet() {}
};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Curvature_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Curvature_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>
class Curvature_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		Curvature_Vertex() {}

        //attributs
        //MinCurvature
        double Kmin;
        //MaxCurvature
        double Kmax;

        //vecteur courbures principales
        Norm VKmin;
        Norm VKmax;
};


template <class kernel, class items>
class Curvature_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Curvature_Polyhedron() {}
};

#endif

#endif // Curvature_ITEMS_H
