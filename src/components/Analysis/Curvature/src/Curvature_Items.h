#ifndef Curvature_ITEMS_H
#define Curvature_ITEMS_H

/*!
 * \file Curvature_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author \author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2008
 */
 
#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

/**
 \class	Curvature_Facet

 \brief	Enriches the Facets of a Polyhedra 

 */
template <class Refs, class T, class P, class Norm, class Plane>
class Curvature_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		Curvature_Facet() {}
};


/*!
 * \class Curvature_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Curvature_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Curvature_Halfedge() {}
};


/*!
 * \class Curvature_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class Curvature_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		Curvature_Vertex() {}

        
              
        double KminCurv;  ///< The minimum curvature
        double KmaxCurv; ///< The minimum curvature

        
        Norm VKminCurv; ///< The minimum curvature direction
        Norm VKmaxCurv; ///< The minimum curvature direction
};

/*!
 * \class Curvature_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class Curvature_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Curvature_Polyhedron() {}
};

#endif

#endif // Curvature_ITEMS_H
