#ifndef MSDM_ITEMS_H
#define MSDM_ITEMS_H

/*!
 * \file MSDM_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
 */
 
#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

/*!
 * \class MSDM_Facet
 * \brief Enriches the Facets of a Polyhedra
 */
template <class Refs, class T, class P, class Norm, class Plane>
class MSDM_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		MSDM_Facet() {}
};

/*!
 * \class MSDM_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MSDM_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		MSDM_Halfedge() {}
};

/*!
 * \class MSDM_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class MSDM_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		MSDM_Vertex() {}
		
		float CourbureCoVariance; ///< Covariance of curvature over the neighborhood.
		float CourbureVariance; ///< Variance of curvature over the neighborhood.
		float CourbureMoyenne; ///< Mean of curvature over the neighborhood.
		
		float MSDM_Local; ///< Local MSDM value.


};

/*!
 * \class MSDM_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class MSDM_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		MSDM_Polyhedron() {}
};

#endif

#endif // MSDM_ITEMS_H
