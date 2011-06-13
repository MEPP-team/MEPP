/*!
 * \file MSDM_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
 */
 
#ifndef MSDM2_ITEMS_H
#define MSDM2_ITEMS_H

#include <mepp_config.h>
 #ifdef BUILD_component_MSDM2

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

/*!
 * \class MSDM_Facet
 * \brief Enriches the Facets of a Polyhedra
 */
template <class Refs, class T, class P, class Norm, class Plane>
class MSDM2_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		MSDM2_Facet() {}
};

/*!
 * \class MSDM_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MSDM2_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		MSDM2_Halfedge() {}
};

/*!
 * \class MSDM_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class MSDM2_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		MSDM2_Vertex() {}
		
		float MSDM2_Local; ///< Local MSDM2 value.
		P match; ///< 3D position of the projection of the vertex on the reference mesh
		float curvmatch; ///< interpolated curvature value of the projection of the vertex on the reference mesh


};

/*!
 * \class MSDM_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class MSDM2_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		MSDM2_Polyhedron() {}
		bool IsDistanceComputed;///< true if the distance has been computed
};

#endif

#endif // MSDM2_ITEMS_H
