#ifndef VSA_ITEMS_H
#define VSA_ITEMS_H

/*!
 * \file VSA_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
 */
 
#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

/*!
 * \class VSA_Facet
 * \brief Enriches the Facets of a Polyhedra
 */
 template <class Refs, class T, class P, class Norm, class Plane>
class VSA_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		VSA_Facet() {}

		
		int LabelVSA;	///< The label from the segmentation
		
		};


/*!
 * \class VSA_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class VSA_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		VSA_Halfedge() {}
};


/*!
 * \class VSA_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class VSA_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		VSA_Vertex() {}

        
};


/*!
 * \class VSA_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class VSA_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		VSA_Polyhedron() {}

		
		int NbFaceLabel; ///< The number of proxies
};

#endif

#endif // VSA_ITEMS_H
