#ifndef Boolean_Operations_ITEMS_H
#define Boolean_Operations_ITEMS_H

/*!
 * \file Boolean_Operations_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author Cyril Leconte
 */

#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"
#include "Boolean_Operations_Polyhedron.h"

/*!
 * \class Boolean_Operations_Facet
 * \brief Enriches the Facets of a Polyhedra
 */
template <class Refs, class T, class Norm>
class Boolean_Operations_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		Boolean_Operations_Facet() {}

		//subdivision :
		bool Issub; /*!< true if the facet has been subdivided*/
		//operations booleennes :
		bool IsExt; /*!< true if the facet belongs to the result*/
		bool IsOK; /*!< true if the facet has been processed*/
		FacetId Label; /*!< an Id for the facet*/
};

/*!
 * \class Boolean_Operations_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Boolean_Operations_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Boolean_Operations_Halfedge() {}
		
		//subdivision :
		bool Isnew; /*!< true if the halfedge has been created or subdivided*/
		//operations booleennes :
		HalfedgeId Label; /*!< an Id for the halfedge*/
};

/*!
 * \class Boolean_Operations_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class Boolean_Operations_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

		Boolean_Operations_Vertex() {}

		//subdivision :
		bool Isnew; /*!< true if the vertex has been created*/
		//operations booleennes :
		VertexId Label; /*!< an Id for the vertex*/
};

/*!
 * \class Boolean_Operations_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class Boolean_Operations_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Boolean_Operations_Polyhedron() {}
};

#endif // BOOLEAN_OPERATIONS

#endif // Boolean_Operations_ITEMS_H
