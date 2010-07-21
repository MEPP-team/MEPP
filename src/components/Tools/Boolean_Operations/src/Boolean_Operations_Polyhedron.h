#ifndef Boolean_Operations_POLYHEDRON_H
#define Boolean_Operations_POLYHEDRON_H

/*!
 * \file Boolean_Operations_Polyhedron.h
 * \brief Defines four types used in polyhedra 
 * \author Cyril Leconte
 */

#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

#include "../../../../mepp/Polyhedron/polyhedron.h"

/*!
 * \typedef VertexId
 * \brief Vertex Id
 */
typedef unsigned long							VertexId;

/*!
 * \typedef HalfedgeId
 * \brief Halfedge Id
 */
typedef unsigned long							HalfedgeId;

/*!
 * \typedef FacetId
 * \brief Facet Id
 */
typedef unsigned long							FacetId;

/*!
 * \typedef InterId
 * \brief Intersection Id
 */
typedef unsigned long							InterId;

#endif // BUILD_component_Boolean_Operations

#endif // Boolean_Operations_POLYHEDRON_H
