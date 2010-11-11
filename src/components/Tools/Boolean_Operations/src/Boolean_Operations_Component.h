#ifndef Boolean_Operations_COMPONENT_H
#define Boolean_Operations_COMPONENT_H

///////////////////////////////////////////////////////////////////////////
// Author: Cyril Leconte
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
//
//  According to :
// 	"Exact and Efficient Booleans for Polyhedra", C. Leconte, H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018, 2010
//
///////////////////////////////////////////////////////////////////////////


#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

#include "../../../../mepp/mepp_component.h"
#include "Boolean_Operations_Definitions.h"

class Viewer;

/*!
 * \class Boolean_Operations_Component
 * \brief Creates a component to use Boolean operations on polyhedra
 */

class Boolean_Operations_Component : public mepp_component
{
public:
	/*!
	 * \brief Constructor
	 */
	Boolean_Operations_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
	{
		cpt_U = cpt_I = cpt_M = 0;
	}
	
	/*!
	 * \brief Destructor
	 */
	~Boolean_Operations_Component() {}
	
	/*!
	 * \brief This method subdivide each facet of a polyhedron
	 * into four facets (by adding a vertex in the middle of each edge)
	 * \param pMesh : The polyhedron to subdivide
	 */
	void SubdiviserPolyedre(PolyhedronPtr pMesh);
	
	/*!
	 * \brief Computes the union of two polyhedra
	 * \param pMin1 : First input polyhedra
	 * \param pMin2 : Second input polyhedra
	 * \param pMout : an empty polyhedron to build the result (pMin1 union pMin2)
	 */
	void Boolean_Union(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout);
	
	/*!
	 * \brief Computes the intersection of two polyhedra
	 * \param pMin1 : First input polyhedra
	 * \param pMin2 : Second input polyhedra
	 * \param pMout : an empty polyhedron to build the result (pMin1 inter pMin2)
	 */
	void Boolean_Inter(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout);
	
	/*!
	 * \brief Computes the subtraction between two polyhedra
	 * \param pMin1 : First input polyhedra
	 * \param pMin2 : Second input polyhedra
	 * \param pMout : an empty polyhedron to build the result (pMin1 minus pMin2)
	 */
	void Boolean_Minus(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout);

	// MT
	int cpt_U, cpt_I, cpt_M;
};

#endif // BOOLEAN_OPERATIONS

#endif // Boolean_Operations_COMPONENT_H
