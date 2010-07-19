#ifndef Boolean_Operations_COMPONENT_H
#define Boolean_Operations_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

#include "../../../../mepp/mepp_component.h"
#include "Boolean_Operations_Definitions.h"

class Viewer;

class Boolean_Operations_Component : public mepp_component
{
public:
	Boolean_Operations_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p) {}
	~Boolean_Operations_Component() {}
	
	void SubdiviserPolyedre(PolyhedronPtr pMesh);
	
	void Boolean_Union(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout);
	void Boolean_Inter(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout);
	void Boolean_Minus(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout);
};

#endif // BOOLEAN_OPERATIONS

#endif // Boolean_Operations_COMPONENT_H
