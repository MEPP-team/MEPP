///////////////////////////////////////////////////////////////////////////
// Author: 
// Year: 
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef Various_Tools_COMPONENT_H
#define Various_Tools_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Tools

#include "../../../../mepp/mepp_component.h"
#include "Various_Tools_Polyhedron.h"

class Viewer;

class Various_Tools_Component : 
  public mepp_component
{
	public:
		Various_Tools_Component(Viewer* v, PolyhedronPtr p);
		~Various_Tools_Component() {}

		bool subdivide_sqrt3(PolyhedronPtr pMesh);
		bool subdivide_sqrt3Twice(PolyhedronPtr pMesh);
		bool subdivide_quad(PolyhedronPtr pMesh);
		bool subdivide_doo(PolyhedronPtr pMesh);
		bool subdivide_loop(PolyhedronPtr pMesh);
		bool subdivide_catmull(PolyhedronPtr pMesh);
		bool triangulate(PolyhedronPtr pMesh);
};

#endif

#endif // Various_Tools_COMPONENT_H
