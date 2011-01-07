///////////////////////////////////////////////////////////////////////////
// Author: 
// Year: 
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include <mepp_config.h>
#ifdef BUILD_component_Various_Tools

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "Various_Tools_Component.h"

// subdivision
#include "../../../../mepp/Tools/Tools_sqrt3.h"
#include "../../../../mepp/Tools/Tools_quad-triangle.h"
#include "../../../../mepp/Tools/Tools_Polyhedron_subdivision.h"

Various_Tools_Component::Various_Tools_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	// MEPP 2
	componentName = "Various_Tools_Component";
	init = 1;
}

bool Various_Tools_Component::subdivide_sqrt3(PolyhedronPtr pMesh)
{
	CSubdivider_sqrt3<Polyhedron,Enriched_kernel> subdivider;
	subdivider.subdivide(*pMesh,1); // one iteration
	pMesh->compute_normals();
	pMesh->compute_type();

	return true;
}

bool Various_Tools_Component::subdivide_sqrt3Twice(PolyhedronPtr pMesh)
{
	CSubdivider_sqrt3<Polyhedron,Enriched_kernel> subdivider;
	subdivider.subdivide(*pMesh,2); // two iterations
	pMesh->compute_normals();
	pMesh->compute_type();

	return true;
}

bool Various_Tools_Component::subdivide_quad(PolyhedronPtr pMesh)
{
	CSubdivider_quad_triangle<Polyhedron,Enriched_kernel> subdivider;
	Polyhedron* new_mesh = new Polyhedron;
	subdivider.subdivide(*pMesh, *new_mesh, true);

	// copy bounding box (approximate, but fast)
	//new_mesh->copy_bounding_box(polyhedron);

	pMesh->copy_from(new_mesh);
	delete new_mesh;

	/*polyhedron->compute_normals();
	polyhedron->compute_type();*/ // dejà fait lors de la copie

	return true;
}

bool Various_Tools_Component::subdivide_doo(PolyhedronPtr pMesh)
{
	Polyhedron_subdivision<Polyhedron>::DooSabin_subdivision(*pMesh,1);
	pMesh->compute_normals();
	pMesh->compute_type();

	return true;
}

bool Various_Tools_Component::subdivide_loop(PolyhedronPtr pMesh)
{
	Polyhedron_subdivision<Polyhedron>::Loop_subdivision(*pMesh,1);
	pMesh->compute_normals();
	pMesh->compute_type();

	return true;
}

bool Various_Tools_Component::subdivide_catmull(PolyhedronPtr pMesh)
{
	Polyhedron_subdivision<Polyhedron>::CatmullClark_subdivision(*pMesh,1);
	pMesh->compute_normals();
	pMesh->compute_type();

	return true;
}

bool Various_Tools_Component::triangulate(PolyhedronPtr pMesh)
{
    pMesh->triangulate();

    return true;
}
#endif
