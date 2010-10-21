#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

///////////////////////////////////////////////////////////////////////////
// Author: Cyril Leconte
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
//
//  According to :
// 	"Exact and Efficient Booleans for Polyhedra", C. Leconte, H. Barki, F. Dupont, Rapport de recherche RR-LIRIS-2010-018, 2010
//
///////////////////////////////////////////////////////////////////////////

#include "BoolPolyhedra.h"
#include "Boolean_Operations_Component.h"

void Boolean_Operations_Component::SubdiviserPolyedre(PolyhedronPtr pMesh)
{
	//Each facet must be triangular
	if(!pMesh->is_pure_triangle())
	{
		pMesh->triangulate();
		return;
	}
	
	Facet_iterator pFacet;
	Vector Vcenter;

	//Initialization of the tags
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		Halfedge_around_facet_circulator pHEcirc = pFacet->facet_begin();
		pFacet->Issub = false;
		pHEcirc->Isnew = false;
		pHEcirc->vertex()->Isnew = false;
		pHEcirc++;
		pHEcirc->Isnew = false;
		pHEcirc->vertex()->Isnew = false;
		pHEcirc++;
		pHEcirc->Isnew = false;
		pHEcirc->vertex()->Isnew = false;
	}
	//For each facet of the polyhedron
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		//We subdivide the facet if it is not already done
		if(!(pFacet->Issub))
		{
			Halfedge_handle pHE = pFacet->facet_begin();
			for(unsigned int i = 0;i!=5;i++)
			{
				if(!pHE->Isnew)
				{
					//each edge is splited in its center
					Vcenter = Vector(0.0, 0.0, 0.0);
					Vcenter = ( (pHE->vertex()->point() - CGAL::ORIGIN) + (pHE->opposite()->vertex()->point() - CGAL::ORIGIN) ) / 2;
					pHE = pMesh->split_edge(pHE);
					pHE->vertex()->point() = CGAL::ORIGIN + Vcenter;
					//update of the tags (the new vertex and the four new halfedges
					pHE->vertex()->Isnew = true;
					pHE->Isnew = true;
					pHE->opposite()->Isnew = true;
					pHE->next()->Isnew = true;
					pHE->next()->opposite()->Isnew = true;
				}
				pHE = pHE->next();
			}
			//Three new edges are build between the three new vertices, and the tags of the facets are updated
			if(!pHE->vertex()->Isnew) pHE = pHE->next();
			pHE = pMesh->split_facet(pHE, pHE->next()->next());
			pHE->opposite()->facet()->Issub = true;
			pHE = pMesh->split_facet(pHE, pHE->next()->next());
			pHE->opposite()->facet()->Issub = true;
			pHE = pMesh->split_facet(pHE, pHE->next()->next());
			pHE->opposite()->facet()->Issub = true;
			pHE->facet()->Issub = true;
		}
	}
}


void Boolean_Operations_Component::Boolean_Union(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout)
{
	BoolPolyhedra(pMin1, pMin2, pMout, UNION);
}

void Boolean_Operations_Component::Boolean_Inter(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout)
{
	//*
	BoolPolyhedra(pMin1, pMin2, pMout, INTER);
	/*/
	pMin1->inside_out();
	pMin2->inside_out();
	BoolPolyhedra(pMin1, pMin2, pMout, UNION);
	pMout->inside_out();
	//*/
}

void Boolean_Operations_Component::Boolean_Minus(PolyhedronPtr pMin1, PolyhedronPtr pMin2, PolyhedronPtr pMout)
{
	//*
	BoolPolyhedra(pMin1, pMin2, pMout, MINUS);
	/*/
	pMin1->inside_out();
	BoolPolyhedra(pMin1, pMin2, pMout, UNION);
	pMout->inside_out();
	//*/
}

#endif // BUILD_component_Boolean_Operations
