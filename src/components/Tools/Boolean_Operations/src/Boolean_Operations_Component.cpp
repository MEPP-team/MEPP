#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

///////////////////////////////////////////////////////////////////////////
// Author: Cyril Leconte
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

#include "BoolPolyhedra.h"
#include "Boolean_Operations_Component.h"

void Boolean_Operations_Component::SubdiviserPolyedre(PolyhedronPtr pMesh)
{
	//on s'assure que chaque face du polyedre est triangulaire :
	if(!pMesh->is_pure_triangle())
	{
		pMesh->triangulate();
		return;
	}
	//declaration des iterateurs et variables :
	Facet_iterator pFacet;
	Vector Vcenter;
	
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
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		if(!(pFacet->Issub))
		{
			Halfedge_handle pHE = pFacet->facet_begin();
			for(unsigned int i = 0;i!=5;i++)
			{
				if(!pHE->Isnew)
				{
					Vcenter = Vector(0.0, 0.0, 0.0);
					Vcenter = ( (pHE->vertex()->point() - CGAL::ORIGIN) + (pHE->opposite()->vertex()->point() - CGAL::ORIGIN) ) / 2;
					pHE = pMesh->split_edge(pHE);
					pHE->vertex()->point() = CGAL::ORIGIN + Vcenter;
					pHE->vertex()->Isnew = true;
					pHE->Isnew = true;
					pHE->opposite()->Isnew = true;
					pHE->next()->Isnew = true;
					pHE->next()->opposite()->Isnew = true;
				}
				pHE = pHE->next();
			}
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

#endif
