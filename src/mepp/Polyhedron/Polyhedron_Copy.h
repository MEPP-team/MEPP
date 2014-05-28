/***************************************************************************
copy_poly.h  -
----------------------------------------------------------------------------
begin                : march 2006
copyright            : (C) 2006 by Celine Roudet - Liris
email                : croudet@liris.cnrs.fr
***************************************************************************/

#ifndef COPY_POLY
#define COPY_POLY

#include "Polyhedron_Builder.h"
#include <CGAL/circulator.h>
#include <CGAL/basic.h>
#include <iostream>

template <class HDS,class Polyhedron,class kernel>
class CModifierCopyPoly : public CGAL::Modifier_base<HDS>
{
private:
	typedef typename HDS::Vertex          Vertex;
	typedef typename Vertex::Point        Point;
	typedef typename HDS::Face_handle     Face_handle;
	typedef typename HDS::Halfedge_handle Halfedge_handle;
	typedef typename CGAL::Enriched_polyhedron_incremental_builder_3<HDS> builder;

	typedef typename kernel::FT FT;
	typedef typename kernel::Vector_3 Vector;  //ajout Céline
	typedef typename kernel::Plane_3 Plane;
	typedef typename CGAL::Triangle_3<kernel> Triangle_3;
	typedef typename Polyhedron::Vertex_handle Vertex_handle;
	typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
	typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
	typedef typename Polyhedron::Edge_iterator Edge_iterator;
	typedef typename Polyhedron::Facet_handle Facet_handle; // MT
	typedef typename Polyhedron::Facet_iterator Facet_iterator;

	typedef typename Polyhedron::Halfedge_around_vertex_circulator
											Halfedge_around_vertex_circulator;  // fin ajout Céline
	typedef typename Polyhedron::Halfedge_around_facet_circulator
											Halfedge_around_facet_circulator;  // fin ajout Céline

	//private fileds :
	Polyhedron *m_pMesh;

public:

	// life cycle
	CModifierCopyPoly(Polyhedron *pMesh)
	{
		CGAL_assertion(pMesh != NULL);
		m_pMesh = pMesh;
	}

	~CModifierCopyPoly() {}

 //////////////////////////////////////////////// BUILDER /////////////////////////////////////////////

	void operator()( HDS& hds)
	{
		builder B(hds,true);        //builder init
		B.begin_surface(3,1,6);     //initialisation of the new polyhedron
			add_vertices(B);            //create all vertices for the new polyhedron
			add_facets(B);
		B.end_surface();
	}

	// add vertices
	void add_vertices(builder &B)
	{
		int index = 0;
		Vertex_iterator pVertex;
		for (pVertex = m_pMesh->vertices_begin() ; pVertex != m_pMesh->vertices_end() ; pVertex++)
		{
			pVertex->tag(index); // tag each original vertex
                        Vertex_handle vertex = B.add_vertex(pVertex->point());  // add original vertices to the new poly

			vertex->color(pVertex->color(0), pVertex->color(1), pVertex->color(2)); // MT: add color
			vertex->texture_coordinates(pVertex->texture_coordinates(0), pVertex->texture_coordinates(1)); // texture coordinates

			index++;
		}
	}


	void add_facets(builder &B)
	{
		Facet_iterator pFacet;
		for (pFacet = m_pMesh->facets_begin() ; pFacet != m_pMesh->facets_end() ; pFacet++)
		{
                    unsigned int degree = 0;
                    degree = Polyhedron::degree(pFacet);  // facet degree
                    CGAL_assertion(degree >= 3);
                    degree = degree; // just for warning with gcc 4.6

                    Halfedge_handle pHalfedge = pFacet->halfedge();
                    //int tag = 0; // MT

                    Facet_handle facet = B.begin_facet();
                    do
                    {
                        B.add_vertex_to_facet(pHalfedge->vertex()->tag());
                        pHalfedge = pHalfedge->next();
                    } while (pHalfedge != pFacet->halfedge());

                    B.end_facet();

					 Halfedge_handle new_h = facet->halfedge();
					 Halfedge_handle old_h = pFacet->halfedge();
					 do
                    {
						new_h->texture_coordinates(old_h->texture_coordinates(0), old_h->texture_coordinates(1));						
						new_h = new_h->next();
						old_h = old_h->next();
					
                    } while (old_h != pFacet->halfedge());

                    facet->color(pFacet->color(0), pFacet->color(1), pFacet->color(2)); // MT: add color
		}
		CGAL_assertion(!B.check_unconnected_vertices());
	}
};

template <class Polyhedron, class kernel>
class CCopyPoly
{
       typedef typename kernel::Vector_3 Vector;
       typedef typename Polyhedron::HalfedgeDS HalfedgeDS;

public:

	CCopyPoly() {}
	~CCopyPoly() {}

public:
   int copy(Polyhedron *OriginalMesh, Polyhedron *NewMesh)
   {
	CModifierCopyPoly<HalfedgeDS, Polyhedron, kernel> builder(OriginalMesh);
        NewMesh->delegate(builder);                     //calls the `operator()' of the `modifier'
        return(0);
   }
};

#endif // COPY_POLY
