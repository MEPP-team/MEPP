///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include "../../../../mepp/mepp_config.h"
#ifdef BUILD_component_CGAL_Example

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "CGAL_Example_Component.h"

// subdivision
#include "../../../../mepp/Tools/Tools_sqrt3.h"

CGAL_Example_Component::CGAL_Example_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	// from IHM
	m_CurrentColor[0] = 1.;
	m_CurrentColor[1] = 0.;
	m_CurrentColor[2] = 0.;

	// MEPP 2
	componentName = "CGAL_Example_Component";
	init = 1;
}

void CGAL_Example_Component::TriangulateAndRandomColorFacets(PolyhedronPtr pMesh)
{
	pMesh->triangulate();

	srand((unsigned)time(NULL));

	Facet_iterator pFacet =	NULL;
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		float rand255r = (float) rand() / (float) RAND_MAX;
		float rand255g = (float) rand() / (float) RAND_MAX;
		float rand255b = (float) rand() / (float) RAND_MAX;

		pFacet->color(rand255r, rand255g, rand255b);
	}

	/*Vertex_iterator pVertex = NULL;
	for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
	{
		float rand255r = (float) rand() / (float) RAND_MAX;
		float rand255g = (float) rand() / (float) RAND_MAX;
		float rand255b = (float) rand() / (float) RAND_MAX;

		pVertex->color(rand255r, rand255g, rand255b);
	}*/
}

double CGAL_Example_Component::round(double num)
{
	return floor(num + 0.5);
}

long CGAL_Example_Component::ColourDistance(unsigned char r1, unsigned char g1, unsigned char b1, unsigned char r2, unsigned char g2, unsigned char b2)
{
  long r,g,b;
  long rmean;

  rmean = ( (int)r1 + (int)r2 ) / 2;
  r = (int)r1 - (int)r2;
  g = (int)g1 - (int)g2;
  b = (int)b1 - (int)b2;

  return (((512+rmean)*r*r)>>8) + 4*g*g + (((767-rmean)*b*b)>>8);
}

void CGAL_Example_Component::CreateCenterVertex(PolyhedronPtr pMesh, bool save)
{
	CSubdivider_sqrt3<Polyhedron,Enriched_kernel> subdivider;

	long dist, maxdist=0, mindist=ColourDistance(0, 0, 0, 255, 255, 255), mindist_pct;
	Facet_iterator pFacet =	NULL;
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		pFacet->tag(0);

		dist = ColourDistance((unsigned char)round(color(0)*255.), (unsigned char)round(color(1)*255.), (unsigned char)round(color(2)*255.),
								(unsigned char)round(pFacet->color(0)*255.), (unsigned char)round(pFacet->color(1)*255.), (unsigned char)round(pFacet->color(2)*255.));

		if (dist < mindist)
			mindist = dist;
		if (dist > maxdist)
			maxdist = dist;
	}

	int pct=3;
	mindist_pct=((maxdist-mindist)*pct/100)+mindist;

	int cpt = 1;
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		if (pFacet->tag()==0)
		{
			dist = ColourDistance((unsigned char)round(color(0)*255.), (unsigned char)round(color(1)*255.), (unsigned char)round(color(2)*255.),
								(unsigned char)round(pFacet->color(0)*255.), (unsigned char)round(pFacet->color(1)*255.), (unsigned char)round(pFacet->color(2)*255.));

			if ((dist >= mindist) && (dist <= mindist_pct))
			{
				pFacet->tag(1);
				Vertex_handle hVertex = subdivider.create_center_vertex(*pMesh, pFacet)->vertex();

				float rand255r = (float) rand() / (float) RAND_MAX;
				float rand255g = (float) rand() / (float) RAND_MAX;
				float rand255b = (float) rand() / (float) RAND_MAX;
				hVertex->color(rand255r, rand255g, rand255b);

				if (save)
				{
					if (cpt <= 5) // example : save only 5 first positives iterations
					{
						time_t ltime;
						time(&ltime);
						/*char filename[4096];

						sprintf(filename, "/tmp/CreateCenterVertex_%03d_%ld.off", cpt, ltime);
						pMesh->write_off(filename, true, false); // save colors but not normals

						sprintf(filename, "/tmp/CreateCenterVertex_%03d_%ld.obj", cpt, ltime);
						pMesh->write_obj(filename);

						sprintf(filename, "/tmp/CreateCenterVertex_%03d_%ld.wrl", cpt, ltime);
						pMesh->write_vrml(filename);*/

						cpt++;
					}
				}
			}
		}
	}
}

void CGAL_Example_Component::ShowBlackAndWhiteFacets(PolyhedronPtr pMesh)
{
	Facet_iterator pFacet =	NULL;
	for (pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
	{
		int all_tags_are_one;

		if (pFacet->tag()==1)
		{
			// circulate around facet
			Halfedge_around_facet_circulator he, end;
			he = end = pFacet->facet_begin();
			all_tags_are_one=1;

			CGAL_For_all(he, end)
			{
				if (he->opposite()->face() != NULL)
				{
					if (he->opposite()->face()->tag()==0)
					{
						he->opposite()->face()->color(1., 1., 1.);
						all_tags_are_one=0;
					}
				}
			}
			if (all_tags_are_one)
				pFacet->color(0., 0., 0.);
		}
	}
}

void CGAL_Example_Component::CreateTetrahedron(PolyhedronPtr pMesh)
{
	Point3d p1( -0.5, -0.5, -0.5);
	Point3d q1( 0.5, -0.5, -0.5);
	Point3d r1( 0.0, -0.5, 0.5);
	Point3d s1( 0.0, 0.5, 0.0);

	Halfedge_handle h = pMesh->make_tetrahedron(p1, q1, r1, s1);
	if (pMesh->is_tetrahedron(h))
	{
		pMesh->compute_bounding_box();

		pMesh->compute_normals();
		pMesh->compute_type();
	}
}

void CGAL_Example_Component::GetClickedVertices(PolyhedronPtr pMesh, double x, double y, int tolerance)
{
	//wxString statusString;
	GLdouble *model ;  GLdouble *proj ;  GLint *view;

	view=new int[4096];
	model=new double[4096];
	proj=new double[4096];

	glGetIntegerv (GL_VIEWPORT, view);
	glGetDoublev (GL_MODELVIEW_MATRIX, model);
	glGetDoublev (GL_PROJECTION_MATRIX, proj);

	y=view[3]-y;

	GLdouble wx ; GLdouble wy ; GLdouble wz;

	int vertexID=0;
	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
		gluProject (pVertex->point().x(),pVertex->point().y(),pVertex->point().z(), model, proj, view, &wx, &wy, &wz);  // on simule la projection du sommet dans l'espace window
		if (wz>0. && wz<1)
		if (x>floor(wx)-tolerance && x<floor(wx)+tolerance)
		if (y>floor(wy)-tolerance && y<floor(wy)+tolerance)  // on fixe une petite tolérance (de 2*5 pixels par exemple) pour la sélection
		{
			pVertex->color(1., 0., 0.);

			/*statusString.Printf(_T("Vertex: %u  -  (%f, %f, %f)"), vertexID, pVertex->point().x(), pVertex->point().y(), pVertex->point().z());
			pFrame->set_status_message(statusString);*/
		}
		vertexID++;
	}

	delete[]view;
	delete[]model;
	delete[]proj;
}
#endif
