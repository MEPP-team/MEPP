///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef CGAL_Example_COMPONENT_H
#define CGAL_Example_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

#include "../../../../mepp/mepp_component.h"
#include "CGAL_Example_Polyhedron.h"

class Viewer;

class CGAL_Example_Component : 
  public mepp_component
{
	public:
		CGAL_Example_Component(Viewer* v, PolyhedronPtr p);
		~CGAL_Example_Component() {}

		void TriangulateAndRandomColorFacets(PolyhedronPtr pMesh);
		void CreateCenterVertex(PolyhedronPtr pMesh, bool save);
		void ShowBlackAndWhiteFacets(PolyhedronPtr pMesh);
		void CreateTetrahedron(PolyhedronPtr pMesh);

		float color(int index) { return m_CurrentColor[index]; };
		void color(float r, float g, float b) { m_CurrentColor[0] = r; m_CurrentColor[1] = g; m_CurrentColor[2] = b; };

		double round(double num);
		long ColourDistance(unsigned char r1, unsigned char g1, unsigned char b1, unsigned char r2, unsigned char g2, unsigned char b2);

		string GetClickedVertices(PolyhedronPtr pMesh, double x, double y, int tolerance);

	// from IHM
	private:
		float m_CurrentColor[3];	
};

#endif

#endif // CGAL_Example_COMPONENT_H
