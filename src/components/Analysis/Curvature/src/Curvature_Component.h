///////////////////////////////////////////////////////////////////////////
// Author: Guillaume Lavoué
// Year: 2008
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////

#ifndef Curvature_COMPONENT_H
#define Curvature_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "../../../../mepp/mepp_component.h"

#include "Curvature_Polyhedron.h"
#include "extract_Vpropres.h"

class Curvature_Component : 
  public mepp_component
{
	public:

		Curvature_Component(Viewer* v, PolyhedronPtr p);
		~Curvature_Component() {}

		void principal_curvature(PolyhedronPtr pMesh,bool IsGeod,double radius);
		void ConstructColorMap(PolyhedronPtr pMesh,int ColorField);//0 rien, 1 min, 2 max.

		void set_displayMinDirections(bool min) { displayMinDirections = min; }
		bool get_displayMinDirections() { return displayMinDirections; }

		void set_displayMaxDirections(bool max) { displayMaxDirections = max; }
		bool get_displayMaxDirections() { return displayMaxDirections; }

	private:

		//les min et max respectif des courbures min et max (utile pour l'affichage)
		double MinNrmMinCurvature,MaxNrmMinCurvature,MinNrmMaxCurvature,MaxNrmMaxCurvature;
		double LUT_CourbureClust[3*256]; //table d'affichage couleur du bleu au rouge, en dégradé

	// from IHM
	private:
		bool displayMinDirections; //est-ce qu'on affiche les directions de courbure minimum
		bool displayMaxDirections; //est-ce qu'on affiche les directions de courbure maximum
};

#endif

#endif // Curvature_COMPONENT_H
