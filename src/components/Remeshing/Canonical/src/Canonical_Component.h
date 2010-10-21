///////////////////////////////////////////////////////////////////////////
// Author: Ho Lee
// Year: 2009
// LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Canonical_COMPONENT_H
#define Canonical_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Canonical

#include "../../../../mepp/mepp_component.h"
#include "Canonical_Polyhedron.h"

#include <queue>

class Viewer;

class Canonical_Component : 
  public mepp_component
{
	public:
		Canonical_Component(Viewer* v, PolyhedronPtr p);
		~Canonical_Component() {}

		// This function decimates vertices which have 3-6 valences and retriangulate the hole left.
		void Decimation_Conquest(PolyhedronPtr pMesh,const bool Normal_flipping,const bool Use_metric,
											const float metric_threshold,const bool Use_forget_metric,const unsigned int &Forget_value);

		// This function decimates vertices of valence 3.
		void Regulation(PolyhedronPtr pMesh,const bool Normal_flipping,const bool Use_metric,
											const float metric_threshold,const bool Use_forget_metric,const unsigned int &Forget_value);

	private:

};

#endif

#endif // Canonical_COMPONENT_H
