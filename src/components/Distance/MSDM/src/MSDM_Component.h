///////////////////////////////////////////////////////////////////////////
// Author: Guillaume LAVOUE
// Year: 2010
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
//According to: Perceptually driven 3d distance metrics with application to watermarking
//			G. Lavoué, E. Drelie Gelasca, F. Dupont, A. Baskurt, and T. Ebrahimi 
//			In SPIE Applications of Digital Image Processing XXIX, vol. 6312, 2006.
//
///////////////////////////////////////////////////////////////////////////

#ifndef MSDM_COMPONENT_H
#define MSDM_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#include "../../../../mepp/mepp_component.h"
#include "MSDM_Polyhedron.h"

class MSDM_Component : 
  public mepp_component
{
	public:
		MSDM_Component(Viewer* v, PolyhedronPtr p);
		~MSDM_Component() {}

		double Processroughness_curve_Dual(PolyhedronPtr P1, PolyhedronPtr p2, double radius, double maxdim,bool IsGauss = true);
		void ComputeDistanceEcartNormalRoughnessPonderate(PolyhedronPtr P1, PolyhedronPtr p2, double Param,double &L);
		double Processroughness_per_vertex_curve(PolyhedronPtr P,Vertex* pVertex,double radius,std::vector<double> &TabDistance,std::vector<Point3d> &TabPoint,double & moyenneRet,double dim,bool IsGauss = false);
		double ProcessCovariance(Vertex* pVertex,double moyenne,double moyenneDeg,std::vector<double> TabDistance,std::vector<Point3d> &TabPoint,std::vector<double> &TabDistanceDeg,std::vector<Point3d> &TabPointDeg, double dim,bool IsGauss);
		void KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef);
		double getMaxDim(PolyhedronPtr polyhedron_ptr);
		void ComputeMaxMin(PolyhedronPtr polyhedron_ptr);
		void ConstructColorMap(PolyhedronPtr pMesh);


private:
		double LUT_CourbureClust[3*256]; //table d'affichage couleur du bleu au rouge, en dégradé
		//les min et max respectif des courbures min et max (utile pour l'affichage)
		double MinMSDM,MaxMSDM;
		bool IsDistanceComputed;
	
			

	
};

#endif

#endif // MSDM_COMPONENT_H
