///////////////////////////////////////////////////////////////////////////
// Author: Guillaume Lavoué
// Year: 2008
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
//
//According to: Variational Shape Approximation
//		David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004. 
///////////////////////////////////////////////////////////////////////////

#ifndef VSA_COMPONENT_H
#define VSA_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "../../../../mepp/mepp_component.h"
#include "VSA_Polyhedron.h"

class Viewer;

class VSA_Component : 
  public mepp_component
{
	public:

		typedef struct
		  {
			//Point Center;
			Vector Normal;
			Facet_iterator Seed;
			double Area;
			std::vector<int> TabAdj;
			double TotalDistorsion;
			Facet_iterator MostDistordedFacet;
			int Label;
		  }
		  Proxy;

		VSA_Component(Viewer* v, PolyhedronPtr p);
		~VSA_Component() {}

		void Variational_Segmentation(PolyhedronPtr pMesh,int NbRegions, int NbIter);//algorithme classique
		//on "seme" des germes et on régularise itérativement les régions
		void Variational_SegmentationIncr(PolyhedronPtr pMesh,int NbRegions, int NbIter);//algorithme icrémental
		//ajout incrémental de régions, plus long mais plus précis

		///pour le détails des fonctions suivantes, voir l'article original:
		//Variational Shape Approximation
		//David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004. 

		void Init(int NbProxy);
		void Flooding();
		void ProxyFitting();
		double DistorsionError(Facet_iterator f,Proxy p);
		void ProxyInsertion();
		void EvalInsertion(int & NumProxMax,double & DistorsionMax);
		void CreateNewProxy(int NumProxMax,double DistorsionMax);

		void ConstructFaceColorMap(PolyhedronPtr pMesh);//creation de la table couleur pour l'affichage des régions

	private:

		int m_NbProxy;//Le nombre de régions
		PolyhedronPtr m_Poly;
		std::vector<Proxy> m_Table_Proxy;//tableau des proxy (i.e. regions)
		double LUT_Seg[3*256]; //table d'affichage de couleurs aléatoires

	// from IHM
	private:
		bool displayFaceLabel;//est-ce qu'on affiche les direction de courbure minimum
};

#endif

#endif // VSA_COMPONENT_H
