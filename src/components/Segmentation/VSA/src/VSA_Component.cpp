///////////////////////////////////////////////////////////////////////////
// Author: Guillaume Lavoué
// Year: 2008
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
//
//According to: Variational Shape Approximation
//		David Cohen-Steiner, Pierre Alliez and Mathieu Desbrun, SIGGRAPH '2004. 
///////////////////////////////////////////////////////////////////////////
#include <mepp_config.h>
#ifdef BUILD_component_VSA

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "VSA_Component.h"

#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>

 typedef struct
  {
	Facet_iterator Facet;
	//distance pour l'algo de LLoyd
	double DistanceLLoyd;
	//cluster visé pour l'algo de LLoyd
	double PossibleCluster;

  }
  FacetToIntegrate;

   typedef struct
  {
	Vertex_iterator Vertex;
	std::vector<int> TabAdjAnch;


  }
  AnchorVertex;

  double AreaFacetTriangleSeg(Facet_iterator &f)
	{
		Halfedge_around_facet_circulator pHalfedge = f->facet_begin();
		Point3d P = pHalfedge->vertex()->point();
		Point3d Q = pHalfedge->next()->vertex()->point();
		Point3d R = pHalfedge->next()->next()->vertex()->point();

		Vector PQ=Q-P;
		Vector PR=R-P;
		Vector QR=R-Q;


		Vector normal	=	CGAL::cross_product(PQ,QR);
		double area=0.5*sqrt(normal*normal);

		return area;

	}

  void VSA_Component::Init(int NbProxy)
	{
		////creation des proxy initiaux par tirage aléatoire de NbProxy triangles
		m_Table_Proxy.clear();

		m_NbProxy=NbProxy;
		int NbFacet=m_Poly->size_of_facets();
		int offset=NbFacet/NbProxy;
		int i=0; //int NumProx=1;
		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{

			if(i%offset==0)///on choisi le triangle
			{
				//on crée un proxy
				Proxy NewProxy;
				NewProxy.Normal=pface->normal();

				NewProxy.Seed=pface;
				///on ajoute le proxy à la liste
				m_Table_Proxy.push_back(NewProxy);


			}

			pface->LabelVSA=-1;
			i++;
		}

		m_Poly->NbFaceLabel=m_NbProxy;

	}

	struct CompFacet
	{
	bool operator()(FacetToIntegrate f1, FacetToIntegrate f2) const
	{
		return(f1.DistanceLLoyd<=f2.DistanceLLoyd);

	}
	};


	double VSA_Component::DistorsionError(Facet_iterator f,Proxy p)
	{
		Vector v=f->normal()-p.Normal;
		double nrm=v*v;
		double area=AreaFacetTriangleSeg(f);
		return nrm*area;
	}


VSA_Component::VSA_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	int i=0;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.515600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.531300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.546900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.562500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.578100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.593800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.609400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.625000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.640600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.656300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.671900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.687500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.703100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.718800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.734400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.750000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.765600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.781300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.796900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.812500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.828100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.843800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.859400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.875000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.890600;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.906300;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.921900;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.937500;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.953100;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.968800;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.984400;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.015600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.031300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.046900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.062500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.078100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.093800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.109400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.125000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.140600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.156300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.171900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.187500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.203100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.218800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.234400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.250000;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.265600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.281300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.296900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.312500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.328100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.343800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.359400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.375000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.390600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.406300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.421900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.437500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.453100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.468800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.484400;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.500000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	1.000000;
	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.015600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;		LUT_Seg[i++]=	0.031300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.984400;		LUT_Seg[i++]=	0.046900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.968800;		LUT_Seg[i++]=	0.062500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.953100;		LUT_Seg[i++]=	0.078100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.937500;		LUT_Seg[i++]=	0.093800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.921900;		LUT_Seg[i++]=	0.109400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.906300;		LUT_Seg[i++]=	0.125000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.890600;		LUT_Seg[i++]=	0.140600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.875000;		LUT_Seg[i++]=	0.156300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.859400;		LUT_Seg[i++]=	0.171900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.843800;		LUT_Seg[i++]=	0.187500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.828100;		LUT_Seg[i++]=	0.203100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.812500;		LUT_Seg[i++]=	0.218800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.796900;		LUT_Seg[i++]=	0.234400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.781300;
	LUT_Seg[i++]=	0.250000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.765600;		LUT_Seg[i++]=	0.265600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.750000;		LUT_Seg[i++]=	0.281300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.734400;		LUT_Seg[i++]=	0.296900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.718800;		LUT_Seg[i++]=	0.312500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.703100;		LUT_Seg[i++]=	0.328100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.687500;		LUT_Seg[i++]=	0.343800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.671900;		LUT_Seg[i++]=	0.359400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.656300;		LUT_Seg[i++]=	0.375000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.640600;		LUT_Seg[i++]=	0.390600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.625000;		LUT_Seg[i++]=	0.406300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.609400;		LUT_Seg[i++]=	0.421900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.593800;		LUT_Seg[i++]=	0.437500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.578100;		LUT_Seg[i++]=	0.453100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.562500;		LUT_Seg[i++]=	0.468800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.546900;		LUT_Seg[i++]=	0.484400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.531300;
	LUT_Seg[i++]=	0.500000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.515600;		LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.500000;		LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.484400;		LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.468800;		LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.453100;		LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.437500;		LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.421900;		LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.406300;		LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.390600;		LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.375000;		LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.359400;		LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.343800;		LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.328100;		LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.312500;		LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.296900;		LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.281300;
	LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.265600;		LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.250000;		LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.234400;		LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.218800;		LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.203100;		LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.187500;		LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.171900;		LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.156300;		LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.140600;		LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.125000;		LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.109400;		LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.093800;		LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.078100;		LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.062500;		LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.046900;		LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.031300;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.015600;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.500000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.484400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.468800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.453100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.437500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.421900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.406300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.390600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.375000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.359400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.343800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.328100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.312500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.296900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.281300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.265600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.250000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.234400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.218800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.203100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.187500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.171900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.156300;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.140600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.125000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.109400;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.093800;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.078100;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.062500;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.046900;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.031300;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.015600;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	1.000000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.984400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.968800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.953100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.937500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.921900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.906300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.890600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.875000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.859400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.843800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.828100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.812500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.796900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.781300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	0.765600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.750000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.734400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.718800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.703100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.687500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.671900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.656300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.640600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.625000;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.609400;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.593800;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.578100;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.562500;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.546900;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;		LUT_Seg[i++]=	0.531300;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;
	LUT_Seg[i++]=	0.515600;	LUT_Seg[i++]=	0.000000;	LUT_Seg[i++]=	0.000000;

	// from IHM
	displayFaceLabel=false;

	// MEPP 2
	componentName = "VSA_Component";
	init = 1;
}

void VSA_Component::Flooding()///repartition des triangles dans leur proxy respectif
	{
		m_Poly->NbFaceLabel=m_NbProxy;
		typedef std::set<FacetToIntegrate,CompFacet> ListFacet_model;
		typedef std::set<FacetToIntegrate,CompFacet>::iterator ListFacet_iterator;

		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{
			pface->LabelVSA=-1;
		}


		ListFacet_model ListFacet;
		for(int i=0;i<m_NbProxy;i++)//pour chaque proxy on initialise les triangle qui pourrait s'étendre
		{
			m_Table_Proxy[i].TabAdj.clear();
			Facet_iterator f=m_Table_Proxy[i].Seed;
			f->LabelVSA=i;


			//on extrait les trois triangle
			Facet_iterator ff1, ff2, ff3;
			FacetToIntegrate f1, f2, f3;

			Halfedge_around_facet_circulator pHalfedge = f->facet_begin();

			if(!pHalfedge->opposite()->is_border())
			{
				ff1 = pHalfedge->opposite()->facet();
				f1.Facet=ff1;
				f1.PossibleCluster=i;
				f1.DistanceLLoyd=DistorsionError(f1.Facet,m_Table_Proxy[i]);
				ListFacet.insert(f1);

			}
			if(!pHalfedge->next()->opposite()->is_border())
			{
				ff2 = pHalfedge->next()->opposite()->facet();
				f2.Facet=ff2;
				f2.PossibleCluster=i;
				f2.DistanceLLoyd=DistorsionError(f2.Facet,m_Table_Proxy[i]);
				ListFacet.insert(f2);
			}


			if(!pHalfedge->next()->next()->opposite()->is_border())
			{
				ff3 = pHalfedge->next()->next()->opposite()->facet();
				f3.Facet=ff3;
				f3.PossibleCluster=i;
				f3.DistanceLLoyd=DistorsionError(f3.Facet,m_Table_Proxy[i]);
				ListFacet.insert(f3);
			}
		}

		ListFacet_iterator it;
		for(it=ListFacet.begin();it!=ListFacet.end();)
		{
			if(it->Facet->LabelVSA==-1)
			{
				it->Facet->LabelVSA=it->PossibleCluster;
				//ensuite on met ses triangles adjacents dans la queue


				//on extrait les trois triangle
				Facet_iterator ff1, ff2, ff3;

				Halfedge_around_facet_circulator pHalfedge = it->Facet->facet_begin();


				FacetToIntegrate f1, f2, f3;
				if(!pHalfedge->opposite()->is_border())
				{
					ff1 = pHalfedge->opposite()->facet();
					if(ff1->LabelVSA==-1 )
					{
						f1.Facet=ff1;
						f1.PossibleCluster=it->PossibleCluster;
						f1.DistanceLLoyd=DistorsionError(f1.Facet,m_Table_Proxy[it->PossibleCluster]);
						ListFacet.insert(f1);
					}
				}
				if(!pHalfedge->next()->opposite()->is_border())
				{
					ff2 = pHalfedge->next()->opposite()->facet();
					if(ff2->LabelVSA==-1)
					{
						f2.Facet=ff2;
						f2.PossibleCluster=it->PossibleCluster;
						f2.DistanceLLoyd=DistorsionError(f2.Facet,m_Table_Proxy[it->PossibleCluster]);
						ListFacet.insert(f2);
					}
				}
				if(!pHalfedge->next()->next()->opposite()->is_border())
				{
					ff3 = pHalfedge->next()->next()->opposite()->facet();
					if(ff3->LabelVSA==-1)
					{
						f3.Facet=ff3;
						f3.PossibleCluster=it->PossibleCluster;
						f3.DistanceLLoyd=DistorsionError(f3.Facet,m_Table_Proxy[it->PossibleCluster]);
						ListFacet.insert(f3);
					}
				}
			}
			ListFacet.erase(it);
			it=ListFacet.begin();


		}

		///on intégre les ionformation de connectivité des proxy
		for(Halfedge_iterator pHalfedge	=	m_Poly->halfedges_begin();
				pHalfedge	!= m_Poly->halfedges_end();
				pHalfedge++)
		{
			if(pHalfedge->is_border()||pHalfedge->opposite()->is_border())
				continue;

			int Label1=pHalfedge->facet()->LabelVSA;
			int Label2=pHalfedge->opposite()->facet()->LabelVSA;
			if(Label1!=Label2)
			{

				bool IsFound=false;
				for(unsigned int i=0;i<m_Table_Proxy[Label1].TabAdj.size();i++)
					if(m_Table_Proxy[Label1].TabAdj[i]==Label2)
						IsFound=true;
				if(IsFound==false)
				{
					m_Table_Proxy[Label1].TabAdj.push_back(Label2);
					m_Table_Proxy[Label2].TabAdj.push_back(Label1);

				}


			}


		}



	}

	void VSA_Component::ProxyFitting()
	{
		Vector * TabNormal=new Vector[m_NbProxy];
		double * TabArea=new double[m_NbProxy];

		double * DistanceMin=new double[m_NbProxy];
		double * DistanceMax=new double[m_NbProxy];

		for (int i=0;i<m_NbProxy;i++)
		{
			TabArea[i]=0;
			TabNormal[i]=CGAL::NULL_VECTOR;
			DistanceMin[i]=100000000;
			DistanceMax[i]=0;



		}
		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{
			double area=AreaFacetTriangleSeg(pface);
			TabArea[pface->LabelVSA]+=area;
			TabNormal[pface->LabelVSA]=TabNormal[pface->LabelVSA]+pface->normal()*area;

		}

		for (int i=0;i<m_NbProxy;i++)
		{

			m_Table_Proxy[i].Normal=TabNormal[i]/TabArea[i];
			m_Table_Proxy[i].Area=TabArea[i];
			m_Table_Proxy[i].TotalDistorsion=0;
			//on se fout du center pour le moment
		}



		// ensuite on assigne une nouvelle seed, à chaque proxy
		for(Facet_iterator	pface	=	m_Poly->facets_begin();
				pface	!= m_Poly->facets_end();
				pface++)
		{
			double distance=DistorsionError(pface,m_Table_Proxy[pface->LabelVSA]);
			m_Table_Proxy[pface->LabelVSA].TotalDistorsion+=distance;
			if(distance<DistanceMin[pface->LabelVSA])
			{

				DistanceMin[pface->LabelVSA]=distance;
				m_Table_Proxy[pface->LabelVSA].Seed=pface;
			}

			//on repère la MostDistordedFacet
			if(distance>DistanceMax[pface->LabelVSA])
			{

				DistanceMax[pface->LabelVSA]=distance;
				m_Table_Proxy[pface->LabelVSA].MostDistordedFacet=pface;
			}


		}

		delete []DistanceMin;
		delete []TabNormal;
		delete []TabArea;
		delete []DistanceMax;



	}


	void VSA_Component::ProxyInsertion()
	{
		int NumProxMax=0;
		double DistorsionMax=0;

		EvalInsertion(NumProxMax,DistorsionMax);
		CreateNewProxy( NumProxMax, DistorsionMax);

	}

	void VSA_Component::EvalInsertion(int & NumProxMax,double & DistorsionMax)
	{
		for (int i=0;i<m_NbProxy;i++)
		{
			//int a=	m_Table_Proxy[i].TotalDistorsion;
			//int areass=m_Table_Proxy[i].Area;
			if(	m_Table_Proxy[i].TotalDistorsion>DistorsionMax)
			{
				NumProxMax=i;
				DistorsionMax=m_Table_Proxy[i].TotalDistorsion;

			}
		}

	}

	void VSA_Component::CreateNewProxy(int NumProxMax,double DistorsionMax)
	{
		Proxy NewProxy;
		if(m_Table_Proxy[NumProxMax].MostDistordedFacet!=m_Table_Proxy[NumProxMax].Seed)
			NewProxy.Seed=m_Table_Proxy[NumProxMax].MostDistordedFacet;
		else
		{
			m_Table_Proxy[NumProxMax].MostDistordedFacet++;
			NewProxy.Seed=m_Table_Proxy[NumProxMax].MostDistordedFacet;
		}


		NewProxy.Normal=m_Table_Proxy[NumProxMax].MostDistordedFacet->normal();
		m_NbProxy++;
		m_Table_Proxy.push_back(NewProxy);
		NewProxy.Seed->tag(15);
	}


	void VSA_Component::Variational_SegmentationIncr(PolyhedronPtr pMesh, int NbRegions, int NbIter)
	{
		m_Poly=pMesh;
		Init(2);
		Flooding();

		for(int i=0;i<NbRegions-2;i++)
		{
			for(int j=0;j<NbIter;j++)
			{
				ProxyFitting();
				Flooding();
			}
			ProxyFitting();
			ProxyInsertion();
			Flooding();

		}
		for(int j=0;j<NbIter;j++)
		{
			ProxyFitting();
			Flooding();
		}

	}

	void VSA_Component::Variational_Segmentation(PolyhedronPtr pMesh, int NbRegions, int NbIter)
	{
		m_Poly=pMesh;
		Init(NbRegions);
		Flooding();


		for(int i=0;i<NbIter;i++)
		{
			ProxyFitting();
			Flooding();
		}

	}

	void VSA_Component::ConstructFaceColorMap(PolyhedronPtr pMesh)
{
		//double R;
		//int indiceLut;
		Vertex_iterator pVertex = NULL;

		Facet_iterator pFacet	=	pMesh->facets_begin();
		for(;pFacet	!= pMesh->facets_end();pFacet++)
		{

			double R=(double)(pFacet->LabelVSA)/(double)pMesh->NbFaceLabel*255.;
			int indiceLut=floor(R);

			pFacet->color(LUT_Seg[3*indiceLut],LUT_Seg[3*indiceLut+1],LUT_Seg[3*indiceLut+2]);

		}
}
#endif
