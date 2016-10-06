/*!
	\file MSDM2_Component.h
	\brief MSDM2 Perceptual distance calculation
	\brief According to: Multiscale Metric for 3D Mesh Visual Quality Assessment
			G. Lavoue
			In Symposium on Geometry Processing (SGP) 2011
	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2011
 */
  
#include <mepp_config.h>
//#ifdef BUILD_component_MSDM2

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "MSDM2_Component.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>

#include <QDir>
#include <QFileDialog>
#include <QMessageBox>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4267)
#include <CGAL/IO/Polyhedron_iostream.h>
#pragma warning(pop)
#else
#include <CGAL/IO/Polyhedron_iostream.h>
#endif

#if defined (_MSC_VER) || ((defined (__linux__) || defined (__APPLE__)) && (CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(3,8,0)))

typedef CGAL::Simple_cartesian<double>					AABB_Kernel;
typedef CGAL::AABB_polyhedron_triangle_primitive<AABB_Kernel,Polyhedron> AABB_Primitive;
typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive>								AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits>												AABB_Tree;
typedef AABB_Tree::Object_and_primitive_id Object_and_primitive_id;
typedef AABB_Tree::Point_and_primitive_id Point_and_primitive_id;

typedef AABB_Kernel::Segment_3 Segment;
typedef AABB_Kernel::Ray_3 Ray;



void MSDM2_Component::Matching_Multires_Init(PolyhedronPtr m_PolyDegrad, PolyhedronPtr m_PolyOriginal , Facet * _TabMatchedFacet)
	{
		
		 // constructs AABB tree
		 AABB_Tree tree(m_PolyOriginal->facets_begin(),m_PolyOriginal->facets_end());
		 tree.accelerate_distance_queries();

		//Searching for the closest point and facet for each vertex

		int ind=0;
		for(Vertex_iterator	pVertex	= m_PolyDegrad->vertices_begin();
					pVertex	!= m_PolyDegrad->vertices_end();
					pVertex++)
		{
			pVertex->MSDM2_Local=0;
			 // computes closest point and primitive id
			Point_and_primitive_id pp = tree.closest_point_and_primitive(pVertex->point());
			Point3d Nearest=pp.first;
			Facet_iterator f_Nearest = pp.second; // closest primitive id

			pVertex->match=Nearest;
			_TabMatchedFacet[ind]=*f_Nearest;

			ind++;

			
			
		}
		
}

void MSDM2_Component::Matching_Multires_Update(PolyhedronPtr m_PolyDegrad, Facet * _TabMatchedFacet)
	{
		
		
		int ind=0;
		for(Vertex_iterator	pVertex	= m_PolyDegrad->vertices_begin();
					pVertex	!= m_PolyDegrad->vertices_end();
					pVertex++)
		{
                        //Point3d Nearest=pVertex->match; // MT
			Facet* f_Nearest=&_TabMatchedFacet[ind];

			pVertex->tag(ind);

			ind++;

			//for debug
			//pVertex->point()=Nearest;

			///calculation of the nearest point curvature value using vertices of the Nearest triangle	
			//we use linear interpolation using barycentric coordinates
			Point3d x1=f_Nearest->halfedge()->vertex()->point();
			Point3d x2=f_Nearest->halfedge()->next()->vertex()->point();
			Point3d x3=f_Nearest->halfedge()->next()->next()->vertex()->point();

			double l1=sqrt((x3-x2)*(x3-x2));
			double l2=sqrt((x1-x3)*(x1-x3));
			double l3=sqrt((x1-x2)*(x1-x2));

			Vector v1=f_Nearest->halfedge()->vertex()->point()-pVertex->point();
			Vector v2=f_Nearest->halfedge()->next()->vertex()->point()-pVertex->point();
			Vector v3=f_Nearest->halfedge()->next()->next()->vertex()->point()-pVertex->point();

			double t1=sqrt(v1*v1);
			double t2=sqrt(v2*v2);
			double t3=sqrt(v3*v3);
			
			double p1=(l1+t2+t3)/2;
			double p2=(t1+l2+t3)/2;
			double p3=(t1+t2+l3)/2;

			double A1=(p1*(p1-l1)*(p1-t3)*(p1-t2));
			double A2=(p2*(p2-l2)*(p2-t3)*(p2-t1));
			double A3=(p3*(p3-l3)*(p3-t1)*(p3-t2));

			if(A1>0) A1=sqrt(A1); else	A1=0;
			if(A2>0) A2=sqrt(A2); else	A2=0;
			if(A3>0) A3=sqrt(A3); else	A3=0;
			
			double c1=f_Nearest->halfedge()->vertex()->KmaxCurv;
			double c2=f_Nearest->halfedge()->next()->vertex()->KmaxCurv;
			double c3=f_Nearest->halfedge()->next()->next()->vertex()->KmaxCurv;

			if((A1+A2+A3)>0)
				pVertex->curvmatch=(A1*c1+A2*c2+A3*c3)/(A1+A2+A3);
			else
				pVertex->curvmatch=(c1+c2+c3)/3;		


			

			
		}
		
}



void MSDM2_Component::ComputeStatistics(Vertex* pVertex, double Param,
	std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2,double radius, double dim)
{
	
	double moyenne1=0;
	double variance1=0;

	double moyenne2=0;
	double variance2=0;

	double covariance=0;

	////gaussian normalisation

	double variance = radius/2;
	double *tab_wi1=new double[TabPoint1.size()];
	double *tab_wi2=new double[TabPoint1.size()];
	double SommeDistance1=0;
	double SommeDistance2=0;
	double SommeWi1=0;
	double SommeWi2=0;
	for(unsigned int i=0;i<TabPoint1.size();i++) // MT
	{
		Vector DistancePt1=TabPoint1[i]-pVertex->point();
		double distPt1=sqrt(DistancePt1*DistancePt1);
		double wi1=1/variance/sqrt(2*3.141592)*exp(-(distPt1*distPt1)/2/variance/variance);
		
		tab_wi1[i]=wi1;
		SommeWi1+=wi1;
		SommeDistance1+=TabDistance1[i]*wi1;

		Vector DistancePt2=TabPoint2[i]-pVertex->match;
		double distPt2=sqrt(DistancePt2*DistancePt2);
		
		double wi2=1/variance/sqrt(2*3.141592)*exp(-(distPt2*distPt2)/2/variance/variance);
		tab_wi2[i]=wi2;
		SommeWi2+=wi2;
		SommeDistance2+=TabDistance2[i]*wi2;
	}
	moyenne1=SommeDistance1/(double)SommeWi1;
	moyenne2=SommeDistance2/(double)SommeWi2;

	
	for(unsigned int i=0;i<TabPoint1.size();i++) // MT
	{
		variance1+=tab_wi1[i]*pow(TabDistance1[i]-moyenne1,2);
		variance2+=tab_wi2[i]*pow(TabDistance2[i]-moyenne2,2);
		covariance+=tab_wi2[i]*(TabDistance1[i]-moyenne1)*(TabDistance2[i]-moyenne2);
	}
		
	variance1=variance1/SommeWi1;
	variance2=variance2/SommeWi2;
	variance1=sqrt(variance1);
	variance2=sqrt(variance2);

	covariance=covariance/SommeWi1;

	///we then compute the MSDM2_Local value
        //double C1=1; // MT
        //double C2=1; // MT
        //double C3=0.5; // MT
			double fact1=(fabs(moyenne1-moyenne2))/(std::max(moyenne1,moyenne2)+1);
			double fact2=(fabs(variance1-variance2))/(std::max(variance1,variance2)+1);
			double fact3=(fabs(variance1*variance2-covariance))/(variance1*variance2+1);
			
			pVertex->MSDM2_Local+=pow((pow(fact1,1)+pow(fact2,1)+Param*pow(fact3,1))/(2.+Param),1./1.);


	delete[] tab_wi1;
	delete[] tab_wi2;
}






	double MSDM2_Component::ProcessMSDM2_Multires(PolyhedronPtr m_PolyOriginal, PolyhedronPtr m_PolyDegrad,int NbLevel, double maxdim,double & FastMSDM, Curvature_ComponentPtr component_ptr_curvature)
{
	
	 
	double somme_MSDM2=0;
	int NbVertex=0;	

	Facet * TabMatchedFacet=new Facet[m_PolyDegrad->size_of_vertices()];

	Matching_Multires_Init(m_PolyDegrad, m_PolyOriginal,TabMatchedFacet);

	
	double RadiusCurvature=0.002;
	for (int i=0;i<NbLevel;i++)
	{
		m_PolyDegrad->set_index_vertices();

		component_ptr_curvature->principal_curvature(m_PolyOriginal,true,RadiusCurvature*maxdim);
		component_ptr_curvature->principal_curvature(m_PolyDegrad,true,RadiusCurvature*maxdim);
				
		KmaxKmean(m_PolyOriginal,maxdim);
		KmaxKmean(m_PolyDegrad,maxdim);

		Matching_Multires_Update(m_PolyDegrad,TabMatchedFacet);

		for(Vertex_iterator	pVertex	=	m_PolyDegrad->vertices_begin();
					pVertex	!= m_PolyDegrad->vertices_end();
					pVertex++)
		{
			

				std::vector<double>  TabDistance1;std::vector<double> TabDistance2;

				TabPoint1.clear();
				TabPoint2.clear();


				ProcessMSDM2_per_vertex(pVertex,RadiusCurvature*5*maxdim,TabDistance1,TabDistance2,TabPoint1,TabPoint2);
				ComputeStatistics((&(*pVertex)), 0.5,TabDistance1,TabDistance2,TabPoint1,TabPoint2,RadiusCurvature*5*maxdim,maxdim);	
				

		}


		RadiusCurvature+=0.001;

	}

	for(Vertex_iterator	pVertex	=	m_PolyDegrad->vertices_begin();
					pVertex	!= m_PolyDegrad->vertices_end();
					pVertex++)
		{
				somme_MSDM2+=pow(pVertex->MSDM2_Local/NbLevel,3);	
				NbVertex++;
	}
			

	FastMSDM=pow(somme_MSDM2/(double)NbVertex,0.33333);

	delete [] TabMatchedFacet;
	m_PolyDegrad->IsDistanceComputed=true;
	return 0;
}

	void MSDM2_Component::KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef)
	{
		
		for(Vertex_iterator	pVertex	=	polyhedron_ptr->vertices_begin();pVertex!= polyhedron_ptr->vertices_end();pVertex++)
		{
			double kmax=pVertex->KmaxCurv*coef;
			double kmin=fabs(pVertex->KminCurv*coef);

			pVertex->KmaxCurv=(kmax+kmin)/2.;

		}

	}

	void MSDM2_Component::ConstructColorMap(PolyhedronPtr P, int MetricOrHausdorff)
	{
		if(MetricOrHausdorff==1)
		{
			if(P->IsDistanceComputed==true)
			{

				double R;
				int indiceLut;
				Vertex_iterator pVertex = NULL;

				if(MaxMSDM2>MinMSDM2)
				{
				  for (pVertex = P->vertices_begin();
					  pVertex != P->vertices_end();
					  pVertex++)
				  {

					
						R=(pVertex->MSDM2_Local-MinMSDM2)/(MaxMSDM2-MinMSDM2)*255;
					
						indiceLut=floor(R);
						pVertex->color(LUT_CourbureClust[3*indiceLut],LUT_CourbureClust[3*indiceLut+1],LUT_CourbureClust[3*indiceLut+2]);

				  }
				}
			}
		}
		else
		{
			if(P->IsDistanceComputed==true)
			{

				double R;
				int indiceLut;
				Vertex_iterator pVertex = NULL;

				if(MaxMSDM2>MinMSDM2)
				{
				  for (pVertex = P->vertices_begin();
					  pVertex != P->vertices_end();
					  pVertex++)
				  {

						float d=sqrt((pVertex->point()-pVertex->match)*(pVertex->point()-pVertex->match));
						R=(d-MinMSDM2)/(MaxMSDM2-MinMSDM2)*255;
					
						indiceLut=floor(R);
						pVertex->color(LUT_CourbureClust[3*indiceLut],LUT_CourbureClust[3*indiceLut+1],LUT_CourbureClust[3*indiceLut+2]);

				  }
				}
			}

		}
	}
	


	void MSDM2_Component::ComputeMaxMin(PolyhedronPtr P, int MetricOrHausdorff)
	{
		if(MetricOrHausdorff==1)
		{
			if(P->IsDistanceComputed==true)
			{
				MinMSDM2=1000;
				MaxMSDM2=0;
				for(Vertex_iterator	pVertex	=	P->vertices_begin();pVertex!= P->vertices_end();pVertex++)
				{
					if(pVertex->MSDM2_Local>MaxMSDM2)
						MaxMSDM2=pVertex->MSDM2_Local;
					if(pVertex->MSDM2_Local<MinMSDM2)
						MinMSDM2=pVertex->MSDM2_Local;

				}
			}
		}
		else
		{
			if(P->IsDistanceComputed==true)
			{
				MinMSDM2=1000;
				MaxMSDM2=0;
				for(Vertex_iterator	pVertex	=	P->vertices_begin();pVertex!= P->vertices_end();pVertex++)
				{
					float d=sqrt((pVertex->point()-pVertex->match)*(pVertex->point()-pVertex->match));
					if(d>MaxMSDM2)
						MaxMSDM2=d;
					if(d<MinMSDM2)
						MinMSDM2=d;

				}
			}

		}
	

	}

	double MSDM2_Component::getMaxDim(PolyhedronPtr polyhedron_ptr)
	{
		polyhedron_ptr->compute_bounding_box();
		
		double max=polyhedron_ptr->xmax()-polyhedron_ptr->xmin();
		if(polyhedron_ptr->ymax()-polyhedron_ptr->ymin()>max)
			max = polyhedron_ptr->ymax()-polyhedron_ptr->ymin();
		if(polyhedron_ptr->zmax()-polyhedron_ptr->zmin()>max)
			max = polyhedron_ptr->zmax()-polyhedron_ptr->zmin();

		return max;

	}


	bool sphere_clip_vector_MSDM2(Point3d &O, double r,const Point3d &P, Vector &V)
    {

        Vector W = P - O ;
        double a = (V*V);
        double b = 2.0 * V * W ;
        double c = (W*W) - r*r ;
        double delta = b*b - 4*a*c ;

	

		if( a==0)
			return true ;

        if(delta < 0) {
            // Should not happen, but happens sometimes (numerical precision)
			
            return true ;
        }
        double t = (- b + std::sqrt(delta)) / (2.0 * a) ;
        if(t < 0.0) {
			
            // Should not happen, but happens sometimes (numerical precision)
            return true ;
        }
        if(t >= 1.0) {
            // Inside the sphere
            return false ;
        }

		if(t==0)
		{
			
			t=0.01;
		}

        V=V*t;

        return true ;
    }

	
	void MSDM2_Component::ProcessMSDM2_per_vertex( Vertex_iterator pVertex,double radius,std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2)
	{

		std::set<int> vertices ;

        std::stack<Vertex_iterator> S ;

		Point3d O = pVertex->point() ;

        S.push(pVertex) ;	
        vertices.insert(pVertex->tag()) ;
		
	
		TabDistance1.push_back(pVertex->KmaxCurv);
		TabPoint1.push_back(pVertex->point());

		TabDistance2.push_back(pVertex->curvmatch);
		TabPoint2.push_back(pVertex->match);

		int NbSommetInSphere=0;
                //double SommeDistance=0; // MT
	

        while(!S.empty())
		{
			Vertex_iterator v = S.top() ;
            S.pop() ;
            Point3d P = v->point() ;
            Halfedge_around_vertex_circulator h = v->vertex_begin();
			Halfedge_around_vertex_circulator pHalfedgeStart = h;
			CGAL_For_all(h,pHalfedgeStart)
			{
                Point3d p1 = h->vertex()->point();
				Point3d p2 = h->opposite()->vertex()->point();

				Point3d p1m = h->vertex()->match;
				Point3d p2m = h->opposite()->vertex()->match;

				Vector V = (p2-p1);
				Vector Vm = (p2m-p1m);

                if(v==pVertex || V * (P - O) > 0.0) 
				{
					double len_old = std::sqrt(V*V);
					bool isect = sphere_clip_vector_MSDM2(O, radius, P, V) ;
					double len_edge = std::sqrt(V*V);

					NbSommetInSphere++;
					
					
					double WeightedCurv1,WeightedCurv2;
					Point3d WeightedP1,WeightedP2;

					bool IsAlreadyIntegrated=false;
					if(!isect) 
					{
						
						Vertex_iterator w=h->opposite()->vertex();
                       if(vertices.find(w->tag()) == vertices.end())
						{
                            vertices.insert(w->tag()) ;
                            S.push(w) ;
                        }
					   else
						   IsAlreadyIntegrated=true;
                    }

					if (IsAlreadyIntegrated==false)
					{
						if(len_old!=0)
						{
							if(isect)
							{
								WeightedCurv1=(1-len_edge/len_old)*h->vertex()->KmaxCurv+len_edge/len_old*h->opposite()->vertex()->KmaxCurv;
								WeightedP1=p1+V;

								WeightedCurv2=(1-len_edge/len_old)*h->vertex()->curvmatch+len_edge/len_old*h->opposite()->vertex()->curvmatch;
								WeightedP2=p1m+(len_edge/len_old)*Vm;
							}
							else
							{
								WeightedCurv1=h->opposite()->vertex()->KmaxCurv;
								WeightedCurv2=h->opposite()->vertex()->curvmatch;
								WeightedP1=p2;
								WeightedP2=p2m;
							}
						}
						else
						{
							WeightedCurv1=h->opposite()->vertex()->KmaxCurv;
							WeightedCurv2=h->opposite()->vertex()->curvmatch;
							WeightedP1=p2;
							WeightedP2=p2m;
						}

						TabDistance1.push_back(WeightedCurv1);
						TabPoint1.push_back(WeightedP1);

						TabDistance2.push_back(WeightedCurv2);
						TabPoint2.push_back(WeightedP2);
					}

					

					
					
				}
                
			}
			
		}
	}
#endif

MSDM2_Component::MSDM2_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	int i=0;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.515600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.531300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.546900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.562500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.578100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.593800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.609400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.625000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.640600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.656300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.671900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.687500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.703100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.718800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.734400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.750000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.765600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.781300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.796900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.812500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.828100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.843800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.859400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.875000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.890600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.906300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.921900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.937500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.953100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.968800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.984400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.015600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.031300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.046900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.062500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.078100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.093800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.109400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.125000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.140600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.156300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.171900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.187500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.203100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.218800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.234400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.250000;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.265600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.281300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.296900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.312500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.328100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.343800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.359400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.375000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.390600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.406300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.421900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.437500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.453100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.468800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.484400;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.500000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.015600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.031300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.984400;		LUT_CourbureClust[i++]=	0.046900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.968800;		LUT_CourbureClust[i++]=	0.062500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.953100;		LUT_CourbureClust[i++]=	0.078100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.937500;		LUT_CourbureClust[i++]=	0.093800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.921900;		LUT_CourbureClust[i++]=	0.109400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.906300;		LUT_CourbureClust[i++]=	0.125000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.890600;		LUT_CourbureClust[i++]=	0.140600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.875000;		LUT_CourbureClust[i++]=	0.156300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.859400;		LUT_CourbureClust[i++]=	0.171900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.843800;		LUT_CourbureClust[i++]=	0.187500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.828100;		LUT_CourbureClust[i++]=	0.203100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.812500;		LUT_CourbureClust[i++]=	0.218800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.796900;		LUT_CourbureClust[i++]=	0.234400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.781300;
	LUT_CourbureClust[i++]=	0.250000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.765600;		LUT_CourbureClust[i++]=	0.265600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.750000;		LUT_CourbureClust[i++]=	0.281300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.734400;		LUT_CourbureClust[i++]=	0.296900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.718800;		LUT_CourbureClust[i++]=	0.312500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.703100;		LUT_CourbureClust[i++]=	0.328100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.687500;		LUT_CourbureClust[i++]=	0.343800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.671900;		LUT_CourbureClust[i++]=	0.359400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.656300;		LUT_CourbureClust[i++]=	0.375000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.640600;		LUT_CourbureClust[i++]=	0.390600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.625000;		LUT_CourbureClust[i++]=	0.406300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.609400;		LUT_CourbureClust[i++]=	0.421900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.593800;		LUT_CourbureClust[i++]=	0.437500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.578100;		LUT_CourbureClust[i++]=	0.453100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.562500;		LUT_CourbureClust[i++]=	0.468800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.546900;		LUT_CourbureClust[i++]=	0.484400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.531300;
	LUT_CourbureClust[i++]=	0.500000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.515600;		LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.500000;		LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.484400;		LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.468800;		LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.453100;		LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.437500;		LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.421900;		LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.406300;		LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.390600;		LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.375000;		LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.359400;		LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.343800;		LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.328100;		LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.312500;		LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.296900;		LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.281300;
	LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.265600;		LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.250000;		LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.234400;		LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.218800;		LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.203100;		LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.187500;		LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.171900;		LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.156300;		LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.140600;		LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.125000;		LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.109400;		LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.093800;		LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.078100;		LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.062500;		LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.046900;		LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.031300;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.015600;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.500000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.484400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.468800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.453100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.437500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.421900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.406300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.390600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.375000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.359400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.343800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.328100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.312500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.296900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.281300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.265600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.250000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.234400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.218800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.203100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.187500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.171900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.156300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.140600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.125000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.109400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.093800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.078100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.062500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.046900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.031300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.015600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;

	p->IsDistanceComputed=false;

	// MEPP 2
	componentName = "MSDM2_Component";
}
//#endif
