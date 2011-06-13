/*!
	\file MSDM2_Component.h
	\brief MSDM2 Perceptual distance calculation
	\brief According to: Multiscale Metric for 3D Mesh Visual Quality Assessment
			G. Lavoue
			In Symposium on Geometry Processing (SGP) 2011
	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2011
 */

#ifndef MSDM2_COMPONENT_H
#define MSDM2_COMPONENT_H

#include <mepp_config.h>
//#ifdef BUILD_component_MSDM2

#include "../../../../mepp/mepp_component.h"
#include "MSDM2_Polyhedron.h"


#include "../../../Analysis/Curvature/src/Curvature_Component.h"

typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;

class MSDM2_Component : 
  public mepp_component
{
	public:
			/*!
		* \brief Constructor
		*/
		MSDM2_Component(Viewer* v, PolyhedronPtr p);

		/*!
		* \brief Destructor
		*/
		~MSDM2_Component() {}

		/**
		 \fn	double MSDM2_Component::ProcessMSDM2_Multires(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal,int NbLevel,
				double maxdim,double & FastMSDM, Curvature_ComponentPtr component_ptr_curvature);
		
		
		 \brief	Computes the multiscale MSDM2 metric
				
		 \param	m_PolyDegrad	  	The first Polyhedron.
		 \param	m_PolyOriginal	  	The second Polyhedron.
		 \param	NbLevel		Number of scales used
		 \param	maxdim	The max dimension of the Bounding Box Length.
		 \param [out] MSDM2Value	The computed value of the MSDM2 metric
		
		
		 */
		double ProcessMSDM2_Multires(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal,int NbLevel, double maxdim,double & MSDM2Value, Curvature_ComponentPtr component_ptr_curvature);
		
		/**
		 \fn	void MSDM2_Component::Matching_Multires_Init(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal, 
				double Length,Facet * _TabMatchedFacet);
				
		
		 \brief	Initialize the matching process
				
		 \param	m_PolyDegrad	  	The first Polyhedron.
		 \param	m_PolyOriginal	  	The second Polyhedron.
		 \param	[out] _TabMatchedFacet	Facets from m_PolyOriginal on which vertices from m_PolyDegrad are projected
				
		
		 */		 
		void Matching_Multires_Init(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal, Facet * _TabMatchedFacet);
		
		/**
		 \fn	void MSDM2_Component::Matching_Multires_Update(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal, 
				double Length,Facet * _TabMatchedFacet);
				
		
		 \brief Updates the matching process
				
		 \param	m_PolyDegrad	  	The first Polyhedron.
		 \param	_TabMatchedFacet	Facets from m_PolyOriginal on which vertices from m_PolyDegrad are projected
				
		
		 */	
		void Matching_Multires_Update(PolyhedronPtr m_PolyDegrad,Facet * _TabMatchedFacet);

		/**
		 \fn	void MSDM_Component::ProcessMSDM2_per_vertex(Vertex_iterator pVertex,double radius,
				std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2);
	
		 \brief	Computes the local neighborhoods
		
		
		 \param	pVertex	The considered vertex
		 \param radius : radius of the neighborhood
		 \param [out] TabDistance1 : Curvature values from the neighborhood of pVertex regarding the first polyhedron
		 \param [out] TabDistance2 : Curvature values from the neighborhood of pVertex regarding the second polyhedron
		 \param [out] TabPoint1 : 3D points from the neighborhoodof pVertex regarding the first polyhedron
		 \param [out] TabPoint2 : 3D points from the neighborhood of pVertex regarding the second polyhedron
		 
		 */
		void ProcessMSDM2_per_vertex(Vertex_iterator pVertex,double radius,std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2);
	
		/**
		 \fn	double MSDM_Component:: ComputeStatistics(Vertex* pVertex, double Param,std::vector<double> & TabDistance1,
				 std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2,
				 double radius,double dim);
		
		 \brief	Calculates the local curvature statistics per vertex.
		
		
		 \param	pVertex	The considered vertex
		
		 \param TabDistance1 : Curvature values from the neighborhood of pVertex regarding the first polyhedron
		 \param TabDistance2 : Curvature values from the neighborhood of pVertex regarding the second polyhedron
		 \param TabPoint1 : 3D points from the neighborhoodof pVertex regarding the first polyhedron
		 \param TabPoint2 : 3D points from the neighborhood of pVertex regarding the second polyhedron
		 \param	radius				   	The radius of the neighborhood.
		 \param	dim	The max dimension of the Bounding Box Length.
		
		
		 */
		void ComputeStatistics(Vertex* pVertex, double Param,std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2,double radius,double dim);
		
		/**
		 \fn	void MSDM_Component::ComputeMaxMin(PolyhedronPtr polyhedron_ptr);
		
		 \brief	Calculates the maximum and minimum local MSDM values for rendering MSDM color map
		
		
		 \param	polyhedron_ptr	The polyhedron.
		 \param MetricOrHausdorff : 0 Hausdorff scalar field , 1 MSDM2 scalar field
		 */
		void ComputeMaxMin(PolyhedronPtr P, int MetricOrHausdorff);
		
		/*!
		* \brief This method map the local MSDM2 scalar field into vertex colors
		* \param pMesh : input polyhedra
		* \param MetricOrHausdorff : 0 Hausdorff scalar field , 1 MSDM2 scalar field
		*/	
		void ConstructColorMap(PolyhedronPtr P, int MetricOrHausdorff);
		
		/**
		 \fn	void MSDM2_Component::KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef);
		
		 \brief	Computes the mean curvature field (kmin+kmax)/2 and normalize it according to the size of the model
				
		 \param	polyhedron_ptr	The polyhedron.
		 \param	coef		  	The normalization coef.
		 */
		void KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef);
		
		/**
		 \fn	double MSDM2_Component::getMaxDim(PolyhedronPtr polyhedron_ptr);
		
		 \brief	Gets the maximum dimension of the object bounding box length.
				
		 \param	polyhedron_ptr	The polyhedron.
		
		 \return	The maximum dimension.
		 */
		double getMaxDim(PolyhedronPtr polyhedron_ptr);

		
		/*! \brief Sets of 3D points of local windows*/
		std::vector<Point3d> TabPoint1;std::vector<Point3d> TabPoint2;
	


private:
		/*! \brief Look-up table for color rendering (from blue to red)*/
		double LUT_CourbureClust[3*256]; 
		
		/*! \brief Minimum and maximum values of the local MSDM field (usefull for color rendering)*/
		double MinMSDM2,MaxMSDM2;

		
		
		
	
			

	
};

#endif

//#endif // MSDM2_COMPONENT_H
