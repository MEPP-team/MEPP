/*!
	\file MSDM_Component.h
	\brief MSDM Perceptual distance calculation
	\brief According to: According to: Perceptually driven 3d distance metrics with application to watermarking
			G. Lavoue, E. Drelie Gelasca, F. Dupont, A. Baskurt, and T. Ebrahimi 
			In SPIE Applications of Digital Image Processing XXIX, vol. 6312, 2006.
	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2010
 */
 
#ifndef MSDM_COMPONENT_H
#define MSDM_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#include "../../../../mepp/mepp_component.h"
#include "MSDM_Polyhedron.h"

/**
 \class	MSDM_Component

 \brief	MSDM component. 

 */
class MSDM_Component : 
  public mepp_component
{
	public:
			/*!
		* \brief Constructor
		*/
		MSDM_Component(Viewer* v, PolyhedronPtr p);

		/*!
		* \brief Destructor
		*/
		~MSDM_Component() {}

		/**
		 \fn	double MSDM_Component::ComputeLocalCurvatureStatistics(PolyhedronPtr P1, PolyhedronPtr p2,
		 		double radius, double maxdim);
		
		 \brief	Calculates local curvature statistics for each vertex withing a given circular neighborhood
				
		 \param	P1	  	The first Polyhedron.
		 \param	p2	  	The second Polyhedron.
		 \param	radius	The radius of the neighborhood
		 \param	maxdim	The max dimension of the Bounding Box Length.
		
		
		 */
		void ComputeLocalCurvatureStatistics(PolyhedronPtr P1, PolyhedronPtr p2, double radius, double maxdim);

		/**
		 \fn	void MSDM_Component::ComputeMSDM_FromStatistics(PolyhedronPtr P1, PolyhedronPtr p2,
		 		double Param,double &L);
		
		 \brief	Calculates the MSDM from statistics.
		
		
		 \param	P1			 	The first Polyhedron.
		 \param	p2			 	The second Polyhedron.
		 \param [out]	L	Global MSDM Value
		 */
		void ComputeMSDM_FromStatistics(PolyhedronPtr P1, PolyhedronPtr p2, double &L);

		/**
		 \fn	double MSDM_Component::ComputeLocalCurvatureStatistics_PerVertex(PolyhedronPtr P,
		 		Vertex* pVertex,double radius,std::vector<double> &TabDistance,std::vector<Point3d> &TabPoint,
		 		double & moyenneRet,double dim);
		
		 \brief	Calculates the local curvature statistics per vertex.
		
		
		 \param	P					   	The polyhedron.
		 \param pVertex					The vertex.
		 \param	radius				   	The radius of the neighborhood.
		 \param [out]	TabDistance		The set of curvature values of the neighborhood vertices.
		 \param [out]	TabPoint   	The set of 3D positions of the neighborhood vertices.
		 \param [out]	moyenneRet 	Mean of the curvature values over the neighborhood vertices.
		 \param	dim					   	The max dimension of the Bounding Box Length.
		
		
		 */

		void ComputeLocalCurvatureStatistics_PerVertex(PolyhedronPtr P,Vertex* pVertex,double radius,std::vector<double> &TabDistance,std::vector<Point3d> &TabPoint,double & moyenneRet,double dim);

		/**
		 \fn	double MSDM_Component::ProcessCovariance_PerVertex(Vertex* pVertex,double moyenne,
		 		double moyenneDeg,std::vector<double> TabDistance,std::vector<Point3d> &TabPoint,
		 		std::vector<double> &TabDistanceDeg,std::vector<Point3d> &TabPointDeg, double dim);
		
		 \brief	Process the covariance per vertex.
			
		 \param 	pVertex		  The vertex.
		 \param	moyenne					  	The mean curvature over the the neighborhood for the original mesh.
		 \param	moyenneDeg				  	The mean curvature over the the neighborhood for the distorded mesh.
		 \param	TabDistance				  	The set of curvature values of the neighborhood vertices for the original mesh.
		 \param 	TabPoint	  	The set of 3D positions of the neighborhood vertices for the original mesh.
		 \param 	TabDistanceDeg	The set of curvature values of the neighborhood vertices for the distorded mesh.
		 \param 	TabPointDeg   	The set of 3D positions of the neighborhood vertices for the distorded mesh.
		 \param	dim						  		The max dimension of the Bounding Box Length.
		
		 \return	The covariance value over the neighborhood.
		 */
		double ProcessCovariance_PerVertex(Vertex* pVertex,double moyenne,double moyenneDeg,std::vector<double> &TabDistance,std::vector<Point3d> &TabPoint,std::vector<double> &TabDistanceDeg,std::vector<Point3d> &TabPointDeg, double dim);

		/**
		 \fn	void MSDM_Component::KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef);
		
		 \brief	Computes the mean curvature field (kmin+kmax)/2 and normalize it according to the size of the model
				
		 \param	polyhedron_ptr	The polyhedron.
		 \param	coef		  	The normalization coef.
		 */

		void KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef);

		/**
		 \fn	double MSDM_Component::getMaxDim(PolyhedronPtr polyhedron_ptr);
		
		 \brief	Gets the maximum dimension of the object bounding box length.
				
		 \param	polyhedron_ptr	The polyhedron.
		
		 \return	The maximum dimension.
		 */

		double getMaxDim(PolyhedronPtr polyhedron_ptr);

		/**
		 \fn	void MSDM_Component::ComputeMaxMin(PolyhedronPtr polyhedron_ptr);
		
		 \brief	Calculates the maximum and minimum local MSDM values for rendering MSDM color map
		
		
		 \param	polyhedron_ptr	The polyhedron.
		 */

		void ComputeMaxMin(PolyhedronPtr polyhedron_ptr);

		/*!
		* \brief This method map the local MSDM scalar field into vertex colors
		* \param pMesh : input polyhedra
		*/	
		void ConstructColorMap(PolyhedronPtr pMesh);


private:
		/*! \brief Look-up table for color rendering (from blue to red)*/
		double LUT_CourbureClust[3*256]; 
		/*! \brief Minimum and maximum values of the local MSDM field (usefull for color rendering)*/
		double MinMSDM,MaxMSDM;
		
		bool IsDistanceComputed;///< true if the distance has been computed
	
			

	
};

#endif

#endif // MSDM_COMPONENT_H
