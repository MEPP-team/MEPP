/*!
	\file Curvature_Component.h
	\brief Curvature computation on polyhedra
	\brief According to: Restricted Delaunay Triangulations and Normal Cycle, David Cohen-Steiner and J.M. Morvan In Proceedings SoCG'03.

	Most of this implementation was originally made by Pierre Alliez and David Cohen-Steiner (C) 2002.
	Thank also to Bruno Levy for his implementation of the geodesic radius
	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2008
 */


#ifndef Curvature_COMPONENT_H
#define Curvature_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "../../../../mepp/mepp_component.h"

#include "Curvature_Polyhedron.h"
#include "extract_Vpropres.h"

/*!
 * \class Curvature_Component
 * \brief Curvature computation on polyhedra
 */
class Curvature_Component : 
  public mepp_component
{
	public:

		/*!
		* \brief Constructor
		*/
		Curvature_Component(Viewer* v, PolyhedronPtr p);

		/*!
		* \brief Destructor
		*/
		~Curvature_Component() {}

		/*!
		 * \brief This method computes the curvature (max, min and directions)
		 * \param pMesh : input polyhedra
		 * \param IsGeod : Using geodesic (true)  or 1-ring (false) neighborhood
		 * \param radius : radius for the geodedisc neighborhood (% BBL)
		 */		 
		void principal_curvature(PolyhedronPtr pMesh,bool IsGeod,double radius);

		/*!
		* \brief This method map the curvature scalar fields (min or max) into vertex colors
		* \param pMesh : input polyhedra
		* \param ColorField : the scalar field to display (0 nothing, 1 min curvature, 2 max curvature.)
		*/	
		void ConstructColorMap(PolyhedronPtr pMesh,int ColorField);

		/*!
		* \brief This method activate or desactivate the rendering of the minimum curvature directions
		* \param min : activate (true) or desactivate (false) the rendering
		*/	
		void set_displayMinDirections(bool min) { displayMinDirections = min; }

		/*!
		* \brief This method tells us if the rendering is activated or not
		*/	
		bool get_displayMinDirections() { return displayMinDirections; }

		/*!
		* \brief This method activate or desactivate the rendering of the maximum curvature directions
		* \param max : activate (true) or desactivate (false) the rendering
		*/	
		void set_displayMaxDirections(bool max) { displayMaxDirections = max; }

		/*!
		* \brief This method tells us if the rendering is activated or not
		*/	
		bool get_displayMaxDirections() { return displayMaxDirections; }

	private:

		/*! \brief Minimum and maximum values of the minimum and maximum curvature fields (usefull for color rendering)*/
		double MinNrmMinCurvature,MaxNrmMinCurvature,MinNrmMaxCurvature,MaxNrmMaxCurvature;

		/*! \brief Look-up table for color rendering (from blue to red)*/
		double LUT_CourbureClust[3*256]; 

	// from IHM
	private:

		/*! \brief Boolean indicating activation or deactivation of the minimum curvature direction rendering*/
		bool displayMinDirections; 
		/*! \brief Boolean indicating activation or deactivation of the maximum curvature direction rendering*/
		bool displayMaxDirections; 
};

#endif

#endif // Curvature_COMPONENT_H
