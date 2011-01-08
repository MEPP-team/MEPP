///////////////////////////////////////////////////////////////////////////
// Author: Kai Wang
// Year: 2009
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef Various_Processing_COMPONENT_H
#define Various_Processing_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "../../../../mepp/mepp_component.h"
#include "Various_Processing_Polyhedron.h"

#include "Various_Processing_Types.h"

class Viewer;

class Various_Processing_Component : 
  public mepp_component
{
	public:
		Various_Processing_Component(Viewer* v, PolyhedronPtr p);
		~Various_Processing_Component() {}

		// General function for noise addition to vertex coordinates, including uniform noise and Gaussian noise
		// The noise intensity is relative to the average distance from vertices to mesh centre
		void NoiseAddition (PolyhedronPtr pMesh, Noise_type noiseType, double noiseIntensity, bool preserveBoundaries);

		// Uniform noise addition
		void NoiseAdditionUniform (PolyhedronPtr pMesh, double noiseIntensity, bool preserveBoundaries);

		// Gaussian noise addition. The noise intensity stands for the standard deviation of the Gaussian distribution.
		void NoiseAdditionGaussian (PolyhedronPtr pMesh, double noiseIntensity, bool preserveBoundaries);

		// For more information about Laplacian smoothing, please refer to the state of art report of Taubin in Eurographics 2000
		void LaplacianSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries);

		// Vertex coordinates are quantized to 2^(bitDepth) levels
		// TODO: fixing after quantization the potentially introduced degeneracies, such as removing the null surface facet
		void CoordinateQuantization (PolyhedronPtr pMesh, int bitDepth);


		void Translation (PolyhedronPtr pMesh, double xTranslation, double yTranslation, double zTranslation);
		
		// (xAxis,yAxis,zAxis) forms the rotation axis
		// angle stands for the counterclockwise rotation angle, seen from the positive side of the rotation axis
		void Rotation (PolyhedronPtr pMesh, double xAxis, double yAxis, double zAxis, double angle);

		void UniformScaling (PolyhedronPtr pMesh, double scalingFactor);

		// Different subdivision methods. For more information, please refer to the CGAL user manual.
		void Subdivision (PolyhedronPtr pMesh, Subdivision_type subdivisionType, int depth);
		void SubdivisionCatmullClark (PolyhedronPtr pMesh, int depth);
		void SubdivisionLoop (PolyhedronPtr pMesh, int depth);
		void SubdivisionDooSabin (PolyhedronPtr pMesh, int depth);
		void SubdivisionSqrt3 (PolyhedronPtr pMesh, int depth);
		void SubdivisionMidpoint (PolyhedronPtr pMesh, int depth);

		// Mesh simplification based on edge collasp. For more information, please refer to the CGAL user manual.
		void Simplification (Polyhedron *pMesh, int targetEdgeNum);
};

#endif

#endif // Various_Processing_COMPONENT_H
