///////////////////////////////////////////////////////////////////////////
// Author: Kai Wang
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////

// as you know, generally I do not write good programs...
// so please send me any of your comments (or rather critics) in order to improve my programming skill and this component.
// Thank you!

#ifndef Processing_COMPONENT_H
#define Processing_COMPONENT_H

#include "../../../../mepp/mepp_component.h"
#include "Compression_Valence_Polyhedron.h"
#include "Processing_Types.h"


class Processing_Component
{
	public:

		Processing_Component() {}
		~Processing_Component() {}

		// General function for noise addition to vertex coordinates, including uniform noise and Gaussian noise
		// The noise intensity is relative to the average distance from vertices to mesh centre
		void NoiseAddition (Polyhedron *pMesh, Noise_type noiseType, double noiseIntensity, bool preserveBoundaries);

		// Uniform noise addition
		void NoiseAdditionUniform (Polyhedron *pMesh, double noiseIntensity, bool preserveBoundaries);

		// Gaussian noise addition. The noise intensity stands for the standard deviation of the Gaussian distribution.
		void NoiseAdditionGaussian (Polyhedron *pMesh, double noiseIntensity, bool preserveBoundaries);

		// For more information about Laplacian smoothing, please refer to the state of art report of Taubin in Eurographics 2000
		void LaplacianSmoothing (Polyhedron *pMesh, double deformFactor, int iteraNum, bool preserveBoundaries);

		// Vertex coordinates are quantized to 2^(bitDepth) levels
		// TODO: fixing after quantization the potentially introduced degeneracies, such as removing the null surface facet
		void CoordinateQuantization (Polyhedron *pMesh, int bitDepth);


		void Translation (Polyhedron *pMesh, double xTranslation, double yTranslation, double zTranslation);
		
		// (xAxis,yAxis,zAxis) forms the rotation axis
		// angle stands for the counterclockwise rotation angle, seen from the positive side of the rotation axis
		void Rotation (Polyhedron *pMesh, double xAxis, double yAxis, double zAxis, double angle);

		void UniformScaling (Polyhedron *pMesh, double scalingFactor);

		// Different subdivision methods. For more information, please refer to the CGAL user manual.
		void Subdivision (Polyhedron *pMesh, Subdivision_type subdivisionType, int depth);
		void SubdivisionCatmullClark (Polyhedron *pMesh, int depth);
		void SubdivisionLoop (Polyhedron *pMesh, int depth);
		void SubdivisionDooSabin (Polyhedron *pMesh, int depth);
		void SubdivisionSqrt3 (Polyhedron *pMesh, int depth);
		void SubdivisionMidpoint (Polyhedron *pMesh, int depth);

		// Mesh simplification based on edge collasp. For more information, please refer to the CGAL user manual.
		void Simplification (Polyhedron *pMesh, int targetEdgeNum);

	private:


};

#endif // Processing_COMPONENT_H
