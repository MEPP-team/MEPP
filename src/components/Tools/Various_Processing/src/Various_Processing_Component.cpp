///////////////////////////////////////////////////////////////////////////
// Author: Kai Wang
// Year: 2009
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "Various_Processing_Component.h"

#include <CGAL/Subdivision_method_3.h>

#ifndef __linux__
#ifndef __APPLE__
	#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
	#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
	#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#endif
#endif

double PI = 3.1415926535;

Various_Processing_Component::Various_Processing_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	// MEPP 2
	componentName = "Various_Processing_Component";
	init = 1;
}

// see Processing_Component.h for the effects of the following functions

void Various_Processing_Component::NoiseAddition (PolyhedronPtr pMesh, Noise_type noiseType, double noiseIntensity, bool preserveBoundaries)
{
	if (noiseType==UNIFORM)
		NoiseAdditionUniform(pMesh, noiseIntensity, preserveBoundaries);
	else if (noiseType==GAUSSIAN)
		NoiseAdditionGaussian(pMesh, noiseIntensity, preserveBoundaries);
}

void Various_Processing_Component::NoiseAdditionUniform (PolyhedronPtr pMesh, double noiseIntensity, bool preserveBoundaries)
{
	// mesh centre calculation based on discrete "moment".
	//TODO: maybe in the future it will be based on volume moment.
	int numVertex = pMesh->size_of_vertices();;
	Vector centroid = Point3d(0,0,0) - CGAL::ORIGIN;
	double distancetoCentroid = 0.0;
	Vertex_iterator	pVertex;
	for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
	{
		Vector vectemp = pVertex->point() - CGAL::ORIGIN;
		centroid = centroid + vectemp;
	}
	centroid = centroid/numVertex;

	// calculate the average distance from vertices to mesh centre
	for (pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
		Vector vectemp = pVertex->point() - CGAL::ORIGIN;
		distancetoCentroid = distancetoCentroid + (double)std::sqrt((vectemp - centroid) * (vectemp - centroid));
	}
	distancetoCentroid = distancetoCentroid/numVertex;

	// add random uniform-distributed (between [-noiseLevel, +noiseLevel])
	srand((unsigned)time(NULL));
	double noisex, noisey, noisez;
	double noiseLevel = distancetoCentroid * noiseIntensity;
	for (pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
		// keep boundaries untouched if demanded by user
		bool is_border_vertex = false;
		bool stopFlag = false;
		Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
		do
		{
			if (hav->is_border()==true)
			{
				is_border_vertex = true;
				stopFlag = true;
			}
			hav++;
		} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

		if ((preserveBoundaries==true)&&(is_border_vertex==true))
			continue;

		noisex = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
		noisey = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
		noisez = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
		Vector temp = Point3d(noisex, noisey, noisez) - CGAL::ORIGIN;
		pVertex->point() = pVertex->point() + temp;
	}

	// for correct rendering, we need to update the mesh normals
	pMesh->compute_normals();
}

// this time, we add Gaussian-distributed additive noise to the mesh vertex coordinates
// the standard deviation of the Gaussian distribution is "noiseLevel = distancetoCentroid * noiseIntensity"
void Various_Processing_Component::NoiseAdditionGaussian (PolyhedronPtr pMesh, double noiseIntensity, bool preserveBoundaries)
{
	int numVertex = pMesh->size_of_vertices();;
	Vector centroid = Point3d(0,0,0) - CGAL::ORIGIN;
	double distancetoCentroid = 0.0;
	Vertex_iterator	pVertex;
	for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
	{
		Vector vectemp = pVertex->point() - CGAL::ORIGIN;
		centroid = centroid + vectemp;
	}
	centroid = centroid/numVertex;

	for (pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
		Vector vectemp = pVertex->point() - CGAL::ORIGIN;
		distancetoCentroid = distancetoCentroid + (double)std::sqrt((vectemp - centroid) * (vectemp - centroid));
	}
	distancetoCentroid = distancetoCentroid/numVertex;

	srand((unsigned)time(NULL));
	double noisex, noisey, noisez;
	double * gaussNumbers = new double[3];
	double noiseLevel = distancetoCentroid * noiseIntensity;
	for (pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
	{
		bool is_border_vertex = false;
		bool stopFlag = false;
		Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
		do
		{
			if (hav->is_border()==true)
			{
				is_border_vertex = true;
				stopFlag = true;
			}
			hav++;
		} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

		if ((preserveBoundaries==true)&&(is_border_vertex==true))
			continue;

		// pseudo-random Gaussian-distributed numbers generation from uniformly-distributed pseudo-random numbers
		double x, y, r2;
		for (int i=0; i<3; i++)
		{
			do
			{
				x = -1.0 + 2.0 * 1.0*rand()/RAND_MAX;
				y = -1.0 + 2.0 * 1.0*rand()/RAND_MAX;
				r2 = x * x + y * y;
			} while ((r2>1.0)||(r2==0.0));
			gaussNumbers[i] = y * sqrt(-2.0 * log(r2) / r2);
		}

		noisex = noiseLevel * gaussNumbers[0];
		noisey = noiseLevel * gaussNumbers[1];
		noisez = noiseLevel * gaussNumbers[2];
		Vector temp = Point3d(noisex, noisey, noisez) - CGAL::ORIGIN;
		pVertex->point() = pVertex->point() + temp;
	}

	pMesh->compute_normals();

	delete [] gaussNumbers;
	gaussNumbers = 0;
}

void Various_Processing_Component::LaplacianSmoothing (PolyhedronPtr pMesh, double deformFactor, int iteraNum, bool preserveBoundaries)
{
	Vertex_iterator	pVertex;
	int numVertex = pMesh->size_of_vertices();
	Vector * newPositions = new Vector[numVertex];

	for (int i=0; i<iteraNum; i++)
	{
		int n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
			Vector currentVector = pVertex->point() - CGAL::ORIGIN;

			// do not smooth the boundary vertices if demanded by user
			bool is_border_vertex = false;
			bool stopFlag = false;
			Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
			do
			{
				if (hav->is_border()==true)
				{
					is_border_vertex = true;
					stopFlag = true;
				}
				hav++;
			} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

			if ((preserveBoundaries==true)&&(is_border_vertex==true))
			{
				newPositions[n] = currentVector;
				n++;
				continue;
			}

			std::size_t degree = (*pVertex).vertex_degree();
			double alpha = 1.0/degree;
			Vector vectemp = Point3d(0,0,0) - CGAL::ORIGIN;
			Halfedge_around_vertex_circulator h = (*pVertex).vertex_begin();
			do
			{
				vectemp = vectemp+(h->opposite()->vertex()->point()-CGAL::ORIGIN-currentVector)*alpha;
				++h;
			} while (h != (*pVertex).vertex_begin());
			newPositions[n] = currentVector + deformFactor*vectemp;
			n++;
		}

		n = 0;
		for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
		{
			pVertex->point() = Point3d(0,0,0) + newPositions[n];
			n++;
		}
	}

	delete [] newPositions;
	newPositions = 0;

	pMesh->compute_normals();
}

// TODO: fixing after quantization the potentially introduced degeneracies, such as removing the null surface facet
void Various_Processing_Component::CoordinateQuantization (PolyhedronPtr pMesh, int bitDepth)
{
	Vertex_iterator	pVertex;
	double quantizationLevel = std::pow(2.0,bitDepth);
	pVertex = pMesh->vertices_begin();
	Point3d point = pVertex->point();
	double xmax = double(point.x());
	double xmin = xmax;
	double ymax = double(point.y());
	double ymin = ymax;
	double zmax = double(point.z());
	double zmin = zmax;
	pVertex++;

	for (; pVertex != pMesh->vertices_end(); pVertex++)
	{
		point = pVertex->point();
		double x = double(point.x());
		double y = double(point.y());
		double z = double(point.z());
		if (x>xmax)
			xmax = x;
		if (x<xmin)
			xmin = x;
		if (y>ymax)
			ymax = y;
		if (y<ymin)
			ymin = y;
		if (z>zmax)
			zmax = z;
		if (z<zmin)
			zmin = z;
	}

	double xstep = (xmax-xmin)/quantizationLevel;
	double ystep = (ymax-ymin)/quantizationLevel;
	double zstep = (zmax-zmin)/quantizationLevel;

	for (pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end();pVertex++)
	{
		point = pVertex->point();
		double x = double(point.x());
		double y = double(point.y());
		double z = double(point.z());

		double xquantified, yquantified, zquantified;

		double xint = 1.0*std::floor((x-xmin)/xstep)*xstep + xmin;
		double xfrac = x - xint;
		if (xfrac<=(0.5*xstep))
			xquantified = xint;
		else
			xquantified = xint + xstep;

		double yint = 1.0*std::floor((y-ymin)/ystep)*ystep + ymin;
		double yfrac = y - yint;
		if (yfrac<=(0.5*ystep))
			yquantified = yint;
		else
			yquantified = yint +ystep;

		double zint = 1.0*std::floor((z-zmin)/zstep)*zstep + zmin;
		double zfrac = z - zint;
		if (zfrac<=(0.5*zstep))
			zquantified = zint;
		else
			zquantified = zint + zstep;

		pVertex->point() = Point3d(xquantified,yquantified,zquantified);
	}

	pMesh->compute_normals();
}

void Various_Processing_Component::Translation (PolyhedronPtr pMesh, double xTranslation, double yTranslation, double zTranslation)
{
	Vector translationVector(xTranslation,yTranslation, zTranslation);
	Affine_transformation translation(CGAL::TRANSLATION, translationVector);
	std::transform(pMesh->points_begin(), pMesh->points_end(), pMesh->points_begin(), translation);

	pMesh->compute_normals();
}

void Various_Processing_Component::Rotation (PolyhedronPtr pMesh, double xAxis, double yAxis, double zAxis, double angle)
{
	// normalize the translation vector
	double normAxis = sqrt(xAxis*xAxis + yAxis*yAxis + zAxis*zAxis);
	xAxis = xAxis / normAxis;
	yAxis = yAxis / normAxis;
	zAxis = zAxis / normAxis;

	// construction of the rotation matrix
	double c = cos(angle/180.0*PI);
	double s = sin(angle/180.0*PI);
	double m00 = xAxis*xAxis + (1.0-xAxis*xAxis)*c;
	double m01 = xAxis*yAxis*(1.0-c) - zAxis*s;
	double m02 = xAxis*zAxis*(1.0-c) + yAxis*s;
	double m10 = xAxis*yAxis*(1.0-c) + zAxis*s;
	double m11 = yAxis*yAxis + (1.0-yAxis*yAxis)*c;
	double m12 = yAxis*zAxis*(1.0-c) - xAxis*s;
	double m20 = xAxis*zAxis*(1.0-c) - yAxis*s;
	double m21 = yAxis*zAxis*(1.0-c) + xAxis*s;
	double m22 = zAxis*zAxis + (1.0-zAxis*zAxis)*c;

	// realize rotation by applying general affine transformation through matrix multiplication
	Affine_transformation rotation(m00, m01, m02, m10, m11, m12, m20, m21, m22);
	std::transform(pMesh->points_begin(), pMesh->points_end(), pMesh->points_begin(), rotation);

	pMesh->compute_normals();
}

void Various_Processing_Component::UniformScaling (PolyhedronPtr pMesh, double scalingFactor)
{
	Affine_transformation uniformScaling(CGAL::SCALING, scalingFactor);
	std::transform(pMesh->points_begin(), pMesh->points_end(), pMesh->points_begin(), uniformScaling);

	pMesh->compute_normals();
}

void Various_Processing_Component::Subdivision (PolyhedronPtr pMesh, Subdivision_type subdivisionType, int depth)
{
	if (subdivisionType==CATMULLCLARK)
		SubdivisionCatmullClark(pMesh, depth);
	else if (subdivisionType==LOOP)
		SubdivisionLoop(pMesh, depth);
	else if (subdivisionType==DOOSABIN)
		SubdivisionDooSabin(pMesh, depth);
	else if (subdivisionType==SQRT3)
		SubdivisionSqrt3(pMesh, depth);
	else if (subdivisionType==MIDPOINT)
		SubdivisionMidpoint(pMesh, depth);

	// after connectivity modification of the mesh, we need to update the related internal properties such as "m_pure_quad" or "m_pure_triangle"
	// other updates can be necessary depending on application
	pMesh->compute_type();
	pMesh->compute_normals();
}

void Various_Processing_Component::SubdivisionCatmullClark (PolyhedronPtr pMesh, int depth)
{
	CGAL::Subdivision_method_3::CatmullClark_subdivision(*pMesh,depth);
}

void Various_Processing_Component::SubdivisionLoop (PolyhedronPtr pMesh, int depth)
{
	CGAL::Subdivision_method_3::Loop_subdivision(*pMesh,depth);
}

void Various_Processing_Component::SubdivisionDooSabin (PolyhedronPtr pMesh, int depth)
{
	CGAL::Subdivision_method_3::DooSabin_subdivision(*pMesh,depth);
}

void Various_Processing_Component::SubdivisionSqrt3 (PolyhedronPtr pMesh, int depth)
{
	CGAL::Subdivision_method_3::Sqrt3_subdivision(*pMesh,depth);

}

// the midpoint subdivision does not introduce any geometry distortion, it simply adds vertices in the middle of edges
// it is ofter used to test the robustness of a geometry processing algorithm
void Various_Processing_Component::SubdivisionMidpoint (PolyhedronPtr pMesh, int depth)
{
	CGAL::Subdivision_method_3::PTQ(*pMesh, Linear_mask_3<Polyhedron>(), depth);
}

// these typedefs are necessary since the embedded simplification algorithm in CGAL
// seems only working with the base Polyhedron class (or its enriched version)
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel,Enriched_items> Surface;

void Various_Processing_Component::Simplification (Polyhedron *pMesh, int targetEdgeNum)
{
#ifndef __linux__
#ifndef __APPLE__
	CGAL::Surface_mesh_simplification::Count_stop_predicate<Surface> stop(targetEdgeNum);

	// upcast... it is indeed ugly... sorry for that. trying to find new solutions...
	Surface * surface = (Surface*)(pMesh);

	/*int r =*/ CGAL::Surface_mesh_simplification::edge_collapse
            (*surface
            ,stop
#if (CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,5,0))
            ,CGAL::vertex_index_map(get(CGAL::vertex_external_index, *surface))
                  .halfedge_index_map(get(CGAL::halfedge_external_index, *surface))
#else
            ,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index, *surface))
                  .edge_index_map(boost::get(CGAL::edge_external_index, *surface))
#endif
            );

	pMesh->compute_type();
	pMesh->compute_normals();
#endif
#endif
}
#endif
