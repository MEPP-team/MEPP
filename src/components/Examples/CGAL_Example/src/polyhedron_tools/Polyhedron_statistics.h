///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2010
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Polyhedron_statistics_H
#define Polyhedron_statistics_H

#include <math.h>

#include "../../../../../mepp/Polyhedron/polyhedron.h"

#include "../useful_methods/Histogram.h"
#include "./Polyhedron_local_computation.h"

static void build_histogram_cos_dihedral_angle( PolyhedronPtr pMesh,
                                                Histogram<double>& H,
                                                bool cumulative=false,
                                                bool only_feature_edge=false,
                                                bool only_normal_edge=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++
        )
        {
            if(only_feature_edge && pEdge->isFeature())
                H.add(cosDihedralAngle(pEdge));
            else if(only_normal_edge && !pEdge->isFeature())
                H.add(cosDihedralAngle(pEdge));
            else if(!only_feature_edge && !only_normal_edge)
                H.add(cosDihedralAngle(pEdge));
        }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Cosinus angles dihédraux", "Cosinus", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Cosinus angles dihédraux", "Cosinus", "Proportion", true);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
    cout << "First non-zero value = " << H.getBinCenter(H.firstNonZero()) << endl;
    cout << "Maximum non-zero value = " << H.getBinCenter(H.argmax()) << endl;
    cout << "Last non-zero value = " << H.getBinCenter(H.lastNonZero()) << endl;
}

static void build_histogram_dihedral_angle( PolyhedronPtr pMesh,
                                            Histogram<double>& H,
                                            bool cumulative=false,
                                            bool only_feature_edge=false,
                                            bool only_normal_edge=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++
        )
        {
            if(only_feature_edge && pEdge->isFeature())
            {
                H.add(dihedralAngle(pEdge)*180.0/M_PI);
            }
            else if(only_normal_edge && !pEdge->isFeature())
                H.add(dihedralAngle(pEdge)*180.0/M_PI);
            else if(!only_feature_edge && !only_normal_edge)
                H.add(dihedralAngle(pEdge)*180.0/M_PI);
        }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Angles dihédraux", "Angles", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Angles dihédraux", "Angles", "Proportion", false);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
    cout << "First non-zero value = " << H.getBinCenter(H.firstNonZero()) << endl;
    cout << "Maximum non-zero value = " << H.getBinCenter(H.argmax()) << endl;
    cout << "Last non-zero value = " << H.getBinCenter(H.lastNonZero()) << endl;
}
////////////////////////////////////////////////////////////////////////
static void build_histogram_triangle_min_angle( PolyhedronPtr pMesh,
                                                Histogram<double>& H,
                                                bool cumulative=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Facet_iterator pFace = pMesh->facets_begin();
            pFace!=pMesh->facets_end();
            pFace++)
    {
        double minAngleT = computeMinTriangleAngles(   pFace->halfedge()->prev()->vertex()->point(),
                           pFace->halfedge()->vertex()->point(),
                           pFace->halfedge()->next()->vertex()->point());
        H.add(minAngleT*180.0/M_PI);
    }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Angles minimaux", "Angles", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Angles minimaux", "Angles", "Proportion", false);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
}

static void build_histogram_triangle_max_angle( PolyhedronPtr pMesh,
                                                Histogram<double>& H,
                                                bool cumulative=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Facet_iterator pFace = pMesh->facets_begin();
            pFace!=pMesh->facets_end();
            pFace++)
    {
        double maxAngleT = computeMaxTriangleAngles(   pFace->halfedge()->prev()->vertex()->point(),
                           pFace->halfedge()->vertex()->point(),
                           pFace->halfedge()->next()->vertex()->point());
        H.add(maxAngleT*180.0/M_PI);
    }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Angles maximaux", "Angles", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Angles maximaux", "Angles", "Proportion", false);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
}


static void build_histogram_max_curvature(  PolyhedronPtr pMesh,
                                            Histogram<double>& H,
                                            bool cumulative=false,
                                            bool only_feature_edge=false,
                                            bool only_normal_edge=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Vertex_iterator pVertex = pMesh->vertices_begin();
            pVertex!=pMesh->vertices_end();
            pVertex++)
    {
        if(only_feature_edge && !isASharpEdgeOnVertexNeighborhood(pVertex))
            continue;

        if(only_normal_edge && isASharpEdgeOnVertexNeighborhood(pVertex))
            continue;

        H.add(pVertex->closestObs()->Kmax);
    }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Max curvatures", "Max curvature", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Max curvatures", "Max curvature", "Proportion", true);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
}

static void build_histogram_normal_curvature_estimate(  PolyhedronPtr pMesh,
                                                        Histogram<double>& H,
                                                        bool cumulative=false,
                                                        bool only_feature_edge=false,
                                                        bool only_normal_edge=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge!=pMesh->edges_end();
            pEdge++)
    {
        if(only_feature_edge && !pEdge->isFeature()) continue;

        if(only_normal_edge && pEdge->isFeature()) continue;

         H.add(0.5*( normalCurvatureEstimateAtEdge(pEdge)+
            normalCurvatureEstimateAtEdge(pEdge->opposite()))   );
    }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Normal curvatures", "Normal curvature", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Normal curvatures", "Normal curvature", "Proportion", true);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
}

static void build_histogram_mean_curvature( PolyhedronPtr pMesh,
                                            Histogram<double>& H,
                                            bool cumulative=false,
                                            bool only_feature_edge=false,
                                            bool only_normal_edge=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge!=pMesh->edges_end();
            pEdge++)
    {
        if(only_feature_edge && !pEdge->isFeature()) continue;

        if(only_normal_edge && pEdge->isFeature()) continue;

         H.add(meanCurvatureAtEdge(pEdge));
    }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Mean curvatures", "Mean curvature", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Mean curvatures", "Mean curvature", "Proportion", true);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
}

static void build_histogram_similarity( PolyhedronPtr pMesh,
                                        Histogram<double>& H,
                                        bool cumulative=false)
{
    H.clear();
    /////////////////////////////////////////////////////////////////////
    for (   Vertex_iterator pVertex = pMesh->vertices_begin();
            pVertex!=pMesh->vertices_end();
            pVertex++)
    {
         H.add(multipleScaleSimilarity(pVertex, 0.07));
    }

    if(cumulative)
    {
        //H.toCumulativeFromRight();
        H.toCumulative();
        //H.truncFromRightStop(.99);
        H.gnuplotOut("Similarity", "Similarity", "Proportion", false);
    }
    else
    {
        //H.toLogScale(); // less good than the gnuplot logscale
        H.gnuplotOut("Similarity", "Similarity", "Proportion", true);
    }

    cout << "Mean = " << H.getMean() << " Min = " << H.getMin() << " Max = " << H.getMax() << " Variance = " << H.getVariance() << endl;
}
#endif // Polyhedron_statistics_H
