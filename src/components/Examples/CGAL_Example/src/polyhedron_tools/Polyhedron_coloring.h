///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2010
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Polyhedron_coloring_H
#define Polyhedron_coloring_H

#include "../../../../../mepp/Polyhedron/polyhedron.h"

#include "./Polyhedron_local_computation.h" // needed for visualize_feature_edges_facet_coloring

#include "../useful_methods/Training_and_testing_data_SVM.h"

static void visualize_region_facet_coloring(std::vector < Facet_handle >& region_f, Vector& RGB)
{
    normalize_Vector(RGB);
    /////////////////////////////////////////////////////////////////////////////////
    std::vector < Facet_handle >::iterator it(region_f.begin()), ite(region_f.end());
    for(;it!=ite;++it)
    {
        (*it)->color(RGB.x(), RGB.y(), RGB.z());
    }
}

static void visualize_region_vertex_coloring(std::vector < Vertex_handle >& region_v, Vector& RGB)
{
    normalize_Vector(RGB);
    /////////////////////////////////////////////////////////////////////////////////
    std::vector < Vertex_handle >::iterator it(region_v.begin()), ite(region_v.end());
    for(;it!=ite;++it)
    {
        (*it)->color(RGB.x(), RGB.y(), RGB.z());
    }
}


/**
    visualize_dihedral_angles_coloring is just for guessing the
    dihedral angles, it is not for exact visualization since
    the painting is order dependent
**/
static void visualize_dihedral_angles_coloring(PolyhedronPtr pMesh)
{
    double precision_used = 0.3;
    for (Edge_iterator pEdge = pMesh->edges_begin(); pEdge != pMesh->edges_end(); pEdge++)
    {
        double d = dihedralAngle(pEdge);
        if (d > precision_used)
        {
            pEdge->facet()->color(1.0, 0.0, 0.0);
            pEdge->opposite()->facet()->color(1.0, 0.0, 0.0);
        }
        else if (d < -precision_used)
        {
            pEdge->facet()->color(0.0, 0.0, 1.0);
            pEdge->opposite()->facet()->color(0.0, 0.0, 1.0);
        }
        else
        {
            pEdge->facet()->color(0.0, 1.0, 0.0);
            pEdge->opposite()->facet()->color(0.0, 1.0, 0.0);
        }
    }
}

static void visualize_feature_edges_facet_coloring(PolyhedronPtr pMesh)
{
    for (Facet_iterator pFacet = pMesh->facets_begin(); pFacet != pMesh->facets_end(); pFacet++)
    {
        unsigned int cpt_feature=0;
        if (pFacet->halfedge()->isFeature()) cpt_feature++;

        if (pFacet->halfedge()->next()->isFeature()) cpt_feature++;

        if (pFacet->halfedge()->prev()->isFeature()) cpt_feature++;

        if (cpt_feature==3)
            pFacet->color(1.0, 0.0, 0.0);
        else if (cpt_feature==2)
            pFacet->color(1.0, 1.0, 0.0);
        else if (cpt_feature==1)
            pFacet->color(0.0, 1.0, 0.0);
        else if (isInvolvedInGeometricVertexFold(pFacet, false))
        {
            pFacet->color(1.0, 1.0, 1.0);
        }
        else
            pFacet->color(0.0, 0.0, 1.0);
    }
}

static void visualize_crest_point_coloring(PolyhedronPtr pMesh, bool noisy=false)
{
    double maximality=500.0; // should be in [500;9000] to identify less sensitive to noise crest points
    for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
    {
        if( isLocalCrestPoint(pVertex, maximality) &&
            (!noisy || !shouldRemoveNonSignificantLines(pVertex, maximality)) )
            pVertex->color(1.0, 0.0, 0.0);
        else
            pVertex->color(0.0, 0.0, 1.0);
    }
}

static void visualize_curvature_vertex_coloring(PolyhedronPtr pMesh,
                                                unsigned int which=2,
                                                bool laplacian_averaging=false)
{
    double minGC = 1e8, maxGC = -1e8, tmpGC; // used to correct color scaling
    for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
    {
        tmpGC = curvature(pVertex, which, laplacian_averaging);

        if(tmpGC<minGC) minGC = tmpGC;
        if(tmpGC>maxGC) maxGC = tmpGC;
    }

    if(maxGC==minGC)maxGC+=1e-8;

    for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); pVertex++)
    {
        tmpGC = curvature(pVertex, which, laplacian_averaging);

        //cout << "tmpGC = " << tmpGC << endl;
        // [minGC;maxGC] => [0;1]
        tmpGC = (tmpGC-minGC)/(maxGC-minGC);

        //cout << "tmpGC = " << tmpGC << endl;

        if(tmpGC<=0.5)
            pVertex->color(0.0, 2.0*tmpGC, 1.0-2.0*tmpGC);
        else
            pVertex->color(2.0*tmpGC-1.0, 2.0-2.0*tmpGC, 0.0);
    }
}

///////////////////////////////////////////////////////////////////////////////////////
// test and coloring
static void test_local_edge_region_growing(PolyhedronPtr pMesh, double radius)
{
    Vector blue(0.0,0.0,1.0), red(1.0,0.0,0.0);
    unsigned int cpt=0;

    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
    {
        std::vector < Vertex_handle > left_region_v;
        std::vector < Vertex_handle > right_region_v;
        std::vector < Halfedge_handle > left_region_h;
        std::vector < Halfedge_handle > right_region_h;
        std::vector < Facet_handle > left_region_f;
        std::vector < Facet_handle > right_region_f;

        if(cpt%500==0)
        {
            local_edge_region_growing(  pEdge,
                                        radius,
                                        left_region_v, right_region_v, left_region_h,
                                        right_region_h, left_region_f, right_region_f
                                    );
            visualize_region_facet_coloring(left_region_f, blue);
            visualize_region_facet_coloring(right_region_f, red);

            visualize_region_vertex_coloring(left_region_v, blue);
            visualize_region_vertex_coloring(right_region_v, red);
        }
        cpt++;
    }


}
#endif // Polyhedron_coloring_H
