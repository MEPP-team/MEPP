///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Polyhedron_local_modification_H
#define Polyhedron_local_modification_H

#include "../useful_methods/General_purpose.h"
#include "../triangle_tools/Triangle_shape_analysis.h"
//#include "../polyhedron_tools/Polyhedron_coloring.h" cannot be used! dependance in the other direction!

#include"../mutable_priority_queue/mutable_priority_queue.h"

#define FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS // improves the results (e.g. cow tail is better)

#define NO_FEATURE_EDGE_RECOMPUTATION

//#define INIT_CURVATURE_SENSITIVE // to be more sensitive to the initial curvature (recommended for smooth surfaces)

//#define HYSTERESIS_TWO_THRESHOLDS

//#define SMOOTHING

// global variables
static double FEATURE_COS_DIHEDRAL_ANGLE = 0.5f;   // 0.82f< <=> >35 degrees ; 0.5f< <=> >60 degrees ; 0.642f< <=> >50 degrees

static double SMOOTHING_FACTOR = 0.07;              // 0.6666 0.9 0.06 0.07 this value is important since it controls the number of accepted smoothing candidates (especially at the beginning)
                                                    // and so it controls the amount of randomness authorized (if at smoothing candidate is not accepted, we try a uniformly taken and then a random one)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// closest observation uptate for distance computation efficiency
// CARE MUST BE TAKEN: this function is correct only for *small moves*
static inline Vertex_handle searchClosestObs(const Vertex_handle pVertex)
{
    if (pVertex->closestObs()==Vertex_handle())
    {
        cout << "searchClosestObs: NULL Vertex_handle observation!!" << endl;
        return Vertex_handle();
    }

    Vertex_handle v, minV;
    double d, minD;


    Vector gradient;
    /*long double initNRJ, NRJtemp;

    initNRJ = featureFuncDataAtt(pVertex,false,gradient);*/
    //////////////////////////////////////////////////////////////////////////////
    // Initialize to the last estimated nearest initial point of the current label.
    minV = pVertex->closestObs();
    minD = (pVertex->point() - minV->point()).squared_length();

    // Iterate on the observed mesh through the neighbors of the
    // last nearest observation and get the closest new one
    //Halfedge_around_vertex_circulator cHeS = minV->vertex_begin();
    //Halfedge_around_vertex_circulator cHeE = cHeS;
    //CGAL_For_all(cHeS,cHeE)
    //{
    //Halfedge_around_vertex_circulator cHe = cHeS->prev()->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator cHe = minV->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {

        d = (pVertex->point() - cHe->prev()->vertex()->point()).squared_length();
        if (d<minD)
        {
            /*pVertex->SetClosestObs(cHe->prev()->vertex());
            NRJtemp = featureFuncDataAtt(pVertex,false,gradient);

            if(NRJtemp<=initNRJ)
            {
                minD=d;
                minV=cHe->prev()->vertex();
            }*/

            minD=d;
            minV=cHe->prev()->vertex();
        }

    }
    //}

    // here we must check that the min will not change if we use again this function
    // However that process can accept far away candidated close to the surface on another location
    Vertex_handle correctminV;
    bool is_OK = false;
    //unsigned int cpt=0;
    while (!is_OK)
    {
        correctminV = minV;
        cHe = minV->vertex_begin();
        cHe2 = cHe;
        CGAL_For_all(cHe,cHe2)
        {

            d = (pVertex->point() - cHe->prev()->vertex()->point()).squared_length();
            if (d<minD )
            {
                minD=d;
                minV=cHe->prev()->vertex();
            }
        }

        is_OK = (correctminV == minV); // the minimum is stable!!

        //cpt++;
    }
    //if(cpt>1) cout << "cpt = " << cpt << endl;

    assert(minV!=Vertex_handle());

    return minV;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline unsigned int valence(const Vertex_handle pVertex)
{
    return (unsigned int)CGAL::circulator_size(pVertex->vertex_begin());
}

static inline bool is_border(const Halfedge_handle h)
{
    if (h->is_border() || h->opposite()->facet()==Facet_handle())
        return true;

    return false;
}

static inline bool is_border(const Vertex_handle pVertex)
{
    assert(pVertex!=Vertex_handle());
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    if (Halfedge_handle(pHalfEdge) == Halfedge_handle()) // isolated vertex
        return true;
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
    {
        if (is_border(pHalfEdge))
            return true;
    }
    return false;
}

static inline double vertexDegreePotential(const Vertex_handle pVertex)
{
    double pot;
    if (is_border(pVertex))
        pot = valence(pVertex)-4.0f;
    else
        pot = valence(pVertex)-6.0f;

    return pot*pot;
}

/**
    meanEdgeLength computes the mean of mesh edge length
**/
static inline double meanEdgeLength(const PolyhedronPtr pMesh)
{
    double sum=0.0;

    for(Edge_iterator   pEdge=pMesh->edges_begin();
                        pEdge!=pMesh->edges_end();
                        ++pEdge
    )
    {
        sum += std::sqrt(pEdge->length2());
    }

    return sum/CGAL::iterator_distance(pMesh->edges_begin(), pMesh->edges_end());
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline bool isvertexOneRingConvex(const Vertex_handle pVertex)
{
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    unsigned int vertex_valence = valence(pVertex);
    if (vertex_valence<=2) return false;
    if ( vertex_valence==3 ) return true;
    // a border vertex have an half-ring neighborhood and we are not interesting testing this case
    if ( is_border(pVertex) ) return false;
    // this piece of source code makes the simplification based on edge collapse not working
    // that is because the one-ring neighborhood is not necessarily convex
    CGAL_For_all(h,he)
    {
        Point3d h_prev_point = h->prev()->vertex()->point(),
                               h_next_point = h->next()->vertex()->point(),
                                              h_opp_next_point = h->opposite()->next()->vertex()->point();

        if (isTriangleAngleGreaterThanPi(h_next_point, h_prev_point, h_opp_next_point)) return false;
    }

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// normals
static inline void compute_planar_facet_normal(const Facet_handle pFace, Vector& NFace)
{
    // Sometimes a null facet is linked to a non-null halfedge
    // in this case we chose to set the normal to zero
    if (pFace == Facet_handle())
    {
        cout << "compute_planar_facet_normal: WARNING: NULL FACE!!!" << endl;
        NFace = CGAL::NULL_VECTOR;
        return;
    }
    Halfedge_handle h = pFace->halfedge();
    NFace = CGAL::cross_product(Vector(h->next()->vertex()->point()-h->vertex()->point()),
                                Vector(h->next()->next()->vertex()->point()-h->next()->vertex()->point()) );

    normalize_Vector(NFace);
}

static inline void compute_region_normal(   const std::vector < Facet_handle >& region_f,
                                            Vector& NRegion)
{
    NRegion = CGAL::NULL_VECTOR;
    std::vector < Facet_handle >::const_iterator it(region_f.begin()), ite(region_f.end());
    for(;it!=ite;++it)
    {
        Vector NFace;
        compute_planar_facet_normal(*it, NFace);
        NRegion = NRegion + NFace*triangleArea(*it);
    }

    normalize_Vector(NRegion);
}


// definition taken from Konrad Polthier
static inline void compute_edge_normal(const Halfedge_handle h, Vector& N)
{
    Vector FN1 = CGAL::NULL_VECTOR, FN2 = CGAL::NULL_VECTOR;
    // A border edge is a feature edge
    if (h->facet()==Facet_handle())
    {
        compute_planar_facet_normal(Halfedge_handle(h->opposite())->facet(), FN2);
    }
    else if (h->opposite()->facet()==Facet_handle())
    {
        compute_planar_facet_normal(h->facet(), FN1);
    }
    else
    {
        compute_planar_facet_normal(h->facet(), FN1);
        compute_planar_facet_normal(Halfedge_handle(h->opposite())->facet(), FN2);
    }

    N = FN1+FN2;

    normalize_Vector(N);
}

static inline void compute_vertex_normal(const Vertex_handle pVertex, Vector& N)
{
    N = CGAL::NULL_VECTOR;
    Vector NFace;
    double Aire; // area weighting gives better surface normal (compared to the real surface)

    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        NFace = CGAL::cross_product(Vector(cHe->next()->vertex()->point()-cHe->vertex()->point()),
                                    Vector(cHe->next()->next()->vertex()->point()-cHe->next()->vertex()->point()) );
        Aire = std::sqrt(NFace*NFace)*0.5f;
        N = N + NFace*Aire;
    }
    // We do not need to divide NT by SAires or by
    // the vertex valence since we normalize it
    normalize_Vector(N);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature detection methods
static inline double cosDihedralAngle(const Halfedge_handle h)
{
    // A border edge is a feature edge
    if (is_border(h))
    {
        return -1; // cos(M_PI): this choice is imposed by the community (see "Feature detection for surface meshes")
    }
    //////////////////////////////////////////////////////////////////////////////////
    Vector FN1, FN2;

    compute_planar_facet_normal(h->facet(), FN1);
    if (FN1==CGAL::NULL_VECTOR) return FEATURE_COS_DIHEDRAL_ANGLE; // true by default
    compute_planar_facet_normal(Halfedge_handle(h->opposite())->facet(), FN2);
    if (FN2==CGAL::NULL_VECTOR) return FEATURE_COS_DIHEDRAL_ANGLE; // true by default

    return FN1*FN2;
}

// used to detect sharp features in "Geometric snakes for triangular meshes" (2002)
static inline double normal_variation(const Vertex_handle pVertex)
{
    double minCosNvNF= 1.0;
    ////////////////////////////////////////////////////////////////
    Vector Nvector;
    compute_vertex_normal(pVertex, Nvector);
    ////////////////////////////////////////////////////////////////
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        Vector Nface;
        compute_planar_facet_normal(h->facet(),Nface);
        ////////////////////////////////////////////////////////////////
        double cosNvNF=Nvector*Nface;
        if(minCosNvNF>cosNvNF) minCosNvNF=cosNvNF;
    }

    return minCosNvNF;
}

// directional:
static inline double minCosDihedralAngleToTheStartorEnd(const Halfedge_handle pE,
                                                        bool start)
{
    // we want to find the largest dihedral angle (so the min cos)!!
    double minCos = 1.0, /*maxCosF = -1.0, */cosTA;
    Halfedge_around_vertex_circulator cHe = (start)?pE->vertex()->vertex_begin():pE->prev()->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        if(Halfedge_handle(cHe)==pE)
            continue;

        cosTA = cosDihedralAngle(cHe);
        if(cosTA<minCos) minCos=cosTA;

    }

    return minCos;
}

// non-directional:
static inline double minCosDihedralAngle(   const Halfedge_handle pE,
                                            bool start)
{
    return std::min<double>(minCosDihedralAngleToTheStartorEnd(pE,true),minCosDihedralAngleToTheStartorEnd(pE,false));
}


static inline double extendedCosDihedralAngle(const Halfedge_handle h)
{
    // A border edge is a feature edge
    if (is_border(h))
    {
        return 1; // cos(M_PI): this choice is imposed by the community (see "Feature detection for surface meshes")
    }
    //////////////////////////////////////////////////////////////////////////////////
    Vector FN1, FN2;
    compute_vertex_normal(h->next()->vertex(), FN1);
    if (FN1==CGAL::NULL_VECTOR) return FEATURE_COS_DIHEDRAL_ANGLE; // true by default
    compute_vertex_normal(h->opposite()->next()->vertex(), FN2);
    if (FN2==CGAL::NULL_VECTOR) return FEATURE_COS_DIHEDRAL_ANGLE; // true by default

    return FN1*FN2;
}

static inline double cosVertexNormals(const Halfedge_handle h)
{
    Vector FN1, FN2;
    compute_vertex_normal(h->vertex(), FN1);
    if (FN1==CGAL::NULL_VECTOR) return FEATURE_COS_DIHEDRAL_ANGLE; // true by default
    compute_vertex_normal(h->prev()->vertex(), FN2);
    if (FN2==CGAL::NULL_VECTOR) return FEATURE_COS_DIHEDRAL_ANGLE; // true by default

    return FN1*FN2;
}

static inline double dihedralAngle(const Halfedge_handle h)
{
    if (is_border(h))
    {
        return M_PI;
    }

    return acos(cosDihedralAngle(h));
}

static inline double meanDihedralAngle(const Vertex_handle pVertex)
{
    double sum=0.0;
    double cpt=0.0;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        sum+=dihedralAngle(h);
        //sum+=dihedralAngle(h->prev());
        cpt+=1.0;
    }

    return sum/cpt;
}

static inline double maxDihedralAngle(const Vertex_handle pVertex)
{
    double maxd=-M_PI;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        double tmp=dihedralAngle(h);
        if(maxd<tmp) maxd=tmp;
    }

    return maxd;
}


static inline double meanCurvatureAtEdge(const Halfedge_handle h)
{
    Vector V = h->vertex()->point() - h->prev()->vertex()->point();
    double dih_angle = dihedralAngle(h), h_len = std::sqrt(V*V);

    //cout << "= " << 2.0*h_len*cos(0.5*dih_angle) << endl;

    return 2.0*h_len*cos(0.5*dih_angle);
}

static inline double meanCurvatureAtEdgeVector(const Halfedge_handle h, Vector& VMC)
{
    Vector Nedge;
    double meanC = meanCurvatureAtEdge(h);
    ////////////////////////////////////////
    compute_edge_normal(h, Nedge);

    VMC = meanC*Nedge;

    return meanC;
}

// directional measure
static inline double normalCurvatureEstimateAtEdge(const Halfedge_handle h)
{
    Vector V = h->vertex()->point() - h->prev()->vertex()->point();
    Vector N;
    compute_vertex_normal(h->vertex(), N);
    return 2.0/(V*V) * (V*N);
}
/////////////////////////////////////////////////////////////////////////////////////
// the notion of strongness is very sensitive to noise!!!
static inline bool thetaStrong(const Halfedge_handle h, double theta)
{
    return cosDihedralAngle(h)<=cos(theta);
}

static inline bool cosThetaStrong(const Halfedge_handle h, double cosTheta)
{
    return cosDihedralAngle(h)<=cosTheta;
}

static inline double cosTurningEdgeAngle(Vector& V1, Vector& V2)
{
    normalize_Vector(V1);
    normalize_Vector(V2);

    return V1*V2;
}

// tangential angle in between 2 edges E and G
static inline double cosTurningEdgeAngle(const Halfedge_handle pE, const Halfedge_handle pG)
{
    Vector  e = pE->vertex()->point() - pE->prev()->vertex()->point(),
            g = pG->vertex()->point() - pG->prev()->vertex()->point();

    normalize_Vector(e);
    normalize_Vector(g);

    // direction correction for g
    if(pE->vertex() == pG->vertex() || pE->prev()->vertex() == pG->prev()->vertex())
    {
        g=-g;
    }

    return e*g; // cosine of the tangential angle
}

static inline double turningEdgeAngle(const Halfedge_handle pE, const Halfedge_handle pG)
{
    return acos(cosTurningEdgeAngle(pE,pG));
}

// more convenient for identifying feature edges:
static inline double dirOneSidedCosTurningDirKmin(const Halfedge_handle pE)
{
    Vector  e = pE->vertex()->point() - pE->prev()->vertex()->point(),
            g = pE->vertex()->closestObs()->VKmin;

    normalize_Vector(e);
    normalize_Vector(g);

    return abs(g*e); // abs cosine of the tangential angle with min principal curcature direction
}

static inline double meanDirOneSidedCosTurningDirKmin(const Vertex_handle pVertex)
{
    double sum=0.0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        sum+=dirOneSidedCosTurningDirKmin(cHe);
    }

    return sum/valence(pVertex);
}
// non-directional:
static inline double oneCosTurningDirKmin(const Halfedge_handle pE)
{
    return std::min<double>(dirOneSidedCosTurningDirKmin(pE), dirOneSidedCosTurningDirKmin(pE->opposite()));
}

// Find the more aligned adjacent edge and return the turning angle
static inline double cosTurningEdgeAngleToTheStartorEnd(const Halfedge_handle pE,
                                                        bool start,
                                                        Halfedge_handle& pE_out)
{
    double maxCos = -1.0, /*maxCosF = -1.0, */cosTA;
    Halfedge_around_vertex_circulator cHe = (start)?pE->vertex()->vertex_begin():pE->prev()->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        if(Halfedge_handle(cHe)==pE)
            continue;

        cosTA = cosTurningEdgeAngle(pE,cHe);
        if(cosTA>maxCos)
        {
            maxCos=cosTA;
            pE_out = cHe;
        }

        //if(cHe->isFeature() && (cosTA>maxCosF)) maxCosF=cosTA;
    }

    //if(maxCosF>-1.0)
    //    return maxCosF;
    //else
        return maxCos;
}

// this notion of weighted face angle is introduced in "Feature detection for surface meshes" (2002)
static inline double weightedFaceAngle( const Halfedge_handle pMainEdge,
                                        const Halfedge_handle pNeighboringEdge)
{
    return abs(cosTurningEdgeAngle(pNeighboringEdge,pMainEdge))*dihedralAngle(pMainEdge);
}
// the notion of strongness is less sensitive to noise!!!
static inline bool strongWRTNeighbor(   const Halfedge_handle pMainEdge,
                                        const Halfedge_handle pNeighboringEdge,
                                        double theta,   // in [0,pi]
                                        double r        // r>=1
                                        )
{
    assert(r>=1.0);
    assert((0.0<=theta) && (M_PI>=theta));

    //if(dihedralAngle(pMainEdge) < turningEdgeAngle(pMainEdge,pNeighboringEdge)) return false;
    if(cosDihedralAngle(pMainEdge) > cosTurningEdgeAngle(pMainEdge,pNeighboringEdge)) return false;

    double wge = weightedFaceAngle(pMainEdge,pNeighboringEdge);
    Vertex_handle pVertex;
    // we assume that pMainEdge and pNeighboringEdge share a vertex!!
    if( (pMainEdge->vertex()==pNeighboringEdge->vertex()) ||
        (pMainEdge->vertex()==pNeighboringEdge->prev()->vertex()))
    {
        pVertex = pMainEdge->vertex();
    }
    else
        pVertex = pMainEdge->prev()->vertex();

    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    double tmpMax=theta;
    CGAL_For_all(cHe,cHe2)
    {
        if( (Halfedge_handle(cHe)==pNeighboringEdge) ||
            cHe->opposite()==pNeighboringEdge ||
            (Halfedge_handle(cHe)==pMainEdge) ||
            cHe->opposite()==pMainEdge
          )
          continue;

        double tmprwhe = r*weightedFaceAngle(cHe,pNeighboringEdge);
        if(tmpMax<tmprwhe) tmpMax = tmprwhe;
        // strong local maxima property:
        if(wge<tmpMax) return false;
    }

    return true;
}

static inline bool urlStrong(   const Halfedge_handle pMainEdge,
                                const Halfedge_handle pNeighboringEdge,
                                double theta_u,   // in [0,pi]
                                double theta_l,   // in [0,pi]
                                double r          // r>=1
                                        )
{
    return thetaStrong(pMainEdge, theta_u) || strongWRTNeighbor(pMainEdge, pNeighboringEdge, theta_l, r);
}

static inline bool urlStronger( const Halfedge_handle pMainEdge,
                                const Halfedge_handle pNeighboringEdge,
                                const Halfedge_handle pNeighboringEdgeRef,
                                double theta_u,   // in [0,pi]
                                double theta_l,   // in [0,pi]
                                double r          // r>=1
                                        )
{
    if( thetaStrong(pMainEdge, theta_u) &&
        (cosDihedralAngle(pMainEdge)<=cosDihedralAngle(pNeighboringEdge))
      ) return true;

    return strongWRTNeighbor(pMainEdge, pNeighboringEdgeRef, theta_l, r);
}

static inline bool strongInOneRing( const Halfedge_handle pMainEdge,
                                    double theta,
                                    double r=1.0         // r>=1
                                        )
{
    Halfedge_around_vertex_circulator cHe = pMainEdge->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        if(!strongWRTNeighbor(pMainEdge, cHe, theta, r)) return false;
    }

    cHe = pMainEdge->prev()->vertex()->vertex_begin();
    cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        if(!strongWRTNeighbor(pMainEdge, cHe, theta, r)) return false;
    }

    return true;
}

// Angle defect is closely related to the Gaussian curvature of smooth surfaces as both as their
// magnitudes measure how a surface deviates from being flat at a point. For this reason
// angle defect is also called simplicial curvature.
static inline double angleDefect(const Vertex_handle pVertex)
{
    double sum=0.0, b=1.0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        sum+=computeTriangleAngle(cHe->prev()->vertex()->point(),
                                  cHe->vertex()->point(),
                                  cHe->next()->vertex()->point());
    }


    if(!is_border(pVertex)) b+=1.0; // if it is not a boundary vertex, the difference is with 2pi!!

    return b*M_PI - sum;
}

static inline double maxAngleDefect(const Halfedge_handle pEdge)
{
    return std::max<double>(angleDefect(pEdge->vertex()), angleDefect(pEdge->prev()->vertex()));
}

static inline double meanAngleDefect(const Halfedge_handle pEdge)
{
    return 0.5*(angleDefect(pEdge->vertex()) + angleDefect(pEdge->prev()->vertex()));
}

static inline double minAngleDefect(const Halfedge_handle pEdge)
{
    return std::min<double>(angleDefect(pEdge->vertex()), angleDefect(pEdge->prev()->vertex()));
}

static inline bool isFeatureAngleDefect(const Vertex_handle pVertex)
{
    return abs(angleDefect(pVertex)) > 0.0;
}

static inline double gaussianCurvatureApproximation(const Vertex_handle pVertex)
{
    double one_r_area = oneRingArea(pVertex);
    if(abs(one_r_area)>1e-16)
        return 3.0*angleDefect(pVertex)/one_r_area;
    else
        return 3.0*angleDefect(pVertex);
}

static inline double maxGaussianCurvatureApproximation(const Halfedge_handle pEdge)
{
    Vector V = pEdge->vertex()->point() - pEdge->prev()->vertex()->point();
    double h_len = std::sqrt(V*V);
    return std::max<double>(gaussianCurvatureApproximation(pEdge->vertex()),
							gaussianCurvatureApproximation(pEdge->prev()->vertex()))*h_len;
}

static inline double meanGaussianCurvatureApproximation(const Halfedge_handle pEdge)
{
    Vector V = pEdge->vertex()->point() - pEdge->prev()->vertex()->point();
    double h_len = std::sqrt(V*V);
    return 0.5*(gaussianCurvatureApproximation(pEdge->vertex()) +
                gaussianCurvatureApproximation(pEdge->prev()->vertex()))*h_len;
}

static inline double minGaussianCurvatureApproximation(const Halfedge_handle pEdge)
{
    Vector V = pEdge->vertex()->point() - pEdge->prev()->vertex()->point();
    double h_len = std::sqrt(V*V);
    return std::min<double>(gaussianCurvatureApproximation(pEdge->vertex()),
							gaussianCurvatureApproximation(pEdge->prev()->vertex()))*h_len;
}

static inline bool isAGoodFeatureEdge(const Halfedge_handle h)
{
#ifndef HYSTERESIS_TWO_THRESHOLDS
    cout << "isAGoodFeatureEdge cannot be called without the hysteresis!!" << endl;
    assert(false);
#endif
    //return strongInOneRing(h, acos(FEATURE_COS_DIHEDRAL_ANGLE));

    double cosDA = cosDihedralAngle(h);
    // angle beteween the two normals greater than 35 degrees (a quite discriminative threshold)
    return (cosDA <= FEATURE_COS_DIHEDRAL_ANGLE);
}

static inline bool isAGoodFeatureEdgeOnVertexNeighborhood(const Vertex_handle pVertex)
{
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;

    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if (isAGoodFeatureEdge(cHe))
        {
            return true;
        }
    }
    return false;
}

static inline bool isAGoodFeatureEdgeOnVertexTwoRingsNeighborhood(const Vertex_handle pVertex)
{
    Halfedge_around_vertex_circulator cH = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe = cH;

    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cH,cHe)
    {
        Halfedge_around_vertex_circulator cH2 = cH->prev()->vertex()->vertex_begin();
        Halfedge_around_vertex_circulator cHe2 = cH2;
        CGAL_For_all(cH2,cHe2)
        {
            if (isAGoodFeatureEdge(cH2))
            {
                return true;
            }
        }
    }
    return false;
}

static inline bool isFeatureEdge(const Halfedge_handle h)
{
    double cosDA = cosDihedralAngle(h);
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    if (cosDA <= -0.94f) return false; // fold-over case (we do not consider a fold-over as a feature edge)

    if (cosDA <= FEATURE_COS_DIHEDRAL_ANGLE)
        return true;
#ifdef HYSTERESIS_TWO_THRESHOLDS
    else if (cosDA >= (1.0+FEATURE_COS_DIHEDRAL_ANGLE)*0.5) // < ?? degrees
        return false;
    else
    {
        // We can improve this by using an interval on which we do not know the answer: [0.82;0.93]
        // : if we are in this interval we return true if there are 2 good feature edges linked to this edge otherwise we return false
        // this approach implies the use of face normal and not vertex normal since the last one is not enough discriminative
        return (isAGoodFeatureEdgeOnVertexNeighborhood(h->vertex()) && isAGoodFeatureEdgeOnVertexNeighborhood(h->prev()->vertex()));
        //return (isAGoodFeatureEdgeOnVertexTwoRingsNeighborhood(h->vertex()) && isAGoodFeatureEdgeOnVertexTwoRingsNeighborhood(h->prev()->vertex()));
    }
#else
    else return false;
#endif
}

static inline void setFeatureEdgesFromThreshold(PolyhedronPtr pMesh)
{
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
        {
            if ( isFeatureEdge(pEdge) )
            {
                pEdge->SetFeature(true);
                pEdge->opposite()->SetFeature(true);
            }
            else
            {
                pEdge->SetFeature(false);
                pEdge->opposite()->SetFeature(false);
            }
        }
}

// this function check and update feature edges geometrically
static inline void setVertexEdgeNeighborhoodToNoFeature(const Vertex_handle pVertex)
{
    //unsigned int nb = 0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;

    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        cHe->SetFeature(false);
    }
}

static inline bool isASharpEdgeOnVertexNeighborhood(const Vertex_handle pVertex)
{
    //unsigned int nb = 0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;

    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if (isFeatureEdge(cHe))
        {
            return true;
        }
#ifdef FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS
        if(!almostCoplanarTriangles(cHe)) return true;
#endif

    }
    return false;
}

// this function check (and update when asked) feature edges geometrically
static inline bool isFeatureEdgeOnVertexNeighborhood(const Vertex_handle pVertex)
{
    //unsigned int nb = 0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;

    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
#ifdef NO_FEATURE_EDGE_RECOMPUTATION
        if ( cHe->isFeature() )
            return true;
#else
        if (isFeatureEdge(cHe))
        {
            cHe->SetFeature(true);
            return true;
        }
        else if ( cHe->isFeature() ) cHe->SetFeature(false);
#endif
    }
    return false;
}


// a feature vertex is linked to 3 feature edges or is a sharp feature
static bool isFeatureVertex(const Vertex_handle pVertex)
{
    assert( pVertex!=Vertex_handle());
    //if ( is_border(pVertex) ) return true;

    //if(isFeatureAngleDefect(pVertex)) return true; // not very useful

    unsigned int cpt=0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if (cHe->isFeature()) //isFeatureEdge(cHe))
        {
            cpt++;
        }
    }
    // this constraint allow to keep unchanged feature edges
    return ( (valence(pVertex)==2 && cpt==2) || (valence(pVertex)>2 && cpt>2) );
}

static unsigned int nbOfAdjacentFeatureEdges(const Vertex_handle pVertex)
{
    unsigned int cpt=0;
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        if (cHe->isFeature()) //isFeatureEdge(cHe))
        {
            cpt++;
        }
    }

    return cpt;
}

static inline bool isFeatureVertexOnVertexNeighborhood(const Vertex_handle pVertex)
{
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;

    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if (isFeatureVertex(cHe->prev()->vertex()))
        {
            cHe->prev()->vertex()->SetFeature(true);
            return true;
        }
        else if ( cHe->prev()->vertex()->isFeature() ) cHe->prev()->vertex()->SetFeature(false);
    }

    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// methods for listing mesh vertex neighboring halfedges addresses (used before some topological modifications)
// that can be usefull for updating some value after that vertex has been removed
static inline void neighboring_halfedges(const Vertex_handle pVertex, set< void* >& nh)
{
    Halfedge_around_vertex_circulator cHe = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    CGAL_For_all(cHe,cHe2)
    {
        nh.insert((void*)( &(*cHe)));
        nh.insert((void*)( &(*cHe->opposite()))); // we take into account the 2 halfedges because when using Edges, we don't know in advance which halfedge is used
    }
}

static inline void neighboring_halfedges(const Halfedge_handle h, set< void* >& nh)
{
    neighboring_halfedges(h->opposite()->vertex(), nh);
    neighboring_halfedges(h->vertex(), nh);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// methods for detecting mesh vertex neighbors (used before some topological modifications)
static inline bool areNeighbors(const Halfedge_handle h1, const Halfedge_handle h2)
{
    // independence for edge collpase!!

    return (h1->vertex()==h2->vertex() ||  h1->opposite()->vertex()==h2->vertex() ||
            h1->vertex()==h2->opposite()->vertex() || h1->opposite()->vertex()==h2->opposite()->vertex() ||
            h1->next()->vertex()==h2->vertex() || h1->next()->vertex()==h2->opposite()->vertex() ||
            h1->opposite()->next()->vertex()==h2->vertex() || h1->opposite()->next()->vertex()==h2->opposite()->vertex()  );
}

static bool areNeighbors2Rings(const Halfedge_handle h1, const Halfedge_handle h2)
{
    Halfedge_around_vertex_circulator h = h1->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        if (areNeighbors(h, h2)) return true;
    }

    h = h1->opposite()->vertex()->vertex_begin();
    he = h;
    CGAL_For_all(h,he)
    {
        if (areNeighbors(h, h2)) return true;
    }

    return false;
}


static bool areNeighbors(const Vertex_handle pVertex1, const Vertex_handle pVertex2)
{
    Halfedge_around_vertex_circulator h1 = pVertex1->vertex_begin();
    Halfedge_around_vertex_circulator he1 = h1;

    CGAL_For_all(h1,he1)
    {
        Halfedge_around_vertex_circulator h2 = pVertex2->vertex_begin();
        Halfedge_around_vertex_circulator he2 = h2;
        CGAL_For_all(h2,he2)
        {
            if (h1->opposite()==h2)
                return true;
        }
    }

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static double min_edge_length_around(const Vertex_handle pVertexCenter)
{
    Halfedge_around_vertex_circulator pHalfEdge = pVertexCenter->vertex_begin();
    Halfedge_around_vertex_circulator end = pHalfEdge;

    Vector vec = pHalfEdge->prev()->vertex()->point()-pHalfEdge->vertex()->point();
    double Minl = vec*vec;

    if (CGAL::circulator_distance(pHalfEdge,end)==1) return std::sqrt(Minl);

    CGAL_For_all(pHalfEdge,end)
    {
        vec = pHalfEdge->prev()->vertex()->point()-pHalfEdge->vertex()->point();
        double Minl_temp = vec*vec;
        if (Minl_temp<Minl) Minl = Minl_temp;
    }
    return std::sqrt(Minl);
}

// find out the minimal geodesic distance
static double min_geodesic_length_around(const Vertex_handle pVertexCenter)
{
    // we must use the closest observation of the label
    Halfedge_around_vertex_circulator pHalfEdge = pVertexCenter->vertex_begin();
    Halfedge_around_vertex_circulator end = pHalfEdge;

    Vector Hyp;
    double Mingl = 1e8;

    // initialization
    if (CGAL::circulator_distance(pHalfEdge,end)==0)
        return (double)0;
    else if (CGAL::circulator_distance(pHalfEdge,end)==1)
    {
        Hyp = pVertexCenter->point() - pHalfEdge->prev()->vertex()->point();
        Mingl = Hyp*Hyp;
        return std::sqrt(Mingl);
    }

    Mingl = 1e32;//min_edge_length_around(pVertexCenter);
    //Mingl = Mingl*Mingl;
    // circulator clockwise order
    CGAL_For_all(pHalfEdge,end)
    {
        Vector adj, v1_v2;//, inter;
        //             center point                        opposite point
        Hyp = pHalfEdge->vertex()->point() - pHalfEdge->prev()->vertex()->point();
        //             next point                          opposite point
        v1_v2 = pHalfEdge->next()->vertex()->point() - pHalfEdge->prev()->vertex()->point(); // v2 - v1
        // we normalize it
        double v1_v2_l = normalize_Vector(v1_v2);
        double hyp_scal_v1_v2 = v1_v2 * Hyp, Mingl_temp;

        //Point3d newPt = pHalfEdge->prev()->vertex()->point()+hyp_scal_v1_v2 * v1_v2;

        if (hyp_scal_v1_v2<0.0)
        {
            // rentre à l'intérieur par la gauche
            Mingl_temp = Hyp*Hyp;
        }
        else if (hyp_scal_v1_v2>v1_v2_l)
        {
            // rentre à l'intérieur par la droite
            Mingl_temp =  (pHalfEdge->vertex()->point() - pHalfEdge->next()->vertex()->point()) *
                          (pHalfEdge->vertex()->point() - pHalfEdge->next()->vertex()->point());
        }
        else
        {
            adj = hyp_scal_v1_v2 * v1_v2;
            Mingl_temp = Hyp*Hyp-adj*adj; // I have checked that we get the same result by direct distance between center and adjacent points
        }

        if (Mingl_temp<1e-15) // in case of numerical imprecision
            Mingl_temp = (double)0;

        if (Mingl_temp<Mingl) Mingl = Mingl_temp;
    }
    //printf("min geodesic length: %f\n",std::sqrt(Mingl));
    return std::sqrt(Mingl);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inverted element detection
static bool fold(const Halfedge_handle h, const Halfedge_handle g)
{
    // Variables declaration
    Vector next_h, Nnext, LeftDir = CGAL::NULL_VECTOR, Vg, next_g, RightDir = CGAL::NULL_VECTOR, Vh;
    //////////////////////////////////////////////////////////////////
    next_h = h->vertex()->point() - h->next()->vertex()->point();
    compute_vertex_normal(h->next()->vertex(),Nnext);
    if (next_h!=CGAL::NULL_VECTOR)
        LeftDir = CGAL::cross_product(next_h,Nnext);
    else
        return true;
    //////////////////////////////////////////////////////////////////
    next_g = g->vertex()->point() - g->next()->vertex()->point();
    compute_vertex_normal(g->next()->vertex(),Nnext);
    if (next_g!=CGAL::NULL_VECTOR)
        RightDir = CGAL::cross_product(next_g,Nnext);
    else
        return true;
    Vg = g->vertex()->point() - h->next()->vertex()->point();
    Vh = h->vertex()->point() - g->next()->vertex()->point();
    // fold to the right     fold to the left
    return (LeftDir*Vg<=0.0f     || RightDir*Vh<=0.0f);
}

static inline bool isASurroundingGeometricVertexFold(const Vertex_handle pVertex, bool is_already_computed)
{
    // Barycentric coordinates cannot be used for detecting fold-over (otherwise
    // true-false (If the center vertex is below or above its one-ring vertices)
    // or false-true(in case of non-convex neighborhood without special care for
    // the local retriangulation) ) can be returned.
    Vector N1, N2, tmp;
    Vector D, left;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    if(!is_already_computed)
    {
        /*
        // this piece of code source do not improve the results while increasing computation time
        unsigned int vertex_valence = valence(pVertex);

        if (vertex_valence<=2) return false;

        if ((vertex_valence==3) &&
                isInTriangle(   h->prev()->vertex()->point(),
                                h->opposite()->next()->vertex()->point(),
                                h->next()->vertex()->point(), pVertex->point())) return false;
        */



        /*
            if(isvertexOneRingConvex(pVertex))
            {
                // this piece of source code makes the simplification based on edge collapse not working
                // that is because the one-ring neighborhood is not necessarily convex
                Point3d refPointTriangles=h->prev()->vertex()->point(); // for fold-over detection
                h++;
                he--;// take care of that!!
                CGAL_For_all(h,he)
                {
                    Point3d h_prev_point = h->prev()->vertex()->point(),
                            h_suiv_point = h->next()->vertex()->point();
                    Vector T_N = triangleNormal(h_prev_point, refPointTriangles, h_suiv_point);

                    Point3d projected_pos = CGAL::ORIGIN +  (pVertex->point()-CGAL::ORIGIN)-
                                                            ((pVertex->point()-CGAL::ORIGIN)*T_N)*T_N;
                    // case in which we are sure there is no fold-over
                    if(isInTriangle(h_prev_point,
                                    refPointTriangles,
                                    h_suiv_point,
                                    projected_pos
                                    )
                        )
                        return false;
                }
            }
        */

        // Now we are in the case where the point normal projection is outside its one-ring neighborhood.
        // and therefore we may have a fold-over present
        // this test alone find more vertex fold-over than the ground true, but it
        // experimentally does not miss one true vertex fold-over
        he = h = pVertex->vertex_begin();
        CGAL_For_all(h,he)
        {
            compute_planar_facet_normal(h->opposite()->facet(), N1);
            compute_planar_facet_normal(h->facet(), N2);
            // to avoid creating too many fold-overs during simplification: we use -0.6 instead of -0.94
            if (N1*N2< -0.6f) // inverted neighboring normals -0.94 (take -0.3 to avoid creating very bad configuration if not care is taken on the volume error)
                return true;
        }
    }
    else
    {
        CGAL_For_all(h,he)
        {
            if (h->prev()->vertex()->isFoldover()) // inverted neighboring normals -0.94 (take -0.3 to avoid creating very bad configuration if not care is taken on the volume error)
                return true;
        }
    }

    return false;
}

static inline bool isInvolvedInGeometricVertexFold(const Halfedge_handle h, bool is_already_computed)
{
    if(h==Halfedge_handle() || is_border(h)) return false;

    if(!is_already_computed)
    {
        Vector N1, N2;
        compute_planar_facet_normal(h->opposite()->facet(), N1);
        compute_planar_facet_normal(h->facet(), N2);
        return (N1*N2< -0.94f); // inverted neighboring normals -0.94 (take -0.3 to avoid creating very bad configuration if not care is taken on the volume error)
    }
    else
    {
        return (h->vertex()->isFoldover() || h->prev()->vertex()->isFoldover()); // inverted neighboring normals -0.94
    }
}

static inline bool isInvolvedInGeometricVertexFold( const vector< Halfedge_handle >& vh,
                                                    bool is_already_computed)
{
    vector< Halfedge_handle >::const_iterator it(vh.begin()), ite(vh.end());

    for(;it!=ite;++it)
    {
        if(isInvolvedInGeometricVertexFold(*it, is_already_computed)) return true;
    }

    return false;
}

// for triangular facet!!
static inline bool isInvolvedInGeometricVertexFold(const Facet_handle pFacet, bool is_already_computed)
{
    if( pFacet==Facet_handle() ) return false;

    return  isInvolvedInGeometricVertexFold(pFacet->halfedge(), is_already_computed) ||
            isInvolvedInGeometricVertexFold(pFacet->halfedge()->prev(), is_already_computed) ||
            isInvolvedInGeometricVertexFold(pFacet->halfedge()->next(), is_already_computed);
}

static inline bool isGeometricVertexFold(const Vertex_handle pVertex)
{
    // a border/feature vertex cannot be responsible for a vertex fold (because it is considered to be a feature)
    if ( //is_border(pVertex) || // a border vertex will automatically have a surrounding feature edge
            isFeatureEdgeOnVertexNeighborhood(pVertex)
       ) return false;

    unsigned int cpt_triangle_geom_fold=0;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        if ( isInvolvedInGeometricVertexFold(h->facet(), false) ) cpt_triangle_geom_fold++;
    }

    return (cpt_triangle_geom_fold>=2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// smoothing operators applied directly
static void vertexOneRingCentroid(Vertex_handle pVertex) // Laplacian smoothing
{
    Vector vec = CGAL::NULL_VECTOR;
    double sum=0.0, wj;
    Halfedge_around_vertex_circulator ha = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator hae = ha;
    CGAL_For_all(ha,hae)
    {
        double areaAi, areaAi_1;
        areaAi = triangleArea(ha->prev()->vertex()->point(), ha->opposite()->next()->vertex()->point(), ha->vertex()->point());
        areaAi_1 = triangleArea(ha->vertex()->point(), ha->next()->vertex()->point(), ha->prev()->vertex()->point());
        // local area dispersion (not divided by the area of the one-ring):
        wj = areaAi * cotangent(ha->prev()->vertex()->point(), ha->opposite()->next()->vertex()->point(), ha->vertex()->point())+
             areaAi_1 * cotangent(ha->vertex()->point(), ha->next()->vertex()->point(), ha->prev()->vertex()->point());

        // to get barycentric coordinates:
        vec = vec + (ha->prev()->vertex()->point() - pVertex->point())*wj;

        sum += wj;
    }
    if (sum) vec = vec*(1.0/sum);

    Vector N;
#ifndef SMOOTHING
#ifdef INIT_CURVATURE_SENSITIVE
    // MORE FIDELITY!!!!!
    // we do that to not changed the NRJ being minimized (here we are proposing a position for the current candidate and the closest observation of the new candidate will be updated when computing the distance)
    N = (pVertex->closestObs())->normal(); // to be more sensitive to the initial curvature
    //N = CGAL::cross_product(pVertex->closestObs()->VKmin, pVertex->closestObs()->VKmax); // too much geometric error
#else
    // more geometric error, since we use the intial surface as reference
    compute_vertex_normal(pVertex, N); // to avoid introducing noise in flat region
#endif
    vec = vec - (vec*N)*N;
#endif
    // Update the vertex position consequently
    if (pVertex->isFoldover())
    {
        if(!isASharpEdgeOnVertexNeighborhood(pVertex))
        {
            pVertex->point() = pVertex->point() + vec; // try to correct the fold-over (without any guarantee)
            if (!isGeometricVertexFold(pVertex) && !isASharpEdgeOnVertexNeighborhood(pVertex) )
                pVertex->SetFoldover(false);
            else
                pVertex->point() = pVertex->point() - vec; // we let the vertex position unchanged
        }
    }
    else
    {
        pVertex->point() = pVertex->point() + SMOOTHING_FACTOR*vec; // Pold + lamda*umbrella(Pold)
    }
}

// Increases the min angle and decreases the max angle
static void vertexOneRingAnglesBasedCentroid(Vertex_handle pVertex) // Angle based smoothing
{
    Vector vec = CGAL::NULL_VECTOR;
    double sumIangles = 0.0f, alpha1, alpha2, Iangle;
    Halfedge_around_vertex_circulator ha = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator hae = ha;
    CGAL_For_all(ha,hae)
    {
        // The angle must be divided in two to take convex and concave cases into account
        alpha1 = computeTriangleAngle(ha->vertex()->point(), ha->prev()->vertex()->point(), ha->opposite()->next()->vertex()->point());
        alpha2 = computeTriangleAngle(ha->next()->vertex()->point(), ha->prev()->vertex()->point(), ha->vertex()->point());
        //cout << "alpha1 = " << alpha1 << " alpha2 = " << alpha2 << endl;
        if (alpha1<0.0 || alpha2<0.0) assert(false);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Vector d1,d2,c,dc;
        d1 = ha->opposite()->next()->vertex()->point() - ha->prev()->vertex()->point();
        normalize_Vector(d1);
        d2 = ha->next()->vertex()->point() - ha->prev()->vertex()->point();
        normalize_Vector(d2);
        c = 0.5*(d1+d2);
        // locally non-convex!
        if ( fold(ha->opposite()->next(),ha->next()) ) continue; //c = -c; (NOT GOOD AT ALL)
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if ( c != CGAL::NULL_VECTOR )
        {
            c = c*(1.0/std::sqrt(c*c));
            Point3d New_ideal_point = ha->prev()->vertex()->point() + c * std::sqrt((ha->vertex()->point() - ha->prev()->vertex()->point())*(ha->vertex()->point() - ha->prev()->vertex()->point()));
            // Calculate the difference between two adjacent angles...
            // Small angles (alpha1+alpha2) must be favorized to avoid foldover
            Iangle = alpha1+alpha2;
            Iangle = 1.0/(Iangle*Iangle);
            dc = (New_ideal_point - pVertex->point());
            //normalize_Vector(dc); // unit direction towards the bissector of that angle
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // the current best results:
            vec = vec + dc * Iangle;
            sumIangles += Iangle;
        }
    }
    if (sumIangles) vec = vec*(1.0/sumIangles);

    Vector N;

#ifndef SMOOTHING
#ifdef INIT_CURVATURE_SENSITIVE
    // MORE FIDELITY!!!!!
    // we do that to not changed the NRJ being minimized
    N = (pVertex->closestObs())->normal(); // to be more sensitive to the initial curvature
    //N = CGAL::cross_product(pVertex->closestObs()->VKmin, pVertex->closestObs()->VKmax); // too much geometric error
#else
    // more geometric error, since we use the intial surface as reference
    compute_vertex_normal(pVertex, N); // to avoid introducing noise in flat region
#endif
    vec = vec - (vec*N)*N;
#endif
    // Update the vertex position consequently
    if (pVertex->isFoldover())
    {
        if(!isASharpEdgeOnVertexNeighborhood(pVertex))
        {
            pVertex->point() = pVertex->point() + vec; // try to correct the fold-over (without any guarantee)
            if (!isGeometricVertexFold(pVertex) && !isASharpEdgeOnVertexNeighborhood(pVertex) )
                pVertex->SetFoldover(false);
            else
                pVertex->point() = pVertex->point() - vec; // we let the vertex position unchanged
        }
    }
    else
    {
        pVertex->point() = pVertex->point() + SMOOTHING_FACTOR*vec; // Pold + lamda*umbrella(Pold)
    }
}

// in the same spirit, but it create a center vertex
static bool createCenterVertex(PolyhedronPtr pMesh, Facet_handle fh, Halfedge_handle& NewOne)  // return true if the center vertex has been added without pb
{
    if ( is_border(fh->halfedge()) ) return false;

    Vertex_handle vOb = fh->halfedge()->vertex()->closestObs();
    Halfedge_around_facet_circulator ha = fh->facet_begin();
    Halfedge_around_facet_circulator hae = ha;
    unsigned int nb=0;
    Vector vec = CGAL::NULL_VECTOR;
    CGAL_For_all(ha,hae)
    {
        vec = vec + (ha->vertex()->point() - CGAL::ORIGIN);
        nb++;
    }
    if (nb) vec = vec*(1.0/double(nb));

    NewOne = pMesh->create_center_vertex(fh->halfedge());
    NewOne->vertex()->point() = CGAL::ORIGIN + vec;

    // we do not want to create either a fold-over or a feature edge
    if (    isGeometricVertexFold( NewOne->vertex() )
            || isASharpEdgeOnVertexNeighborhood(NewOne->vertex())
            || (computeMinTriangleAngles(NewOne->vertex()) < 25.0/180.*M_PI)
            || (computeMaxTriangleAngles(NewOne->vertex()) > 95.0/180.*M_PI)
       )
    {
        pMesh->erase_center_vertex (NewOne);
        return false;
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    setVertexEdgeNeighborhoodToNoFeature(NewOne->vertex());
    ///////////////////////////////////////////////////////////////////////////////////////
    // the data structure must be updated consequently:
    //AddvertexProp(NewOne->vertex(), vOb);
    //SetvertexProp(NewOne->vertex(), searchNearestObsForLabel(NewOne->vertex()));
    NewOne->vertex()->SetClosestObs(vOb);
    NewOne->vertex()->SetClosestObs(searchClosestObs(NewOne->vertex()));
    ///////////////////////////////////////////////////////////////////////////////////////

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// topological operators

static inline bool isAflipableEdge(PolyhedronPtr pMesh, Halfedge_handle h)
{
    if (h->isFeature() //isFeatureEdge(h)
            ||
            is_border(h) // border halfedge (isFeature may be false because of fold-over detection
       ) return false;

#ifdef ASSERTION
    // The edge is between two triangles
    assert(h == h->next()->next()->next());
    assert(h->opposite() == h->opposite()->next()->next()->next());
#endif

    // we cannot flip if the new edge will be out of the four vertices polygon
    if (areNeighbors(h->next()->vertex(), h->opposite()->next()->vertex() ) // necessary but not sufficient for non-flipable edge
            ||
            fold(h->next(), h->opposite()->next())
       ) return false;

    // We do not flip a bad local configuration
    if ( isInvolvedInGeometricVertexFold(h, true) ||
            isInvolvedInGeometricVertexFold(h->next(), true) ||
            isInvolvedInGeometricVertexFold(h->prev(), true) ||
            isInvolvedInGeometricVertexFold(h->opposite()->next(), true) ||
            isInvolvedInGeometricVertexFold(h->opposite()->prev(), true) ) return false;

    // we cannot flip an edge if we may change the feature edge detection result (for instance,
    // by changing the edge direction alignement)
    /*if ( isFeatureEdgeOnVertexNeighborhood(h->vertex()) ||
            isFeatureEdgeOnVertexNeighborhood(h->next()->vertex()) ||
            isFeatureEdgeOnVertexNeighborhood(h->opposite()->vertex()) ||
            isFeatureEdgeOnVertexNeighborhood(h->opposite()->next()->vertex()) ) return false;*/

    // we flip the initial position to make sure we will not introduce to much geometric error

    h = pMesh->flip_edge (h);
    if(isFeatureEdge(h))
    {
        h = pMesh->flip_edge (h);
        return false;
    }
    else
    {
        h = pMesh->flip_edge (h);
//#ifdef FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS
//        if(!almostCoplanarTriangles(h)) return false;
//#endif
        return true;
    }

}


static void edgesSurroundingEdgeCollpase(   const Halfedge_handle h,
                                            vector<Halfedge_handle> & edges,
                                            bool inlude_direct_edges
                                        )
{
    edges.clear(); // we empty the initial vector

    Halfedge_around_vertex_circulator ha = h->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator hae = ha;
    if(!is_border(h))
    {
        CGAL_For_all(ha,hae)
        {
            if ( (Halfedge_handle(ha)!=h)
                    &&(inlude_direct_edges ||(Halfedge_handle(ha)!=h->opposite()->prev())
                        &&(Halfedge_handle(ha)!=h->next()->opposite())  )
               )
            {
                edges.push_back(ha);
            }
        }

        ha = h->opposite()->vertex()->vertex_begin();
        hae = ha;
        CGAL_For_all(ha,hae)
        {
            if ( (Halfedge_handle(ha)!=h->opposite())
                    &&(inlude_direct_edges|| (Halfedge_handle(ha)!=h->prev())
                        &&(Halfedge_handle(ha)!=h->opposite()->next()->opposite()))
               )
            {
                edges.push_back(ha);
            }
        }
    }
    else if(h->facet()!=Facet_handle())
    {
        CGAL_For_all(ha,hae)
        {
            if ( (Halfedge_handle(ha)!=h)
                    &&(inlude_direct_edges ||Halfedge_handle(ha)!=h->next()->opposite())
               )
            {
                edges.push_back(ha);
            }
        }

        ha = h->opposite()->vertex()->vertex_begin();
        hae = ha;
        CGAL_For_all(ha,hae)
        {
            if ( (Halfedge_handle(ha)!=h->opposite())
                    &&(inlude_direct_edges || Halfedge_handle(ha)!=h->prev())
               )
            {
                edges.push_back(ha);
            }
        }
    }
    else
    {
        CGAL_For_all(ha,hae)
        {
            if ( (Halfedge_handle(ha)!=h)
                    &&(inlude_direct_edges ||Halfedge_handle(ha)!=h->opposite()->prev())
               )
            {
                edges.push_back(ha);
            }
        }

        ha = h->opposite()->vertex()->vertex_begin();
        hae = ha;
        CGAL_For_all(ha,hae)
        {
            if ( (Halfedge_handle(ha)!=h->opposite())
                    &&(inlude_direct_edges ||Halfedge_handle(ha)!=h->opposite()->next()->opposite())
               )
            {
                edges.push_back(ha);
            }
        }
    }
}

// This method does not take into account any triangle for one of the 2 h vertices
// if the vertex valence is <=3
static void trianglesSurroundingEdgeCollpase(const Halfedge_handle h, vector<Halfedge_handle> & triangles)
{
    triangles.clear(); // we empty the initial vector

    Halfedge_around_vertex_circulator ha = h->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator hae = ha;

    if (!is_border(h))
    {
        CGAL_For_all(ha,hae)
        {
            if ( Halfedge_handle(ha)!=h && ha->next()!=h->opposite()
               )
            {
                triangles.push_back(ha);
            }
        }

        ha = h->opposite()->vertex()->vertex_begin();
        hae = ha;
        CGAL_For_all(ha,hae)
        {
            if ( ha->opposite()!=h && ha->next()!=h  )
            {
                triangles.push_back(ha);
            }
        }
    }
    else if(h->facet()!=Facet_handle())
    {
        CGAL_For_all(ha,hae)
        {
            if ( Halfedge_handle(ha)!=h )
            {
                triangles.push_back(ha);
            }
        }

        ha = h->opposite()->vertex()->vertex_begin();
        hae = ha;
        CGAL_For_all(ha,hae)
        {
            if ( ha->opposite()!=h && ha->next()!=h  )
            {
                triangles.push_back(ha);
            }
        }
    }
    else
    {
        CGAL_For_all(ha,hae)
        {
            if ( Halfedge_handle(ha)!=h && ha->next()!=h->opposite()
               )
            {
                triangles.push_back(ha);
            }
        }

        ha = h->opposite()->vertex()->vertex_begin();
        hae = ha;
        CGAL_For_all(ha,hae)
        {
            if ( ha->opposite()!=h )
            {
                triangles.push_back(ha);
            }
        }
    }
}

static inline bool isOneDifferentFeatureEdge(const Halfedge_handle h)
{
#ifdef ASSERTION
    assert(pVertex != Vertex_handle());
#endif

    if (is_border(h)) return false;

    Halfedge_around_vertex_circulator cHe = h->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if ((Halfedge_handle(cHe)!=h) && cHe->isFeature())
        {
            return true;
        }

#ifdef FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS
        if((Halfedge_handle(cHe)!=h) && !almostCoplanarTriangles(cHe)) return true;
#endif
    }

    cHe = h->opposite()->vertex()->vertex_begin();
    cHe2 = cHe;
    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if ((cHe->opposite()!=h) && cHe->isFeature())
        {
            return true;
        }

#ifdef FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS
        if((cHe->opposite()!=h) && !almostCoplanarTriangles(cHe)) return true;
#endif
    }

    return false;
}



static inline bool isAcollapsableEdge(PolyhedronPtr pMesh, Halfedge_handle h)
{
    // we collapse short edges only
    //if (!isEasyShortEdge(pMesh, h)) return false; // we want to collapse not too much!!!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( h->vertex()->isFeature() && h->prev()->vertex()->isFeature() ) return false;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS
    if(!almostCoplanarTriangles(h)) return false;
#endif
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // We do not collapse a bad local configuration
    vector< Halfedge_handle > surrounding_edges;
    edgesSurroundingEdgeCollpase(h, surrounding_edges, true);
    if( isInvolvedInGeometricVertexFold(surrounding_edges, true) ) return false;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // case of a neighboring feature edge different from h (we do not want to make them move)
    // + we cannot collapse an edge if we may change the feature edge detection result (for instance,
    // by changing the edge direction alignement)
    //if(isOneDifferentFeatureEdge(h)) return false; // this check is now done in edge_collapse

    return true;
}

static inline bool noFeatureEdgeInOtherDir( const vector<Halfedge_handle> & edges, // do not contain ref_h
                                            const Halfedge_handle ref_h)
{
    vector<Halfedge_handle>::const_iterator iter = edges.begin(), itere=edges.end();
    for (;iter!=itere;++iter)
    {
        if ( (*iter)->isFeature()// || isFeatureEdge(*iter)
            )
        {
            Point3d a, b, c;
            b = (*iter)->vertex()->point();
            c = (*iter)->prev()->vertex()->point();

            if( (ref_h->vertex()==(*iter)->vertex()) ||
                (ref_h->vertex()==(*iter)->prev()->vertex())
              )
            {
                a = ref_h->prev()->vertex()->point();
            }
            else if((ref_h->prev()->vertex()==(*iter)->vertex()) ||
                    (ref_h->prev()->vertex()==(*iter)->prev()->vertex())
              )
            {
                a = ref_h->vertex()->point();
            }
            else
            {
                assert(false);
                return false; // pb detected
            }

            if( !CGAL::collinear(a, b, c))
                return false;
            //cout << "cosTurningEdgeAngle(*iter, ref_h) = " << cosTurningEdgeAngle(*iter, ref_h) << endl;
            //if(abs(cosTurningEdgeAngle(*iter, ref_h))<0.98) return false;
        }
    }

    return true;
}

static inline double deltaTriangleShapePotentialCollapse(   PolyhedronPtr pMesh,
                                                            Halfedge_handle h,
                                                            Point3d& collapsePos // output collapse pos
                                                            )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    collapsePos = h->vertex()->point();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    unsigned int nbaf_h=nbOfAdjacentFeatureEdges(h->vertex()), nbaf_hp=nbOfAdjacentFeatureEdges(h->prev()->vertex());

    // to make sure that after an edge collapse the feature edge direction is let unchanged!!
    if( (nbaf_h>=3)  && (nbaf_hp>=2) || (nbaf_h>=2)  && (nbaf_hp>=3)
        ) return PROHIB_ENERGY;

    if( (nbaf_h==2) && (nbaf_hp==2) )
    {
        vector<Halfedge_handle> edgesForDir;
        edgesSurroundingEdgeCollpase(h, edgesForDir, true);

        if(!noFeatureEdgeInOtherDir(edgesForDir, h)) return PROHIB_ENERGY;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Vector vec = (h->opposite()->vertex()->point() - h->vertex()->point());
    Point3d hPos = h->vertex()->point(), hopPos = h->opposite()->vertex()->point(), NewPos;
    vector<Halfedge_handle> triangles_sur_EC;
    trianglesSurroundingEdgeCollpase(h, triangles_sur_EC);

    double initShapePot, ShapePot;
    initShapePot = ShapePot = triangleShapePotential(triangles_sur_EC);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // initial bad configuration
    if(initShapePot>=PROHIB_ENERGY) return PROHIB_ENERGY;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FORBID_TOO_GEOMETRICALY_DEGRADING_TOPOLOGICAL_OPERATIONS // here it degrades the results
    if(isOneDifferentFeatureEdge(h)) return PROHIB_ENERGY;
#endif
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( (nbaf_h>0) && (nbaf_hp==0) )
    {
        h->vertex()->point() = h->opposite()->vertex()->point() = collapsePos;
        ShapePot = triangleShapePotential(triangles_sur_EC);
    }
    else if( (nbaf_hp>0) && (nbaf_h==0) )
    {
       collapsePos = h->prev()->vertex()->point();
       h->vertex()->point() = h->opposite()->vertex()->point() = collapsePos;
       ShapePot = triangleShapePotential(triangles_sur_EC);
    }
    else
    {
        unsigned short nbCandidates = 5, i = 0;
        while (i<=nbCandidates)
        {
            NewPos = hPos+vec*double(i)/nbCandidates;

            h->vertex()->point() = h->opposite()->vertex()->point() = NewPos;
            double ShapePotAfterCollapse = triangleShapePotential(triangles_sur_EC);
            if ( (ShapePotAfterCollapse < ShapePot)
               )
            {
                // if we find out a better configuration, we have to check whether there is no foldover
                h->vertex()->point() = hPos;
                if (!isASurroundingGeometricVertexFold(h->opposite()->vertex(), false) &&
                        !isASurroundingGeometricVertexFold(h->next()->vertex(), false) &&
                        !isASurroundingGeometricVertexFold(h->opposite()->next()->vertex(), false)
                   )
                {
                    h->vertex()->point() = NewPos;
                    h->opposite()->vertex()->point() = hopPos;
                    if (!isASurroundingGeometricVertexFold(h->vertex(), false) /*&&
                        !isASurroundingGeometricVertexFold(h->next()->vertex(), false) &&
                        !isASurroundingGeometricVertexFold(h->opposite()->next()->vertex(), false)*/)
                    {
                        collapsePos = NewPos;
                        ShapePot = ShapePotAfterCollapse; // we update the new best configuration score
                    }
                }
            }
            ++i;
        }
    }

    h->vertex()->point() = hPos;
    h->opposite()->vertex()->point() = hopPos;

    return ShapePot - (initShapePot+triangleShapePotential(h));
}


// compute mixed area in the one-ring neighborhood
static double mixedArea(const Vertex_handle pVertex)
{
    double mArea = 0.0;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
		  if(isObtuseTriangle(  h->prev()->vertex()->point(),
                                h->vertex()->point(),
                                h->next()->vertex()->point())   )
          {
              if(isObtuseAngle( h->prev()->vertex()->point(),
                                h->vertex()->point(),
                                h->next()->vertex()->point()) )
              {
                    mArea+=0.5*triangleArea(h->prev()->vertex()->point(),
                                            h->vertex()->point(),
                                            h->next()->vertex()->point());
              }
              else
              {
                    mArea+=0.25*triangleArea(   h->prev()->vertex()->point(),
                                                h->vertex()->point(),
                                                h->next()->vertex()->point());
              }
          }
          else
          { // for non-obtuse triangle
              mArea+=voronoiAreaTriangle(   h->prev()->vertex()->point(),
                                            h->vertex()->point(),
                                            h->next()->vertex()->point());
          }
	}

	return mArea;
}

static double meanCurvatureAtVertex(const Vertex_handle pVertex)
{
    Vector V(0,0,0);
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        double  cotaij=cotangent(   h->vertex()->point(),
                                    h->next()->vertex()->point(),
                                    h->prev()->vertex()->point()),
                cotbij=cotangent(   h->prev()->vertex()->point(),
                                    h->opposite()->next()->vertex()->point(),
                                    h->vertex()->point());

        V=V+(cotaij+cotbij)*(h->vertex()->point()-h->prev()->vertex()->point());
	}

    //return 0.25*std::sqrt(V*V); //unormalized
	return 0.25*std::sqrt(V*V)/mixedArea(pVertex);
}

/*

    sphere center <=> edge midpoint
    sphere radius <=> radius
    ax+by+cz+d=0 <=> plane equation such that above plane (positive signed distance)
                     is to the right and bellow plane is to the left (negative signed distance)
*/
static void local_edge_region_growing(  const Halfedge_handle pEdge,
                                        double radius,
                                        std::vector < Vertex_handle > &left_region_v,
                                        std::vector < Vertex_handle > &right_region_v,
                                        std::vector < Halfedge_handle > &left_region_h,
                                        std::vector < Halfedge_handle > &right_region_h,
                                        std::vector < Facet_handle > &left_region_f,
                                        std::vector < Facet_handle > &right_region_f
                                        )
{
    //////////////////////////////////////////////////////////////////////////////////////////
    left_region_v.clear();
    right_region_v.clear();
    left_region_h.clear();
    right_region_h.clear();
    left_region_f.clear();
    right_region_f.clear();
    //////////////////////////////////////////////////////////////////////////////////////////
    if(is_border(pEdge)) return; // for the time being we do not process border edges
    //////////////////////////////////////////////////////////////////////////////////////////
    // clipping sphere
    Point3d separating_plane_p = midpoint(pEdge->vertex()->point(), pEdge->prev()->vertex()->point());

    radius = std::max<double>(radius,2.0*std::sqrt( std::min<double>(pEdge->length2(),
													std::min<double>(   pEdge->prev()->length2(),
																		pEdge->next()->length2()    ) )));
    //////////////////////////////////////////////////////////////////////////////////////////
    // separating plane equation
    Vector Ne, separating_plane_N;
    compute_edge_normal(pEdge, Ne);
    separating_plane_N = CGAL::cross_product(pEdge->vertex()->point()-pEdge->prev()->vertex()->point(),Ne);
    //////////////////////////////////////////////////////////////////////////////////////////
    //std::stack< Vertex_handle > pile_v;
    //pile_v.push(h->vertex());
    //pile_v.push(h->prev()->vertex());
    ///////////////////////////////////////////////////
    // there is one key/priority associated to each element (which is represented by an id)
    // several elements can have the same priority
    //                  key(==priority)  val(==id)     key_comp           val_comp
    mutable_priority_queue<double, void*, std::less<double>, std::greater<void*> > mut_p_q;
    mutable_priority_queue<double, void*, std::less<double>, std::greater<void*> >::iterator iter_mut_p_q;
    map< void*, Vertex_handle, std::greater<void*> > selected_vectors;
    mut_p_q.insert((void*)( &(*pEdge->vertex()) ), 0.0);
    selected_vectors[(void*)&(*pEdge->vertex())]=pEdge->vertex();
    mut_p_q.insert((void*)( &(*pEdge->prev()->vertex()) ), 0.0);
    selected_vectors[(void*)&(*pEdge->prev()->vertex())]=pEdge->prev()->vertex();
    //////////////////////////////////////////////////////////////////////////////////////////
    // void* is the generic address type!
    std::set< void* > treated_Vertex;
    std::set< void* > treated_Halfedge;
    std::set< void* > treated_facet;

    treated_Vertex.insert(&(*pEdge->vertex()));
    treated_Vertex.insert(&(*pEdge->prev()->vertex()));

    treated_Halfedge.insert(&(*pEdge));
    treated_Halfedge.insert(&(*pEdge->opposite()));

    treated_facet.insert(&(*pEdge->facet()));
    treated_facet.insert(&(*pEdge->opposite()->facet()));

    Vector N_proxy_right, N_proxy_left, N_tmp;

    // we will keep trace of rejected facets during the region growing process
    // in order to test again them with the last proxy representation
    std::vector < Facet_handle > left_region_f_rejected, right_region_f_rejected;

    double thr_accept = 0.86;// above 0.7, it is very sensitive to the scan of triangles
                             // if higher threshold is needed, we recommend to use
                             // a priority queue (priority = distance to separating plane)
                             // instead of a stack
    if(point_to_plane_signed_distance(separating_plane_N, separating_plane_p, pEdge->prev()->prev()->vertex()->point()) >= 0.0)
    {
        right_region_f.push_back(pEdge->facet());
        compute_planar_facet_normal(pEdge->facet(), N_proxy_right);
    }
    else
    {
        left_region_f.push_back(pEdge->facet());
        compute_planar_facet_normal(pEdge->facet(), N_proxy_left);
    }

    if(point_to_plane_signed_distance(separating_plane_N, separating_plane_p, pEdge->opposite()->prev()->prev()->vertex()->point()) >= 0.0)
    {
        right_region_f.push_back(pEdge->opposite()->facet());
        compute_planar_facet_normal(pEdge->opposite()->facet(), N_proxy_right);
    }
    else
    {
        left_region_f.push_back(pEdge->opposite()->facet());
        compute_planar_facet_normal(pEdge->opposite()->facet(), N_proxy_left);
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    unsigned int cpt_v=0, cpt_h=0, cpt_f=0;
    //////////////////////////////////////////////////////////////////////////////////////////
    while(!mut_p_q.empty())
    {
        //std::cout << "priority: " << mut_p_q.begin()->first << endl;
        Vertex_handle v = selected_vectors[mut_p_q.begin()->second];
        //pile_v.top() ;
        //pile_v.pop() ;
        mut_p_q.erase(mut_p_q.begin()->second);
        //cout << "Stack size = " << pile_v.size() << endl;
        /////////////////////////////////////////////////////////////
        bool atLeastOneGoodTriangle=false;
        std::vector< Vertex_handle > treated_Vertex_tmp;
        Halfedge_around_vertex_circulator h = v->vertex_begin();
        Halfedge_around_vertex_circulator he = h;
        CGAL_For_all(h,he)
        {
            if( isInSphere(separating_plane_p, radius, h->prev()->vertex()) &&
               (treated_Vertex.find(&(*h->prev()->vertex()))==treated_Vertex.end())
               )
            {
                treated_Vertex_tmp.push_back(h->prev()->vertex());

/*
                mut_p_q.insert( (void*)( &(*h->prev()->vertex())),
                                std::abs(point_to_plane_signed_distance(separating_plane_N,
                                                                        separating_plane_p,
                                                                        h->prev()->vertex()->point())));
                selected_vectors[(void*)&(*h->prev()->vertex())]=h->prev()->vertex();
                treated_Vertex.insert( &(*h->prev()->vertex()) );

                cpt_v++;
*/

                // the facet is a new one and it is present at least in part in the
                // local sphere
                if(point_to_plane_signed_distance(separating_plane_N, separating_plane_p, h->prev()->vertex()->point()) >= 0.0)
                {
                    right_region_v.push_back(h->prev()->vertex());
                }
                else
                {
                    left_region_v.push_back(h->prev()->vertex());
                }
            }

            if( isInSphere(separating_plane_p, radius, h) &&
               (treated_Halfedge.find(&(*h))==treated_Halfedge.end())
               )
            {
                treated_Halfedge.insert(&(*h->prev()));
                treated_Halfedge.insert(&(*h->prev()->opposite()));
                cpt_h++;
                // the facet is a new one and it is present at least in part in the
                // local sphere
                if(point_to_plane_signed_distance(separating_plane_N, separating_plane_p, midpoint(h->vertex()->point(),h->prev()->vertex()->point())) >= 0.0)
                {
                    right_region_h.push_back(h);
                }
                else
                {
                    left_region_h.push_back(h);
                }
            }

            if( (h->facet()!=Facet_handle()) &&
                isInSphere(separating_plane_p, radius, h->facet()) &&
               (treated_facet.find(&(*h->facet()))==treated_facet.end())
               )
            {
                treated_facet.insert(&(*h->facet()));
                cpt_f++;
                compute_planar_facet_normal(h->facet(), N_tmp);
                // the facet is a new one and it is present at least in part in the
                // local sphere
                if(point_to_plane_signed_distance(separating_plane_N, separating_plane_p,
                    barycenter( h->vertex()->point(),
                                h->prev()->vertex()->point(),
                                h->next()->vertex()->point())) >= 0.0)
                {
                    if(N_proxy_right*N_tmp>thr_accept)
                    {
                        atLeastOneGoodTriangle=true;
                        right_region_f.push_back(h->facet());
                        compute_region_normal(right_region_f, N_proxy_right);
                    }
                    else
                    {
                        right_region_f_rejected.push_back(h->facet());
                    }
                }
                else
                {
                    if(N_proxy_left*N_tmp>thr_accept)
                    {
                        atLeastOneGoodTriangle=true;
                        left_region_f.push_back(h->facet());
                        compute_region_normal(left_region_f, N_proxy_left);
                    }
                    else
                    {
                        left_region_f_rejected.push_back(h->facet());
                    }
                }
            }
        }

        if(atLeastOneGoodTriangle) // to have a continuous region to the right and to the left
        {
            std::vector< Vertex_handle >::iterator it(treated_Vertex_tmp.begin()), ite(treated_Vertex_tmp.end());

            for(;it!=ite;++it)
            {
                mut_p_q.insert( (void*)( &(**it) ),
                                std::abs(point_to_plane_signed_distance(separating_plane_N,
                                                                        separating_plane_p,
                                                                        (*it)->point())));
                selected_vectors[(void*)&(**it)]=*it;
                treated_Vertex.insert( &(**it) );

                cpt_v++;
            }
        }

    }

    // we will do a second pass on rejected facets
    std::vector < Facet_handle >::iterator  it(left_region_f_rejected.begin()),
                                            ite(left_region_f_rejected.end());
    for(; it!=ite; ++it)
    {
        compute_planar_facet_normal(*it, N_tmp);

        if(N_proxy_left*N_tmp>thr_accept)
        {
            left_region_f.push_back(*it);
            compute_region_normal(left_region_f, N_proxy_left);
        }
    }

    it = right_region_f_rejected.begin();
    ite = right_region_f_rejected.end();
    for(; it!=ite; ++it)
    {
        compute_planar_facet_normal(*it, N_tmp);

        if(N_proxy_right*N_tmp>thr_accept)
        {
            right_region_f.push_back(*it);
            compute_region_normal(right_region_f, N_proxy_right);
        }
    }

    //cout << "Nb vertices = " << cpt_v << " Nb edges = " << cpt_h << " Nb triangles = " << cpt_f << endl;
}

static void setFeatureEdgesFromSurroundingProxies(PolyhedronPtr pMesh, double radius)
{
    for(   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
    {
        std::vector < Vertex_handle > left_region_v;
        std::vector < Vertex_handle > right_region_v;
        std::vector < Halfedge_handle > left_region_h;
        std::vector < Halfedge_handle > right_region_h;
        std::vector < Facet_handle > left_region_f;
        std::vector < Facet_handle > right_region_f;

        local_edge_region_growing(  pEdge,
                                    radius,
                                    left_region_v, right_region_v, left_region_h,
                                    right_region_h, left_region_f, right_region_f
                                );

        Vector N_right, N_left;
        compute_region_normal(right_region_f, N_right);
        compute_region_normal(left_region_f, N_left);

        if(N_right*N_left < FEATURE_COS_DIHEDRAL_ANGLE)
        {
            pEdge->SetFeature(true);
            pEdge->opposite()->SetFeature(true);
        }
        else
        {
            pEdge->SetFeature(false);
            pEdge->opposite()->SetFeature(false);
        }
    }
}

static inline double regionBasedCosDihedralAngle(   const Halfedge_handle h,
                                                    double radius)
{
    // A border edge is a feature edge
    if (is_border(h))
    {
        return -1; // cos(M_PI): this choice is imposed by the community (see "Feature detection for surface meshes")
    }
    //////////////////////////////////////////////////////////////////////////////////
    Vector FN1, FN2;
    std::vector < Vertex_handle > left_region_v;
    std::vector < Vertex_handle > right_region_v;
    std::vector < Halfedge_handle > left_region_h;
    std::vector < Halfedge_handle > right_region_h;
    std::vector < Facet_handle > left_region_f;
    std::vector < Facet_handle > right_region_f;

    local_edge_region_growing(  h,
                                radius,
                                left_region_v, right_region_v, left_region_h,
                                right_region_h, left_region_f, right_region_f
                            );

    compute_region_normal(right_region_f, FN1);
    compute_region_normal(left_region_f, FN2);

    return FN1*FN2;
}


static inline bool levelOfMaximalityReached(const Vertex_handle pVertex,
                                            const Vertex_handle pVertexinVKmax,
                                            const Vertex_handle pVertexinMinusVKmax,
                                            double maximality
                                            )
{
    double kvmax = std::max<double>( abs(pVertex->Kmax), abs(pVertex->Kmin) );
    double kvmaxinVKmax = std::max<double>( abs(pVertexinVKmax->Kmax), abs(pVertexinVKmax->Kmin) );
    double kvmaxinMinusVKmax = std::max<double>( abs(pVertexinMinusVKmax->Kmax), abs(pVertexinMinusVKmax->Kmin) );

    kvmax = kvmax*kvmax;
    kvmaxinVKmax = kvmaxinVKmax*kvmaxinVKmax;
    kvmaxinMinusVKmax = kvmaxinMinusVKmax*kvmaxinMinusVKmax;

    return (kvmax-kvmaxinVKmax > maximality) && (kvmax-kvmaxinMinusVKmax > maximality);
}

/**
    This method tries to identify local crest points according to the definition
    in the paper "Crest lines extraction from 3D triangular meshes"
**/
static inline bool isLocalCrestPoint(const Vertex_handle pVertex, double maximality)
{
    if( valence(pVertex) <= 2 ) return true;

    Vertex_handle pVertexinVKmax, pVertexinMinusVKmax;

    Vector DirK = (abs(pVertex->Kmax)>=abs(pVertex->Kmin))?pVertex->VKmax:pVertex->VKmin;

    double maxg = -1.0, ming = 1.0;

    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        Vector Dirtmp = h->prev()->vertex()->point() - pVertex->point();
        normalize_Vector(Dirtmp);

        if(Dirtmp*DirK>maxg)
        {
            maxg = Dirtmp*DirK;
            pVertexinVKmax = h->prev()->vertex();
        }

        if(Dirtmp*DirK<ming)
        {
            ming = Dirtmp*DirK;
            pVertexinMinusVKmax = h->prev()->vertex();
        }
    }

    return levelOfMaximalityReached(pVertex, pVertexinVKmax, pVertexinMinusVKmax, maximality);
}

static inline bool significantEdge(const Halfedge_handle pEdge)
{
    if( !pEdge->isFeature() ||
        (   !isLocalCrestPoint(pEdge->vertex(), 500.0) &&
            !isLocalCrestPoint(pEdge->prev()->vertex(), 500.0)
        )
       ) return false;

    return true;
}

static inline bool shouldRemoveNonSignificantLines( const Vertex_handle pVertex,
                                                    double maximality)
{
    if( !isLocalCrestPoint(pVertex, maximality) ) return false;

    unsigned int nbpointsmin=3;

    unsigned int cpt=1;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        if( isLocalCrestPoint(h->prev()->vertex(), maximality) )
        {
            cpt++;

            Halfedge_around_vertex_circulator h2 = h->prev()->vertex()->vertex_begin();
            Halfedge_around_vertex_circulator he2 = h2;
            CGAL_For_all(h2,he2)
            {
                if( (h2->prev()->vertex()!=pVertex) &&
                    isLocalCrestPoint(h2->prev()->vertex(), maximality) ) cpt++;
            }

            if(cpt >= nbpointsmin) return false;
        }
    }

    return true;
}

/**
    DOES NOT WORK: UNDER REVISION
**/
static inline bool belongsToSignificantLines(const Halfedge_handle pEdge)
{
    if(!pEdge->isFeature()) return false;

    unsigned int cpt=1;

    /*
    vector < Halfedge_handle > edgestofollow;
    edgestofollow.push_back(pEdge);
    while(cpt<=2)
    {
        Halfedge_handle currenth = edgestofollow.back(); // take the last element
        edgestofollow.pop_back(); // remove the last element
        ///////////////////////////////////////////////////////////////////////////
        Halfedge_around_vertex_circulator h = currenth->vertex()->vertex_begin();
        Halfedge_around_vertex_circulator he = h;
        CGAL_For_all(h,he)
        {
            if( (Halfedge_handle(h)!=currenth) &&
                (h->opposite()!=currenth) &&
                    h->isFeature()
              )
            {
                cpt++;
                edgestofollow.push_back(h);
            }
        }
    }
    */

    Halfedge_around_vertex_circulator h = pEdge->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        if( (Halfedge_handle(h)!=pEdge) &&
                h->isFeature()
          )
        {
            cpt++;

            Halfedge_around_vertex_circulator h2 = h->prev()->vertex()->vertex_begin();
            Halfedge_around_vertex_circulator he2 = h2;
            CGAL_For_all(h2,he2)
            {
                if( (Halfedge_handle(h2)!=h->opposite()) &&
                        h2->isFeature()
                  )
                {
                    cpt++;
                }
            }
        }
    }

    h = pEdge->prev()->vertex()->vertex_begin();
    he = h;
    CGAL_For_all(h,he)
    {
        if((h->opposite()!=pEdge) &&
                h->isFeature()
          )
        {
            cpt++;

            Halfedge_around_vertex_circulator h2 = h->prev()->vertex()->vertex_begin();
            Halfedge_around_vertex_circulator he2 = h2;
            CGAL_For_all(h2,he2)
            {
                if( (Halfedge_handle(h2)!=h->opposite()) &&
                        h2->isFeature()
                  )
                {
                    cpt++;
                }
            }
        }
    }

    return (cpt>2);
}

/**
    This method tries to remove non-local maxima edges in the regionBasedCosDihedralAngle sense.

**/
static inline bool isLocalMaximaAmongOthers(const Halfedge_handle pEdge,
                                            double radius)
{

    if(!pEdge->isFeature()) return false;
    else if(is_border(pEdge)) return true;
    else if(!isOneDifferentFeatureEdge(pEdge) )
    {   // we remove isolated feature edges
        pEdge->SetFeature(false);
        pEdge->opposite()->SetFeature(false);
        return false;
    }
    else if(isFeatureVertex(pEdge->vertex()) ||
            isFeatureVertex(pEdge->prev()->vertex()) )
    { // we let feature edges around a corner unchanged
        return true;
    }
    //else if(significantEdge(pEdge)) return true; // significant edges are considered as local maxima
    else if(!belongsToSignificantLines(pEdge))
    {   // we remove isolated feature edges
        pEdge->SetFeature(false);
        pEdge->opposite()->SetFeature(false);
        return false;
    }


    double refAD = regionBasedCosDihedralAngle(pEdge, radius);
    double costhreshold = 0.6;
    Halfedge_around_vertex_circulator h = pEdge->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        if( (Halfedge_handle(h)!=pEdge) && h->isFeature()

            && (cosTurningEdgeAngle(pEdge,h)<costhreshold) // non-aligned!
          )
        {
            double tmpAD = regionBasedCosDihedralAngle(h, radius);
            if((tmpAD<refAD))
            { // pEdge is not a local maxima!
                pEdge->SetFeature(false);
                pEdge->opposite()->SetFeature(false);
                return false; // we remove isolated feature edges!!!
            }
        }
    }

    he = h = pEdge->opposite()->vertex()->vertex_begin();
    CGAL_For_all(h,he)
    {
        if( (Halfedge_handle(h)!=pEdge->opposite()) && h->isFeature()

            && (cosTurningEdgeAngle(pEdge->opposite(),h)<costhreshold) // non-aligned!
          )
        {
            double tmpAD = regionBasedCosDihedralAngle(h, radius);
            if(tmpAD<refAD)
            {
                pEdge->SetFeature(false);
                pEdge->opposite()->SetFeature(false);
                return false; // we remove isolated feature edges!!!
            }
        }
    }

    return true;
}

static inline void removeNonMaximaFeatureEdges(PolyhedronPtr pMesh, double radius)
{
    unsigned int cpt=0, lastcpt=1;
    while( (cpt!=lastcpt) )
    {
        cpt = lastcpt = 0;
        for(   Edge_iterator pEdge = pMesh->edges_begin();
                pEdge != pMesh->edges_end();
                pEdge++)
        {
            if(!isLocalMaximaAmongOthers(pEdge, radius)) lastcpt++;
        }

        for(   Edge_iterator pEdge = pMesh->edges_begin();
                pEdge != pMesh->edges_end();
                pEdge++)
        {
            if(!isLocalMaximaAmongOthers(pEdge, radius)) cpt++;
        }
    }
}

/**
    isEndFeatureLineEdge returns true if the current edge is a feature edge
    located at the end of a feature line. In that case DirLookFor is a unit direction
    vector which goes towards the direction closest to the "ideal next feature edge".
**/
static inline bool isEndFeatureLineEdge(const Halfedge_handle pEdge, Vector &DirLookFor)
{
    DirLookFor = CGAL::NULL_VECTOR;
    if(!pEdge->isFeature()) return false;

    /////////////////////////////////////////////////////////////////////////////
    Halfedge_handle TheOtherFeatureEdge = Halfedge_handle();
    bool F=false;
    /////////////////////////////////////////////////////////////////////////////
    unsigned int cpt=0;
    Halfedge_around_vertex_circulator cHe = pEdge->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator cHe2 = cHe;
    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if ((Halfedge_handle(cHe)!=pEdge) && cHe->isFeature())
        {
            TheOtherFeatureEdge=Halfedge_handle(cHe);
            F=true;
            cpt++;
        }
    }

    cHe = pEdge->prev()->vertex()->vertex_begin();
    cHe2 = cHe;
    // circulator clockwise oriented (because the halfedges are
    // counter-clockwise oriented inside a facet)
    CGAL_For_all(cHe,cHe2)
    {
        if ((cHe->opposite()!=pEdge) && cHe->isFeature())
        {
            TheOtherFeatureEdge=Halfedge_handle(cHe);
            cpt++;
        }
    }

    if(cpt==1)
    {
        if(F)
            DirLookFor = pEdge->prev()->vertex()->point() - TheOtherFeatureEdge->prev()->vertex()->point();
        else
            DirLookFor = pEdge->vertex()->point() - TheOtherFeatureEdge->prev()->vertex()->point();

        normalize_Vector(DirLookFor);
    }

    return (cpt==1);
}

/**
    findClosestSurroundingNonFeatureEdgeInDir is called after isEndFeatureLineEdge
    to select the "ideal next feature edge" when isEndFeatureLineEdge returns true.
    That can be used for feature lines gaps filling.
**/
static inline void findClosestSurroundingNonFeatureEdgeInDir(   const Halfedge_handle pEdge,
                                                                const Vector &Dir,
                                                                Halfedge_handle& outinDir)
{
    assert(pEdge->isFeature()); // check up
    //////////////////////////////////////////////////////////////////////////////////////
    Vector DirCurrentEdge = pEdge->vertex()->point() - pEdge->prev()->vertex()->point();
    double cosTurning, cosTurningTmp;
    Halfedge_around_vertex_circulator cHe,cHe2;
    if(DirCurrentEdge*Dir>=0.0)
    {
        cHe2 = cHe = pEdge->vertex()->vertex_begin();
        outinDir = Halfedge_handle(cHe);
        cosTurning = cosTurningEdgeAngle(pEdge,outinDir);
        // circulator clockwise oriented (because the halfedges are
        // counter-clockwise oriented inside a facet)
        CGAL_For_all(cHe,cHe2)
        {
            if ((Halfedge_handle(cHe)!=pEdge) && !cHe->isFeature()
                )
            {
                cosTurningTmp = cosTurningEdgeAngle(pEdge,Halfedge_handle(cHe));
                if(cosTurningTmp>cosTurning)
                {
                     cosTurning=cosTurningTmp;
                     outinDir = Halfedge_handle(cHe);
                }
            }
        }
    }
    else
    {
        cHe2 = cHe = pEdge->opposite()->vertex()->vertex_begin();
        outinDir = Halfedge_handle(cHe);
        cosTurning = cosTurningEdgeAngle(pEdge,outinDir);
        // circulator clockwise oriented (because the halfedges are
        // counter-clockwise oriented inside a facet)
        CGAL_For_all(cHe,cHe2)
        {
            if ((cHe->opposite()!=pEdge) && !cHe->isFeature())
            {
                cosTurningTmp = cosTurningEdgeAngle(pEdge,cHe->opposite());
                 if(cosTurningTmp>cosTurning)
                 {
                     cosTurning=cosTurningTmp;
                     outinDir = Halfedge_handle(cHe);
                 }
            }
        }
    }
}

static inline void votingGapFiling( PolyhedronPtr pMesh,
                                    std::vector < Facet_handle >& region_f
                                    //double radius
                                    )
{
    Vector DirLookFor;
    //   halfedge address + cpt
    map< void*, unsigned int, std::greater<void*> > selected_edges;
    map< void*, Halfedge_handle, std::greater<void*> > edges;
    for(Edge_iterator pEdge=pMesh->edges_begin();pEdge!=pMesh->edges_end();++pEdge)
    {
        selected_edges[(void*)&(*pEdge)]=0;
        edges[(void*)&(*pEdge)]=pEdge;
    }

    for(Edge_iterator pEdge=pMesh->edges_begin();pEdge!=pMesh->edges_end();++pEdge)
    {
        if(isEndFeatureLineEdge(pEdge, DirLookFor))
        {
            Halfedge_handle outinDir;
            findClosestSurroundingNonFeatureEdgeInDir(pEdge,DirLookFor, outinDir);
            ///////////////////////////////////////////////////////////////////////////////
            selected_edges[(void*)&(*outinDir)]++; // vote for that edge
            ///////////////////////////////////////////////////////////////////////////////
            //region_f.push_back(outinDir->facet());
            //region_f.push_back(outinDir->opposite()->facet());
        }
    }

    map< void*, unsigned int, std::greater<void*> >::const_iterator it(selected_edges.begin()), ite(selected_edges.end());
    for(;it!=ite;++it)
    {
        if(it->second>0)
        {
            cout << it->second << endl;
            /*edges[it->first]->SetFeature(true);
            edges[it->first]->opposite()->SetFeature(true);

            region_f.push_back(edges[it->first]->facet());
            region_f.push_back(edges[it->first]->opposite()->facet());*/
        }
    }

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**********************************************
// principal curvature for a vertextriangleArea
// taken from Guillaume Lavoué
//**********************************************
//**********************************************
// compute v x v^t
//**********************************************
static void vector_times_transpose_mult(double pVector[3],
                                               double ppMatrix[3][3],
                                               double coeff)
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      ppMatrix[i][j] = coeff * pVector[i] * pVector[j];
}

//**********************************************
// add two matrices
//**********************************************
static void add(double pMatrix[3][3],
                       double pMatrixSum[3][3])
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      pMatrixSum[i][j] += pMatrix[i][j];
}

//**********************************************
// fix sine
//**********************************************
static double fix_sine(double sine)
{
  if (sine >= 1)
    return M_PI/2;
  else
    if (sine <= -1)
      return -M_PI/2;
    else
      return std::asin(sine);
}

static bool sphere_clip_vector(Point3d &O, double r,const Point3d &P, Vector &V)
{

    Vector W = P - O ;
    double a = (V*V);
    double b = 2.0 * V * W ;
    double c = (W*W) - r*r ;
    double delta = b*b - 4*a*c ;
    if (delta < 0) {
        // Should not happen, but happens sometimes (numerical precision)
        return true ;
    }
    double t = (- b + ::sqrt(delta)) / (2.0 * a) ;
    if (t < 0.0) {
        // Should not happen, but happens sometimes (numerical precision)
        return true ;
    }
    if (t >= 1.0) {
        // Inside the sphere
        return false ;
    }

    V=V*t;

    return true ;
}

static void principal_curvature_per_vert(   Vertex pVertex,
                                            double ppMatrix_sum[3][3])
{

	double area=0;

  // iterate over all edges
  Halfedge_around_vertex_circulator pHalfedge = pVertex.vertex_begin();
  Halfedge_around_vertex_circulator pHalfedgeStart = pHalfedge;
  CGAL_For_all(pHalfedge,pHalfedgeStart)
  {

    // build edge vector and comput its norm
	  Point3d p1 = pHalfedge->vertex()->point();
	  Point3d p2 = pHalfedge->opposite()->vertex()->point();
	  Vector edge = (p1-p2);
		double len_edge = std::sqrt(edge*edge);
		if (len_edge == 0) // avoid divide by zero
		continue;

    // compute (signed) angle between two incident faces, if exists
    Facet_handle pFacet1 = pHalfedge->facet();
    Facet_handle pFacet2 = pHalfedge->opposite()->facet();
    CGAL_assertion(pFacet1 != pFacet2);
    if (pFacet1 == NULL || pFacet2 == NULL)
      continue; // border edge

	//area+=AreaFacetTriangle(pHalfedge->facet());
	area+=triangleArea(pFacet1);

	Vector normal1 = pFacet1->normal();
	Vector normal2 = pFacet2->normal();

    double sine = (CGAL::cross_product(normal1,normal2)*edge)/len_edge;
    double beta = fix_sine(sine);

    // compute edge * edge^t * coeff, and add it to current matrix
    double pVector_edge[3] = {edge.x(),edge.y(),edge.z()};
    double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    vector_times_transpose_mult(pVector_edge,ppMatrix,beta/len_edge);
    add(ppMatrix,ppMatrix_sum);
  }
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
}

static void geodes_principal_curvature_per_vert(Vertex *pVertex,
                                                double ppMatrix_sum[3][3], double radius)
    {

       std::set<Vertex*> vertices ;
        Point3d O = pVertex->point() ;
        std::stack<Vertex*> S ;
        S.push(pVertex) ;
        vertices.insert(pVertex) ;
		int iter=0;
        while(!S.empty())
		{
			Vertex* v = S.top() ;
            S.pop() ;
            Point3d P = v->point() ;
            Halfedge_around_vertex_circulator h = v->vertex_begin();
			Halfedge_around_vertex_circulator pHalfedgeStart = h;
			CGAL_For_all(h,pHalfedgeStart)
			{
                Point3d p1 = h->vertex()->point();
				Point3d p2 = h->opposite()->vertex()->point();
				Vector V = (p2-p1);
                if (v==pVertex || V * (P - O) > 0.0)
				{
					//double len_old = std::sqrt(V*V);
					bool isect = sphere_clip_vector(O, radius, P, V) ;

                    if (!h->is_border_edge())
					{
						double len_edge = std::sqrt(V*V);
                         // compute (signed) angle between two incident faces, if exists
						Facet_handle pFacet1 = h->facet();
						Facet_handle pFacet2 = h->opposite()->facet();
						CGAL_assertion(pFacet1 != pFacet2);
						if (pFacet1 == NULL || pFacet2 == NULL)
							continue; // border edge
						Vector normal1 = pFacet1->normal();
						Vector normal2 = pFacet2->normal();

						double sine = (CGAL::cross_product(normal1,normal2)*V)/len_edge;
						double beta = fix_sine(sine);

                         // compute edge * edge^t * coeff, and add it to current matrix


						double pVector_edge[3] = {V.x(),V.y(),V.z()};
						double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
						vector_times_transpose_mult(pVector_edge,ppMatrix,beta/len_edge);
						add(ppMatrix,ppMatrix_sum);
                    }


                    if (!isect) {

						Vertex_iterator w=h->opposite()->vertex();
                        if (vertices.find(&(*w)) == vertices.end()) {
                            vertices.insert(&(*w)) ;
                            S.push(&(*w)) ;
                        }
                    }
				}

			}
			iter++;
		}

		double area=M_PI*radius*radius;
		for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
    }

static void principal_curvature(PolyhedronPtr pMesh,
                                bool IsGeod,
                                double radius,
                                double &MinNrmMinCurvature,
                                double &MaxNrmMinCurvature,
                                double &MinNrmMaxCurvature,
                                double &MaxNrmMaxCurvature)
{
	MinNrmMinCurvature=100000;
	MaxNrmMinCurvature=-100000;

	MinNrmMaxCurvature=100000;
	MaxNrmMaxCurvature=-100000;


  Vertex_iterator pVertex = NULL;
  for (pVertex = pMesh->vertices_begin();
      pVertex != pMesh->vertices_end();
      pVertex++)
  {
    double ppMatrix_sum[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	double eigenvalues[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	double eigenvectors[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

    if (IsGeod==true)//voisinage geodesique
		geodes_principal_curvature_per_vert((&(*pVertex)),ppMatrix_sum,radius);
	else//voisinage topologique 1-ring
		principal_curvature_per_vert(*pVertex,ppMatrix_sum);

	//valeurs propres
			double **CovMat=(double**)malloc((4)*sizeof(double*));
			double **VectPro=(double**)malloc((4)*sizeof(double*));
			double **Valpro=(double**)malloc((4)*sizeof(double*));
			for (int i=0;i<(4);i++)
			{
				CovMat[i]=(double*)malloc((4)*sizeof(double));
				VectPro[i]=(double*)malloc((4)*sizeof(double));
				Valpro[i]=(double*)malloc((4)*sizeof(double));
			}

			CovMat[1][1]=ppMatrix_sum[0][0];
			CovMat[1][2]=ppMatrix_sum[0][1];
			CovMat[1][3]=ppMatrix_sum[0][2];
			CovMat[2][1]=ppMatrix_sum[1][0];
			CovMat[2][2]=ppMatrix_sum[1][1];
			CovMat[2][3]=ppMatrix_sum[1][2];
			CovMat[3][1]=ppMatrix_sum[2][0];
			CovMat[3][2]=ppMatrix_sum[2][1];
			CovMat[3][3]=ppMatrix_sum[2][2];

			//la matrice n'est elle pas déja diagonale?
			if (ppMatrix_sum[0][1]==0 && ppMatrix_sum[0][2]==0 &&
				ppMatrix_sum[1][0]==0 && ppMatrix_sum[1][2]==0 &&
				ppMatrix_sum[2][1]==0 && ppMatrix_sum[2][0]==0)
			{
				for (int i=1;i<4;i++)
					for (int j=1;j<4;j++)
					{
						Valpro[i][j]=CovMat[i][j];
						if (i==j)
						VectPro[i][j]=1;
						else
						VectPro[i][j]=0;


					}



			}
			else
			{
				//recherche des vecteurs et valeurs propres
				if (ValPro(3,CovMat,1e-15,10000.,VectPro,Valpro)==-1)
				{
					pVertex->VKmax=CGAL::NULL_VECTOR;
					pVertex->VKmin=CGAL::NULL_VECTOR;
					return;
				}
			}
				//  Call the Jacovi subroutine
			for (int u=0;u<4;u++)
				for (int v=0;v<4;v++)
				{
					Valpro[u][v]=fabs(Valpro[u][v]);

				}
			EigSrt(Valpro,VectPro,3);
			Vector VKmax(VectPro[1][2],VectPro[2][2],VectPro[3][2]);
			Vector VKmin(VectPro[1][1],VectPro[2][1],VectPro[3][1]);


			eigenvalues[0][0]=Valpro[1][1];
			eigenvalues[0][1]=Valpro[1][2];
			eigenvalues[0][2]=Valpro[1][3];
			eigenvalues[1][0]=Valpro[2][1];
			eigenvalues[1][1]=Valpro[2][2];
			eigenvalues[1][2]=Valpro[2][3];
			eigenvalues[2][0]=Valpro[3][1];
			eigenvalues[2][1]=Valpro[3][2];
			eigenvalues[2][2]=Valpro[3][3];

			eigenvectors[0][0]=VectPro[1][1];
			eigenvectors[0][1]=VectPro[1][2];
			eigenvectors[0][2]=VectPro[1][3];
			eigenvectors[1][0]=VectPro[2][1];
			eigenvectors[1][1]=VectPro[2][2];
			eigenvectors[1][2]=VectPro[2][3];
			eigenvectors[2][0]=VectPro[3][1];
			eigenvectors[2][1]=VectPro[3][2];
			eigenvectors[2][2]=VectPro[3][3];




			pVertex->VKmax=VKmax;
			pVertex->VKmin=VKmin;

			pVertex->Kmax=Valpro[1][1];
			pVertex->Kmin=Valpro[2][2];


			for (int i=0;i<(3);i++)
			{
				free(CovMat[i]);
				free(VectPro[i]);
				free(Valpro[i]);
			}
			free(CovMat);
			free(VectPro);
			free(Valpro);


#ifdef _MSC_VER
			MinNrmMinCurvature=min(MinNrmMinCurvature,pVertex->Kmin);
			MaxNrmMinCurvature=max(MaxNrmMinCurvature,pVertex->Kmin);

			MinNrmMaxCurvature=min(MinNrmMaxCurvature,pVertex->Kmax);
			MaxNrmMaxCurvature=max(MaxNrmMaxCurvature,pVertex->Kmax);
#else
			MinNrmMinCurvature=CGAL::min(MinNrmMinCurvature,pVertex->Kmin);
			MaxNrmMinCurvature=CGAL::max(MaxNrmMinCurvature,pVertex->Kmin);

			MinNrmMaxCurvature=CGAL::min(MinNrmMaxCurvature,pVertex->Kmax);
			MaxNrmMaxCurvature=CGAL::max(MaxNrmMaxCurvature,pVertex->Kmax);
#endif

  }




}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



#endif // Polyhedron_local_modification_H
