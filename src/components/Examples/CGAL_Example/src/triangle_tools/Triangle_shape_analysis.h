///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Triangle_shape_analysis_H
#define Triangle_shape_analysis_H


#include "../../../../../mepp/Polyhedron/polyhedron.h"

#include "../useful_methods/extract_Vpropres.h" // for principal curvature computation (copied from curvature component)

#include "../useful_methods/Math_methods.h"

#ifndef PROHIB_ENERGY
    #define PROHIB_ENERGY 99999999
#endif

static inline double trianglePerimeter(const Point3d& a, const Point3d& b, const Point3d& c)
{
    Vector ab(b-a), bc(c-b), ca(a-c);
    return (std::sqrt(ab*ab)+std::sqrt(bc*bc)+std::sqrt(ca*ca));
};

static inline Vector triangleNormalUnnormalized(const Point3d& a, const Point3d& b, const Point3d& c)
{
    Vector c_a(a-c), c_b(b-c);
    return CGAL::cross_product(c_a, c_b);
};

static inline Vector triangleNormal(const Point3d& a, const Point3d& b, const Point3d& c)
{
    Vector N = triangleNormalUnnormalized(a, b, c);
    normalize_Vector(N);

    return N;
};

static inline bool almostCoplanarTriangles(const Halfedge_handle h)
{
    if(h->facet() == Facet_handle() || h->opposite()->facet() == Facet_handle()) return false;

    Vector N1 = triangleNormal( h->prev()->vertex()->point(),
                                h->vertex()->point(),
                                h->next()->vertex()->point());

    Vector N2 = triangleNormal( h->prev()->vertex()->point(),
                                h->vertex()->point(),
                                h->next()->vertex()->point());

    return (N1*N2>0.99);
}


static inline double triangleArea(const Point3d& a, const Point3d& b, const Point3d& c)
{
    Vector N = triangleNormalUnnormalized(a, b, c);
    return std::sqrt(N*N)*0.5;
};

static inline double triangleArea(const Facet_handle pFacet)
{
    return triangleArea(pFacet->halfedge()->vertex()->point(),
                        pFacet->halfedge()->prev()->vertex()->point(),
                        pFacet->halfedge()->next()->vertex()->point());
};

// returns d in the plane equation ax+by+cz+d=0 with (a,b,c) the normalized normal
static inline double plane_distance_from_origin(const Vector& N, const Point3d& planePos)
{
    return -(N.x()*planePos.x()+N.y()*planePos.y()+N.z()*planePos.z());
}

static inline double point_to_plane_signed_distance(const Vector& N, const Point3d& planePos, const Point3d& Pos)
{
    return plane_distance_from_origin(N, planePos) + N.x()*Pos.x()+N.y()*Pos.y()+N.z()*Pos.z();
}

// returns d in the plane equation ax+by+cz+d=0 with (a,b,c) the normalized normal
static inline double plane_distance_from_origin(const Point3d& p1, const Point3d& p2, const Point3d& p3)
{
    Vector N = triangleNormal(p1, p2, p3);

    return plane_distance_from_origin(N, p1);
}

static inline bool plane_plane_intersection(const Vector& N1,
                                            const double d1,
                                            const Vector& N2,
                                            const double d2,
                                            Vector& line_dir,
                                            Point3d& P0)
{
    line_dir = CGAL::cross_product(N1, N2);  // it can also be interpreted as the normal of an
                                                    // orthonormal third plane
    // parallel planes case
    if(std::sqrt(line_dir*line_dir) < 1e-8) return false;

    // 2 non-parallel planes intersect!!
    // one must first select a nonzero coordinate of line_dir, and then set that coordinate of P0 = 0.
    /*if(abs(line_dir.z())>1e-8)
    {
        P0.z()=0;
        P0.x()= (N1.y()*d2 - N2.y()*d1) / ();
    }
    else if(abs(line_dir.y())>1e-8)
    {

    }
    if(abs(line_dir.x())>1e-8)
    {

    }
    */

    cout << "plane_plane_intersection: not implemented!!" << endl;

    assert(false);

    return true;
}

static inline double tetrahedronVolume( const Point3d& a,
                                        const Point3d& b,
                                        const Point3d& c,
                                        const Point3d& d)
{
    return abs((a-d)*CGAL::cross_product(b-d, c-d))/6.0;
};

static inline double tetrahedronVolume(const Halfedge_handle h)
{
    // we can compute the tetrahedron volume only if the edge is not a border edge!!
    if(h->facet() == Facet_handle() || h->opposite()->facet() == Facet_handle())
        return PROHIB_ENERGY;

    return tetrahedronVolume(   h->next()->vertex()->point(),
                                h->vertex()->point(),
                                h->prev()->vertex()->point(),
                                h->opposite()->next()->vertex()->point()
                            );
};

static inline double oneRingArea(const Vertex_handle pVertex)
{
    double sum=0.0;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        sum+=triangleArea(  h->prev()->vertex()->point(),
                            h->vertex()->point(),
                            h->next()->vertex()->point());
    }

    return sum;
}

static inline Point3d triangleBarycenter(const Point3d& a, const Point3d& b, const Point3d& c)
{
    double tmp = 1.0/3.0;
    return Point3d((a.x()+b.x()+c.x())*tmp, (a.y()+b.y()+c.y())*tmp, (a.z()+b.z()+c.z())*tmp);
}

static inline long double triangleShapePotential(   const Point3d& a,
                                                    const Point3d& b, // the central point for which we give the gradient
                                                    const Point3d& c)
{
    long double D[3], circumradius, potential = PROHIB_ENERGY, den;

	// Length of the triangle sides
	D[0] = sqrt((a-b).squared_length());
	D[1] = sqrt((a-c).squared_length());
	D[2] = sqrt((b-c).squared_length());

	den = std::min<double>(D[2],std::min<double>(D[0],D[1]));
    if(den>1e-8 && !CGAL::collinear(a, b, c))
    {
        circumradius = D[0]*D[1]*D[2]/(4.0*triangleArea(a,b,c));
        potential = circumradius/den;
    }

    if(potential > PROHIB_ENERGY) potential = PROHIB_ENERGY; // to truncate the shape error to a manageable max value

    return potential;
};

// return the sum of the two adjacent triangle potentials (if the edge is not a border edge)
static inline long double triangleShapePotential(const Halfedge_handle h)
{
    double potential;
    if(h->opposite()->facet()==Facet_handle())
    {
        potential = triangleShapePotential(h->prev()->vertex()->point(),
                                           h->vertex()->point(), // the central point for which we give the gradient
                                           h->next()->vertex()->point());
    }
    else if(h->facet()==Facet_handle())
    {
        potential = triangleShapePotential(h->prev()->vertex()->point(),
                                           h->opposite()->next()->vertex()->point(), // the central point for which we give the gradient
                                           h->vertex()->point());
    }
    else
    {
        potential = triangleShapePotential(h->prev()->vertex()->point(),
                                           h->vertex()->point(), // the central point for which we give the gradient
                                           h->next()->vertex()->point()) +
                    triangleShapePotential(h->prev()->vertex()->point(),
                                           h->opposite()->next()->vertex()->point(), // the central point for which we give the gradient
                                           h->vertex()->point());
    }

    return potential;
}

static inline long double triangleShapePotential(const vector<Halfedge_handle> & triangles)
{
   double sumPot = 0.0;
   unsigned int nbe = triangles.size();
   for(unsigned int i=0; i<nbe; ++i)
   {
		sumPot += triangleShapePotential(   triangles[i]->prev()->vertex()->point(),
                                            triangles[i]->vertex()->point(),
                                            triangles[i]->next()->vertex()->point() );
   }

   return sumPot;
}

static inline long double triangleShapePotentialAfterFlip(  PolyhedronPtr pMesh,
                                                            Halfedge_handle h)
{
    double potential;
    if(h->opposite()->facet()==Facet_handle() || h->facet()==Facet_handle())
    { // border edge case
        potential = PROHIB_ENERGY;
    }
    else
    {
        pMesh->flip_edge(h);
        potential = triangleShapePotential(h);
        pMesh->flip_edge(h);
    }

    return potential;
}


static inline long double triangleShapePotentialAfterSplit(const Halfedge_handle h)
{
	Point3d NewPos = h->vertex()->point()+(h->prev()->vertex()->point() - h->vertex()->point())*0.5f;
    if( h->opposite()->facet()==Facet_handle() ) // border edge
	{
        return  triangleShapePotential(NewPos, h->vertex()->point(), h->next()->vertex()->point())+
                triangleShapePotential(h->next()->vertex()->point(), h->prev()->vertex()->point(), NewPos);
	}
    else if( h->facet()==Facet_handle() ) // border edge
	{
        return  triangleShapePotential(NewPos, h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point()),
                triangleShapePotential(h->opposite()->next()->vertex()->point(), h->vertex()->point(), NewPos);
	}
	else
	{ // normal edge
        return triangleShapePotential(NewPos, h->vertex()->point(), h->next()->vertex()->point())+
               triangleShapePotential(h->next()->vertex()->point(), h->prev()->vertex()->point(), NewPos)+
               triangleShapePotential(h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point(), NewPos)+
               triangleShapePotential(NewPos, h->opposite()->next()->vertex()->point(), h->vertex()->point());
	}
}

// return the difference of shape potential after and before an edge flip
static inline long double deltaTriangleShapePotentialFlip(  PolyhedronPtr pMesh,
                                                            Halfedge_handle h)
{
    double init_pot = triangleShapePotential(h);
    if(init_pot>=PROHIB_ENERGY) return PROHIB_ENERGY; // bad init configuration

    return triangleShapePotentialAfterFlip(pMesh, h)-init_pot;
}

// return the difference of shape potential after and before an edge split
static inline long double deltaTriangleShapePotentialSplit(const Halfedge_handle h)
{
    double init_pot = triangleShapePotential(h);
    if(init_pot>=PROHIB_ENERGY) return PROHIB_ENERGY; // bad init configuration

    return triangleShapePotentialAfterSplit(h)-triangleShapePotential(h);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// triangle quality measures => voir artilce de SHEWCHUK "What is a good linear element?"

// q > 0.6 the triangle is of acceptable quality. q = 1 when h1 = h2 = h3.
static inline double triangle_PDE_quality_measure(  const Point3d& a,
                                                    const Point3d& b,
                                                    const Point3d& c)
{
    Vector ab(b-a), bc(c-b), ca(a-c);
    return (4.0*triangleArea(a,b,c)*std::sqrt(3.0))/(ab*ab+bc*bc+ca*ca);
}


// Interpolation quality measure (greater the result better the triangle shape and size)
// which has scale-invariance property.
static inline double triangleSmoothSizeAndShape(const Point3d& a, const Point3d& b, const Point3d& c)
{
    double D[3], Area = triangleArea(a,b,c);
	// Squared Length of the triangle sides
	D[0] = (a-b).squared_length();
	D[1] = (a-c).squared_length();
	D[2] = (b-c).squared_length();

    return Area/pow(D[0]*D[1]*D[2],1.0/3.0);
};

static inline double triangleSmoothSizeAndShape(Halfedge_handle h)
{
    if(h->is_border() || h->opposite()->facet()==Facet_handle())
        return triangleSmoothSizeAndShape(h->prev()->vertex()->point(),h->vertex()->point(), h->next()->vertex()->point());

    return std::min<double>(triangleSmoothSizeAndShape(h->prev()->vertex()->point(),h->vertex()->point(), h->next()->vertex()->point()),
                    triangleSmoothSizeAndShape(h->prev()->vertex()->point(),h->opposite()->next()->vertex()->point(), h->vertex()->point()));
};

static inline double triangleAspectRatio(const Point3d& a, const Point3d& b, const Point3d& c)
{
    double perimeter = trianglePerimeter(a,b,c);
    if(perimeter>0.0)
        return triangleArea(a,b,c)/(perimeter*perimeter);
    else
        return PROHIB_ENERGY;
};

static inline double triangleAspectRatio(Halfedge_handle h)
{
    return std::min<double>(triangleAspectRatio(h->prev()->vertex()->point(), h->vertex()->point(), h->next()->vertex()->point()),
                    triangleAspectRatio(h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point(), h->vertex()->point()));
};
// return the mean of the three edges square length
static inline double triangleMeanSquare(const Point3d& a, const Point3d& b, const Point3d& c)
{
    Vector ab(b-a), bc(c-b), ca(a-c);
    return 1.0/3.0*(ab*ab+bc*bc+ca*ca);
};

static inline double triangleRootMeanSquare(const Point3d& a, const Point3d& b, const Point3d& c)
{
    return std::sqrt(triangleMeanSquare(a, b, c));
};

// Conditioning quality measure (greater the result better the triangle conditioning)
// be aware that the creation of smaller elements can ocasionally decrease the min eigenvalue
// and worsen the conditioning of the stiffness matrix
// To refine without worsening the conditioning, one option is to compare the current element with the
// worst new element that will appear if the refinement occurs: the worst new local element must be better
// than the original element by a margin large enough to justify the smaller elements
static inline double triangleScaleInvariantConditioning(const Point3d& a, const Point3d& b, const Point3d& c)
{
    double Area = triangleArea(a,b,c), l=3.0*triangleMeanSquare(a, b, c);

    return Area/(l+std::sqrt(l*l-48.0*Area*Area));
};

static inline double triangleSmoothScaleInvariantConditioning(const Point3d& a, const Point3d& b, const Point3d& c)
{
    double Area = triangleArea(a,b,c), l=triangleMeanSquare(a, b, c);

    return Area/l;
};

// harmonic mean is better for combining quality values than arithmetic mean
// greater the result better the triangle quality
static inline double triangleQuality(const Point3d& a, const Point3d& b, const Point3d& c)
{
    //return triangleAspectRatio(a,b,c);
    return 1.0/triangleShapePotential(a,b,c);
    //return triangle_PDE_quality_measure(a,b,c);
    //return triangleSmoothSizeAndShape(a,b,c); // decreases large angles
    //return triangleSmoothScaleInvariantConditioning(a,b,c); // increases small angles
    //return harmonicMean(triangleSmoothSizeAndShape(a,b,c), triangleScaleInvariantConditioning(a,b,c), 0.5);
};

// greater the result better the triangle quality
static inline double triangleQuality(const Halfedge_handle h)
{
    if( h->opposite()->facet()==Facet_handle() )
    {
        return triangleQuality(h->prev()->vertex()->point(),h->vertex()->point(), h->next()->vertex()->point());
    }
    else if( h->facet()==Facet_handle() )
    {
        return triangleQuality(h->prev()->vertex()->point(),h->opposite()->next()->vertex()->point(), h->vertex()->point());
    }
    else
    {
        return std::min<double>(triangleQuality(h->prev()->vertex()->point(),h->vertex()->point(), h->next()->vertex()->point()),
                        triangleQuality(h->prev()->vertex()->point(),h->opposite()->next()->vertex()->point(), h->vertex()->point()) );
    }
};
// greater the result better the triangle quality
static inline double triangleQualityAfterFlip(const Halfedge_handle h)

{
    if( h->opposite()->facet()==Facet_handle() || h->opposite()->facet()==Facet_handle())
    {
        return -PROHIB_ENERGY; // it is forbidden to flip such an edge
    }
    else
    {
        return std::min<double>(triangleQuality(h->vertex()->point(),h->next()->vertex()->point(),h->opposite()->next()->vertex()->point()),
                        triangleQuality(h->next()->vertex()->point(),h->prev()->vertex()->point(),h->opposite()->next()->vertex()->point()));
    }
};

static inline double triangleQuality(const std::vector< Halfedge_handle >& edges)
{
    double minQuality = PROHIB_ENERGY;
    std::vector< Halfedge_handle >::const_iterator it(edges.begin()), ite(edges.end());
    for(;it!=ite;++it)
    {
        double tq = triangleQuality(*it);
        if(minQuality > tq) minQuality = tq;
    }

    return minQuality;
}



static inline double triangleQualityAfterSplit(const Halfedge_handle h)

{
	Point3d NewPos = h->vertex()->point()+(h->prev()->vertex()->point() - h->vertex()->point())*0.5f;
    if( h->opposite()->facet()==Facet_handle() )
	{
        return std::min<double>(triangleQuality(NewPos, h->vertex()->point(), h->next()->vertex()->point()),
                        triangleQuality(h->next()->vertex()->point(), h->prev()->vertex()->point(), NewPos));
	}
	else if( h->facet()==Facet_handle() )
	{
        return std::min<double>(triangleQuality(NewPos, h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point()),
                        triangleQuality(h->opposite()->next()->vertex()->point(), h->vertex()->point(), NewPos));
	}
	else
	{
        return std::min<double>(    std::min<double>(triangleQuality(NewPos, h->vertex()->point(), h->next()->vertex()->point()),
                                     triangleQuality(h->next()->vertex()->point(), h->prev()->vertex()->point(), NewPos)),
                            std::min<double>(triangleQuality(h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point(), NewPos),
                                     triangleQuality(NewPos, h->opposite()->next()->vertex()->point(), h->vertex()->point())));
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// is the projection of point Pos on the triangle plane in the triangle (a,b,c)?
static inline bool isInTriangle(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& Pos)
{ // Pos =
	 // we are going to compute the barycentric coordinates of Pos   (from numerical recipes)
	 Vector p1p = a-c, p2p = b-c, Posp = Pos-c;
	 double p1p2 = p1p*p1p, p2p2 = p2p*p2p, p1pp2p = p1p*p2p, p1pPosp = p1p*Posp, p2pPosp = p2p*Posp;
	 /////////////////////////////////////////////////////////////////////////////////////////////////
	 double den = p1p2*p2p2-p1pp2p*p1pp2p;
	 double alpha = (p2p2*p1pPosp-p1pp2p*p2pPosp)/den, beta = (p1p2*p2pPosp-p1pp2p*p1pPosp)/den, gamma;
	 gamma = 1-(alpha+beta);
	 /////////////////////////////////////////////////////////////////////////////////////////////////
	 return !((alpha<0.0)||(beta<0.0)||(gamma<0.0));
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline double bilinear_interpolation_triangle(   const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& Pos,
                                                        double vala, double valb, double valc)
{
 	 Vector p1p = a-c, p2p = b-c, Posp = Pos-c;
	 double p1p2 = p1p*p1p, p2p2 = p2p*p2p, p1pp2p = p1p*p2p, p1pPosp = p1p*Posp, p2pPosp = p2p*Posp;
	 /////////////////////////////////////////////////////////////////////////////////////////////////
	 double den = p1p2*p2p2-p1pp2p*p1pp2p;
	 double alpha = (p2p2*p1pPosp-p1pp2p*p2pPosp)/den, beta = (p1p2*p2pPosp-p1pp2p*p1pPosp)/den, gamma;
	 gamma = 1-(alpha+beta);
     /////////////////////////////////////////////////////////////////////////////////////////////////
     if(alpha<0) alpha=0.0;
     if(alpha>1) alpha=1.0;

     if(beta<0) beta=0.0;
     if(beta>1) beta=1.0;

     if(gamma<0) gamma=0.0;
     if(gamma>1) gamma=1.0;
     /////////////////////////////////////////////////////////////////////////////////////////////////
     return alpha*vala + beta*valb + gamma*valc;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// triangle angle methods
static inline double computeTriangleAngle(   const Point3d& a,
                                             const Point3d& b, // the central point on which we compute the angle
                                             const Point3d& c)
{
    Vector b_a(a-b), b_c(c-b);
    double lba = std::sqrt(b_a*b_a), lbc = std::sqrt(b_c*b_c);
    if(lba>1e-8)
        b_a = b_a/lba;
    else
        b_a = CGAL::NULL_VECTOR;
    if(lbc>1e-8)
        b_c = b_c/lbc;
    else
        b_c = CGAL::NULL_VECTOR;

    double cosabc = b_a*b_c;

    return acos(cosabc);
}

static inline bool isTriangleAngleGreaterThanPi( const Point3d& a,
                                                 const Point3d& b, // the central point on which we compute the angle
                                                 const Point3d& c)
{
    return (computeTriangleAngle(a, b, c) > M_PI);
}

static inline double computeMinTriangleAngles(const Point3d& a, const Point3d& b, const Point3d& c)
{
    return std::min<double>(computeTriangleAngle(a,b,c), std::min<double>(computeTriangleAngle(b,c,a), computeTriangleAngle(c,a,b)));
}

static inline double computeMaxTriangleAngles(const Point3d& a, const Point3d& b, const Point3d& c)
{
    return std::max<double>(computeTriangleAngle(a,b,c), std::max<double>(computeTriangleAngle(b,c,a), computeTriangleAngle(c,a,b)));
}
// compute the min angle in the two adjacent triangles
static inline double computeMinTriangleAngles(const Halfedge_handle h)
{
	 return std::min<double>(   computeMinTriangleAngles(h->opposite()->next()->vertex()->point(), h->vertex()->point(), h->prev()->vertex()->point()),
                                computeMinTriangleAngles(h->prev()->vertex()->point(), h->vertex()->point(), h->next()->vertex()->point()));
}
// compute the max angle in the two adjacent triangles
static inline double computeMaxTriangleAngles(const Halfedge_handle h)
{
	 return std::max<double>(   computeMaxTriangleAngles(h->opposite()->next()->vertex()->point(), h->vertex()->point(), h->prev()->vertex()->point()),
                                computeMaxTriangleAngles(h->prev()->vertex()->point(), h->vertex()->point(), h->next()->vertex()->point()));
}

// compute the min angle for all the triangles in the one-ring neighborhood
static double computeMinTriangleAngles(const Vertex_handle pVertex)
{
    double minA = M_PI;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
		  minA = std::min<double>(minA,	computeMinTriangleAngles( h->prev()->vertex()->point(),
                                                                  h->vertex()->point(),
                                                                  h->next()->vertex()->point())   );
	}

	return minA;
}
// compute the max angle for all the triangles in the one-ring neighborhood
static double computeMaxTriangleAngles(const Vertex_handle pVertex)
{
    double maxA = 0.0;
    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
		  maxA = std::max<double>(maxA,	computeMaxTriangleAngles( h->prev()->vertex()->point(),
                                                                  h->vertex()->point(),
                                                                  h->next()->vertex()->point())   );
	}

	return maxA;
}
// compute the min angle for all the triangles in the halfedge list (we use the 3 points obtained
// by the halfedge and the next-halfedge )
static double computeMinTriangleAngles(const vector<Halfedge_handle> & triangles)
{
   double minA = M_PI;
   unsigned int nbe = triangles.size();
   for(unsigned int i=0; i<nbe; ++i)
   {
		minA = std::min<double>(minA, computeMinTriangleAngles(   triangles[i]->prev()->vertex()->point(),
                                                                  triangles[i]->vertex()->point(),
                                                                  triangles[i]->next()->vertex()->point())   );

   }

   return minA;
}
// compute the max angle for all the triangles in the halfedge list (we use the 3 points obtained
// by the halfedge and the next-halfedge )
static double computeMaxTriangleAngles(const vector<Halfedge_handle> & triangles)
{
   double maxA = 0.0;
   unsigned int nbe = triangles.size();
   for(unsigned int i=0; i<nbe; ++i)
   {
		maxA = std::max<double>(maxA, computeMaxTriangleAngles(   triangles[i]->prev()->vertex()->point(),
                                                                  triangles[i]->vertex()->point(),
                                                                  triangles[i]->next()->vertex()->point())   );

   }

   return maxA;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// triangle angle methods used for connectivity modification evaluation
static inline double computeMinTriangleAnglesAfterFlip(const Halfedge_handle h)
{
	 return std::min<double>(   computeMinTriangleAngles(h->opposite()->next()->vertex()->point(), h->next()->vertex()->point(), h->prev()->vertex()->point()),
                                computeMinTriangleAngles(h->opposite()->next()->vertex()->point(), h->vertex()->point(), h->next()->vertex()->point()));
}

static inline double computeMaxTriangleAnglesAfterFlip(const Halfedge_handle h)
{
	 return std::max<double>(   computeMaxTriangleAngles(h->opposite()->next()->vertex()->point(), h->next()->vertex()->point(), h->prev()->vertex()->point()),
                                computeMaxTriangleAngles(h->opposite()->next()->vertex()->point(), h->vertex()->point(), h->next()->vertex()->point()));
}
// compute the minimum angle obtained for the 4 triangles resulting from an edge split
static inline double computeMinTriangleAnglesAfterSplit(const Halfedge_handle h)
{
	Point3d NewPos = h->vertex()->point()+(h->prev()->vertex()->point() - h->vertex()->point())*0.5f;
    if( h->opposite()->facet()==Facet_handle() ) // border edge
	{
        return std::min<double>(computeMinTriangleAngles(NewPos, h->vertex()->point(), h->next()->vertex()->point()),
                                computeMinTriangleAngles(h->next()->vertex()->point(), h->prev()->vertex()->point(), NewPos));
	}
    else if( h->facet()==Facet_handle() ) // border edge
	{
        return std::min<double>(computeMinTriangleAngles(NewPos, h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point()),
                                computeMinTriangleAngles(h->opposite()->next()->vertex()->point(), h->vertex()->point(), NewPos));
	}
	else
	{ // normal edge
        return std::min<double>(computeMinTriangleAngles(NewPos, h->vertex()->point(), h->next()->vertex()->point()),
                                std::min<double>(	computeMinTriangleAngles(h->next()->vertex()->point(), h->prev()->vertex()->point(), NewPos),
                                                    std::min<double>(	computeMinTriangleAngles(h->prev()->vertex()->point(), h->opposite()->next()->vertex()->point(), NewPos),
                                                                        computeMinTriangleAngles(NewPos, h->opposite()->next()->vertex()->point(), h->vertex()->point()) ) ) );
	}
}

// the greater, the better
static inline double flipQuality(PolyhedronPtr pMesh, Halfedge_handle h)// to use for triangle mesh
{
    double minAfter = computeMinTriangleAnglesAfterFlip(h);
    double angleimprovingMin = minAfter - computeMinTriangleAngles(h);

    // to do a flip that is not so bad
    if (angleimprovingMin<0.0 && minAfter>0.35) angleimprovingMin = 0.01*(minAfter-0.35);

    if (0.0<angleimprovingMin)
    {
        double inita = triangleAspectRatio(h);
        h = pMesh->flip_edge (h); // we flip the initial position since we try to maximize the minimum angle
        double aftera = triangleAspectRatio(h);

        h = pMesh->flip_edge (h); // we go back to the initial configuration

        // We do not want to decrease the triangle aspect ratio (or to create a feature edge => but it is forbidden to use a threshold here!!)
        if ( inita>aftera )
        {
            return 0.0f;
        }

        return angleimprovingMin;
    }

    return 0.0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// methods used for mixed area computation (see "Discrete differential-geometry operators
// for triangulated 2-manifolds")

static inline bool isObtuseAngle(const Point3d& a, const Point3d& b, const Point3d& c)
{
    return (computeTriangleAngle(a, b, c)>M_PI/3.0);
}

static inline bool isObtuseTriangle(const Point3d& a, const Point3d& b, const Point3d& c)
{
    return (computeMaxTriangleAngles(a, b, c)>M_PI/3.0);
}

static inline double voronoiAreaTriangle(const Point3d& a, const Point3d& b, const Point3d& c)
{ // for non-obtuse triangle!!!
    assert(!isObtuseTriangle(a, b, c));

    Vector ac(c-a), ab(b-a);

    return 0.125*(ac*ac*cotangent(a,b,c) + ab*ab*cotangent(b,c,a));
}

#endif // Triangle_shape_analysis_H
