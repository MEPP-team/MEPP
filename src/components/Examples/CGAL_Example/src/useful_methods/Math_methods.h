///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Math_methods_H
#define Math_methods_H

#include "../../../../../mepp/Polyhedron/polyhedron.h"

#include "../triangle_manager/CGAL_includes_and_types.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
normalize_Vector efficiently normalize a Vector and return its L2 norm
**/
static inline double normalize_Vector(Vector& V)
{
    double nV2 = V*V;

    if (nV2>1e-16)
    {
        nV2 = std::sqrt(nV2);
        V = V*(1.0/nV2);
    }
    else
    {
        //cout << "normalize_Vector: WARNING NULL VECTOR" << endl;
        nV2=0.0;
    }

    return nV2;
}

static inline void clamp(double& x, double minx, double maxx)
{
    if (x<minx) x=minx;
    else if (x>maxx) x=maxx;
};

static inline void clamp(unsigned int& x, unsigned int minx, unsigned int maxx)
{
    if (x<minx) x=minx;
    else if (x>maxx) x=maxx;
};

static inline void clamp(double& x, double& y, double minx, double maxx, double miny, double maxy)
{
    if (x<minx) x=minx;
    else if (x>maxx) x=maxx;

    if (y<miny) y=miny;
    else if (y>maxy) y=maxy;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static inline double harmonicMean(double Q1, double Q2, double wQ1)
{
    if (wQ1<0.0) wQ1=0.0;
    else if (wQ1>1.0) wQ1=1.0;

    double Qinv = wQ1/std::max<double>(Q1,1e-8)+(1.0-wQ1)/std::max<double>(Q2,1e-8);

    return 1.0/Qinv;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static inline double cotangent(const Point3d& a, const Point3d& b, const Point3d& c)
{
    Vector ba = a-b, bc = c-b, cross;
    cross = CGAL::cross_product(bc, ba);
    if (cross!=CGAL::NULL_VECTOR)
        return (ba*bc)/std::sqrt(cross*cross);
    else
        return 1e16;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Noise generation

// Get two unit normal randoms
static inline void gaussianNoise(double& N1, double& N2, bool initsrand = true)
{
    double V1, V2, S;
    if (initsrand) srand((unsigned)time(NULL));
    // we use the polar method to turn two uniformly distributed randoms, U1 and U2, into two unit normal randoms, X and Y
    do
    {
        //U1=rand()/double(RAND_MAX);	// U1=[0,1]
        //U2=rand()/double(RAND_MAX);	// U2=[0,1]
        //V1=2.0 * U1 -1.0;			// V1=[-1,1]
        //V2=2.0 * U2 - 1.0;			// V2=[-1,1]
        V1=2.0 * rand()/double(RAND_MAX) -1.0;
        V2=2.0 * rand()/double(RAND_MAX) - 1.0;
        S=V1 * V1 + V2 * V2;
    }
    while (S >=1.0 || S==0.0);

    N1=std::sqrt(-2.0 * log(S) / S) * V1;
    N2=std::sqrt(-2.0 * log(S) / S) * V2;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

// the four first variables are input variables and the two last ones are output variables (valid if the function return true)
static inline bool two_2Dcercles_intersection(	const Point3d& c1,
												const Point3d& c2,
												const double rad_1,
												const double rad_2,
												Point3d& inter1,
												Point3d& inter2)
{
    double a = 2.0*(c2.x()-c1.x()), b = 2.0*(c2.y()-c1.y()), c = (c2.x()-c1.x())*(c2.x()-c1.x())+(c2.y()-c1.y())*(c2.y()-c1.y())-rad_2*rad_2+rad_1*rad_1;
    double aa = a*a, bb = b*b, cc = c*c, _2ac = 2.0*a*c;
    double delta = _2ac*_2ac - 4.0*(aa+bb)*(cc-bb*rad_1*rad_1);
    if (delta>=0)
    {
        double inter1x, inter2x, inter1y, inter2y;
        if (b!=0.0)
        {
            // delta>0
            double dd = 2.0*(aa+bb);
            inter1x = c1.x()+(_2ac-std::sqrt(delta))/dd;
            inter2x = c1.x()+(_2ac+std::sqrt(delta))/dd;
            inter1y = c1.y() + (c-a*(inter1x-c1.x()))/b;
            inter2y = c1.y() + (c-a*(inter2x-c1.x()))/b;
        }
        else
        {
            // delta==0
            inter1x = inter2x = c1.x()+(a*c)/(aa);
            double dd=(2*c-a*a)/(2*a);
            double rr = std::sqrt(rad_2*rad_2-(dd*dd));
            inter1y = c1.y() - rr;
            inter2y = c1.y() + rr;
        }
        inter1 = Point3d(inter1x, inter1y, 0);
        inter2 = Point3d(inter2x, inter2y, 0);
        return true;
    }
    else
        return false;
};


static inline bool two_2Dcercles_intersection(const Point3d& c1,
											  const Point3d& c2,
											  const double rad_1,
											  const double rad_2,
											  Point3d& c3,
											  double& rad_3)
{
    Point3d inter1, inter2;
    if (two_2Dcercles_intersection(c1, c2, rad_1, rad_2, inter1, inter2))
    {
        c3 = inter1+0.5*(inter2-inter1); // midpoint of intersections
        rad_3 = std::sqrt((inter1-c3).squared_length()); // the radius is zero only if inter2 == inter1
        return true;
    }
    else return false;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline double ACCR(double X1, double X2, double Y1, double Y2){ return (Y2-Y1)/(X2-X1); }
static inline double ACCR_P2D(const Point3d& P1, const Point3d& P2){ return ACCR(P1.x(),P2.x(),P1.y(),P2.y()); }
static inline Point3d MIN2(const Point3d& P1, const Point3d& P2)
{
    return Point3d(std::min<double>(P1.x(),P2.x()),
				   std::min<double>(P1.y(),P2.y()),
				   std::min<double>(P1.z(),P2.z()) );
}
static inline Point3d MIN3(const Point3d& P1, const Point3d& P2, const Point3d& P3)
{
    return MIN2(P1,MIN2(P2,P3));
}
static inline Point3d MAX2(const Point3d& P1, const Point3d& P2)
{
    return Point3d(std::max<double>(P1.x(),P2.x()),std::max<double>(P1.y(),P2.y()),std::max<double>(P1.z(),P2.z()));
}
static inline Point3d MAX3(const Point3d& P1, const Point3d& P2, const Point3d& P3)
{
    return MAX2(P1,MAX2(P2,P3));
}

static inline Point MINnef(const Point& P1, const Point& P2)
{
    return Point(min(P1.x(),P2.x()), min(P1.y(),P2.y()));
}
static inline Point MIN3nef(const Point& P1, const Point& P2, const Point& P3)
{
    return MINnef(P1,MINnef(P2,P3));
}
static inline Point MAXnef(const Point& P1, const Point& P2)
{
    return Point(max(P1.x(),P2.x()),max(P1.y(),P2.y()));
}
static inline Point MAX3nef(const Point& P1, const Point& P2, const Point& P3)
{
    return MAXnef(P1,MAXnef(P2,P3));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// other
static inline Point3d midpoint(const Point3d& P1, const Point3d& P2)
{
    return P1+0.5*(P2-P1);
}

static inline Point3d barycenter(const Point3d& P1, const Point3d& P2, const Point3d& P3)
{
    return Point3d((P1.x()+P2.x()+P3.x())/3.0, (P1.y()+P2.y()+P3.y())/3.0,(P1.z()+P2.z()+P3.z())/3.0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif // Math_methods_H
