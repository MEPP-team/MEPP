///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Triangle_H
#define Triangle_H

#include "../triangle_manager/CGAL_includes_and_types.h"

#include "../triangle_tools/Triangle_shape_analysis.h"

#ifdef USE_SEGMENT
    #include "Unique_segment.h"
#endif

/*
    STATICS METHODS USEFULL FOR TRIANGLES
*/
static inline bool isValidTriangle(const Point_dt& p1, const Point_dt& p2, const Point_dt& p3)
{
    if( CGAL::collinear(p1, p2, p3) || Triangle_dt(p1, p2, p3).is_degenerate() ) return false;

    return true;
}

static inline bool isValidTriangle(const Point3d& p1, const Point3d& p2, const Point3d& p3)
{
    Point_dt    p1_dt(p1.x(),p1.y()),
                p2_dt(p2.x(),p2.y()),
                p3_dt(p3.x(),p3.y());

    if( CGAL::collinear(p1, p2, p3) || !isValidTriangle(p1_dt, p2_dt, p3_dt) ) return false;
//(triangleShapePotential(p3d1,p3d2,p3d3)<TH
    return true;
}

/*
    END OF STATICS METHODS USEFULL FOR TRIANGLES
*/



// Now we start our algorithm with the set of all possible triangles taking care that

// for the time being the comparison between two triangles is inefficient because
// I do not use a key (for exemple I could use ijk where i, j an k are the 3 points indices
// composing the triangle)
class triangle
{
    public:
    typedef CGAL::Vector_3< CGAL::Cartesian<double> > vector;

    Point m_p1, m_p2, m_p3; // 2D nef points
    double m_potential; // quality measure of the triangle
    // the following variables will be used to speed up the intersection computation (approximation for easier computation):
    Point3d m_circumc; // circumcenter
    double m_squared_radius;
    Point3d m_min_AABB, m_max_AABB; // approximated axis-aligned bounding box

    public:
    triangle():m_p1(Point(0,0)),m_p2(Point(1,0)),m_p3(Point(1,1))
    {
        m_circumc = circumcenter(Point3d(0,0,0),Point3d(1,0,0),Point3d(1,1,0));
        m_squared_radius = (m_circumc-Point3d(0,0,0)).squared_length();
        m_potential = 0;
    };

    triangle(const Point& p1, const Point& p2, const Point& p3, double potential, const Point3d& circumc, double squared_radius, const Point3d& min_AABB, const Point3d& max_AABB):m_p1(p1),m_p2(p2),m_p3(p3),m_potential(potential),m_circumc(circumc),m_squared_radius(squared_radius),m_min_AABB(min_AABB), m_max_AABB(max_AABB)
    {
    };

    triangle(const Point3d& p3d1, const Point3d& p3d2, const Point3d& p3d3) // build directly a counterclockwise oriented triangle!!
    {
        //cout << "triangle: IN..." << endl;
        assert(isValidTriangle(p3d1,p3d2,p3d3));

        if(p3d1.x() <= p3d2.x() && p3d1.x() <= p3d3.x())
        {
            m_p1 = Point(p3d1.x(),p3d1.y());
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // for counter-clockwise orientation, the second point must have the min growth rate
            if(p3d1.x()==p3d2.x())
            {
                if(p3d2.y() <= p3d1.y())
                {
                    m_p2 = Point(p3d2.x(),p3d2.y()); m_p3 = Point(p3d3.x(),p3d3.y());
                }
                else
                {
                    m_p2 = Point(p3d3.x(),p3d3.y()); m_p3 = Point(p3d2.x(),p3d2.y());
                }
            }
            else if(p3d1.x()==p3d3.x())
            {
                if(p3d3.y() <= p3d1.y())
                {
                    m_p2 = Point(p3d3.x(),p3d3.y()); m_p3 = Point(p3d2.x(),p3d2.y());
                }
                else
                {
                    m_p2 = Point(p3d2.x(),p3d2.y()); m_p3 = Point(p3d3.x(),p3d3.y());
                }
            }
            else
            { // here we can compute the growth rate to choose the right order
                if( ACCR_P2D(p3d1,p3d2) <= ACCR_P2D(p3d1,p3d3) )
                {
                    m_p2 = Point(p3d2.x(),p3d2.y()); m_p3 = Point(p3d3.x(),p3d3.y());
                }
                else
                {
                    m_p2 = Point(p3d3.x(),p3d3.y()); m_p3 = Point(p3d2.x(),p3d2.y());
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        else if(p3d2.x()<= p3d1.x() && p3d2.x()<= p3d3.x())
        {
            m_p1 = Point(p3d2.x(),p3d2.y());
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // for counter-clockwise orientation, the second point must have the min growth rate
            if(p3d2.x()==p3d1.x())
            {
                if(p3d1.y() <= p3d2.y())
                {
                    m_p2 = Point(p3d1.x(),p3d1.y()); m_p3 = Point(p3d3.x(),p3d3.y());
                }
                else
                {
                    m_p2 = Point(p3d3.x(),p3d3.y()); m_p3 = Point(p3d1.x(),p3d1.y());
                }
            }
            else if(p3d2.x()==p3d3.x())
            {
                if(p3d3.y() <= p3d2.y())
                {
                    m_p2 = Point(p3d3.x(),p3d3.y()); m_p3 = Point(p3d1.x(),p3d1.y());
                }
                else
                {
                    m_p2 = Point(p3d1.x(),p3d1.y()); m_p3 = Point(p3d3.x(),p3d3.y());
                }
            }
            else
            { // here we can compute the growth rate to choose the right order
                if( ACCR_P2D(p3d2,p3d1) <= ACCR_P2D(p3d2,p3d3) )
                {
                    m_p2 = Point(p3d1.x(),p3d1.y()); m_p3 = Point(p3d3.x(),p3d3.y());
                }
                else
                {
                    m_p2 = Point(p3d3.x(),p3d3.y()); m_p3 = Point(p3d1.x(),p3d1.y());
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        else
        {
            m_p1 = Point(p3d3.x(),p3d3.y());
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // for counter-clockwise orientation, the second point must have the min growth rate
            if(p3d3.x()==p3d1.x())
            {
                if(p3d1.y() <= p3d3.y())
                {
                    m_p2 = Point(p3d1.x(),p3d1.y()); m_p3 = Point(p3d2.x(),p3d2.y());
                }
                else
                {
                    m_p2 = Point(p3d2.x(),p3d2.y()); m_p3 = Point(p3d1.x(),p3d1.y());
                }
            }
            else if(p3d3.x()==p3d2.x())
            {
                if(p3d2.y() <= p3d3.y())
                {
                    m_p2 = Point(p3d2.x(),p3d2.y()); m_p3 = Point(p3d1.x(),p3d1.y());
                }
                else
                {
                    m_p2 = Point(p3d1.x(),p3d1.y()); m_p3 = Point(p3d2.x(),p3d2.y());
                }
            }
            else
            { // here we can compute the growth rate to choose the right order
                if( ACCR_P2D(p3d3,p3d2) <= ACCR_P2D(p3d3,p3d1) )
                {
                    m_p2 = Point(p3d2.x(),p3d2.y()); m_p3 = Point(p3d1.x(),p3d1.y());
                }
                else
                {
                    m_p2 = Point(p3d1.x(),p3d1.y()); m_p3 = Point(p3d2.x(),p3d2.y());
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }


        //cout << "triangle: data storing..." << endl;
        m_potential = triangleShapePotential(p3d1, p3d2, p3d3);
        if(m_potential<PROHIB_ENERGY)
        {
            m_circumc = CGAL::circumcenter(p3d1, p3d2, p3d3);
            m_squared_radius = (p3d1-m_circumc).squared_length() + 1e-6;
        }
        else
        {
            m_squared_radius = flatTriangleBarycenter(m_circumc);
        }

        m_min_AABB = MIN3(p3d1, p3d2, p3d3);
        m_max_AABB = MAX3(p3d1, p3d2, p3d3);

        //print();

        //cout << "triangle: OUT..." << endl;
    }

     // using delaunay points:
    triangle(const Point_dt& pd1, const Point_dt& pd2, const Point_dt& pd3) // build directly a counterclockwise oriented triangle using delaunay points!!
    {
        assert(isValidTriangle(pd1,pd2,pd3));
        //Point3d p3d1(pd1.x().to_double(),pd1.y().to_double(),0), p3d2(pd2.x().to_double(),pd2.y().to_double(),0), p3d3(pd3.x().to_double(),pd3.y().to_double(),0);
        Point3d p3d1(CGAL::to_double(pd1.x()),CGAL::to_double(pd1.y()),0), p3d2(CGAL::to_double(pd2.x()),CGAL::to_double(pd2.y()),0), p3d3(CGAL::to_double(pd3.x()),CGAL::to_double(pd3.y()),0);
        //Point3d p3d1(pd1.x(),pd1.y(),0), p3d2(pd2.x(),pd2.y(),0), p3d3(pd3.x(),pd3.y(),0);
        *this = triangle(p3d1, p3d2, p3d3); // we call the constructor above to set all the values
    }


    ~triangle(){};
    /////////////////////////////////////////////////////////////////////////////////////
    inline void operator= (const triangle& t)
    {
        m_p1=t.m_p1; m_p2=t.m_p2; m_p3=t.m_p3;
        m_potential = t.m_potential;
        m_circumc = t.m_circumc;
        m_squared_radius = t.m_squared_radius;
        m_min_AABB = t.m_min_AABB;
        m_max_AABB = t.m_max_AABB;
    }


    inline bool operator== (const triangle& t) const
    {
        return (( (m_p1==t.m_p1) && ( ((m_p2==t.m_p2) && (m_p3==t.m_p3)) || ((m_p2==t.m_p3) && (m_p3==t.m_p2)) ) ) ||
                ( (m_p1==t.m_p2) && ( ((m_p2==t.m_p3) && (m_p3==t.m_p1)) || ((m_p2==t.m_p1) && (m_p3==t.m_p3)) ) ) ||
                ( (m_p1==t.m_p3) && ( ((m_p2==t.m_p1) && (m_p3==t.m_p2)) || ((m_p2==t.m_p2) && (m_p3==t.m_p1)) ) )
                );
    }

    inline bool operator!= (const triangle& t) const
    {
        return !(*this==t);
    }

    inline bool operator< (const triangle& t) const
    { // strict weak ordering: if a is less than b then b is not less than a, if a is less than b and b is less than c then a is less than c
        //return (*this!=t) && (m_p1<t.m_p1 || (m_p1 <= t.m_p1 && m_p2<t.m_p2) || (m_p1 <= t.m_p1 && m_p2 <= t.m_p2 && m_p3<t.m_p3));
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /*Triangle_dt T1(Point_dt(m_p1.x().to_double(),m_p1.y().to_double()),Point_dt(m_p2.x().to_double(),m_p2.y().to_double()),Point_dt(m_p3.x().to_double(),m_p3.y().to_double()));
        Triangle_dt T2(Point_dt(t.m_p1.x().to_double(),t.m_p1.y().to_double()),Point_dt(t.m_p2.x().to_double(),t.m_p2.y().to_double()),Point_dt(t.m_p3.x().to_double(),t.m_p3.y().to_double()));
        return T1<T2;*/
        if((*this==t)) return false; // we do not compare identical triangles!!
        // to speed up the computation we use the triangle AABB
        ///////////////////////////////////////////////////////////////////////////////////////
        // the following code make the triangle comparison false
        /*double areathis=area(), areat=t.area();
        if(areathis<areat) return true; // put the bigger triangle firstly
        else if(areathis>areat) return false;*/
        ///////////////////////////////////////////////////////////////////////////////////////
        //if(is_strictly_to_the_left(t) || is_strictly_below(t)) return true;
        //if(t.is_strictly_to_the_left(*this) || t.is_strictly_below(*this)) return false;
        if(weak_less(t)) return true;
        if(t.weak_less(*this)) return false;
        ///////////////////////////////////////////////////////////////////////////////////////
        Point minthis = MIN3nef(m_p1,m_p2,m_p3), mint = MIN3nef(t.m_p1,t.m_p2,t.m_p3), maxthis = MAX3nef(m_p1,m_p2,m_p3), maxt = MAX3nef(t.m_p1,t.m_p2,t.m_p3);

        // we can get a weak ordering with x < t.x || x==x and y<t.y
        // we first compare the y coordinate and then the x coordinate
        if( is_first_less_than_two(minthis, maxthis, mint, maxt) ) return true;
        //else if ( is_first_less_than_two(mint, maxt, minthis, maxthis) ) return false;
        else if((minthis.x() == mint.x()) && (minthis.y() == mint.y()) &&
                (maxthis.x() == maxt.x()) && (maxthis.y() == maxt.y()))
        {
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // here minthis==mint && maxthis==maxt
            Point newminthis, newmaxthis, newmint, newmaxt, pointaloneMinthis, pointaloneMint, pointaloneMaxthis, pointaloneMaxt;
            second_min(newminthis,pointaloneMinthis);
            t.second_min(newmint,pointaloneMint);
            second_max(newmaxthis,pointaloneMaxthis);
            t.second_max(newmaxt,pointaloneMaxt);

            if( newminthis.x()==newmint.x() && newminthis.y()==newmint.y() && newmaxthis.x()==newmaxt.x() && newmaxthis.y()==newmaxt.y() )
            {
                // [newminthis,newmaxthis] and [newmint,newmaxt] describe the same segment (not necessarily a segment of a triangle),
                // so we can only compare the points which do not belong to that segment:

                if( pointaloneMinthis.x()==pointaloneMint.x() && pointaloneMinthis.y()==pointaloneMint.y() &&
                    pointaloneMaxthis.x()==pointaloneMaxt.x() && pointaloneMaxthis.y()==pointaloneMaxt.y()      )
                {
                    if(m_circumc.x()==t.m_circumc.x() && m_circumc.y()==t.m_circumc.y())
                    {
                        return critical(t);
                        /*cout << "BIG PB FOUND!!" << endl;
                        print();
                        t.print();
                        exit(-1);*/
                    }
                    else
                    {
                        // for such close triangles, we compare their circumcenter in the other direction:
                        return ((m_circumc.x()<t.m_circumc.x()) || ((m_circumc.x()==t.m_circumc.x()) && (m_circumc.y()<t.m_circumc.y())));
                    }
                }

                return is_first_less_than_two(pointaloneMinthis, pointaloneMaxthis, pointaloneMint, pointaloneMaxt);
            }

            return is_first_less_than_two(newminthis, newmaxthis, newmint, newmaxt);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        else return false;
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }

    /////////////////////////////////////////////////////////////////////////////////////
    inline bool is_strictly_to_the_left(const triangle& t) const
    {
        return (m_max_AABB.x()<t.m_min_AABB.x());
    }

    inline bool is_strictly_below(const triangle& t) const
    {
        return m_max_AABB.y()<t.m_min_AABB.y();
    }

    inline bool weak_less(const triangle& t) const
    {
        return ((m_max_AABB.x()<t.m_min_AABB.x() ||
                (m_max_AABB.x()==t.m_min_AABB.x() && m_max_AABB.y()<t.m_min_AABB.y())) ||

                ((m_min_AABB.x()<t.m_min_AABB.x())&&(m_max_AABB.x()<t.m_max_AABB.x())) ||
                (   ((m_min_AABB.x()==t.m_min_AABB.x())&&(m_max_AABB.x()==t.m_max_AABB.x()))&&
                    ((m_min_AABB.y()<t.m_min_AABB.y())&&(m_max_AABB.y()<t.m_max_AABB.y())))
                );
    }
#ifdef USE_SEGMENT
    inline bool intersect(const triangle& t) const
    {
        if(m_max_AABB.x()<=t.m_min_AABB.x()) return false; // t is to the right

        if(t.m_max_AABB.x()<=m_min_AABB.x()) return false; // t is to the left

        if(m_max_AABB.y()<=t.m_min_AABB.y()) return false; // t is above

        if(t.m_max_AABB.y()<=m_min_AABB.y()) return false; // t is below
        //return (m_circumc - t.m_circumc).squared_length()<2*(m_squared_radius+t.m_squared_radius); // in exact arithmetic we can only use that
        if( (m_circumc - t.m_circumc).squared_length()>=m_squared_radius+t.m_squared_radius+2.0*std::sqrt(m_squared_radius*t.m_squared_radius) ) return false;

        /*if( is_inside(t.m_p1) || is_inside(t.m_p2) || is_inside(t.m_p3) ||
            t.is_inside(m_p1) || t.is_inside(m_p2) || t.is_inside(m_p3)     ) return true;*/

        Triangle_dt T1(Point_dt(m_p1.x().to_double(),m_p1.y().to_double()),Point_dt(m_p2.x().to_double(),m_p2.y().to_double()),Point_dt(m_p3.x().to_double(),m_p3.y().to_double()));
        Triangle_dt T2(Point_dt(t.m_p1.x().to_double(),t.m_p1.y().to_double()),Point_dt(t.m_p2.x().to_double(),t.m_p2.y().to_double()),Point_dt(t.m_p3.x().to_double(),t.m_p3.y().to_double()));

        CGAL::Object result = CGAL::intersection(T1, T2);
        if (const CGAL::Point_2<k_delaunay> *ipoint = CGAL::object_cast<CGAL::Point_2<k_delaunay> >(&result))
        { // handle the point intersection case with *ipoint.

            return false; // strict inter
        }
        else if (const Segment_dt *iseg = CGAL::object_cast< Segment_dt >(&result))
        { // handle the segment intersection case with *iseg.

            return false; // strict inter

        }
        else if (const CGAL::Triangle_2<k_delaunay> *itri = CGAL::object_cast<CGAL::Triangle_2<k_delaunay> >(&result))
        { // handle the triangle intersection case

            return true;
        }
        else if (const std::vector< CGAL::Point_2<k_delaunay> >  *ivec = CGAL::object_cast< std::vector< CGAL::Point_2<k_delaunay> > >(&result))
        { // handle the polygon intersection case

            return true;
        }
        else
        { // handle the no intersection case.
            return false;
        }
    }

    inline bool intersect(const Segment_dt& s, std::set<unique_segment>& sub_segments) const
    {   // compute an intersection between the current triangle and the segment s, and insert in sub_segments
        // the intersection result if it is a segment (at the end we have the uniqueness)
        Triangle_dt T(Point_dt(m_p1.x().to_double(),m_p1.y().to_double()),Point_dt(m_p2.x().to_double(),m_p2.y().to_double()),Point_dt(m_p3.x().to_double(),m_p3.y().to_double()));
        CGAL::Object result = CGAL::intersection(T, s);
        if (const CGAL::Point_2<k_delaunay> *ipoint = CGAL::object_cast<CGAL::Point_2<k_delaunay> >(&result))
        { // handle the point intersection case with *ipoint.
            return false;
        }
        else if (const Segment_dt *iseg = CGAL::object_cast< Segment_dt >(&result))
        { // handle the segment intersection case with *iseg.
            if(unique_segment::isValid(iseg->point(0), iseg->point(1)))
            {
                sub_segments.insert(unique_segment(iseg->point(0), iseg->point(1)));
                return true;
            }
            return false;
        }
        else
        { // handle the no intersection case.
            return false;
        }
    }
#endif

    inline vector normal(double z1=0, double z2=0, double z3=0) const
    {
        vector Vab(m_p2.x().to_double()-m_p1.x().to_double(), m_p2.y().to_double()-m_p1.y().to_double(),z2-z1);
        vector Vac(m_p3.x().to_double()-m_p1.x().to_double(), m_p3.y().to_double()-m_p1.y().to_double(),z3-z1);

        return CGAL::cross_product(Vab, Vac);
    }

    inline bool is_3D_triangle_possible(double z1, double z2, double z3) const
    {
        return (normal(z1, z2, z3).x()>0.42);
    }

    inline double area() const
    {
        Point3d a(m_p1.x().to_double(),m_p1.y().to_double(),0), b(m_p2.x().to_double(),m_p2.y().to_double(),0), c(m_p3.x().to_double(),m_p3.y().to_double(),0);
        return triangleArea(a,b,c);
    }

    inline Point3d barycenter() const
    {
        Point3d a(m_p1.x().to_double(),m_p1.y().to_double(),0), b(m_p2.x().to_double(),m_p2.y().to_double(),0), c(m_p3.x().to_double(),m_p3.y().to_double(),0);
        return Point3d((a.x()+b.x()+c.x())/3.0, (a.y()+b.y()+c.y())/3.0, 0.0);
    }

    inline Point barycenterExact() const // throw error
    {
        return Point((m_p1.x()+m_p2.x()+m_p3.x())/3.0, (m_p1.y()+m_p2.y()+m_p3.y())/3.0, 0.0);
    }

    inline Point3d middle_side1() const
    {
        Point3d a(m_p1.x().to_double(),m_p1.y().to_double(),0), b(m_p2.x().to_double(),m_p2.y().to_double(),0);
        return Point3d((a.x()+b.x())/2.0, (a.y()+b.y())/2.0, 0.0);
    }

    inline Point3d middle_side2() const
    {
        Point3d a(m_p2.x().to_double(),m_p2.y().to_double(),0), b(m_p3.x().to_double(),m_p3.y().to_double(),0);
        return Point3d((a.x()+b.x())/2.0, (a.y()+b.y())/2.0, 0.0);
    }

    inline Point3d middle_side3() const
    {
        Point3d a(m_p3.x().to_double(),m_p3.y().to_double(),0), b(m_p1.x().to_double(),m_p1.y().to_double(),0);
        return Point3d((a.x()+b.x())/2.0, (a.y()+b.y())/2.0, 0.0);
    }

    inline bool is_inside(const Point& P) const // with a nef point
    {
        assert((P.x().denominator()!=0) && (P.y().denominator()!=0));
        //CGAL::Gmpq px(P.x().numerator(), P.x().denominator());
        //CGAL::Gmpq py(P.y().numerator(), P.y().denominator());

        Point_dt P_dt(P.x().to_double(), P.y().to_double());
        return is_inside_dt(P_dt);
    };

    inline bool is_inside(const Point_ar& P) const
    {
        //Point_dt P_dt(CGAL::to_double(P.x()),CGAL::to_double(P.y()));

        return is_inside_dt(P);
    };

    inline bool is_inside(const Point3d& P) const
    {
        Point_dt P_dt(P.x(),P.y());
        return is_inside_dt(P_dt);
    };

    inline bool is_inside_dt(const Point_dt& P) const
    {
        Triangle_dt T(Point_dt(m_p1.x().to_double(),m_p1.y().to_double()),Point_dt(m_p2.x().to_double(),m_p2.y().to_double()),Point_dt(m_p3.x().to_double(),m_p3.y().to_double()));
        /*CGAL::Gmpq p1x(m_p1.x().numerator(),m_p1.x().denominator());
        CGAL::Gmpq p1y(m_p1.y().numerator(),m_p1.y().denominator());
        CGAL::Gmpq p2x(m_p2.x().numerator(),m_p2.x().denominator());
        CGAL::Gmpq p2y(m_p2.y().numerator(),m_p2.y().denominator());
        CGAL::Gmpq p3x(m_p3.x().numerator(),m_p3.x().denominator());
        CGAL::Gmpq p3y(m_p3.y().numerator(),m_p3.y().denominator());
        Triangle_dt T(Point_dt(p1x,p1y),Point_dt(p2x,p2y),Point_dt(p3x,p3y));*/

        if(!T.is_degenerate())
        {
            //Bounded_side B = T.bounded_side ( T.vertex(0) );
            //Bounded_side B = T.bounded_side ( T.vertex(0)+0.5*(T.vertex(1)-T.vertex(0))+0.25*(T.vertex(2)-T.vertex(1)) );
            Bounded_side B = T.bounded_side ( P );
            switch(B)
            {
                    case ON_BOUNDARY:
                        //cout << "On boundary!!" << endl;
                        return false;
                        break;
                    case ON_BOUNDED_SIDE:
                        //cout << "Inside the triangle!!" << endl;
                        return true;
                        break;
                    case ON_UNBOUNDED_SIDE:
                        //cout << "Outside the triangle!!" << endl;
                        return false;
                        break;
                    default:
                        cout << "What is wrong in is_inside_cgal_triangle????" << endl;
                        assert(false);
                        return false;
                    break;
            }
        }

        return false;
    };
    /////////////////////////////////////////////////////////////////////////////////////

    inline void print() const
    {
        cout << "********************************************" << endl;
        cout << "Triangle points: P1("<< m_p1.x().to_double() << "," << m_p1.y().to_double() << ")" << "; P2("<< m_p2.x().to_double() << "," << m_p2.y().to_double() << ")" << "; P3("<< m_p3.x().to_double() << "," << m_p3.y().to_double() << ")" << endl;
        cout << "Triangle potential: " << m_potential << endl;
        cout << "Triangle circumcircle: center (" << m_circumc.x() << "," << m_circumc.y() << ") with R*R = " << m_squared_radius << endl;
        cout << "Triangle AABB: min (" << m_min_AABB.x() << "," << m_min_AABB.y() << ") and max (" << m_max_AABB.x() << "," << m_max_AABB.y() << ")" << endl;
        cout << "********************************************" << endl;
    }

    private:

    inline double flatTriangleBarycenter(Point3d & bar) const
    {
        Point3d p3d1(CGAL::to_double(m_p1.x()),CGAL::to_double(m_p1.y()),0), p3d2(CGAL::to_double(m_p2.x()),CGAL::to_double(m_p2.y()),0), p3d3(CGAL::to_double(m_p3.x()),CGAL::to_double(m_p3.y()),0);

        Vector v12=p3d2-p3d1, v23=p3d3-p3d2, v31=p3d1-p3d3;

        double d12 = v12*v12, d23 = v23*v23, d31 = v31*v31, circumRad2=0.0;

        if(d12>=d23 && d12>=d31)
        {
            bar = p3d1+0.5*v12;
            circumRad2 = d12*0.25;
        }
        else if(d23>=d12 && d23>=d31)
        {
            bar = p3d2+0.5*v23;
            circumRad2 = d23*0.25;
        }
        else
        {
            bar = p3d3+0.5*v31;
            circumRad2 = d31*0.25;
        }

        return circumRad2;
    }

    inline void second_min(Point& newmin, Point& pointalone) const
    { // the minimum without taking into account the point the more on the left (and the more in the south if needed)
      // computation is based on the fact that the triangle is counterclockwise oriented
        if(m_p1.x() <= m_p2.x() && m_p1.x() <= m_p3.x()) // one or two points are the more on the left
        {
            if(m_p1.x() == m_p2.x())
            { newmin = MINnef(m_p1,m_p3); pointalone = m_p2; }
            else
            { newmin = MINnef(m_p2,m_p3); pointalone = m_p1; }
        }
        else if(m_p2.x() <= m_p1.x() && m_p2.x() <= m_p3.x())
        {
            if(m_p2.x() == m_p3.x()) { newmin = MINnef(m_p1,m_p2); pointalone = m_p3; }
            else { newmin = MINnef(m_p1,m_p3); pointalone = m_p2; }
        }
        else
        {
            if(m_p3.x() == m_p1.x()) { newmin = MINnef(m_p2,m_p3); pointalone = m_p1; }
            else { newmin = MINnef(m_p1,m_p2); pointalone = m_p3; }
        }
    }

    inline void second_max(Point& newmax, Point& pointalone) const
    {
        if(m_p1.x() >= m_p2.x() && m_p1.x() >= m_p3.x())
        {
            if(m_p1.x() == m_p2.x()) { newmax = MAXnef(m_p1,m_p3); pointalone = m_p2; }
            else { newmax = MAXnef(m_p2,m_p3); pointalone = m_p1; }
        }
        else if(m_p2.x() >= m_p1.x() && m_p2.x() >= m_p3.x())
        {
            if(m_p2.x() == m_p3.x()) { newmax = MAXnef(m_p1,m_p2); pointalone = m_p3; }
            else { newmax = MAXnef(m_p1,m_p3); pointalone = m_p2; }
        }
        else
        {
            if(m_p3.x() == m_p1.x()) { newmax = MAXnef(m_p2,m_p3); pointalone = m_p1; }
            else { newmax = MAXnef(m_p1,m_p2); pointalone = m_p3; }
        }
    }

    inline bool is_first_less_than_two(const Point& mint1, const Point& maxt1, const Point& mint2, const Point& maxt2) const
    {
        return ((mint1.x() < mint2.x()) || (mint1.x() == mint2.x() && mint1.y() < mint2.y()) ||
                (mint1.x() == mint2.x() && mint1.y() == mint2.y() && maxt1.x() < maxt2.x()) ||
                (mint1.x() == mint2.x() && mint1.y() == mint2.y() && maxt1.x() == maxt2.x() && maxt1.y() < maxt2.y())
                );
    }

    static inline bool is_first_point_same_second(const Point& P1, const Point& P2)
    {
        return ( (P1.x()==P2.x()) && (P1.y()==P2.y()) );
    }

    static inline bool is_first_point_less_than_second(const Point& P1, const Point& P2)
    {
        return ((P1.x()<P2.x()) || ( (P1.x()==P2.x()) && (P1.y()<P2.y()) ));
    }

    inline bool critical(const triangle& t) const
    {
        bool p1issame = false, p2issame = false, p3issame = false;
        bool p1isp2 = false, p1isp3 = false;
        bool p2isp1 = false, p2isp3 = false;
        bool p3isp1 = false, p3isp2 = false;
        // we look for 2 points shared between the 2 triangles then we compared the remaining point
        if( is_first_point_same_second(m_p1, t.m_p1) )
        {
            p1issame = true;
        }
        else if( is_first_point_same_second(m_p1, t.m_p2) )
        {
            p1isp2 = true;
        }
        else if( is_first_point_same_second(m_p1, t.m_p3) )
        {
            p1isp3 = true;
        }

        if( is_first_point_same_second(m_p2, t.m_p2) )
        {
            p2issame = true;
        }
        else if( is_first_point_same_second(m_p2, t.m_p1) )
        {
            p2isp1 = true;
        }
        else if( is_first_point_same_second(m_p2, t.m_p3) )
        {
            p2isp3 = true;
        }

        if( is_first_point_same_second(m_p3, t.m_p3) )
        {
            p3issame = true;
        }
        else if( is_first_point_same_second(m_p3, t.m_p1) )
        {
            p3isp1 = true;
        }
        else if( is_first_point_same_second(m_p3, t.m_p2) )
        {
            p3isp2 = true;
        }
        ////////////////////////////////////////////////////////////
        if( p1issame && p2isp3 )
        {
            return is_first_point_less_than_second(m_p3, t.m_p2);
        }
        else if( p1isp2 && p2isp1 )
        {
            return is_first_point_less_than_second(m_p3, t.m_p3);
        }
        else if( p1isp3 && p2issame)
        {
            return is_first_point_less_than_second(m_p3, t.m_p1);
        }
        /////
        else if( p2isp1 && p3issame )
        {
            return is_first_point_less_than_second(m_p1, t.m_p2);
        }
        else if( p2issame && p3isp1 )
        {
            return is_first_point_less_than_second(m_p1, t.m_p3);
        }
        else if( p2isp3 && p3isp2)
        {
            return is_first_point_less_than_second(m_p1, t.m_p1);
        }
        /////
        else if( p3isp1 && p1isp3 )
        {
            return is_first_point_less_than_second(m_p2, t.m_p2);
        }
        else if( p3isp2 && p1issame )
        {
            return is_first_point_less_than_second(m_p2, t.m_p3);
        }
        else if( p3issame && p1isp2)
        {
            return is_first_point_less_than_second(m_p2, t.m_p1);
        }
        /////
        else if( p1issame && p2issame )
        {
            return is_first_point_less_than_second(m_p3, t.m_p3);
        }
        else if( p2issame && p3issame )
        {
            return is_first_point_less_than_second(m_p1, t.m_p1);
        }
        else if( p3issame && p1issame)
        {
            return is_first_point_less_than_second(m_p2, t.m_p2);
        }
        else
        {
            assert(false);
            return false;
        }
    }
};

#endif
