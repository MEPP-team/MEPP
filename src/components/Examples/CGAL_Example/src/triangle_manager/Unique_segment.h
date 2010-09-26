///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Unique_segment_H
#define Unique_segment_H

#include "CGAL_includes_and_types.h"
#define EPS double(1e-8)

class unique_segment
{
    public:
    CGAL::Point_2< k_delaunay > m_p1, m_p2;
    public:
    unique_segment() {}
    unique_segment(const CGAL::Point_2<k_delaunay>& p1, const CGAL::Point_2<k_delaunay>& p2)
    {
        assert(isValid(p1, p2)); // we do not want segment restricted to a unique point
        if( (p1.x()<p2.x()) || (abs(p1.x()-p2.x()) < EPS && p1.y()<p2.y()) )
        {
            m_p1 = p1;
            m_p2 = p2;
        }
        else
        {
            m_p1 = p2;
            m_p2 = p1;
        }
    }

    ~unique_segment() {}

    inline bool operator== (const unique_segment& seg) const
    {
        return (
                    (  (abs(m_p1.x()-seg.m_p1.x()) < EPS) && (abs(m_p1.y()-seg.m_p1.y()) < EPS) &&
                       (abs(m_p2.x()-seg.m_p2.x()) < EPS) && (abs(m_p2.y()-seg.m_p2.y()) < EPS)  )
                       ||
                    (  (abs(m_p2.x()-seg.m_p1.x()) < EPS) && (abs(m_p2.y()-seg.m_p1.y()) < EPS) &&
                       (abs(m_p1.x()-seg.m_p2.x()) < EPS) && (abs(m_p1.y()-seg.m_p2.y()) < EPS)  )
                );
        /*return ((m_p1.x()==seg.m_p1.x()) && (m_p1.y()==seg.m_p1.y()) &&
                (m_p2.x()==seg.m_p2.x()) && (m_p2.y()==seg.m_p2.y())) ||
                ((m_p2.x()==seg.m_p1.x()) && (m_p2.y()==seg.m_p1.y()) &&
                (m_p1.x()==seg.m_p2.x()) && (m_p1.y()==seg.m_p2.y()));*/
    }

    inline bool operator!= (const unique_segment& seg) const
    {
        return !(*this==seg);
    }

    inline bool operator< (const unique_segment& seg) const
    {
        return ((m_p1.x()<seg.m_p1.x()) || (m_p1.x()==seg.m_p1.x()&&m_p1.y()<seg.m_p1.y()) ||
                (m_p1.x()==seg.m_p1.x()&&m_p1.y()==seg.m_p1.y() && m_p2.x()<seg.m_p2.x()) ||
                (m_p1.x()==seg.m_p1.x()&&m_p1.y()==seg.m_p1.y() && m_p2.x()==seg.m_p2.x() && m_p2.y()<seg.m_p2.y()) );
    }

    static inline bool isValid(const CGAL::Point_2<k_delaunay>& p1, const CGAL::Point_2<k_delaunay>& p2)
    {
        return ((abs(p1.x()-p2.x()) >= EPS) || (abs(p1.y()-p2.y())>= EPS));
    }

    inline void print() const
    {
        cout << "Segment made of: p1(" << m_p1.x() << "," << m_p1.y() << ") and p2(" << m_p2.x() << "," << m_p2.y() << ")" << endl;
    }
};





#endif // Unique_segment_H
