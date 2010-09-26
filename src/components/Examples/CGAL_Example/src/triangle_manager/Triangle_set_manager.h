///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Triangle_set_manager_H
#define Triangle_set_manager_H


#include "Triangle.h"
#include "Index_set.h"
#include "Unique_segment.h"
#include "../useful_methods/General_purpose.h"
#include "../point_tools/Point_reading_generation_saving.h"
#include "../triangle_tools/Triangle_reading_generation_saving.h"

// STL
#include <vector>
//////////////////////////////////////////////

//#define USE_OPEN_CV
// OpenCV
#ifdef USE_OPEN_CV
//#include <opencv/highgui.h>
#include <opencv/cv.h>
#endif
//////////////////////////////////////////////


// KD-tree node for geometrical partitioning
class node
{
    public:
    std::vector< unsigned int > m_triangle_indices;
    node* m_left, *m_right; // left and right nodes
    Point3d m_min_AABB, m_max_AABB; // approximated axis-aligned bounding box
    unsigned int m_d;
    public:
    node():m_triangle_indices(),m_left(NULL),m_right(NULL),m_min_AABB(0,0,0),m_max_AABB(0,0,0),m_d(0){}

    node(const std::vector< unsigned int >& triangle_indices, node* left, node* right, const Point3d& min_AABB, const Point3d& max_AABB,unsigned int d):m_triangle_indices(triangle_indices),m_left(left),m_right(right),m_min_AABB(min_AABB),m_max_AABB(max_AABB),m_d(d){}

    ~node()
    {
        if(m_left) delete m_left;
        if(m_right) delete m_right;
    }

    inline unsigned int size() const { return m_triangle_indices.size(); }
    inline bool isLeaf() const { return !m_left&&!m_right; }
    inline unsigned int max_depth() const { if(isLeaf()) return m_d; return std::max<int>(m_left->max_depth(),m_right->max_depth()); }

    inline double w() const { return m_max_AABB.x()-m_min_AABB.x(); }
    inline double h() const { return m_max_AABB.y()-m_min_AABB.y(); }

    inline bool contain(const Point3d& p) const // contain p
    {

        if( (m_min_AABB.x() <= p.x()) &&
            (m_min_AABB.y() <= p.y()) &&
            (m_min_AABB.z() <= p.z()) &&

            (p.x() <= m_max_AABB.x()) &&
            (p.y() <= m_max_AABB.y()) &&
            (p.z() <= m_max_AABB.z()) )
            {
                return true;
            }

        return false;
    }

    // a node contains a whole triangle iff it contains its AABB
    inline bool contain(const triangle& t) const
    {
        return contain(t.m_min_AABB) && contain(t.m_max_AABB);
    }

    // contain the segment [p1,p2]
    inline bool contain(const Point3d& p1, const Point3d& p2) const { return contain(p1) && contain(p2); }
#ifdef USE_SEGMENT
    inline bool contain(const Segment_dt& s) const // contain the segment [p1,p2]
    {
        //Point3d p1(s.point(0).x().to_double(), s.point(0).y().to_double(), 0.0), p2(s.point(1).x().to_double(), s.point(1).y().to_double(), 0.0);
        Point3d p1(CGAL::to_double(s.point(0).x()), CGAL::to_double(s.point(0).y()), 0.0), p2(CGAL::to_double(s.point(1).x()), CGAL::to_double(s.point(1).y()), 0.0);
        //Point3d p1(s.point(0).x(), s.point(0).y(), 0.0), p2(s.point(1).x(), s.point(1).y(), 0.0);

        return contain(p1, p2);
    }
#endif
    //////////////////////////////////////////////////////////////////////////////////
    inline void print() const
    {
        cout << "*********************************************" << endl;
        cout << "Node with depth == " << m_d << endl;
        cout << "There are " << m_triangle_indices.size() << " triangles in that node." << endl;
        cout << "Node AABB: min (" << m_min_AABB.x() << "," << m_min_AABB.y() << ") and max (" << m_max_AABB.x() << "," << m_max_AABB.y() << ")" << endl;
        if(!isLeaf())
        {
            cout << "left son:" << endl;
            m_left->print();
            cout << "right son:" << endl;
            m_right->print();
        }
        cout << "*********************************************" << endl;
    }
};


class triangle_manager // to mainly do geometrical queries based on triangle
{
    private:
    std::vector< triangle > m_Alltriangles;
    node* m_root;

    double m_known_feasible_triangulation_cost;

    private:

    void split_node_list(node* n, long int d) const;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // index management methods
    static void printIntSetFlux(ofstream& fileflux, const std::set< unsigned int >& ind);
    static inline void printIntSet(const std::set< unsigned int >& ind);

    static bool is_constraint_present(  const std::set< unsigned int >& constraint,
                                               const std::vector< std::set< unsigned int > >& triangle_constraints);
    static bool are_two_sets_identical(const std::set< unsigned int >& ind1, const std::set< unsigned int >& ind2);

    static inline void intersection(const std::set< unsigned int >& ind1, const std::set< unsigned int >& ind2, std::set< unsigned int >& indinter);
    static inline void union_(const std::set< unsigned int >& ind1, const std::set< unsigned int >& ind2, std::set< unsigned int >& indunion);

    void compute_AABB_given_triangle_indices(const std::vector < unsigned int >& triangle_indices, Point3d& min_AABB, Point3d& max_AABB) const;

    inline bool all_triangles_contain_p(const Point3d& P, const std::set< unsigned int >& triangle_indices) const;
    bool all_triangles_share_a_common_point(unsigned int indexNewT, const std::set< unsigned int >& triangle_indices) const;

    inline void add_new_constraint(std::set< unsigned int >& ind1, std::set< unsigned int >& ind2, std::set< unsigned int >& indadd) const;
    unsigned int validate_constraint_one_step(std::vector< std::set< unsigned int > >& in_triangle_constraints, std::vector< std::set< unsigned int > >& out_triangle_constraints) const;

    static void print_ccb (Arrangement::Ccb_halfedge_const_circulator circ);
    static void print_face (Arrangement::Face_const_handle f);
    static void print_arrangement (const Arrangement& arr);

    static Point_ar center_of_mass(Arrangement::Face_const_handle f);

    public:
    // constructors
    triangle_manager() { m_root = NULL; m_known_feasible_triangulation_cost = DBL_MAX; }
    triangle_manager(const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_triangles);
    // destructor
    ~triangle_manager() { if(m_root) delete m_root; }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // make sure you use a valid triangulation without T-junction for that!!
    void set_known_feasible_triangulation_cost(double feasible_triangulation_cost);
    inline double get_triangles_total_energy() const;

    static inline void printConstraints(const std::vector< std::set< unsigned int > >& triangle_constraints);

    inline void triangles_containing_p(const Point3d& P, const std::set< unsigned int >& triangle_indices, std::set< unsigned int >& out_triangle_indices) const;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static void generate_Arrangement_2D_from_triangles( const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& triangles,
                                                        Arrangement & arr // out
                                                        );
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void generate_2D_region_center_of_mass( Arrangement & arr, // in
                                            std::vector <Point_ar>& out_centers // out (must be not empty!!)
                                                    ) const;

    void generate_region_constraints_from_region_center_of_mass(const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& triangles) const;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef USE_SEGMENT
    static void generate_segment_set_from_triangles(const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& triangles,
                                                    std::set< unique_segment,  std::less<unique_segment> > & vseg // out
                                                    )
    {
        vseg.clear(); // because it is a set we do not need to shrink the capacity

        std::set<triangle, std::less<triangle>, std::allocator<triangle> >::const_iterator its(triangles.begin()), ite(triangles.end());
        for(;its!=ite;++its)
        {
            CGAL::Point_2<k_delaunay>  p1(its->m_p1.x().to_double(), its->m_p1.y().to_double()),
                                p2(its->m_p2.x().to_double(), its->m_p2.y().to_double()),
                                p3(its->m_p3.x().to_double(), its->m_p3.y().to_double());
            /////////////////////////////////////////
            if(unique_segment::isValid(p1, p2)) vseg.insert(unique_segment(p1, p2));
            if(unique_segment::isValid(p2, p3)) vseg.insert(unique_segment(p2, p3));
            if(unique_segment::isValid(p3, p1)) vseg.insert(unique_segment(p3, p1));
        }
    }


    static inline bool are_intersecting_segments(   const unique_segment& seg1,
                                                    const unique_segment& seg2,
                                                    std::vector <Point3d>& out_point_inter,
                                                    std::set <unique_segment>& out_unique_segment_inter)
    {
        assert((seg1.m_p1 != seg1.m_p2) && (seg2.m_p1 != seg2.m_p2));
        Segment_dt _seg1(seg1.m_p1, seg1.m_p2), _seg2(seg2.m_p1, seg2.m_p2);
        CGAL::Object result = CGAL::intersection(_seg1, _seg2);

        if (const CGAL::Point_2<k_delaunay> *ipoint = CGAL::object_cast<CGAL::Point_2<k_delaunay> >(&result))
        { // handle the point intersection case with *ipoint.
            //out_point_inter.push_back(Point3d(ipoint->x().to_double(), ipoint->y().to_double(), 0.0));
            out_point_inter.push_back(Point3d(CGAL::to_double(ipoint->x()), CGAL::to_double(ipoint->y()), 0.0));
            //out_point_inter.push_back(Point3d(ipoint->x(), ipoint->y(), 0.0));

            return true;
        }
        else if (const Segment_dt *iseg = CGAL::object_cast< Segment_dt >(&result))
        { // handle the segment intersection case with *iseg.
            if( (abs(iseg->point(0).x()-iseg->point(1).x()) < EPS) &&
                (abs(iseg->point(0).y()-iseg->point(1).y()) < EPS)  )
            { // we do not want segment reduced to a point
                return false;
            }
            out_unique_segment_inter.insert(unique_segment(iseg->point(0),iseg->point(1)));
            return true;
        }
        else
        { // handle the no intersection case.
            return false;
        }
    }


    static void compute_segments_intersection(  const std::set< unique_segment > & vseg,
                                                std::vector <Point3d>& out_point_inter,
                                                std::set <unique_segment>& out_unique_segment_inter    )

    {   // linear method => inefficient
        // this method computes all segment intersection
        std::set< unique_segment >::const_iterator its(vseg.begin()), ite(vseg.end());
        std::set< unique_segment >::const_iterator its2, ite2(vseg.end());
        for(;its!=ite;++its)
        {
            its2=its; // to not do twice the same calculus!!
            ++its2;
            for(;its2!=ite2;++its2)
            {
                //if(its==its2) continue;
                are_intersecting_segments(*its, *its2, out_point_inter, out_unique_segment_inter);
            }
        }
    }

    static void compute_recursive_segment_intersection( const std::set< unique_segment >& vseg_in,
                                                        std::set< unique_segment >& vseg_out    )

    {   // this method compute recursively unique_segment intersection till not subsegment is created
        cout << "compute_recursive_segment_intersection: Number of input segments: " << vseg_in.size() << endl;
        std::vector <unique_segment> file_unique_segment;
        std::set< unique_segment >::const_iterator its(vseg_in.begin()), ite(vseg_in.end());
        for(;its!=ite;++its)
            file_unique_segment.push_back(*its);

        unsigned int nbe = file_unique_segment.size();
        //cout << nbe << " segments!!" << endl;
        while(nbe>0)
        {
            unsigned int cur_seg_id = nbe-1;
            bool change = false;
            for(unsigned int i=0; i<nbe-1; ++i)
            {
                /////////////////////////////////////////////////////////////////////////
                // the 2 segments must be different for the processing
                if(file_unique_segment[cur_seg_id]==file_unique_segment[i])
                {   // it is possible to have that after a few steps, because
                    // different segment intersections can produce the same intersection
                    // in that case we simply remove the last element (without inserting it!!)
                    cout << "The two segment are the same." << endl;
                    file_unique_segment.pop_back(); // remove the last element
                    nbe--;
                    cur_seg_id = nbe-1;
                    continue;
                }
                /////////////////////////////////////////////////////////////////////////
                std::vector <Point3d> out_point_inter;
                std::set <unique_segment> out_unique_segment_inter;
                if(are_intersecting_segments(file_unique_segment[cur_seg_id], file_unique_segment[i], out_point_inter, out_unique_segment_inter))
                {
                    cout << "New segment intersection computation:" << endl;
                    file_unique_segment[cur_seg_id].print();
                    file_unique_segment[i].print();
                    if(!out_unique_segment_inter.empty())
                    {   // case where the intersection is a segment
                        cout << "Case where the intersection is a segment:" << endl;
                        // we have to update file_unique_segment[i] and file_unique_segment[nbe-1]
                        // plus add the intersection!!

                        unique_segment inter = *(out_unique_segment_inter.begin());
                        cout << "Resulting intersection: ";
                        inter.print();

                        if((file_unique_segment[cur_seg_id]==inter) || (file_unique_segment[i]==inter))
                        {   // remind that we cannot have file_unique_segment[nbe-1]==file_unique_segment[i]
                            cout << "Case where one segment is included in the other (intersection is one of the two segments)." << endl;
                            if(file_unique_segment[cur_seg_id]==inter)
                            {
                                if(file_unique_segment[i].m_p1==inter.m_p1)
                                {
                                    if(unique_segment::isValid(inter.m_p2,file_unique_segment[i].m_p2))
                                    {
                                        file_unique_segment[i] = unique_segment(inter.m_p2,file_unique_segment[i].m_p2);
                                        change=true;
                                    }
                                }
                                else if(file_unique_segment[i].m_p2==inter.m_p2)
                                {
                                    if(unique_segment::isValid(file_unique_segment[i].m_p1, inter.m_p1))
                                    {
                                        file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1, inter.m_p1);
                                        change=true;
                                    }
                                }
                                else if(file_unique_segment[i].m_p1==inter.m_p2)
                                {
                                    if(unique_segment::isValid(inter.m_p1, file_unique_segment[i].m_p2))
                                    {
                                        file_unique_segment[i] = unique_segment(inter.m_p1, file_unique_segment[i].m_p2);
                                        change=true;
                                    }
                                }
                                else if(file_unique_segment[i].m_p2==inter.m_p1)
                                {
                                    if(unique_segment::isValid(inter.m_p2, file_unique_segment[i].m_p1))
                                    {
                                        file_unique_segment[i] = unique_segment(inter.m_p2, file_unique_segment[i].m_p1);
                                        change=true;
                                    }
                                }
                                else
                                { // case where the intersection is strictly inside file_unique_segment[i]
                                    if( (file_unique_segment[i].m_p1-inter.m_p1).squared_length() <
                                        (file_unique_segment[i].m_p1-inter.m_p2).squared_length() )
                                    {
                                        if( unique_segment::isValid(inter.m_p2,file_unique_segment[i].m_p2) &&
                                            unique_segment::isValid(file_unique_segment[i].m_p1, inter.m_p1) )
                                        {
                                            file_unique_segment.push_back(unique_segment(inter.m_p2,file_unique_segment[i].m_p2));
                                            file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1, inter.m_p1);
                                            cout << "element ajoute:" << endl;
                                            file_unique_segment[nbe].print();
                                            change=true;
                                            nbe++;
                                        }
                                    }
                                    else
                                    {
                                        if( unique_segment::isValid(inter.m_p1,file_unique_segment[i].m_p2) &&
                                            unique_segment::isValid(file_unique_segment[i].m_p1, inter.m_p2) )
                                        {
                                            file_unique_segment.push_back(unique_segment(inter.m_p1,file_unique_segment[i].m_p2));
                                            file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1, inter.m_p2);
                                            cout << "element ajoute:" << endl;
                                            file_unique_segment[nbe].print();
                                            change=true;
                                            nbe++;
                                        }
                                    }
                                }
                            }
                            else
                            {// file_unique_segment[i]==inter
                                if(file_unique_segment[cur_seg_id].m_p1==inter.m_p1)
                                {
                                    if(unique_segment::isValid(inter.m_p2,file_unique_segment[cur_seg_id].m_p2))
                                    {
                                        file_unique_segment[cur_seg_id] = unique_segment(inter.m_p2,file_unique_segment[cur_seg_id].m_p2);
                                        change=true;
                                    }
                                }
                                else if(file_unique_segment[cur_seg_id].m_p2==inter.m_p2)
                                {
                                    if(unique_segment::isValid(file_unique_segment[cur_seg_id].m_p1, inter.m_p1))
                                    {
                                        file_unique_segment[cur_seg_id] = unique_segment(file_unique_segment[cur_seg_id].m_p1, inter.m_p1);
                                        change=true;
                                    }
                                }
                                else if(file_unique_segment[cur_seg_id].m_p1==inter.m_p2)
                                {
                                    if(unique_segment::isValid(inter.m_p1, file_unique_segment[cur_seg_id].m_p2))
                                    {
                                        file_unique_segment[cur_seg_id] = unique_segment(inter.m_p1, file_unique_segment[cur_seg_id].m_p2);
                                        change=true;
                                    }
                                }
                                else if(file_unique_segment[cur_seg_id].m_p2==inter.m_p1)
                                {
                                    if(unique_segment::isValid(inter.m_p2, file_unique_segment[cur_seg_id].m_p1))
                                    {
                                        file_unique_segment[cur_seg_id] = unique_segment(inter.m_p2, file_unique_segment[cur_seg_id].m_p1);
                                        change=true;
                                    }
                                }
                                else
                                { // case where the intersection is strictly inside file_unique_segment[nbe-1]
                                    if( (file_unique_segment[cur_seg_id].m_p1-inter.m_p1).squared_length() <
                                        (file_unique_segment[cur_seg_id].m_p1-inter.m_p2).squared_length() )
                                    {
                                        if( unique_segment::isValid(inter.m_p2,file_unique_segment[cur_seg_id].m_p2) &&
                                            unique_segment::isValid(file_unique_segment[cur_seg_id].m_p1, inter.m_p1))
                                        {
                                            file_unique_segment.push_back(unique_segment(inter.m_p2,file_unique_segment[cur_seg_id].m_p2));
                                            file_unique_segment[cur_seg_id] = unique_segment(file_unique_segment[cur_seg_id].m_p1, inter.m_p1);
                                            change=true;
                                            nbe++;
                                        }
                                    }
                                    else
                                    {
                                        if( unique_segment::isValid(inter.m_p1,file_unique_segment[cur_seg_id].m_p2) &&
                                            unique_segment::isValid(file_unique_segment[cur_seg_id].m_p1, inter.m_p2))
                                        {
                                            file_unique_segment.push_back(unique_segment(inter.m_p1,file_unique_segment[cur_seg_id].m_p2));
                                            file_unique_segment[cur_seg_id] = unique_segment(file_unique_segment[cur_seg_id].m_p1, inter.m_p2);

                                            change=true;
                                            nbe++;
                                        }
                                    }

                                }
                            }
                        }
                        else
                        { // case where the intersection is a new segment (whose extremity belongs to the initial segments)

                            cout << "Case where intersection is different from the two initial segments." << endl;

                            if(inter.m_p2==file_unique_segment[cur_seg_id].m_p2)
                            {
                                if( unique_segment::isValid(file_unique_segment[cur_seg_id].m_p1,inter.m_p1) &&
                                    unique_segment::isValid(inter.m_p2,file_unique_segment[i].m_p2) &&
                                    unique_segment::isValid(inter.m_p2,file_unique_segment[i].m_p1))
                                {
                                    file_unique_segment[cur_seg_id] = unique_segment(file_unique_segment[cur_seg_id].m_p1,inter.m_p1);

                                    if(inter.m_p1==file_unique_segment[i].m_p1)
                                        file_unique_segment[i] = unique_segment(inter.m_p2,file_unique_segment[i].m_p2);
                                    else
                                        file_unique_segment[i] = unique_segment(inter.m_p2,file_unique_segment[i].m_p1);

                                    change=true;
                                }
                            }
                            else if(inter.m_p1==file_unique_segment[cur_seg_id].m_p1)
                            {
                                if( unique_segment::isValid(inter.m_p2,file_unique_segment[cur_seg_id].m_p2) &&
                                    unique_segment::isValid(file_unique_segment[i].m_p1,inter.m_p1) &&
                                    unique_segment::isValid(file_unique_segment[i].m_p2,inter.m_p1))
                                {
                                    file_unique_segment[cur_seg_id] = unique_segment(inter.m_p2,file_unique_segment[cur_seg_id].m_p2);

                                    if(inter.m_p2==file_unique_segment[i].m_p2)
                                        file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1,inter.m_p1);
                                    else
                                        file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p2,inter.m_p1);

                                    change=true;
                                }
                            }
                            else if(inter.m_p2==file_unique_segment[cur_seg_id].m_p1)
                            {
                                if( unique_segment::isValid(inter.m_p1,file_unique_segment[cur_seg_id].m_p2) &&
                                    unique_segment::isValid(inter.m_p2,file_unique_segment[i].m_p2) &&
                                    unique_segment::isValid(inter.m_p2,file_unique_segment[i].m_p1))
                                {
                                    file_unique_segment[cur_seg_id] = unique_segment(inter.m_p1,file_unique_segment[cur_seg_id].m_p2);

                                    if(inter.m_p1==file_unique_segment[i].m_p1)
                                        file_unique_segment[i] = unique_segment(inter.m_p2,file_unique_segment[i].m_p2);
                                    else
                                        file_unique_segment[i] = unique_segment(inter.m_p2,file_unique_segment[i].m_p1);

                                    change=true;
                                }
                            }
                            else //if(inter.m_p1==file_unique_segment[nbe-1].m_p2)
                            {
                                if( unique_segment::isValid(inter.m_p2,file_unique_segment[cur_seg_id].m_p1) &&
                                    unique_segment::isValid(inter.m_p1,file_unique_segment[i].m_p2) &&
                                    unique_segment::isValid(file_unique_segment[i].m_p1, inter.m_p1))
                                {
                                    file_unique_segment[cur_seg_id] = unique_segment(inter.m_p2,file_unique_segment[cur_seg_id].m_p1);

                                    if(inter.m_p2==file_unique_segment[i].m_p1)
                                        file_unique_segment[i] = unique_segment(inter.m_p1,file_unique_segment[i].m_p2);
                                    else
                                        file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1, inter.m_p1);

                                    change=true;
                                }
                            }

                            if(change)
                            {
                                file_unique_segment.push_back(inter);
                                nbe++;
                            }
                        }

                        cout << "After intersection computation the two initial segments are:" << endl;
                        file_unique_segment[cur_seg_id].print();
                        file_unique_segment[i].print();

                        cout << "and the number of element is: " << nbe << endl;

                        if(change) break;
                    }
                    else if(!out_point_inter.empty())
                    {   // CRITICAL PART!!!!!
                        // case where the intersection is a point
                        cout << "Case where the intersection is a point:" << endl;
                        CGAL::Point_2<k_delaunay> pt((out_point_inter.begin())->x(),(out_point_inter.begin())->y());
                        cout << "Resulting intersection: ";
                        cout << "p(" << pt.x() << "," << pt.y() << ")" << endl;
                        // we simply split the two segments if the intersecting point is not a common extremity
                        if( file_unique_segment[i].m_p1!=file_unique_segment[cur_seg_id].m_p1 &&
                            file_unique_segment[i].m_p1!=file_unique_segment[cur_seg_id].m_p2 &&
                            file_unique_segment[i].m_p2!=file_unique_segment[cur_seg_id].m_p1 &&
                            file_unique_segment[i].m_p2!=file_unique_segment[cur_seg_id].m_p2)
                        {
                            cout << "No common extremity case: " << endl;

                            if( pt==file_unique_segment[i].m_p1 || pt==file_unique_segment[i].m_p2 )
                            {
                                if( unique_segment::isValid(pt,file_unique_segment[cur_seg_id].m_p2) &&
                                    unique_segment::isValid(file_unique_segment[cur_seg_id].m_p1,pt)    )
                                {
                                    file_unique_segment.push_back(unique_segment(pt,file_unique_segment[cur_seg_id].m_p2));
                                    file_unique_segment[cur_seg_id] = unique_segment(file_unique_segment[cur_seg_id].m_p1,pt);
                                    //cout << "element ajoute:" << endl;
                                    //file_unique_segment[nbe].print();
                                    nbe++;

                                    change=true;
                                    break;
                                }
                            }
                            else if( pt==file_unique_segment[cur_seg_id].m_p1 || pt==file_unique_segment[cur_seg_id].m_p2 )
                            {
                                if( unique_segment::isValid(pt,file_unique_segment[i].m_p2) &&
                                    unique_segment::isValid(file_unique_segment[i].m_p1,pt))
                                {
                                    file_unique_segment.push_back(unique_segment(pt,file_unique_segment[i].m_p2));
                                    file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1,pt);
                                    //cout << "element ajoute:" << endl;
                                    //file_unique_segment[nbe].print();
                                    nbe++;

                                    change=true;
                                    break;
                                }
                            }
                            else
                            {

                                if( unique_segment::isValid(pt,file_unique_segment[cur_seg_id].m_p2) &&
                                    unique_segment::isValid(pt,file_unique_segment[i].m_p2) &&
                                    unique_segment::isValid(file_unique_segment[cur_seg_id].m_p1,pt) &&
                                    unique_segment::isValid(file_unique_segment[i].m_p1,pt))
                                {
                                    file_unique_segment.push_back(unique_segment(pt,file_unique_segment[cur_seg_id].m_p2));
                                    file_unique_segment.push_back(unique_segment(pt,file_unique_segment[i].m_p2));

                                    file_unique_segment[cur_seg_id] = unique_segment(file_unique_segment[cur_seg_id].m_p1,pt);
                                    file_unique_segment[i] = unique_segment(file_unique_segment[i].m_p1,pt);
                                    //cout << "elements ajoutes:" << endl;
                                    //file_unique_segment[nbe].print();
                                    //file_unique_segment[nbe+1].print();
                                    nbe+=2;

                                    change=true;
                                    break;
                                }

                            }
                        }
                        else
                        {
                            cout << "Common extremity case: nothing to do." << endl;
                            assert(!change);
                        }
                    }
                    else
                    {
                        cout << "compute_recursive_segment_intersection: unexpected error occurs." << endl;
                        assert(false); // error
                    }
                }
            }

            if(!change)
            {   // if the current segment does not generate a new segment by intersecting other segment we have
                // finished for this segment
                vseg_out.insert(file_unique_segment[cur_seg_id]);
                file_unique_segment.pop_back(); // remove the last element
                nbe--;
                cout << "No new intersection with the current segment; the number of element is now: " << nbe << endl;
            }
        }
        cout << "compute_recursive_segment_intersection: Number of output segments: " << vseg_out.size() << endl;
    }


    /**
    This method computes all the intersections between triangle edges and gets a set of points

    DO NOT compute the Delaunay triangulation: because for the time being we cannot see how to
    efficiently build a triangulation of all the regions which respect the region boundaries
     **/
    void get_region_segment_offset( const std::set<triangle, std::less<triangle> >& triangles, // initial triangle set
                                    std::vector<Point3d>& point_offset) const
    {
        point_offset.clear();
        if(point_offset.capacity())
        {
            std::vector<Point3d> empty1(point_offset);
            point_offset.swap(empty1);
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        std::set< unique_segment,  std::less<unique_segment> > vseg;
        generate_segment_set_from_triangles(triangles, vseg); // we generate a segment set from the triangle set (OK)
        cout << "The initial segment set is generated. " << vseg.size() << " segments."  << endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        compute_minimal_segment_intersection_and_generate_point_offset(vseg, point_offset);
    }

    void compute_minimal_segment_intersection_and_generate_point_offset(const std::set< unique_segment,std::less<unique_segment> >& vseg,
                                                                        std::vector<Point3d>& point_offset) const
    {
        std::set<unique_segment> all_sub_segments; // segments used to compute the point offsets
        std::set< unique_segment,std::less<unique_segment> >::const_iterator its(vseg.begin()), ite(vseg.end());
        for(;its!=ite;++its)
        {
            Segment_dt s(its->m_p1, its->m_p2);
            std::set< unsigned int > out_inside_triangle_indices;
            std::set<unique_segment> all_sub_segments_for_one_s; // subsegment from s
            if(unique_segment::isValid(its->m_p1, its->m_p2))
            {
                ////////////////////////////////////////////////////////////////////////////////////////////////
                get_triangles_which_intersect_s(s, out_inside_triangle_indices, all_sub_segments_for_one_s);
                ////////////////////////////////////////////////////////////////////////////////////////////////
                // Here we have to transform all_sub_segments_for_one_s into successive subsegment (doing intersection)
                //std::set<unique_segment> inter_all_sub_segments_for_one_s;
                compute_recursive_segment_intersection(all_sub_segments_for_one_s, all_sub_segments); // OK
            }
        }
        //////////////////////////////////////////////////////////////////////////////////////////
        // once we have the correct subsegment set, we can compute point_offset
        generate_point_offset(all_sub_segments, point_offset); // OK
    }

    /**
    generate_point_offset generate the vector of 3d points to use to compute
    the constraints (given the segments). 2 points are added for each segment (left and right).
    **/
    inline void generate_point_offset(  const std::set<unique_segment>& all_sub_segments, // correct set!!
                                        std::vector<Point3d>& point_offset) const
    {
        double offset = 1e-6;
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        point_offset.reserve(point_offset.size()+all_sub_segments.size());
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        std::set<unique_segment>::const_iterator its(all_sub_segments.begin()), ite(all_sub_segments.end());
        for(;its!=ite;++its)
        {
            //Vector D(its->m_p2.x().to_double()-its->m_p1.x().to_double(), its->m_p2.y().to_double()-its->m_p1.y().to_double(),0);
            Vector D(CGAL::to_double(its->m_p2.x())-CGAL::to_double(its->m_p1.x()), CGAL::to_double(its->m_p2.y())-CGAL::to_double(its->m_p1.y()),0);
            //Vector D(its->m_p2.x()-its->m_p1.x(), its->m_p2.y()-its->m_p1.y(),0);

            //Point3d start(its->m_p1.x().to_double(), its->m_p1.y().to_double(), 0);
            Point3d start(CGAL::to_double(its->m_p1.x()), CGAL::to_double(its->m_p1.y()), 0);
            //Point3d start(its->m_p1.x(), its->m_p1.y(), 0);

            Vector orthoD(-D.y(),D.x(),0);
            orthoD = orthoD/std::sqrt(orthoD*orthoD);
            // we add two offset points:
            point_offset.push_back(Point3d(start+0.5*D+offset*orthoD));
            point_offset.push_back(Point3d(start+0.5*D-offset*orthoD));
        }
    }
    /**

    **/

    inline void generate_region_constraints_from_segment_offset(const std::set<triangle, std::less<triangle> >& triangles,
                                                                std::vector< std::set< unsigned int >* >& point_constraints) const
    {
        point_constraints.clear();
        /////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<Point3d> point_offset;  // we do not use a set for that (we will ensure latter the uniqueness
                                            // of constraints)!!
        get_region_segment_offset(triangles, point_offset); // we generate the point set
        /////////////////////////////////////////////////////////////////////////////////////////////
        get_triangle_constraints_from_points(point_offset, point_constraints);
    }
#endif
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline void get_triangles_which_contain_p_linearly( const node*n, const Point3d& P,
                                                        std::set<triangle, std::less<triangle>, std::allocator<triangle> >& out_inside_triangles
                                                        ) const
    { // linear search in the number of triangles (sub-optimal)
        if(!n) return; // we do nothing for the null node
        //out_inside_triangles.clear(); // clear the initial triangle list if needed
        std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
        for(;it_index!=it_index_e;++it_index)
        {
            if(m_Alltriangles[*it_index].is_inside(P)) out_inside_triangles.insert(triangle(m_Alltriangles[*it_index])); // if the point is inside we add it to the corresponding set
        }
    }

    inline void get_triangles_which_contain_p_linearly( const node*n, const Point3d& P,
                                                        std::set< unsigned int >& out_inside_triangle_indices // we used set to be sure to have a triangle appearing only once
                                                        ) const
    { // linear search in the number of triangles (sub-optimal)
        if(!n) return; // we do nothing for the null node
        //out_inside_triangles.clear(); // clear the initial triangle list if needed
        std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
        for(;it_index!=it_index_e;++it_index)
        {
            if(m_Alltriangles[*it_index].is_inside(P)) out_inside_triangle_indices.insert(*it_index); // if the point is inside we add it to the corresponding set
        }
    }

    inline void get_triangles_which_contain_p_linearly( const node*n, const Point_ar& P,
                                                        std::set< unsigned int >& out_inside_triangle_indices // we used set to be sure to have a triangle appearing only once
                                                        ) const
    { // linear search in the number of triangles (sub-optimal)
        if(!n) return; // we do nothing for the null node
        //out_inside_triangles.clear(); // clear the initial triangle list if needed
        std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
        for(;it_index!=it_index_e;++it_index)
        {
            if(m_Alltriangles[*it_index].is_inside(P)) out_inside_triangle_indices.insert(*it_index); // if the point is inside we add it to the corresponding set
        }
    }

    // sublinear method:
    void get_triangles_which_contain_p(const Point3d& P, std::set<triangle, std::less<triangle>, std::allocator<triangle> >& out_inside_triangles) const
    {
        double split_val;
        node* cn = m_root;
        while(!cn->isLeaf())
        { // while the current node is not a leaf node, we continue to go deper in the binary tree
            if(cn->w() > cn->h())
            {// splitting along Ox
                split_val = cn->m_min_AABB.x() + 0.5*cn->w();
                if(P.x()<=split_val)
                {
                    if( (P.x()==split_val) && (cn->m_left->size() > cn->m_right->size()) )
                    {
                        cn = cn->m_right;
                    }
                    else cn = cn->m_left;
                }
                else
                    cn = cn->m_right;
            }
            else
            {// splitting along Oy
                split_val = cn->m_min_AABB.y() + 0.5*cn->h();
                if(P.y()<=split_val)
                {
                    if( (P.y()==split_val) && (cn->m_left->size() > cn->m_right->size()) )
                    {
                        cn = cn->m_right;
                    }
                    else cn = cn->m_left;
                }
                else
                    cn = cn->m_right;
            }
        }

        get_triangles_which_contain_p_linearly(cn, P, out_inside_triangles);
    }

    void get_triangles_which_contain_p(const Point3d& P, std::set< unsigned int >& out_inside_triangle_indices) const
    {
        double split_val;
        node* cn = m_root;
        while(!cn->isLeaf())
        { // while the current node is not a leaf node, we continue to go deper in the binary tree
            if(cn->w() > cn->h())
            {// splitting along Ox
                split_val = cn->m_min_AABB.x() + 0.5*cn->w();
                if(P.x()<=split_val)
                {
                    if( (P.x()==split_val) && (cn->m_left->size() > cn->m_right->size()) )
                    {
                        cn = cn->m_right;
                    }
                    else cn = cn->m_left;
                }
                else
                    cn = cn->m_right;
            }
            else
            {// splitting along Oy
                split_val = cn->m_min_AABB.y() + 0.5*cn->h();
                if(P.y()<=split_val)
                {
                    if( (P.y()==split_val) && (cn->m_left->size() > cn->m_right->size()) )
                    {
                        cn = cn->m_right;
                    }
                    else cn = cn->m_left;
                }
                else
                    cn = cn->m_right;
            }
        }

        get_triangles_which_contain_p_linearly(cn, P, out_inside_triangle_indices);
    }

    void get_triangles_which_contain_p(const Point_ar& P, std::set< unsigned int >& out_inside_triangle_indices) const
    {
        double split_val;
        node* cn = m_root;
        while(!cn->isLeaf())
        { // while the current node is not a leaf node, we continue to go deper in the binary tree
            if(cn->w() > cn->h())
            {// splitting along Ox
                split_val = cn->m_min_AABB.x() + 0.5*cn->w();
                if(P.x()<=split_val)
                {
                    if( (P.x()==split_val) && (cn->m_left->size() > cn->m_right->size()) )
                    {
                        cn = cn->m_right;
                    }
                    else cn = cn->m_left;
                }
                else
                    cn = cn->m_right;
            }
            else
            {// splitting along Oy
                split_val = cn->m_min_AABB.y() + 0.5*cn->h();
                if(P.y()<=split_val)
                {
                    if( (P.y()==split_val) && (cn->m_left->size() > cn->m_right->size()) )
                    {
                        cn = cn->m_right;
                    }
                    else cn = cn->m_left;
                }
                else
                    cn = cn->m_right;
            }
        }

        get_triangles_which_contain_p_linearly(cn, P, out_inside_triangle_indices);
    }

    inline void get_triangle_constraints_from_points(   const std::vector<Point3d>& pt,
                                                        std::vector< std::set< unsigned int >* >& point_constraints,
                                                        bool remove_doubloon=true) const
    {   // that is the method which uses a lot of RAM!!!
        point_constraints.clear();
        if(point_constraints.capacity())
        {
            std::vector< std::set< unsigned int >* > empty1(point_constraints);
            point_constraints.swap(empty1); // shrink capacity to zero
        }
        point_constraints.reserve(pt.size()); // enlarge capacity to the desired size
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        cout << "get_triangle_constraints_from_points: Number max of constraints: " << pt.size() << endl;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::set< index_set > point_constraints_inter; // used for removing the doubloon
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        unsigned int nbel=0, nbmaxel=0;

        std::vector<Point3d>::const_iterator its(pt.begin()), ite(pt.end());
        for(;its!=ite;++its)
        {
            std::set< unsigned int > out_inside_triangle_indices;
            //////////////////////////////////////////////////////////////////////////////
            get_triangles_which_contain_p(*its, out_inside_triangle_indices);
            if(!out_inside_triangle_indices.empty())
            {
                if(remove_doubloon)
                    point_constraints_inter.insert(index_set(out_inside_triangle_indices)); // will add only a set of constraint if new
                else
                {
                    point_constraints.push_back(new std::set< unsigned int >(out_inside_triangle_indices));
                    //cout << " Size current constraint: " << out_inside_triangle_indices.size() << endl;
                    nbel+=out_inside_triangle_indices.size();
                    if(out_inside_triangle_indices.size()>nbmaxel) nbmaxel = out_inside_triangle_indices.size();
                }
            }
        }

        if(remove_doubloon)
        {   // the final constraint set will not contain doubloon
            cout << " Number of constraints: " << point_constraints_inter.size() << endl;
            cout << "get_triangle_constraints_from_points: just before the convertion..." << endl;
            index_set::from_set_set_to_vector_set(point_constraints_inter, point_constraints);
        }
        else
        {
            cout << "get_triangle_constraints_from_points: " << nbel*1.0/point_constraints.size() << " elements per constraint on average. " << nbmaxel << " at max." << endl;
            cout << "get_triangle_constraints_from_points: " << nbel*4.0/1048576 << " Mo needed." << endl;
        }
    }

    inline void get_triangles_which_contain_triangle_barycenter(const triangle& t,
                                                                std::set< unsigned int >& out_inside_triangle_indices) const
    {
        out_inside_triangle_indices.clear(); // at the end we only want the constraints related to the
                                             // triangle barycenter
        Point3d b = t.barycenter();
        get_triangles_which_contain_p(b, out_inside_triangle_indices);
    }

    inline void get_triangle_barycenter_constraints(const std::vector<triangle>& vt,
                                                    std::vector< std::set< unsigned int > >& triangle_constraints) const
    {
        triangle_constraints.clear(); // at the end we only want the constraints related to vt barycenters
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::set< index_set > triangle_constraints_inter; // used for removing the doubloon
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<triangle>::const_iterator its(vt.begin()), ite(vt.end());
        for(;its!=ite;++its)
        {
            std::set< unsigned int > out_inside_triangle_indices;
            get_triangles_which_contain_triangle_barycenter(*its, out_inside_triangle_indices);
            index_set Id(out_inside_triangle_indices);
            triangle_constraints_inter.insert(Id); // will add only a set of constraint if new
        }
        // the final constraint set will not contain doubloon:
        index_set::from_set_set_to_vector_set(triangle_constraints_inter, triangle_constraints);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef USE_SEGMENT
    void get_triangles_which_intersect_t_linearly(  node*n, const triangle& t,
                                                    std::set< unsigned int >& out_inside_triangle_indices
                                                    ) const
    {
        if(!n) return; // we do nothing for the null node
        //out_inside_triangles.clear(); // clear the initial triangle list if needed
        std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
        for(;it_index!=it_index_e;++it_index)
        {
            if(m_Alltriangles[*it_index].intersect(t))
            {
                // here and only here the current triangle can intersect t
                out_inside_triangle_indices.insert(*it_index);
            }
        }
    }

    void get_triangles_which_intersect_t(   const triangle& t,
                                            std::set< unsigned int >& out_inside_triangle_indices) const
    { // submodular method
        node* cn = m_root;
        while(!cn->isLeaf() && (cn->m_left->contain(t) || cn->m_right->contain(t)))
        {
            if(cn->m_left->contain(t))
            {
                cn = cn->m_left;
            }
            else
            {
                cn = cn->m_right;
            }
        }
        get_triangles_which_intersect_t_linearly(cn, t, out_inside_triangle_indices);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void get_triangles_which_intersect_s_linearly(  node*n, const Segment_dt& s,
                                                    std::set< unsigned int >& out_inside_triangle_indices,
                                                    std::set<unique_segment>& sub_segments // to get the subsegment of interest
                                                    ) const
    {
        if(!n) return; // we do nothing for the null node
        //out_inside_triangles.clear(); // clear the initial triangle list if needed
        std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
        for(;it_index!=it_index_e;++it_index)
        {
            if(m_Alltriangles[*it_index].intersect(s, sub_segments))
            {
                // here and only here the current triangle can intersect t
                out_inside_triangle_indices.insert(*it_index);
            }
        }
    }

    void get_triangles_which_intersect_s(   const Segment_dt& s,
                                            std::set< unsigned int >& out_inside_triangle_indices,
                                            std::set<unique_segment>& sub_segments) const
    { // submodular method
        node* cn = m_root;

        while( !cn->isLeaf() && (cn->m_left->contain(s) || cn->m_right->contain(s)) )
        {
            if(cn->m_left->contain(s))
            {
                cn = cn->m_left;
            }
            else
            {
                cn = cn->m_right;
            }
        }
        get_triangles_which_intersect_s_linearly(cn, s, out_inside_triangle_indices, sub_segments);
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void get_some_triangles_which_intersect_t(const triangle& t, std::set< unsigned int >& out_inside_triangle_indices) const
    { // give all the triangle which intersect the triangle
        out_inside_triangle_indices.clear();
        Point3d P1(t.m_p1.x().to_double(),t.m_p1.y().to_double(),0), P2(t.m_p2.x().to_double(),t.m_p2.y().to_double(),0), P3(t.m_p3.x().to_double(),t.m_p3.y().to_double(),0);
        Vector V1(P2-P1), V2(P3-P1);
        // the precision should be adapted to the triangle size!!
        double v1prec=0.1, v2prec=0.1;
        for(;v1prec<=0.9;v1prec+=0.2)
            for(v2prec=0.1;v2prec<=0.9;v2prec+=0.2)
            {
                Point3d P = P1+v1prec*V1+v2prec*V2;
                get_triangles_which_contain_p(P, out_inside_triangle_indices);
            }
    }

    void get_max_triangles_which_intersect_t_same_region(const triangle& t, std::set< unsigned int >& out_inside_triangle_indices) const
    { // give all the triangle which intersect the triangle on the same part of t: it gives a max constraint (limited by precision)!!
        out_inside_triangle_indices.clear();
        unsigned int nbtinconstraint = 0;
        Point3d P1(t.m_p1.x().to_double(),t.m_p1.y().to_double(),0), P2(t.m_p2.x().to_double(),t.m_p2.y().to_double(),0), P3(t.m_p3.x().to_double(),t.m_p3.y().to_double(),0);
        Vector V1(P2-P1), V2(P3-P1);
        // the precision should be adapted to the triangle size!!
        double v1prec=0.1, v2prec=0.1;
        for(;v1prec<=0.9;v1prec+=0.2)
            for(v2prec=0.1;v2prec<=0.9;v2prec+=0.2)
            {
                std::set< unsigned int > out_inside_triangle_indices_inter;
                Point3d P = P1+v1prec*V1+v2prec*V2;
                get_triangles_which_contain_p(P, out_inside_triangle_indices_inter);
                unsigned int nbtNewconstraint = out_inside_triangle_indices_inter.size();
                if(nbtNewconstraint>nbtinconstraint)
                {
                    out_inside_triangle_indices = out_inside_triangle_indices_inter;
                    nbtinconstraint = nbtNewconstraint;
                }
            }
    }

    inline bool is_worthwhile_constraint(const std::set< unsigned int >& current_constraint, const std::vector< std::set< unsigned int > >& current_constraint_set) const
    { // linear method in the number of constraints!! not efficient
        unsigned int size_current_constraint = current_constraint.size();
        // particular case:
        if(size_current_constraint <= 1) return false;
        // general case:
        //cout << "current constraint: " << endl;
        //printIntSet(current_constraint);
        std::vector< std::set< unsigned int > >::const_iterator its(current_constraint_set.begin()), ite(current_constraint_set.end());
        for(;its!=ite;++its)
        {
            if( are_two_sets_identical(current_constraint, *its)  )
            {
                 return false; // the constraint is not worthwhile
            }
        }

        return true;
    }

    void get_partial_t_constraints(const triangle& t, set < index_set >& out_triangle_constraints) const
    { // give all the triangle which intersect the triangle on the same part of t: it gives a max constraint (limited by precision)!!
        //unsigned int nbtMaxconstraint = 1; // number of constraints of the longuest constraint set currently present (we do not want constraint of one triangle!!)
        Point3d P1(t.m_p1.x().to_double(),t.m_p1.y().to_double(),0), P2(t.m_p2.x().to_double(),t.m_p2.y().to_double(),0), P3(t.m_p3.x().to_double(),t.m_p3.y().to_double(),0);
        Vector V1(P2-P1), V2(P3-P1);
        // the precision should be adapted to the triangle size!!
        double precision = 0.2;
        double v1prec=0.1, v2prec=0.1;
        for(;v1prec<=0.9;v1prec+=precision)
            for(v2prec=0.1;v2prec<=0.9;v2prec+=precision)
            {
                std::set< unsigned int > out_inside_triangle_indices_inter;
                Point3d P = P1+v1prec*V1+v2prec*V2;
                get_triangles_which_contain_p(P, out_inside_triangle_indices_inter);
                unsigned int nbtNewconstraint = out_inside_triangle_indices_inter.size();
                if(nbtNewconstraint>1)
                {
                    index_set id(out_inside_triangle_indices_inter);
                    out_triangle_constraints.insert(id);
                }

                /*


                // we got a new set of constraints which can be added to the final constraint set if:
                // - its number of constraints is strictly greater than any previous constraint
                // - its intersection size with any previous constraint is < min(nbtNewconstraint,currentconstraintSize)
                //   which means that at least one new element of the new constraint is different and thus worthwhile
                if(nbtNewconstraint>nbtMaxconstraint)
                {
                    out_triangle_constraints.push_back(out_inside_triangle_indices_inter);
                    nbtMaxconstraint = nbtNewconstraint; // we update the current longuest constraint set
                }
                else if(is_worthwhile_constraint(out_inside_triangle_indices_inter, out_triangle_constraints))
                {
                    out_triangle_constraints.push_back(out_inside_triangle_indices_inter);
                }*/
            }
    }

    void get_partial_constraint_set(std::vector< std::set< unsigned int > >& out_triangle_constraints) const
    {
        out_triangle_constraints.clear(); // we do not use old constraints!!
        //std::vector< std::set< unsigned int > > out_triangle_constraints_inter;
        //unsigned int tid=0;
        set < index_set > unique_triangle_constraints;
        std::vector< triangle >::const_iterator its(m_Alltriangles.begin()), ite(m_Alltriangles.end());
        for(;its!=ite;++its)
        {
            /*std::set< unsigned int > out_inside_triangle_indices;
            //get_some_triangles_which_intersect_t(*its, out_inside_triangle_indices); // index of triangles which intersect a part of t
            get_max_triangles_which_intersect_t_same_region(*its, out_inside_triangle_indices);
            out_triangle_constraints.push_back(out_inside_triangle_indices);*/

            get_partial_t_constraints(*its, unique_triangle_constraints);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /*std::set< unsigned int >::const_iterator its2(out_inside_triangle_indices.begin()), ite2(out_inside_triangle_indices.end());
            for(;its2!=ite2;++its2)
            {
                std::set< unsigned int > OK_constraints;
                OK_constraints.insert(tid);
                if(*its2>tid) // instead of using != we use > for avoiding doubloon
                {
                    OK_constraints.insert(*its2); // index of one triangle intersecting t
                    out_triangle_constraints_inter.push_back(OK_constraints);
                }
            }*/
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //tid++;
        }
        index_set::from_set_set_to_vector_set(unique_triangle_constraints, out_triangle_constraints);
/*
        printConstraints(out_triangle_constraints_inter);
        unsigned int nbchanges=validate_constraint_one_step(out_triangle_constraints_inter, out_triangle_constraints);
        while(nbchanges>0)
        {
            cout << "Number of changes = " << nbchanges << endl;
            out_triangle_constraints_inter = out_triangle_constraints; // we restart from the last constraint set
            printConstraints(out_triangle_constraints_inter);
            nbchanges=validate_constraint_one_step(out_triangle_constraints_inter, out_triangle_constraints);
        };
        */
    }
    ///////////////////////////////////////////////////////////////////////////////
    inline void generate_2Dtriangle_file() const
     {
        ofstream flux("./triangles.txt", ios::out | ios::trunc);
        // open the file:
        if(!flux)
        {
            cout << "Cannot open the file: ./triangles.txt!!!" << endl;
            exit(-1);
        }
        //flux << "# 2D triangles list (8 coordinates for the 3 points in the counter-clockwise order+shape potential+triangle area)\n";
        //flux << "# Total area of the convex polygon = " << total_polygon_area << endl;
        //flux << mytriangles.size() << endl;
        std::vector< triangle >::const_iterator sit(m_Alltriangles.begin()), send(m_Alltriangles.end());
        for(;sit!=send;++sit)
        {
            //cout << "Triangle " << cpt << " is composed of P1("<< sit->m_p1.x().to_double() << "," << sit->m_p1.y().to_double() << ")" << "; P2("<< sit->m_p2.x().to_double() << "," << sit->m_p2.y().to_double() << ")" << "; P3("<< sit->m_p3.x().to_double() << "," << sit->m_p3.y().to_double() << ")" << endl;
            flux << sit->m_p1.x().to_double() << " " << sit->m_p1.y().to_double() << " " << sit->m_p2.x().to_double() << " " << sit->m_p2.y().to_double() << " " << sit->m_p3.x().to_double() << " "<< sit->m_p3.y().to_double();
            flux << " " << sit->m_potential << " " << sit->area() << endl;
        }
        flux.close();
    }

     static inline void generate_DEM_3Dtriangle_file(   const vector <unsigned int>& DEM,
                                                        unsigned int nbline,
                                                        unsigned int nbcolumn)
     {
        ofstream flux("./trianglesDEM.txt", ios::out | ios::trunc);
        // open the file:
        if(!flux)
        {
            cout << "Cannot open the file: ./trianglesDEM.txt!!!" << endl;
            exit(-1);
        }
        assert(nbline>0);
        assert(nbcolumn>0);

        double height_scale = 0.2/512.0*std::min<double>(nbline,nbcolumn);
        for(long int line=0; line<(long int)nbline-1; line+=1) // Oy
            for(unsigned int col=0; col<nbcolumn-1; col+=1) // Ox
            { // recall that the DEM frame is different than the cartesian frame!!
                // 2 triangles per pixel
                    Point3d p1(col,(long int)(nbline)-1-line,height_scale*get_image_value(DEM, nbline, nbcolumn, col, line)),
                            p2(col,(long int)(nbline)-2-line,height_scale*get_image_value(DEM, nbline, nbcolumn, col, line+1)),
                            p3(col+1,(long int)(nbline)-2-line,height_scale*get_image_value(DEM, nbline, nbcolumn, col+1, line+1)),
                            p4(col+1,(long int)(nbline)-1-line,height_scale*get_image_value(DEM, nbline, nbcolumn, col+1, line));

                    flux << p1.x() << " " << p1.y() << " " << p1.z() << " "; // p1
                    flux << p2.x() << " " << p2.y() << " " << p2.z() << " "; // p2
                    flux << p3.x() << " " << p3.y() << " " << p3.z() << " "; // p3
                    flux << " " << triangleShapePotential(p1, p2, p3) << " " << triangleArea(p1, p2, p3) << endl;

                    flux << p3.x() << " " << p3.y() << " " << p3.z() << " "; // p3
                    flux << p4.x() << " " << p4.y() << " " << p4.z() << " "; // p4
                    flux << p1.x() << " " << p1.y() << " " << p1.z() << " "; // p1
                    flux << " " << triangleShapePotential(p3, p4, p1) << " " << triangleArea(p3, p4, p1) << endl;
            }
        flux.flush();
        flux.close();
    }

     inline void generate_3Dtriangle_file(   const vector <unsigned int>& DEM, unsigned int nbline, unsigned int nbcolumn) const
     {
        ofstream flux("./triangles3D.txt", ios::out | ios::trunc);
        // open the file:
        if(!flux)
        {
            cout << "Cannot open the file: ./triangles3D.txt!!!" << endl;
            exit(-1);
        }
        flux << m_known_feasible_triangulation_cost << endl;
        double height_scale = 0.2/512.0*std::min<double>(nbline,nbcolumn);
        std::vector< triangle >::const_iterator sit(m_Alltriangles.begin()), send(m_Alltriangles.end());
        for(;sit!=send;++sit)
        { // recall that the DEM frame is different than the cartesian frame!!
            flux << sit->m_p1.x().to_double() << " " <<
                    sit->m_p1.y().to_double() << " " <<
                    height_scale*get_image_value(   DEM, nbline, nbcolumn,
                                                    sit->m_p1.x().to_double(),
                                                    nbline-1-sit->m_p1.y().to_double()) << " ";
            flux << sit->m_p2.x().to_double() << " " <<
                    sit->m_p2.y().to_double() << " " <<
                    height_scale*get_image_value(   DEM, nbline, nbcolumn,
                                                    sit->m_p2.x().to_double(),
                                                    nbline-1-sit->m_p2.y().to_double()) << " ";
            flux << sit->m_p3.x().to_double() << " "<<
                    sit->m_p3.y().to_double() << " " <<
                    height_scale*get_image_value(   DEM, nbline, nbcolumn,
                                                    sit->m_p3.x().to_double(),
                                                    nbline-1-sit->m_p3.y().to_double()) << " ";
            flux << " " << sit->m_potential << " " << sit->area() << endl;
        }
        flux.flush();
        flux.close();
    }
    ///////////////////////////////////////////////////////////////////////////////
#ifdef USE_OPEN_CV    
	 void generate_region_constraints_render_triangles(int minRegionSize, std::vector< std::set< unsigned int > >& triangle_constraints, const std::vector< unsigned int >&  triangle_indices_pb) const // modified Victor method
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(!triangle_indices_pb.empty()) // case where we are looking for specific constraint set
        { // we will zoom on that region to add new constraints and improve the results
            Point3d min_AABB, max_AABB;
            compute_AABB_given_triangle_indices(triangle_indices_pb, min_AABB, max_AABB);
            node* n = new node(triangle_indices_pb, NULL, NULL, min_AABB, max_AABB,0);
            split_node_list(n, 8);
            generate_region_constraints_render_triangles(n, triangle_constraints, minRegionSize);
            delete n;
        }
        else generate_region_constraints_render_triangles(m_root, triangle_constraints, minRegionSize); // initial set of constraints
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
#endif
    inline void generate_region_constraints_from_points(const std::vector<Point_ar>& pt) const
    {
        ofstream file("./regions.txt", ios::out | ios::trunc);
        // open the file:
        if(!file)
        {
            cout << "Cannot open the file: ./regions.txt!!!" << endl;
            exit(-1);
        }
        /////////////////////////////////////////////////////////////////////////////////
        std::vector<Point_ar>::const_iterator its(pt.begin()), ite(pt.end());
        for(;its!=ite;++its)
        {
            std::set< unsigned int > out_inside_triangle_indices;
            //////////////////////////////////////////////////////////////////////////////
            get_triangles_which_contain_p(*its, out_inside_triangle_indices);

            //assert(out_inside_triangle_indices.size()>0); // to check that at least one triangle is present per region
            if(out_inside_triangle_indices.size()>0)
            {
                std::set< unsigned int >::const_iterator its2(out_inside_triangle_indices.begin()), ite2(out_inside_triangle_indices.end());
                for(;its2!=ite2;++its2) file << *its2 << " ";
                file << endl;
            }

        }

        // close the file:
        file.flush();
        file.close();
    }

    static void generate_region_file_given_constraints(const std::vector< std::set< unsigned int > >& triangle_constraints)
    {
        ofstream file("./regions.txt", ios::out | ios::trunc);
        // open the file:
        if(!file)
        {
            cout << "Cannot open the file: ./regions.txt!!!" << endl;
            exit(-1);
        }
        //flux << "# 2D regions list (the indices refer to triangle)\n";
        std::vector< std::set< unsigned int > >::const_iterator its(triangle_constraints.begin()), ite(triangle_constraints.end());
        for(;its!=ite;++its)
        {
            std::set< unsigned int >::const_iterator its2((*its).begin()), ite2((*its).end());
            for(;its2!=ite2;++its2) file << *its2 << " ";
            file << endl;
        }
        // close the file:
        file.flush();
        file.close();
    }

    static void generate_region_file_given_constraints(const std::vector< std::set< unsigned int >* >& triangle_constraints)
    {
        ofstream file("./regions.txt", ios::out | ios::trunc);
        // open the file:
        if(!file)
        {
            cout << "Cannot open the file: ./regions.txt!!!" << endl;
            exit(-1);
        }
        //flux << "# 2D regions list (the indices refer to triangle)\n";
        std::vector< std::set< unsigned int >* >::const_iterator its(triangle_constraints.begin()), ite(triangle_constraints.end());
        for(;its!=ite;++its)
        {
            std::set< unsigned int >::const_iterator its2((**its).begin()), ite2((**its).end());
            for(;its2!=ite2;++its2) file << *its2 << " ";
            file << endl;
        }
        // close the file:
        file.flush();
        file.close();
    }

    void generate_region_file_given_2Dpoints(char* point_file) const
    {
        ifstream file(point_file, ios::in);  // on ouvre le fichier en lecture
        if(!file)
        {
            cerr << "Cannot open file: " << point_file << "!" << endl;
            exit(-1);
        }
        ofstream file2("./regions.txt", ios::out | ios::trunc);
        // open the file:
        if(!file2)
        {
            cout << "Cannot open the file: ./regions.txt!!!" << endl;
            exit(-1);
        }
        string line;
        while( getline(file, line) )
        { // for each point
            istringstream iss(line);
            double d1, d2;
            iss >> d1; iss >> d2;
            //cout << d1 << ", " << d2 << endl;
            Point3d P(d1,d2,0);
            std::set< unsigned int > out_inside_triangle_indices;
            get_triangles_which_contain_p(P, out_inside_triangle_indices);
            printIntSetFlux(file2, out_inside_triangle_indices);
        }
        // close the flux:
        file.close();
        file2.flush();
        file2.close();
    }
    ///////////////////////////////////////////////////////////////////////////////
    void efficient_scalar_field_integration_over_triangles( const vector <unsigned int>& DEM,
                                                            unsigned int nbline,
                                                            unsigned int nbcolumn)
    {   // efficient in case of locallity!!
        // To integrate efficiently over triangles, we will use directly their member m_potential,
        // so we first modify the value of m_potential according to its weight and  we then compute the
        // integral...
        double m_potential_weight=0.25; // if greater than 0.75 => not enough adaptative

        //double max_inter=0.0, min_inter=1e8, max_pot=0.0, min_pot=1e8; // debug

        double height_scale = 0.2/512.0*std::min<double>(nbline,nbcolumn);

        // we firstly initialize to zero every triangle potential
        unsigned int nbt = m_Alltriangles.size();
        for(unsigned int i=0; i<nbt; ++i)
        {
            m_Alltriangles[i].m_potential = 0.0;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        double integral_weight = (1.0-m_potential_weight)/10.0; // we divide by 10 to decrease the value of the total erreur while getting approximately the same results (I mean adaptivity is still good)
        double h=0.25; // precision of integration (1 for one pixel and less for better precision)
        double h2=h*h;
        for(double line=0.0; line<double(nbline); line+=h)
            for(double col=0.0; col<double(nbcolumn); col+=h)
            { // recall that the DEM frame is different than the cartesian frame!!
                Point3d sampleP_DEM_FRAME(col, line, 0), // sample point in the DEM frame
                        sampleP_cartesian_FRAME(col, double(nbline)-1-line, 0); // sample point in the cartesian frame
                // this is the part where we collect the triangles (in the cartesian frame)
                std::set< unsigned int > out_inside_triangle_indices;
                ////////////////////////////////////////////////////////////////////////////////////////
                get_triangles_which_contain_p(sampleP_cartesian_FRAME, out_inside_triangle_indices);
                ////////////////////////////////////////////////////////////////////////////////////////
                // this is the part where we update their integral
                std::set< unsigned int >::const_iterator its(out_inside_triangle_indices.begin()),ite(out_inside_triangle_indices.end());
                for(;its!=ite;++its)
                {
                    Point3d p1(m_Alltriangles[*its].m_p1.x().to_double(),m_Alltriangles[*its].m_p1.y().to_double(),0),
                            p2(m_Alltriangles[*its].m_p2.x().to_double(),m_Alltriangles[*its].m_p2.y().to_double(),0),
                            p3(m_Alltriangles[*its].m_p3.x().to_double(),m_Alltriangles[*its].m_p3.y().to_double(),0);

                    m_Alltriangles[*its].m_potential +=  get_interpolation_triangle_abs_error(p1, p2, p3, DEM, nbline, nbcolumn, sampleP_DEM_FRAME);
                }
            }

        // to take into account triangle shape potential
        for(unsigned int i=0; i<nbt; ++i)
        {
            Point3d p1( m_Alltriangles[i].m_p1.x().to_double(),
                        m_Alltriangles[i].m_p1.y().to_double(),
                        height_scale*get_image_value(DEM, nbline, nbcolumn, m_Alltriangles[i].m_p1.x().to_double(),
                                                                            nbline-1-m_Alltriangles[i].m_p1.y().to_double())),
                        p2( m_Alltriangles[i].m_p2.x().to_double(),
                            m_Alltriangles[i].m_p2.y().to_double(),
                            height_scale*get_image_value(DEM, nbline, nbcolumn, m_Alltriangles[i].m_p2.x().to_double(),
                                                                                nbline-1-m_Alltriangles[i].m_p2.y().to_double())),
                        p3( m_Alltriangles[i].m_p3.x().to_double(),
                            m_Alltriangles[i].m_p3.y().to_double(),
                            height_scale*get_image_value(DEM, nbline, nbcolumn, m_Alltriangles[i].m_p3.x().to_double(),
                                                                                nbline-1-m_Alltriangles[i].m_p3.y().to_double()));

            double t_pot = triangleShapePotential(p1,p2,p3);
            // 1e20 is the limit for my LP solver; current max: 1e8 (to decrease if the current LP pb cannot be solve because of this limit)
            // we try to make the LP pb quicker solvable:
            if(t_pot > 1e8) t_pot = 1e8; // value truncation to make it solvable

            //if(t_pot<min_pot) min_pot=t_pot;// debug
            //if(t_pot>max_pot) max_pot=t_pot;// debug

            //if(std::sqrt(h2 * m_Alltriangles[i].m_potential)<min_inter) min_inter=std::sqrt(h2 * m_Alltriangles[i].m_potential);// debug
            //if(std::sqrt(h2 * m_Alltriangles[i].m_potential)>max_inter) max_inter=std::sqrt(h2 * m_Alltriangles[i].m_potential);// debug


            m_Alltriangles[i].m_potential = m_potential_weight * t_pot + integral_weight * h2 * m_Alltriangles[i].m_potential; // the right part is a volume
        }

            //cout << "triangle interpolation error: MIN = " << min_inter << " MAX = " << max_inter << endl;
            //cout << "triangle shape potential: MIN = " << min_pot << " MAX = " << max_pot << endl;
    }

    ///////////////////////////////////////////////////////////////////////////////
    inline void print() const
    {
        cout << "There are " << m_Alltriangles.size() << " triangles at all." << endl;
        if(m_root) m_root->print();
    }

    ///////////////////////////////////////////////////////////////////////////////
    // TEST
#ifdef USE_SEGMENT
    static void test_compute_recursive_segment_intersection()
    {
        cout << "TEST: compute_recursive_segment_intersection" << endl;
        std::set< unique_segment > vseg_in, vseg_out;
        CGAL::Point_2<k_delaunay> p1(0,0), p2(10,0), p3(4,0), p5(2,0), p6(6,0), p7(9,0), p8(0,-1), p9(2,1), p10(0.5,0.5), p11(0.666667,0.666667);
        unique_segment s1(p1,p2), s2(p1,p3), s3(p5,p2), s4(p6,p7), s5(p8,p9), s6(p10,p11), s7(p1,p11);
        vseg_in.insert(s1); vseg_in.insert(s2); vseg_in.insert(s3); vseg_in.insert(s4); vseg_in.insert(s5);
        vseg_in.insert(s6); vseg_in.insert(s7);
        cout << "Initial segment set: " << endl;
        std::set< unique_segment >::const_iterator its(vseg_in.begin()), ite(vseg_in.end());
        for(;its!=ite;++its) its->print();
        ////////////////////////////////////////////////////////////
        compute_recursive_segment_intersection(vseg_in, vseg_out);
        ////////////////////////////////////////////////////////////
        cout << "Final segment set: " << endl;
        its = vseg_out.begin();
        ite = vseg_out.end();
        for(;its!=ite;++its) its->print();
        cout << "END TEST: compute_recursive_segment_intersection" << endl;
    }
#endif

    void test(const std::vector< Point3d >& cloud_point) const
    {
        cout << "TEST of triangle_manager" << endl;
        cout << "Number of triangles: " << m_Alltriangles.size() << endl;
        cout << "Built binary tree max depth = " << m_root->max_depth() << endl;
        //print();
        std::set<triangle, std::less<triangle>, std::allocator<triangle> > out_inside_triangles1, out_inside_triangles2;
        std::set< unsigned int > out_inside_triangle_indices;
        long start_time_s = clock()/CLOCKS_PER_SEC;
        //Point3d P = m_root->m_min_AABB+0.25*(m_root->m_max_AABB-m_root->m_min_AABB);
        std::vector< Point3d >::const_iterator its(cloud_point.begin()), ite(cloud_point.end());
        for(;its!=ite;++its)
            get_triangles_which_contain_p(*its, out_inside_triangles1);

        print_elapsed_time(start_time_s);
        cout << "The point P is contained in " << out_inside_triangles1.size() << " triangles: " << endl;
        //std::set<triangle, std::less<triangle>, std::allocator<triangle> >::const_iterator sit(out_inside_triangles1.begin()), send(out_inside_triangles1.end());
        //for(;sit!=send;++sit) sit->print();

        start_time_s = clock()/CLOCKS_PER_SEC;
        std::vector< std::set< unsigned int > > out_triangle_constraints;
        get_partial_constraint_set(out_triangle_constraints);

        cout << "triangle partial constraint generation:";
        print_elapsed_time(start_time_s);

        cout << "TEST of triangle_manager successfully terminated" << endl;
    }

    private:
#ifdef USE_OPEN_CV
    void generate_region_constraints_render_triangles(node* n, std::vector< std::set< unsigned int > >& triangle_constraints, int minRegionSize) const // modified Victor method
    { // add new region constraints given the current node without clearing the initial constraint set
        if(!n) return;

        unsigned int lim=0;
        // this approach is inefficient in terms of computation time and in terms of tightness of solution of the Linear Program
        if(n->m_d<std::min(n->max_depth(),lim) && !n->isLeaf())
        {
            generate_region_constraints_render_triangles(n->m_left, triangle_constraints, minRegionSize);
            generate_region_constraints_render_triangles(n->m_right, triangle_constraints, minRegionSize);
        }

        // here n is a leaf node
        double minx = n->m_min_AABB.x(), miny = n->m_min_AABB.y(), maxx = n->m_max_AABB.x(), maxy = n->m_max_AABB.y();
        assert(maxx>=minx && maxy>=miny);
        unsigned int nTriangles = n->size(); // number of triangles present in that node
        double *vertices = new double[6*nTriangles];
        // the bounding box is slightly increased
        minx -= 1;
        miny -= 1;
        maxx += 1;
        maxy += 1;
        // set the maximum scale available:
        double precision = 32768.0/(std::min(n->max_depth(),lim)+1); // with 16384 precision the solver solution is not enough often tight
        int scale = int(precision/std::max((maxx-minx),(maxy-miny))); // give reasonable computational time

        for(unsigned int i = 0; i < nTriangles; i++)
        {
            vertices[6*i] = m_Alltriangles[n->m_triangle_indices[i]].m_p1.x().to_double();
            vertices[6*i+1] = m_Alltriangles[n->m_triangle_indices[i]].m_p1.y().to_double();
            vertices[6*i+2] = m_Alltriangles[n->m_triangle_indices[i]].m_p2.x().to_double();
            vertices[6*i+3] = m_Alltriangles[n->m_triangle_indices[i]].m_p2.y().to_double();
            vertices[6*i+4] = m_Alltriangles[n->m_triangle_indices[i]].m_p3.x().to_double();
            vertices[6*i+5] = m_Alltriangles[n->m_triangle_indices[i]].m_p3.y().to_double();
        }
        // the bounding box min becomes the new origin (+the scale can be changed):
        for(unsigned int i = 0; i < nTriangles*6; i+=2)
        {
            vertices[i] = (vertices[i]-minx)*scale;
            vertices[i+1] = (vertices[i+1]-miny)*scale;
        }

        // image creation:
        cout << "minx = " << minx << " miny = " << miny << " maxx = " << maxx << " maxy = " << maxy << endl;
        cout << "scale = " << scale << " (maxx-minx)*scale = " << (maxx-minx)*scale << " (maxy-miny)*scale = " << (maxy-miny)*scale << endl;
        int imWidth = (unsigned int)((maxx-minx)*scale);
        int imHeight = (unsigned int)((maxy-miny)*scale);

        cout << "Create an image with width = " << imWidth << " and imHeight = " << imHeight << endl;
        IplImage *buf = cvCreateImage(cvSize(imWidth, imHeight), IPL_DEPTH_8U, 1);

        cout << "Image created." << endl;
        cout << "Creating lines..." << endl;
        cvSet(buf, cvScalar(0));
        for(unsigned int i = 0; i < nTriangles; i++)
        {
            cvLine(buf, cvPoint(vertices[i*6], vertices[i*6+1]), cvPoint(vertices[i*6+2], vertices[i*6+3]), cvScalar(255), 1, 4);
            cvLine(buf, cvPoint(vertices[i*6], vertices[i*6+1]), cvPoint(vertices[i*6+4], vertices[i*6+5]), cvScalar(255), 1, 4);
            cvLine(buf, cvPoint(vertices[i*6+2], vertices[i*6+3]), cvPoint(vertices[i*6+4], vertices[i*6+5]), cvScalar(255), 1, 4);
        }
        cout << "Lines created." << endl;
        delete [] vertices;
        cout << "Flood filling..." << endl;
        //          im   starting point  new val for painting
        cvFloodFill(buf, cvPoint(0,0),    cvScalar(128),        cvScalarAll(0), cvScalarAll(0), NULL, 8);
        // since the image is slightly bigger than the number_type_delaunayresented geometry domain (we enlarge the bounding box)
        // we will firstly fill the contour
        CvConnectedComp c;
        // for all pixels:
        for(int x = 0; x < imWidth; x++)
            for(int y = 0; y < imHeight; y++)
                if( ((uchar*)(buf->imageData + buf->widthStep*y))[x] == 0 ) // if a pixel is black: it has not been visited yet
                {
                    cvFloodFill(buf, cvPoint(x, y), cvScalar(128), cvScalarAll(0), cvScalarAll(0), &c, 8);

                    if(c.area < minRegionSize) // we are not interested in too small region
                        continue;

                    std::set< unsigned int > out_inside_triangle_indices;
                    Point3d P(double(x+0.5)/scale + minx, double(y+0.5)/scale + miny,0);
                    //printf("%lf %lf\n", double(x+0.5)/scale + minx, double(y+0.5)/scale + miny);
                    get_triangles_which_contain_p(P, out_inside_triangle_indices);
                    // we only add non empty constraint set (and only new constraint!!)
                    if( !out_inside_triangle_indices.empty() &&
                        !is_constraint_present(out_inside_triangle_indices, triangle_constraints)) triangle_constraints.push_back(out_inside_triangle_indices);
                    //fprintf(out, "%lf %lf\n", double(x+0.5)/scale + minx, double(y+0.5)/scale + miny);

                    ((uchar*)(buf->imageData + buf->widthStep*y))[x] = 1;
                }
        cvReleaseImage(&buf);
        cout << "Flood fill done." << endl;
    }
#endif
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a region represent either a triangle or an intersection of several triangles (it is used to set constraints on recovering triangles)
class region
{
    public:
        Nef_polyhedron m_r; // the geometrical representation of the region
        std::vector< unsigned int > m_indices; // the indices of all triangles involved in the intersection which has produced m_r
        // for the region, to not get bored by approximation error, I will use approximated center and radius (with some care)
        Point3d m_circumc; // approximated circumcenter
        double m_squared_radius; // approximated circumradius
        Point3d m_min_AABB, m_max_AABB; // approximated axis-aligned bounding box
    public:

        region():m_r(Nef_polyhedron::EMPTY), m_indices(),m_circumc(), m_squared_radius(), m_min_AABB(), m_max_AABB() {}

        region(const Nef_polyhedron& r, const std::vector< unsigned int >& indices, const Point3d& circumc, double squared_radius, const Point3d& min_AABB, const Point3d& max_AABB):m_r(r), m_indices(indices), m_circumc(circumc), m_squared_radius(squared_radius), m_min_AABB(min_AABB), m_max_AABB(max_AABB) {}

        ~region(){}

        inline bool may_intersect(const region& r)
        {
            if(m_max_AABB.x()<=r.m_min_AABB.x()) return false; // r is to the right

            if(r.m_max_AABB.x()<=m_min_AABB.x()) return false; // r is to the left

            if(m_max_AABB.y()<=r.m_min_AABB.y()) return false; // r is above

            if(r.m_max_AABB.y()<=m_min_AABB.y()) return false; // r is below

            //return (m_circumc - r.m_circumc).squared_length()<2*(m_squared_radius+r.m_squared_radius);
            return (m_circumc - r.m_circumc).squared_length()<m_squared_radius+r.m_squared_radius+2.0*std::sqrt(m_squared_radius*r.m_squared_radius); // better test if sqrt is usable
        }

        inline void print() const
        {
        cout << "********************************************" << endl;
        cout << "Region circumcircle: center (" << m_circumc.x() << "," << m_circumc.y() << ") with R*R = " << m_squared_radius << endl;
        cout << "Region AABB: min (" << m_min_AABB.x() << "," << m_min_AABB.y() << ") and max (" << m_max_AABB.x() << "," << m_max_AABB.y() << ")" << endl;
        cout << "********************************************" << endl;
        }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////



#endif
