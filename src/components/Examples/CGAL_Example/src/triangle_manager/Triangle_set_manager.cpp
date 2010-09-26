///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#include "Triangle_set_manager.h"


triangle_manager::triangle_manager(const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_triangles)
{
        m_known_feasible_triangulation_cost = DBL_MAX;
        ////////////////////////////////////////////////////////////////////
        double maximum = 1e8;
        Point3d min_AABB(maximum,maximum,0), max_AABB(-maximum,-maximum,0);
        unsigned int cpt = 0;
        std::vector < unsigned int > indices;
        std::set<triangle, std::less<triangle>, std::allocator<triangle> >::const_iterator sit(set_of_triangles.begin()), send(set_of_triangles.end());
        for(;sit!=send;++sit)
        {
            min_AABB = MIN2(min_AABB,sit->m_min_AABB);
            max_AABB = MAX2(max_AABB,sit->m_max_AABB);
            //cout << "min_AABB.x() = " << min_AABB.x() << " min_AABB.y() = " << min_AABB.y() << " min_AABB.z() = " << min_AABB.z() << endl;
            m_Alltriangles.push_back(*sit); // add a new triangle
            indices.push_back(cpt);
            cpt++;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        unsigned int max_d = 0; // maximum allowed depth
        if(!set_of_triangles.size())
        {
            cout << "WARNING: no triangle have been put int the triangle_manager!!!!" << endl;
            min_AABB = max_AABB = Point3d(0,0,0);
        }
        else
        {
            max_d = 3*int((double)log((double)set_of_triangles.size()));//2*int(log(set_of_triangles.size())); // maximum allowed depth
        }

        m_root = new node(indices, NULL, NULL, min_AABB, max_AABB,0);

        if(0<max_d) split_node_list(m_root,max_d-1);
}

void triangle_manager::split_node_list(node* n, long int d) const
{
        if(!n || n->m_left || n->m_right) return; // for null node we do nothing or if one child does not exist

        std::vector< unsigned int > triangle_indices_left, triangle_indices_right;
        Point3d maxAABB_left, minAABB_right;
        // we split along the longest side
        double split_val;
        if(n->w() > n->h())
        { // splitting along Ox
            split_val = n->m_min_AABB.x() + 0.5*n->w();
            maxAABB_left = Point3d(split_val, n->m_max_AABB.y(), 0);
            minAABB_right = Point3d(split_val, n->m_min_AABB.y(), 0);
            std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
            for(;it_index!=it_index_e;++it_index)
            {
                if(m_Alltriangles[*it_index].m_max_AABB.x() <= split_val)
                {
                    triangle_indices_left.push_back(*it_index);
                }
                else if(m_Alltriangles[*it_index].m_min_AABB.x() >= split_val)
                {
                    triangle_indices_right.push_back(*it_index);
                }
                else
                { // push the triangle in both list
                    triangle_indices_left.push_back(*it_index);
                    triangle_indices_right.push_back(*it_index);
                }
            }
        }
        else
        { // splitting along Oy
            split_val = n->m_min_AABB.y() + 0.5*n->h();
            maxAABB_left = Point3d(n->m_max_AABB.x(),split_val, 0);
            minAABB_right = Point3d(n->m_min_AABB.x(),split_val, 0);
            std::vector< unsigned int >::const_iterator it_index(n->m_triangle_indices.begin()), it_index_e(n->m_triangle_indices.end());
            for(;it_index!=it_index_e;++it_index)
            {
                if(m_Alltriangles[*it_index].m_max_AABB.y() <= split_val)
                {
                    triangle_indices_left.push_back(*it_index);
                }
                else if(m_Alltriangles[*it_index].m_min_AABB.y() >= split_val)
                {
                    triangle_indices_right.push_back(*it_index);
                }
                else
                { // push the triangle in both list
                    triangle_indices_left.push_back(*it_index);
                    triangle_indices_right.push_back(*it_index);
                }
            }
        }

        // put the children of n:
        n->m_left = new node(triangle_indices_left, NULL, NULL, n->m_min_AABB, maxAABB_left,n->m_d + 1);
        n->m_right = new node(triangle_indices_right, NULL, NULL, minAABB_right, n->m_max_AABB,n->m_d + 1);

        if(d>0)
        {
            if(n->size()!=triangle_indices_left.size() && triangle_indices_left.size()>20) split_node_list(n->m_left, d-1);
            if(n->size()!=triangle_indices_right.size() && triangle_indices_right.size()>20) split_node_list(n->m_right, d-1);
        }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// index management methods
inline void triangle_manager::printIntSetFlux(ofstream& fileflux, const std::set< unsigned int >& ind)
{
    std::set< unsigned int >::const_iterator its(ind.begin()),ite(ind.end());
    for(;its!=ite;++its)
    {
        fileflux << *its << " ";
    }
    fileflux << endl;
}

inline void triangle_manager::printIntSet(const std::set< unsigned int >& ind)
{
    cout << "The set contains: ";
    std::set< unsigned int >::const_iterator its(ind.begin()),ite(ind.end());
    for(;its!=ite;++its)
    {
        cout << *its << " ";
    }
    cout << endl;
}

inline bool triangle_manager::is_constraint_present(  const std::set< unsigned int >& constraint,
                                           const std::vector< std::set< unsigned int > >& triangle_constraints)
{   // suboptimal linear method
    if(triangle_constraints.empty()) return false;
    std::vector< std::set< unsigned int > >::const_iterator its(triangle_constraints.begin()),ite(triangle_constraints.end());
    for(;its!=ite;++its)
        if(*its == constraint) return true; // it works since the global operator == is well defined for set

    return false;
}

inline bool triangle_manager::are_two_sets_identical(const std::set< unsigned int >& ind1, const std::set< unsigned int >& ind2)
{ // sublinear method!!!
    std::set< unsigned int >::const_iterator its(ind1.begin()),ite(ind1.end());
    for(;its!=ite;++its)
    {
        if(ind2.find(*its)==ind2.end()) // if we did not find out the value => the set are not the same
            return false;
    }

    return true;
}

inline void triangle_manager::intersection(const std::set< unsigned int >& ind1, const std::set< unsigned int >& ind2, std::set< unsigned int >& indinter)
{
    indinter.clear();
    std::set< unsigned int >::const_iterator its(ind1.begin()),ite(ind1.end());
    for(;its!=ite;++its)
    {
        if(ind2.find(*its)!=ind2.end()) // if we find out the value in both ind1 and ind2 we add it to the indinter
            indinter.insert(*its);
    }
}


inline void triangle_manager::union_(const std::set< unsigned int >& ind1, const std::set< unsigned int >& ind2, std::set< unsigned int >& indunion)
{
    indunion.clear();
    indunion.insert(ind1.begin(), ind1.end());
    indunion.insert(ind2.begin(), ind2.end());
}

inline void triangle_manager::compute_AABB_given_triangle_indices(const std::vector < unsigned int >& triangle_indices, Point3d& min_AABB, Point3d& max_AABB) const
{
    double maximum(1e8);
    if(triangle_indices.empty()) maximum = 0;
    min_AABB = Point3d(maximum,maximum,0);
    max_AABB = Point3d(-maximum,-maximum,0);
    std::vector < unsigned int >::const_iterator its(triangle_indices.begin()), ite(triangle_indices.end());
    for(;its!=ite;++its)
    {
        min_AABB = MIN2(min_AABB,m_Alltriangles[*its].m_min_AABB);
        max_AABB = MAX2(max_AABB,m_Alltriangles[*its].m_max_AABB);
    }
}

inline bool triangle_manager::all_triangles_contain_p(const Point3d& P, const std::set< unsigned int >& triangle_indices) const
{
    std::set< unsigned int >::const_iterator its(triangle_indices.begin()), ite(triangle_indices.end());
    for(;its!=ite;++its)
    {
        if( !(m_Alltriangles[*its].is_inside(P)) ) return false; // if the point is not inside we return false
    }
    return true;
}

bool triangle_manager::all_triangles_share_a_common_point(unsigned int indexNewT, const std::set< unsigned int >& triangle_indices) const
{ // unusefull method
    assert(indexNewT<m_Alltriangles.size());
    Point3d P1(m_Alltriangles[indexNewT].m_p1.x().to_double(),m_Alltriangles[indexNewT].m_p1.y().to_double(),0), P2(m_Alltriangles[indexNewT].m_p2.x().to_double(),m_Alltriangles[indexNewT].m_p2.y().to_double(),0), P3(m_Alltriangles[indexNewT].m_p3.x().to_double(),m_Alltriangles[indexNewT].m_p3.y().to_double(),0);
    Vector V1(P2-P1), V2(P3-P1);
    // the precision should be adapted to the triangle size!!
    double v1prec=0.1, v2prec=0.1;
    for(;v1prec<=0.9;v1prec+=0.1)
        for(v2prec=0.1;v2prec<=0.9;v2prec+=0.1)
        {
            Point3d P = P1+v1prec*V1+v2prec*V2;
            if(all_triangles_contain_p(P,triangle_indices)) return true;
        }

    return false;
}

inline void triangle_manager::add_new_constraint(std::set< unsigned int >& ind1, std::set< unsigned int >& ind2, std::set< unsigned int >& indadd) const
{ // add a new constraint and update the two input variables consequently (because these two input variables should not add twice the same constraint)
    std::set< unsigned int > indinter;
    indadd.clear();
    intersection(ind1, ind2, indinter);
    unsigned int nbcomElt = indinter.size();
    // no constrainst to add if the number of common indices is strictly less than 2
    if(nbcomElt<2) return; // no intersection

    std::set< unsigned int >::const_iterator its(indinter.begin()),ite(indinter.end());
    // at least these two triangles are in the new constraint:
    indadd.insert(*its); ++its;
    indadd.insert(*its); ++its;
    //if(nbcomElt==2) return; // one common triangle
    for(;its!=ite;++its)
    { // for each remaining common triangle, we add it to the current constraint set
      // if it shares a common point with all the current triangles in the constraint set
        if( all_triangles_share_a_common_point(*its, indadd) ) indadd.insert(*its);
    }

    // we have to remove the indices from indadd present in both entries: false!!! because a triangle of
    // the intersection can also be still present in part in one triangle!!
    /*its = indadd.begin(); ite = indadd.end();
    for(;its!=ite;++its)
    {
        ind1.erase(*its);
        ind2.erase(*its);
    }
    */
}

unsigned int triangle_manager::validate_constraint_one_step(std::vector< std::set< unsigned int > >& in_triangle_constraints, std::vector< std::set< unsigned int > >& out_triangle_constraints) const
{
    unsigned int nbchanges=0;
    out_triangle_constraints.clear(); // we clear it because we built a new version of it!!
    std::vector< std::set< unsigned int > >::iterator its(in_triangle_constraints.begin()), ite(in_triangle_constraints.end());
    for(;its!=ite;++its)
    {
        std::vector< std::set< unsigned int > >::iterator its2(its+1), ite2(in_triangle_constraints.end());
        for(;its2!=ite2;++its2)
        {
            std::set< unsigned int > indadd;
            add_new_constraint(*its, *its2, indadd); // find out a common intersection and update entries consequently
            unsigned int nbEltComRegion = indadd.size();
            if(nbEltComRegion>1)
            {
                out_triangle_constraints.push_back(indadd);
                ++nbchanges; // we notice one change!!
            }
        }
        // once we did all possible intersection, if *its is not empty we add it:
        if(!(*its).empty()) out_triangle_constraints.push_back(*its);
    }

    return nbchanges;
}

void triangle_manager::print_ccb (Arrangement::Ccb_halfedge_const_circulator circ)
{
  Arrangement::Ccb_halfedge_const_circulator curr = circ;
  std::cout << "(" << curr->source()->point() << ")";
  do {
    Arrangement::Halfedge_const_handle he = curr;
    std::cout << "   [" << he->curve() << "]   "
              << "(" << he->target()->point() << ")";
  } while (++curr != circ);
  std::cout << std::endl;
}


void triangle_manager::print_face (Arrangement::Face_const_handle f)
{
  // Print the outer boundary.
  if (f->is_unbounded())
    std::cout << "Unbounded face. " << std::endl;
  else {
    std::cout << "Outer boundary: ";
    print_ccb (f->outer_ccb());
  }

  // Print the boundary of each of the holes.
  Arrangement::Hole_const_iterator hi;
  int                                 index = 1;
  for (hi = f->holes_begin(); hi != f->holes_end(); ++hi, ++index) {
    std::cout << "    Hole #" << index << ": ";
    print_ccb (*hi);
  }

  // Print the isolated vertices.
  Arrangement::Isolated_vertex_const_iterator iv;
  for (iv = f->isolated_vertices_begin(), index = 1;
       iv != f->isolated_vertices_end(); ++iv, ++index)
  {
    std::cout << "    Isolated vertex #" << index << ": "
              << "(" << iv->point() << ")" << std::endl;
  }
}

void triangle_manager::print_arrangement (const Arrangement& arr)
{
  // Print the arrangement vertices.
  Arrangement::Vertex_const_iterator vit;
  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "(" << vit->point() << ")";
    if (vit->is_isolated())
      std::cout << " - Isolated." << std::endl;
    else
      std::cout << " - degree " << vit->degree() << std::endl;
  }

  // Print the arrangement edges.
  Arrangement::Edge_const_iterator eit;
  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
    std::cout << "[" << eit->curve() << "]" << std::endl;

  // Print the arrangement faces.
  Arrangement::Face_const_iterator fit;
  std::cout << arr.number_of_faces() << " faces:" << std::endl;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    print_face (fit);
}

Point_ar triangle_manager::center_of_mass(Arrangement::Face_const_handle f)
{
    //Point3d c(0,0,0), ori(0,0,0);
    Point_ar c(0,0), ori(0,0);
    if(f->is_fictitious())
    {
        cout << "center_of_mass: The face does not exist." << endl;
        exit(-1);
        return c;
    }

    Arrangement::Ccb_halfedge_const_circulator circ_cur = f->outer_ccb(), circ_end;
    circ_end = circ_cur;
    Lazy_exact_nt_type div = 0;
    do {
        //cout << "************************************" << endl;
        //cout << "init point: " << circ_cur->target()->point() << endl;
        //cout << "proposed value: (" << CGAL::to_double(circ_cur->target()->point().x());
        //cout << "," << CGAL::to_double(circ_cur->target()->point().y()) << ")" << endl;
        /*c = c + (Point3d( CGAL::to_double(circ_cur->target()->point().x()),
                           CGAL::to_double(circ_cur->target()->point().y()),0) -ori);*/
          c = c + (circ_cur->target()->point() - ori);

        div+=1.0;
        //cout << "current center: (" << c.x()/div << "," << c.y()/div << ")" << endl;
        //cout << "************************************" << endl;

    } while (++circ_cur != circ_end);

    if(div!=0) c=ori+(c-ori)/div;//Point3d(c.x()/div,c.y()/div,0);


    //c3d = Point3d(CGAL::to_double(c.x()), CGAL::to_double(c.y()),0);

    return c;
}






void triangle_manager::set_known_feasible_triangulation_cost(double feasible_triangulation_cost)
{
        if(1.1*feasible_triangulation_cost < m_known_feasible_triangulation_cost)
        {
            m_known_feasible_triangulation_cost = 1.1*feasible_triangulation_cost;
        }
}

inline double triangle_manager::get_triangles_total_energy() const
{ // linear search in the number of triangles (sub-optimal)
        if(!m_root) return 0.0; // we do nothing for the null node
        double sum=0;
        std::vector< unsigned int >::const_iterator it_index(m_root->m_triangle_indices.begin()), it_index_e(m_root->m_triangle_indices.end());
        for(;it_index!=it_index_e;++it_index)
        {
            sum+=m_Alltriangles[*it_index].m_potential;
        }

        return sum;
}

inline void triangle_manager::printConstraints(const std::vector< std::set< unsigned int > >& triangle_constraints)
{
    cout << "******************************************" << endl;
    cout << "The region constraints: ";
    std::vector< std::set< unsigned int > >::const_iterator its(triangle_constraints.begin()),ite(triangle_constraints.end());
    for(;its!=ite;++its)
    {
        printIntSet(*its);
    }
    cout << "******************************************" << endl;
}

inline void triangle_manager::triangles_containing_p(const Point3d& P, const std::set< unsigned int >& triangle_indices, std::set< unsigned int >& out_triangle_indices) const
{
    out_triangle_indices.clear();
    std::set< unsigned int >::const_iterator its(triangle_indices.begin()), ite(triangle_indices.end());
    for(;its!=ite;++its)
    { // for each triangle of the input list
        if( m_Alltriangles[*its].is_inside(P) ) out_triangle_indices.insert(*its); // if the point is inside the triangle we add the triangle index to the out list
    }
}

void triangle_manager::generate_Arrangement_2D_from_triangles( const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& triangles,
                                                    Arrangement & arr // out
                                                    )
{
    arr.clear();  	 // clears the arrangement
    ////////////////////////////////////////////////
    std::vector< Segment_ar > segs;
    segs.reserve(3*triangles.size());
    std::set<triangle, std::less<triangle>, std::allocator<triangle> >::const_iterator its(triangles.begin()), ite(triangles.end());
    for(;its!=ite;++its)
    {
        Point_ar    p1(its->m_p1.x().to_double(), its->m_p1.y().to_double()),
                    p2(its->m_p2.x().to_double(), its->m_p2.y().to_double()),
                    p3(its->m_p3.x().to_double(), its->m_p3.y().to_double());

        /*Segment_ar tri[3];
        tri[0] = Segment_ar(p1,p2);
        tri[1] = Segment_ar(p2,p3);
        tri[2] = Segment_ar(p3,p1);*/
        // we insert the triangle in the 2D arrangement:

        segs.push_back(Segment_ar(p1,p2));
        segs.push_back(Segment_ar(p2,p3));
        segs.push_back(Segment_ar(p3,p1));

        //CGAL::insert_x_monotone_curves (arr, &tri[0], &tri[2]);
        //print_arrangement (arr);
    }
    /*
      Insert a range of x-monotone curves into the arrangement (aggregated
      insertion). The inserted x-monotone curves may intersect one another and
      may also intersect the existing arrangement.
     */
    CGAL::insert_x_monotone_curves (arr, segs.begin(), segs.end()); // seems to be correct!! (check done for 4 points)

    assert(arr.is_valid());
    //print_arrangement (arr);
}

void triangle_manager::generate_2D_region_center_of_mass(   Arrangement & arr, // in
                                                            std::vector <Point_ar>& out_centers // out (must be not empty!!)
                                                                    ) const
{
    //print_arrangement (arr);
    if(out_centers.empty())
    {
        std::cout << "WARNING: generate_2D_region_center_of_mass: out_centers is empty." << std::endl;
        return;
    }
    //cout << "out_centers must be fillet with " << out_centers.size() << " points." << endl;

    //out_centers.reserve(arr.number_of_faces()-1);
    ////////////////////////////////////////////////////////////////////////////////////
    std::cout << "In generate_2D_region_center_of_mass:" << std::endl;
    std::cout << arr.number_of_vertices() << " vertices; ";
    std::cout << arr.number_of_edges() << " edges; ";
    std::cout << arr.number_of_faces() << " faces; ";
    std::cout << "whose " << arr.number_of_unbounded_faces() << " unbounded faces." << std::endl;
    ////////////////////////////////////////////////////////////////////////////////////
    Arr_accessor<Arrangement> acc (arr); // constructs an accessor attached to the given arrangement arr.

    unsigned int i=0;
    Arrangement::Face_const_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    {

        if(!fit->is_unbounded() )// && !fit->is_fictitious())
        {
            assert(!fit->is_fictitious());
            //print_face (fit);
            out_centers[i++] = center_of_mass(fit);
            //////////////////////////////////////////////////////////////////////////////
            // to make sure everything is correct:
            // c la merde : fit->halfedge() est private!
            //assert( acc.point_is_in (center_of_mass(fit), (Halfedge_ar)(fit->halfedge()) ) );
            //////////////////////////////////////////////////////////////////////////////
            std::set< unsigned int > out_inside_triangle_indices;
            get_triangles_which_contain_p(out_centers[i-1], out_inside_triangle_indices);
            //////////////////////////////////////////////////////////////////////////////
            if(out_inside_triangle_indices.size()==0)
            {
                cout << "generate_2D_region_center_of_mass: Error: no triangle contains a given center point. Aborded." << endl;
                print_face (fit);
                exit(-1);
            }
        }
        //else cout << "Unbounded face found." << endl;
    }
}

void triangle_manager::generate_region_constraints_from_region_center_of_mass(const std::set<triangle, std::less<triangle>, std::allocator<triangle> >& triangles) const
{
    //////////////////////////////////////////////////////////////////////
    long start_time_s = clock()/CLOCKS_PER_SEC;
    Arrangement arr;
    generate_Arrangement_2D_from_triangles(triangles, arr);
    cout << "The 2D arrangement generation is terminated.";
    print_elapsed_time(start_time_s);
    //////////////////////////////////////////////////////////////////////
    std::vector<Point_ar> out_centers(arr.number_of_faces()-1);
    start_time_s = clock()/CLOCKS_PER_SEC;
    generate_2D_region_center_of_mass(arr, out_centers);
    cout << "The 2D region center of mass generation is terminated.";
    print_elapsed_time(start_time_s);
    //////////////////////////////////////////////////////////////////////
    arr.clear(); // we do not need the arrangement anymore
    //////////////////////////////////////////////////////////////////////
    generate_region_constraints_from_points(out_centers);
}
