///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Triangle_reading_generating_saving_H
#define Triangle_reading_generating_saving_H

#include "../triangle_manager/Triangle.h"

static bool is_second_order_neighbor(const Delaunay& dt, const vh v1, const vh v2)// for finite points
{
    if(dt.is_infinite(v1) || dt.is_infinite(v2)) return false;

    if(v1==v2) return true; // 0 order
    v_cir vc = dt.incident_vertices(v1); // see its one-ring neighbors
    v_cir ve = vc;
    CGAL_For_all(vc,ve)
    {
        if(dt.is_infinite(vc)) continue;

        if(vc==v2) return true; // 1 order

        v_cir vc2 = dt.incident_vertices(vc); // see its one-ring neighbors
        v_cir ve2 = vc2;
        CGAL_For_all(vc2,ve2)
        {
            if(dt.is_infinite(vc2)) continue;

            if(vc2==v2) return true;  // 2 order
        }
    }

    return false;
};



static void build_2D_delaunay_triangulation(const std::vector< Point3d >& cloud_point,
                                            Delaunay& dt,
                                            std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_authorized_triangle)
{
    //dt.clear();
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "Building Delaunay triangulation... ";
    std::vector< Point3d >::const_iterator its(cloud_point.begin()), ite(cloud_point.end());
    for(;its!=ite;++its)
        dt.insert(Point_dt(its->x(), its->y())); // we incrementally build the Delaunay triangulation
    cout << "Done. ";
    assert(dt.is_valid()); // we check if the delaunay triangulation is valid
    cout << "The triangulation is valid regarding the Delaunay property. ";
    //cout << dt.number_of_faces() << " triangles." << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    f_it pFacet = dt.finite_faces_begin(), pFacetEnd = dt.finite_faces_end(); // to go through finite triangles
    double maxs=0;
    for(;pFacet!=pFacetEnd;++pFacet)
    {
        Triangle_dt T = dt.triangle(pFacet);

        if( isValidTriangle(T[0], T[1], T[2]) )
        { // to make sure that the triangle is not flat
            triangle tri(T[0], T[1], T[2]);
            if(tri.m_potential>maxs) maxs=tri.m_potential;
            set_of_authorized_triangle.insert(tri);
        }
    }
    cout << "build_2D_delaunay_triangulation: max triangle potential = " << maxs << endl;
}

static void build_2D_delaunay_triangulation_2_rings(Delaunay& dt,
                                                    std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_authorized_triangle
                                                    //, Number_type_delaunay& totalArea
                                                    )
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    assert(dt.is_valid()); // we check if the delaunay triangulation is valid
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //totalArea = 0.0;
    // here we will compute the total area of the convex domain triangulated by Delaunay:
    /*f_it pFacet = dt.finite_faces_begin(), pFacetEnd = dt.finite_faces_end();
    for(;pFacet!=pFacetEnd;++pFacet)
    {
        Triangle_dt T = dt.triangle(pFacet);
        //totalArea += T.area();
    }*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // now we are looking for a particular subset of triangles (which intersect locally)
    //double TH = 1e8; // it is useful to avoid vertical triangles at the border
    double maxs=0;
    //set_of_authorized_triangle.clear();
    v_it pVertex = dt.finite_vertices_begin(), pVertexEnd = dt.finite_vertices_end();
    for(;pVertex!=pVertexEnd;++pVertex) // for each finite vertex of the Delaunay triangulation...
    {
        Point_dt p2; // 2D delaunay points
        p2 = pVertex->point(); // we are looking for triangle which contain pVertex->point() as vertex
        //cout << "new vertex of the Delaunay triangulation..." << endl;
        //cout << "point p2 = (" << p2.x() << "," << p2.y() << ")" << endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        v_cir vc = dt.incident_vertices(pVertex); // see its one-ring neighbors
        v_cir ve = vc;
        CGAL_For_all(vc,ve)
        {   // we know that the points here are first order neighbor so second order neighbor (and !=pVertex)

            // first of all, we must add all the possible triangles in the one-ring:
            if(dt.is_infinite(vc)) continue; // we do not touch to infinite vertices!!

            Point_dt p1 = vc->point(), p3;
            //cout << "point p1 = (" << p1.x() << "," << p1.y() << ")" << endl;
            //////////////////////////////////////////////////////////////////////////////////
            // the goal here is to find out all the possible triangles in the one-ring
            v_cir vc2 = dt.incident_vertices(pVertex); // see its one-ring neighbors again
            v_cir ve2 = vc2;
            CGAL_For_all(vc2,ve2)
            {
                if(vc2==vc || dt.is_infinite(vc2)) continue;

                p3 = vc2->point(); // another point than vc->point();
                //cout << "point p3 = (" << p3.x() << "," << p3.y() << ")" << endl;
                // add a triangle in counterclockwise order:
                //Point3d p3d1(p1.x().to_double(),p1.y().to_double(),0), p3d2(p2.x().to_double(),p2.y().to_double(),0), p3d3(p3.x().to_double(),p3.y().to_double(),0);
                Point3d p3d1(CGAL::to_double(p1.x()),CGAL::to_double(p1.y()),0), p3d2(CGAL::to_double(p2.x()),CGAL::to_double(p2.y()),0), p3d3(CGAL::to_double(p3.x()),CGAL::to_double(p3.y()),0);
                //Point3d p3d1(p1.x(),p1.y(),0), p3d2(p2.x(),p2.y(),0), p3d3(p3.x(),p3.y(),0);

                if( isValidTriangle(p3d1, p3d2, p3d3) )
                {
                    assert((p1.x()!=p2.x()) || (p1.x()!=p3.x()));
                    assert((p1.y()!=p2.y()) || (p1.y()!=p3.y()));
                    triangle tri(p3d1, p3d2, p3d3);
                    set_of_authorized_triangle.insert(tri);
                    if(tri.m_potential>maxs) maxs=tri.m_potential;
                }
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // secondly, we must add the possible triangles which share a vertex on the one-ring and another on the second-ring
            vc2 = dt.incident_vertices(vc); // see its one-ring neighbors
            ve2 = vc2;
            CGAL_For_all(vc2,ve2)
            { // we know that the points here are second order neighbor
                if(dt.is_infinite(vc2) || vc2==pVertex) continue; // we just skip this case

                //if(is_second_order_neighbor(dt, vc, vc2)) // we are sure that vc!=vc2
                //{
                // that part is not relevant (the triangles will be added by the next part
                   /* p1 = vc->point();
                    p3 = vc2->point();
                    // add a triangle in counterclockwise order:
                    Point3d p3d1(p1.x(),p1.y(),0), p3d2(p2.x(),p2.y(),0), p3d3(p3.x(),p3.y(),0);
                    if(!CGAL::collinear(p1, p2, p3) && !CGAL::collinear(p3d1, p3d2, p3d3) && (triangleShapePotential(p3d1,p3d2,p3d3)<TH))
                    {
                        set_of_authorized_triangle.insert(triangle(p3d1, p3d2, p3d3));
                    }
                    */
                //}


                // finally, we must add the possible triangles which share two vertices on the second-ring
                v_cir vc3 = dt.incident_vertices(vc); // see its one-ring neighbors
                v_cir ve3 = vc3;
                CGAL_For_all(vc3,ve3)
                {
                    if(dt.is_infinite(vc3) || vc3==pVertex || vc3==vc2) continue;

                    //if(is_second_order_neighbor(dt, vc2, vc3))
                    //{
                        p1 = vc2->point();
                        p3 = vc3->point();
                        // add a triangle in counterclockwise order:
                        //Point3d p3d1(p1.x().to_double(),p1.y().to_double(),0), p3d2(p2.x().to_double(),p2.y().to_double(),0), p3d3(p3.x().to_double(),p3.y().to_double(),0);
                        Point3d p3d1(CGAL::to_double(p1.x()),CGAL::to_double(p1.y()),0), p3d2(CGAL::to_double(p2.x()),CGAL::to_double(p2.y()),0), p3d3(CGAL::to_double(p3.x()),CGAL::to_double(p3.y()),0);
                        //Point3d p3d1(p1.x(),p1.y(),0), p3d2(p2.x(),p2.y(),0), p3d3(p3.x(),p3.y(),0);

                        if( isValidTriangle(p3d1, p3d2, p3d3) )
                        {
                            assert((p1.x()!=p2.x()) || (p1.x()!=p3.x()));
                            assert((p1.y()!=p2.y()) || (p1.y()!=p3.y()));
                            triangle tri(p3d1, p3d2, p3d3);
                            set_of_authorized_triangle.insert(tri);
                            if(tri.m_potential>maxs) maxs=tri.m_potential;
                        }
                    //}
                 }

             }


        }

    }
    cout << "build_2D_delaunay_triangulation_2_rings: max triangle potential = " << maxs << endl;
}

static void generate_triangles_2D_delaunay_triangulation_refinement_sheme(  std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_authorized_triangle,
                                                                            std::vector< Point3d >& cloud_point)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //assert(dt.is_valid()); // we check if the delaunay triangulation is valid
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::vector< Point3d > new_points;
    //f_it pFacet = dt.finite_faces_begin(), pFacetEnd = dt.finite_faces_end(); // iterator over finite faces
    //for(;pFacet!=pFacetEnd;++pFacet)
    std::set<triangle, std::less<triangle>, std::allocator<triangle> > inter_triangles;
    std::set<triangle, std::less<triangle>, std::allocator<triangle> >::const_iterator its(set_of_authorized_triangle.begin()), ite(set_of_authorized_triangle.end());
    for(;its!=ite;++its)
    {
        bool add_p12=false, add_p23=false, add_p31=false;
        //Triangle_dt T = dt.triangle(pFacet);
        // do something with T...
        //triangle tri(T[0], T[1], T[2]);
        triangle tri(*its);

        Point3d p12 = tri.middle_side1(), p23 = tri.middle_side2(), p31 = tri.middle_side3();
        Point3d p3d1(tri.m_p1.x().to_double(),tri.m_p1.y().to_double(),0), p3d2(tri.m_p2.x().to_double(),tri.m_p2.y().to_double(),0), p3d3(tri.m_p3.x().to_double(),tri.m_p3.y().to_double(),0);
        ///////////////////////////////////////////////////////////////
        // 10 TRIANGLES REFINEMENT SCHEME
        // insert the 4 basic triangles
        //unsigned int nbt=0;
        if( isValidTriangle(p12,p23,p31) &&
            isValidTriangle(p3d1,p12,p31) &&
            isValidTriangle(p12,p3d2,p23) &&
            isValidTriangle(p23,p3d3,p31) )
            {
                inter_triangles.insert(triangle(p12,p23,p31));
                inter_triangles.insert(triangle(p3d1,p12,p31));
                inter_triangles.insert(triangle(p12,p3d2,p23));
                inter_triangles.insert(triangle(p23,p3d3,p31));
                //nbt+=4;
                add_p12=true;
                add_p23=true;
                add_p31=true;
            }

        // insert the 6 bisector triangles
        if( isValidTriangle(p3d1,p3d2,p23) &&
            isValidTriangle(p23,p3d3,p3d1))
            {
                inter_triangles.insert(triangle(p3d1,p3d2,p23));
                inter_triangles.insert(triangle(p23,p3d3,p3d1));
                //nbt+=2;
                add_p23=true;
            }

        if( isValidTriangle(p3d2,p3d3,p31) &&
            isValidTriangle(p31,p3d1,p3d2))
            {
                inter_triangles.insert(triangle(p3d2,p3d3,p31));
                inter_triangles.insert(triangle(p31,p3d1,p3d2));
                //nbt+=2;
                add_p31=true;
            }

        if( isValidTriangle(p3d3,p3d1,p12) &&
            isValidTriangle(p12,p3d2,p3d3))
            {
                inter_triangles.insert(triangle(p3d3,p3d1,p12));
                inter_triangles.insert(triangle(p12,p3d2,p3d3));
                //nbt+=2;
                add_p12=true;
            }

        //if(nbt<10) cout  << "nb triangles added for one step: " << nbt << endl;
        /*Point3d bar = tri.barycenter();
        if( isValidTriangle(p3d1,p3d2,bar) &&
            isValidTriangle(p3d2,p3d3,bar) &&
            isValidTriangle(p3d3,p3d1,bar) )
            {
                inter_triangles.insert(triangle(p3d1,p3d2,bar));
                inter_triangles.insert(triangle(p3d2,p3d3,bar));
                inter_triangles.insert(triangle(p3d3,p3d1,bar));

                cloud_point.push_back(bar);
            }
            */
        ///////////////////////////////////////////////////////////////
        // insert the point in the vector to build another Delaunay triangulation later
        if(add_p12) cloud_point.push_back(p12);
        if(add_p23) cloud_point.push_back(p23);
        if(add_p31) cloud_point.push_back(p31);
    }
    set_of_authorized_triangle.insert(inter_triangles.begin(), inter_triangles.end());
}


static void generate_triangles_2D_barycenter_refinement_sheme(  std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_authorized_triangle,
                                                                std::vector< Point3d >& cloud_point)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //assert(dt.is_valid()); // we check if the delaunay triangulation is valid
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::vector< Point3d > new_points;
    //f_it pFacet = dt.finite_faces_begin(), pFacetEnd = dt.finite_faces_end(); // iterator over finite faces
    //for(;pFacet!=pFacetEnd;++pFacet)
    std::set<triangle, std::less<triangle>, std::allocator<triangle> > inter_triangles;

    std::set<triangle, std::less<triangle>, std::allocator<triangle> >::const_iterator its(set_of_authorized_triangle.begin()), ite(set_of_authorized_triangle.end());
    for(;its!=ite;++its)
    {
        triangle tri(*its);
        Point3d bar = tri.barycenter();
        Point3d p3d1(tri.m_p1.x().to_double(),tri.m_p1.y().to_double(),0), p3d2(tri.m_p2.x().to_double(),tri.m_p2.y().to_double(),0), p3d3(tri.m_p3.x().to_double(),tri.m_p3.y().to_double(),0);
        ///////////////////////////////////////////////////////////////
        if( isValidTriangle(p3d1,p3d2,bar) &&
            isValidTriangle(p3d2,p3d3,bar) &&
            isValidTriangle(p3d3,p3d1,bar) )
            {
                inter_triangles.insert(triangle(p3d1,p3d2,bar));
                inter_triangles.insert(triangle(p3d2,p3d3,bar));
                inter_triangles.insert(triangle(p3d3,p3d1,bar));

                cloud_point.push_back(bar);
            }
    }
    set_of_authorized_triangle.insert(inter_triangles.begin(), inter_triangles.end());
}


static void all_possible_non_flat_triangles(const std::vector< Point3d >& cloud_point,
                                            std::set<triangle, std::less<triangle>, std::allocator<triangle> >& set_of_authorized_triangle)
{
    cout << "all_possible_non_flat_triangles: IN" << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    unsigned long nbp = cloud_point.size(), nb_point_alignement=0, nb_added_tri=0;
    if(nbp>100) cout << "WARNING: all_possible_non_flat_triangles: intractable number of points." << endl;
    for(unsigned long i=0; i<(nbp-2); ++i)
    {
        for(unsigned long j=i+1; j<(nbp-1); ++j)
        {
            for(unsigned long k=j+1; k<nbp; ++k)
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // we cannot generate a flat triangle, since its region without edges does not exist!!! So if the 3 points are aligned, we just skiped
                if( !isValidTriangle(cloud_point[i], cloud_point[j], cloud_point[k]) )
                {
                    nb_point_alignement++; // to monitor what's happened
                    continue;
                }
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                set_of_authorized_triangle.insert(triangle(cloud_point[i], cloud_point[j], cloud_point[k])); // counter-clockwise oriented (I have checked it on several examples)
                nb_added_tri++;
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "all_possible_non_flat_triangles: results checking... ";
    cout << "We must have added at maximum " << nbp*(nbp-1)*(nbp-2)/6 << " triangles and we have added ";
    cout << nb_added_tri << " triangles and " << nb_point_alignement << " triangles have not been added because they were flat." << endl;
    assert(nbp*(nbp-1)*(nbp-2)==6*(nb_point_alignement+nb_added_tri));
    cout << "all_possible_non_flat_triangles: OUT" << endl;
}


static void read_2Dtriangle_file(const char* filename, std::set<triangle, std::less<triangle>, std::allocator<triangle> >& mytriangles)
{
    mytriangles.clear();
    ifstream flux(filename, ios::in);
    // open the file:
    if(!flux)
    {
        cout << "read_2Dtriangle_file: Cannot open the file: "<< filename << "!!!" << endl;
        exit(-1);
    }
    string line;
    while ( getline(flux, line) )
    { // for each triangle
        istringstream iss(line);
        double d1, d2, d3, d4, d5, d6, d7, d8;
        iss >> d1; iss >> d2; iss >> d3; iss >> d4; iss >> d5; iss >> d6; iss >> d7; iss >> d8;
        //cout << d1 << ", " << d2 << ", " << d3 << ", " << d4 << ", " << d5 << ", " << d6 << ", " << d7 << endl;
        Point3d p3d1(d1,d2,0), p3d2(d3,d4,0), p3d3(d5,d6,0);
        mytriangles.insert(triangle(p3d1,p3d2,p3d3));
    }
    flux.close();
}

static void read_triangle_indices(const char* filename, vector< unsigned int >& triangle_indices)
{
    triangle_indices.clear();
    ifstream flux(filename, ios::in);
    // open the file:
    if(!flux)
    {
        cout << "read_triangle_indices: Cannot open the file: "<< filename << "!!!" << endl;
        //exit(-1); // we do not exit here
        return;
    }
    string line;
    while ( getline(flux, line) )
    { // for each triangle
        istringstream iss(line);
        unsigned int i;
        iss >> i;
        triangle_indices.push_back(i);
    }
    flux.close();
}

#endif // Triangle_reading_generating_saving_H
