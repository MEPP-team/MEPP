///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_POLYHEDRON
#define HEADER_POLYHEDRON

#include <boost/shared_ptr.hpp>

// for typeFuncOpenSave
#include <QString>
class Viewer;
// for typeFuncOpenSave

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "polyhedron_enriched_polyhedron.h"

typedef double number_type;
typedef CGAL::Simple_cartesian<number_type> Enriched_kernel;

typedef Enriched_polyhedron<Enriched_kernel, Enriched_items>	Polyhedron;

//----

typedef Polyhedron::Point_3 Point3d;

typedef Polyhedron::Vertex Vertex;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Vertex_iterator Vertex_iterator;

typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::HalfedgeDS HalfedgeDS;

typedef Polyhedron::Facet Facet;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;

//----

typedef Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

typedef Facet::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

//----

typedef CGAL::Vector_3<Enriched_kernel> Vector;

//----

typedef boost::shared_ptr<Polyhedron> PolyhedronPtr;

typedef int (*typeFuncOpenSave)(PolyhedronPtr, QString, Viewer*);

#endif