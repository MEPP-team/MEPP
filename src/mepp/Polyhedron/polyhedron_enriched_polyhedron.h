#ifndef HEADER_ENRICHED_POLYHEDRON
#define HEADER_ENRICHED_POLYHEDRON

#include <mepp_config.h>

#include "polyhedron_shared_items.h"

#include <polyhedron_enrichment_include.h>

template <class Refs, class T, class P, class Norm, class Plane>
class Enriched_facet :
	/*************** HERITAGE FACETTE ***************/
	#include <polyhedron_enrichment_facet.h>
	/*************************************************/
	virtual public MEPP_Common_Facet<Refs, T, Norm>
{
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Enriched_halfedge :
	/*************** HERITAGE HALFEDGE ***************/
	#include <polyhedron_enrichment_halfedge.h>
	/*************************************************/
	virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
};

// a refined vertex with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Enriched_vertex :
	/*************** HERITAGE VERTEX ***************/
	#include <polyhedron_enrichment_vertex.h>
	/*************************************************/
	virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:
		Enriched_vertex()
		{
		}

		// Le constructeur PAS par défaut appelle les constructeurs PAR DEFAUT des classes ancetres
		// en plus de celles appelles appellées explicitement.
		// On a besoin d'appeller explicitement le constructeur(pt) de base pour la création de polyhedre
		Enriched_vertex(const P& pt) : MEPP_Common_Vertex<Refs, T, P, Norm>(pt)
		{
			this->point() = pt;
		}

		// La creation du polyhedre implique un appel au constructeur par copie,
		// qui appelle le constructeur par défaut de base et ne copie pas le point.
		// Il faut donc avoir un constructeur par copie explicite qui s'occupe du point.
		Enriched_vertex(const Enriched_vertex& v)
		{
			this->point() = v.point();
		}
};

struct Enriched_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template <class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_vertex<Refs,
                          CGAL::Tag_true,
                          Point,
                          Normal> Vertex;
    };

    // wrap face
    template <class Refs, class Traits>
    struct Face_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
		typedef typename Traits::Plane_3 Plane;
        typedef Enriched_facet<Refs,
                         CGAL::Tag_true,
                         Point,
                         Normal,
						 Plane> Face;
    };

    // wrap halfedge
    template <class Refs, class Traits>
    struct Halfedge_wrapper
    {
        typedef typename Traits::Vector_3 Normal;
        typedef Enriched_halfedge<Refs,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            Normal> Halfedge;
    };
};

template <class kernel, class items>
class Enriched_polyhedron :
	/*************** HERITAGE POLYHEDRON ***************/
	#include <polyhedron_enrichment_polyhedron.h>
	/*************************************************/
	virtual public MEPP_Common_Polyhedron<kernel,items>
{
};

#endif
