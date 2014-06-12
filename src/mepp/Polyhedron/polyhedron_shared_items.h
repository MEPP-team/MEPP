#ifndef HEADER_POLYHEDRON_SHARED_ITEMS
#define HEADER_POLYHEDRON_SHARED_ITEMS

#ifdef _MSC_VER
	#if !defined(NOMINMAX)
		#define NOMINMAX
	#endif
#endif

#ifdef _MSC_VER
#include <windows.h>
#endif

//#include <GL/glew.h>
#ifdef __APPLE__
#  include <OpenGL/glu.h>
#else
#  include <GL/glu.h>
#endif

#include <QGLViewer/vec.h>
#include <QGLViewer/quaternion.h>

#include <fstream>
#include <list>

//#include <CGAL/IO/Polyhedron_VRML_2_ostream.h> // for vrml writing

// For OFF Loading
#include "Polyhedron_OFF_CGALImporter.h"

// For OBJ/SMF Loading
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#include "Polyhedron_OBJ_CGALImporter.h"
#include "Polyhedron_SMF_CGALImporter.h"
#pragma warning(pop)
#else
#include "Polyhedron_OBJ_CGALImporter.h"
#include "Polyhedron_SMF_CGALImporter.h"
#endif

// For PLY Loading
#include "Polyhedron_PLY_CGALImporter.h"

// For X3D Loading
#include "Polyhedron_X3D_CGALImporter.h"
#include "X3D_old.h"

#ifdef WITH_ASSIMP
	// For ASSIMP Loading
	#include "Polyhedron_ASSIMP_CGALImporter.h"
#endif

// For Polyhedron copy
#include "Polyhedron_Copy.h"

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

using namespace std;
using namespace CGAL;

struct Texture
{
	QString m_name;
	QImage m_data;
	GLuint m_id;
};

// compute facet normal
struct Facet_normal // (functor)
{
	template <class Facet>
	void operator()(Facet& f)
	{
		typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
		typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
		do
		{
			typename Facet::Normal_3 normal = CGAL::cross_product(
				h->next()->vertex()->point() - h->vertex()->point(),
				h->next()->next()->vertex()->point() - h->next()->vertex()->point());
			double sqnorm = to_double(normal * normal);
			if (sqnorm != 0)
				normal = normal / (float)std::sqrt(sqnorm);
			sum = sum + normal;
		}
		while(++h != f.facet_begin());
		double sqnorm = to_double(sum * sum);
		if (sqnorm != 0.0)
			f.normal() = sum / std::sqrt(sqnorm);
		else
		{
			f.normal() = CGAL::NULL_VECTOR;
			//TRACE("degenerate face\n");
		}
	}
};
// from Hichem
/*struct Facet_normal
{
	template <class Facet>
	void operator()(Facet& f)
	{
		typedef typename Facet::Normal_3		Vector_3;
		Vector_3								facet_normal;
		typename Facet::Halfedge_around_facet_const_circulator h = f.facet_begin();
		do
		{
			facet_normal = CGAL::cross_product(h->next()->vertex()->point() - h->vertex()->point(),
				h->next()->next()->vertex()->point() - h->next()->vertex()->point());
		}
		while (facet_normal != CGAL::NULL_VECTOR && ++h	!= f.facet_begin());
		if(facet_normal != CGAL::NULL_VECTOR)
		{
			f.normal() = facet_normal;
		}
		else // All consecutive facet edges are collinear --> degenerate (0 area) facet
		{
			std::cerr <<std::endl << "Degenerate polyhedron facet" << std::endl << std::flush;
			f.normal() = CGAL::NULL_VECTOR;
			assert(facet_normal != CGAL::NULL_VECTOR);			
		}
	}
};*/

// compute vertex normal
struct Vertex_normal // (functor)
{
    template <class Vertex>
    void operator()(Vertex& v)
    {
        typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
        typename Vertex::Halfedge_around_vertex_const_circulator pHalfedge = v.vertex_begin();
        typename Vertex::Halfedge_around_vertex_const_circulator begin = pHalfedge;
        CGAL_For_all(pHalfedge,begin)
          if (!pHalfedge->is_border())
            normal = normal + pHalfedge->facet()->normal();
        double sqnorm = to_double(normal * normal);
        if (sqnorm != 0.0f)
          v.normal() = normal / (float)std::sqrt(sqnorm);
        else
          v.normal() = CGAL::NULL_VECTOR;
    }
};

#if (0)
// compute halfedge normal : ajout Céline
template <class kernel>
struct Halfedge_normal // (functor)
{
	typedef typename kernel::FT FT;

    template <class Halfedge>
    void operator()(Halfedge& h)
    {
        typename Halfedge::Normal_3 normal = CGAL::NULL_VECTOR;

        if (!h.is_border())
            normal = normal + (h.facet()->normal());
		if (!(h.opposite()->is_border()))
            normal = normal + (h.opposite()->facet()->normal());

        FT sqnorm = normal * normal;
        if (sqnorm != 0.0f)
          h.normal() = normal / (FT)std::sqrt(sqnorm);
        else
          h.normal() = CGAL::NULL_VECTOR;
    }
};
#endif

template <class Refs, class T, class Norm>
class MEPP_Common_Facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
	protected:
		// tag
		int m_tag;

		// normal
		Norm m_normal;

		// color
		float m_color[3];

	public:
		// life cycle
		MEPP_Common_Facet()
		{
			color(0.5f, 0.5f, 0.5f);
		}

		// tag
		const int& tag() const { return m_tag;  }
		int& tag() { return m_tag;  }
		void tag(const int& t)  { m_tag = t; }

		// normal
		typedef Norm Normal_3;
		Normal_3& normal() { return m_normal; }
		const Normal_3& normal() const { return m_normal; }

		// color
		float color(int index) { return m_color[index]; }
		void color(float r, float g, float b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MEPP_Common_Halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
	protected:
		// tag
		int m_tag;

#if (0)
		// normal
		Norm m_normal;  // AJOUT Céline : halfedge normal (sum of 2 incident facet normals)
#endif

		// option for edge superimposing
		bool m_control_edge;

		// texture coordinates : AJOUT Laurent Chevalier
		float m_texture_coordinates[2];

	public:
		// life cycle
		MEPP_Common_Halfedge()
		{
			m_control_edge = true;

			// texture coordinates : AJOUT Laurent Chevalier
			texture_coordinates(0.0f, 0.0f);
		}

		// tag
		const int& tag() const { return m_tag;  }
		int& tag() { return m_tag;  }
		void tag(const int& t)  { m_tag = t; }

#if (0)
		// normal : AJOUT Céline
		typedef Norm Normal_3;
		Normal_3& normal() { return m_normal; }
		const Normal_3& normal() const { return m_normal; }
#endif

		// texture coordinates : AJOUT Laurent Chevalier
		float texture_coordinates(int index) { return m_texture_coordinates[index]; }
		void texture_coordinates(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }

		// control edge
		bool& control_edge()  { return m_control_edge; }
		const bool& control_edge()  const { return m_control_edge; }
		void control_edge(const bool& flag) { m_control_edge  = flag; }
};

template <class Refs, class T, class P, class Norm>
class MEPP_Common_Vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
	protected:
		// tag
		int m_tag;

		// normal
		Norm m_normal;

		// color
		float m_color[3];

		// texture coordinates
		float m_texture_coordinates[2];

	public:
		// life cycle
		MEPP_Common_Vertex()
		{
			color(0.5f, 0.5f, 0.5f);
			texture_coordinates(0.0f, 0.0f);
		}
		// repeat mandatory constructors
		MEPP_Common_Vertex(const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
		{
			color(0.5f, 0.5f, 0.5f);
			texture_coordinates(0.0f, 0.0f);
		}

		// color
		float color(int index) { return m_color[index]; }
		void color(float r, float g, float b) { m_color[0] = r; m_color[1] = g; m_color[2] = b; }

		// texture coordinates
		float texture_coordinates(int index) { return m_texture_coordinates[index]; }
		void texture_coordinates(float u, float v) { m_texture_coordinates[0] = u; m_texture_coordinates[1] = v; }

		// normal
		typedef Norm Normal_3;
  		//typedef Norm Vector;

		Normal_3& normal() { return m_normal; }
		const Normal_3& normal() const { return m_normal; }

		// tag
		int& tag() {  return m_tag; }
		const int& tag() const {  return m_tag; }
		void tag(const int& t)  { m_tag = t; }
};

template <class kernel, class items>
class MEPP_Common_Polyhedron : public CGAL::Polyhedron_3<kernel,items>
{
	public:
		typedef typename kernel::FT FT;
		typedef typename kernel::Point_3 Point;
		typedef typename kernel::Vector_3 Vector;
		typedef typename kernel::Iso_cuboid_3 Iso_cuboid;
		typedef typename MEPP_Common_Polyhedron::Facet_handle Facet_handle;
		typedef typename MEPP_Common_Polyhedron::Vertex_handle Vertex_handle;
		typedef typename MEPP_Common_Polyhedron::Halfedge_handle Halfedge_handle;
		typedef typename MEPP_Common_Polyhedron::Facet_iterator Facet_iterator;
		typedef typename MEPP_Common_Polyhedron::Vertex_iterator Vertex_iterator;
		typedef typename MEPP_Common_Polyhedron::Halfedge_iterator Halfedge_iterator;
		typedef typename MEPP_Common_Polyhedron::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
		typedef typename MEPP_Common_Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
		typedef typename MEPP_Common_Polyhedron::Point_iterator Point_iterator;
		typedef typename MEPP_Common_Polyhedron::Edge_iterator Edge_iterator;
		typedef typename MEPP_Common_Polyhedron::HalfedgeDS HalfedgeDS;
		typedef typename HalfedgeDS::Face Facet;
		typedef typename Facet::Normal_3 Normal;
		typedef Aff_transformation_3<kernel> Affine_transformation;

	protected:
		// tag : AJOUT Céline
		int m_tag;

		// bounding box
		Iso_cuboid m_bbox;

		bool m_pure_quad;
		bool m_pure_triangle;

		bool m_has_color;
		bool m_has_texture_coordinates;

		// MT
		unsigned int m_nb_components;
		unsigned int m_nb_boundaries;
		
		GLuint id_cube;

		vector<Texture> m_texture_array;

		bool m_use_halfedge_texture_coordinates;

	public:
		// life cycle
		MEPP_Common_Polyhedron()
		{
			m_pure_quad = false;
			m_pure_triangle = false;
			m_has_color = false;
			m_has_texture_coordinates = false;
			m_use_halfedge_texture_coordinates = false;

			// MT
			m_nb_components = 0;
			m_nb_boundaries = 0;

			id_cube = 0;
			m_texture_array.clear();
		}

		virtual ~MEPP_Common_Polyhedron()
		{
			// todo: commented because there is a problem when we close a child and that there is another child in rotation

			/*if (id_cube)
				glDeleteLists(id_cube, 1);*/

			for (unsigned int i = 0; i < m_texture_array.size(); i++)
			{
				if ( m_texture_array[i].m_id > 0 )
					glDeleteTextures(1, &m_texture_array[i].m_id);
			}
		}

		// MT
		unsigned int nb_components() {  return m_nb_components; }
		unsigned int nb_boundaries() {  return m_nb_boundaries; }
		string pName;
		bool pShow;
		qglviewer::Vec pInitialCameraPosition;
		qglviewer::Quaternion pInitialCameraOrientation;

		// tag : AJOUT Céline
		int& tag() {  return m_tag; }
		const int& tag() const {  return m_tag; }
		void tag(const int& t)  { m_tag = t; }

		// type
		bool is_pure_triangle() { return m_pure_triangle; }
		bool is_pure_quad() { return m_pure_quad; }

		// normals (per facet, then per vertex)
		void compute_normals_per_facet()
		{
			std::for_each(this->facets_begin(),this->facets_end(),Facet_normal());
		}
		void compute_normals_per_vertex()
		{
			std::for_each(this->vertices_begin(),this->vertices_end(),Vertex_normal());
		}
#if (0)
		// ajout Céline :
	    void compute_normals_per_halfedge()
	    {
			std::for_each(this->halfedges_begin(), this->halfedges_end(), Halfedge_normal<kernel>());
		}
#endif

		void compute_normals()
		{
			compute_normals_per_facet();
			compute_normals_per_vertex();
#if (0)
			compute_normals_per_halfedge();  // ajout Céline
#endif
		}

		// bounding box
		Iso_cuboid& bbox() { return m_bbox; }
		const Iso_cuboid bbox() const { return m_bbox; }

		// compute bounding box
		void compute_bounding_box()
		{
			if (this->size_of_vertices() == 0)
			{
				return;
			}

			FT xmin,xmax,ymin,ymax,zmin,zmax;
			Vertex_iterator pVertex = this->vertices_begin();
			xmin = xmax = pVertex->point().x();
			ymin = ymax = pVertex->point().y();
			zmin = zmax = pVertex->point().z();

			for (;pVertex !=  this->vertices_end();pVertex++)
			{
				const Point& p = pVertex->point();

				xmin = std::min(xmin,p.x());
				ymin = std::min(ymin,p.y());
				zmin = std::min(zmin,p.z());

				xmax = std::max(xmax,p.x());
				ymax = std::max(ymax,p.y());
				zmax = std::max(zmax,p.z());
			}

			m_bbox = Iso_cuboid(xmin,ymin,zmin,xmax,ymax,zmax);
		}

		// bounding box
		FT xmin() { return m_bbox.xmin(); }
		FT xmax() { return m_bbox.xmax(); }
		FT ymin() { return m_bbox.ymin(); }
		FT ymax() { return m_bbox.ymax(); }
		FT zmin() { return m_bbox.zmin(); }
		FT zmax() { return m_bbox.zmax(); }

		// copy bounding box
		void copy_bounding_box(MEPP_Common_Polyhedron<kernel,items> *pMesh)
		{
			m_bbox = pMesh->bbox();
		}

		// degree of a face
		static size_t degree(Facet_handle pFace)
		{
			return CGAL::circulator_size(pFace->facet_begin());
		}

		// valence of a vertex
		static unsigned int valence(Vertex_handle pVertex)
		{
			return CGAL::circulator_size(pVertex->vertex_begin());
		}

		// check wether a vertex is on a boundary or not
		static bool is_border(Vertex_handle pVertex)
		{
			Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
			if (pHalfEdge == NULL) // isolated vertex
			return true;
			Halfedge_around_vertex_circulator d = pHalfEdge;
			CGAL_For_all(pHalfEdge,d)
			if (pHalfEdge->is_border())
				return true;
			return false;
		}

		// get any border halfedge attached to a vertex
		Halfedge_handle get_border_halfedge(Vertex_handle pVertex)
		{
			Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
			Halfedge_around_vertex_circulator d = pHalfEdge;
			CGAL_For_all(pHalfEdge,d)
			if (pHalfEdge->is_border())
				return pHalfEdge;
			return NULL;
		}

		// tag all halfedges
		void tag_halfedges(const int tag)
		{
			for (Halfedge_iterator pHalfedge = this->halfedges_begin();
					pHalfedge != this->halfedges_end(); pHalfedge++)
			pHalfedge->tag(tag);
		}

		// tag all facets
		void tag_facets(const int tag)
		{
			for (Facet_iterator pFacet = this->facets_begin();pFacet != this->facets_end(); pFacet++)
			pFacet->tag(tag);
		}

		// set index for all vertices
		void set_index_vertices()
		{
			int index = 0;
			for (Vertex_iterator pVertex = this->vertices_begin();
				pVertex != this->vertices_end();
				pVertex++)
			pVertex->tag(index++);
		}

		// is pure degree ?
		bool is_pure_degree(unsigned int d)
		{
			for (Facet_iterator pFace  = this->facets_begin();
				pFace != this->facets_end();
				pFace++)
			if (degree(pFace) != d)
				return false;
			return true;
		}

		// compute type
		void compute_type()
		{
			m_pure_quad = is_pure_degree(4);
			m_pure_triangle = is_pure_degree(3);
		}

		// compute facet center
		void compute_facet_center(Facet_handle pFace, Point& center)
		{
			Halfedge_around_facet_circulator pHalfEdge = pFace->facet_begin();
			Halfedge_around_facet_circulator end = pHalfEdge;
			Vector vec(0.0,0.0,0.0);
			int degree = 0;
			CGAL_For_all(pHalfEdge,end)
			{
				vec = vec + (pHalfEdge->vertex()->point()-CGAL::ORIGIN);
				degree++;
			}
			center = CGAL::ORIGIN + (vec/(FT)degree);
		}

		// compute average edge length around a vertex
		FT average_edge_length_around(Vertex_handle pVertex)
		{
			FT sum = 0.0;
			Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
			Halfedge_around_vertex_circulator end = pHalfEdge;
			Vector vec(0.0,0.0,0.0);
			int degree = 0;
			CGAL_For_all(pHalfEdge,end)
			{
				Vector vec = pHalfEdge->vertex()->point()-
							pHalfEdge->opposite()->vertex()->point();
				sum += std::sqrt(to_double(vec*vec));
				degree++;
			}
			return sum / (FT) degree;
		}

		virtual void gl_draw(bool smooth_shading, bool use_normals, bool use_vertex_color, bool use_face_color, bool use_texture)
		{
			// texture
			if (use_texture && has_texture())
			{
				int r = 0;
				//r = rand() % m_texture_array.size();

				GLuint id = m_texture_array[r].m_id;
				QString name = m_texture_array[r].m_name;

				glBindTexture(GL_TEXTURE_2D, id);
			}

			// draw polygons
			Facet_iterator pFacet = this->facets_begin();
			for (;pFacet != this->facets_end();pFacet++)
			{				
				// begin polygon assembly
				::glBegin(GL_POLYGON);
					if (use_face_color)
					{
						float r = pFacet->color(0);
						float g = pFacet->color(1);
						float b = pFacet->color(2);
						::glColor3f(r, g, b);
					}
					gl_draw_facet(pFacet, smooth_shading, use_normals, use_vertex_color, use_face_color, use_texture);
				::glEnd();
				// end polygon assembly
			}
			glFlush();
		}

		virtual void gl_draw_facet(Facet_handle pFacet, bool smooth_shading, bool use_normals,
									bool use_vertex_color, bool use_face_color, bool use_texture)
		{
			use_face_color = use_face_color ; // avoid warning

			// one normal per face
			if (use_normals && !smooth_shading)
			{
				//const Facet::Normal_3& normal = pFacet->normal();
				const Normal& normal = pFacet->normal();
				::glNormal3d(to_double(normal[0]),to_double(normal[1]),to_double(normal[2]));
			}

			// revolve around current face to get vertices
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			do
			{
				// one normal per vertex
				if (use_normals && smooth_shading)
				{
					const Normal& normal = pHalfedge->vertex()->normal();
					::glNormal3d(to_double(normal[0]),to_double(normal[1]),to_double(normal[2]));
				}

				// color
				if (use_vertex_color)
				{
					float r = pHalfedge->vertex()->color(0);
					float g = pHalfedge->vertex()->color(1);
					float b = pHalfedge->vertex()->color(2);
					::glColor3f(r, g, b);
				}

				// texture
				if (use_texture && has_texture())
				{
					float s,t;
					if ( use_halfedge_texture_coordinates() ) {
						s = pHalfedge->texture_coordinates(0);
						t = pHalfedge->texture_coordinates(1);
					} else
					{
						s = pHalfedge->vertex()->texture_coordinates(0);
						t = pHalfedge->vertex()->texture_coordinates(1);
					}

					::glTexCoord2f(s, t);
				}

				// polygon assembly is performed per vertex
				const Point& point  = pHalfedge->vertex()->point();
				::glVertex3d(to_double(point[0]),to_double(point[1]),to_double(point[2]));
			}
			while(++pHalfedge != pFacet->facet_begin());
		}

		// superimpose edges
		virtual void superimpose_edges(bool voronoi_edge = false)
		{
			::glBegin(GL_LINES);
			for (Edge_iterator h = this->edges_begin(); h != this->edges_end(); h++)
			{
				if (voronoi_edge)
				{
					Facet_handle pFace1 = h->facet();
					Facet_handle pFace2 = h->opposite()->facet();
					if (pFace1 == NULL || pFace2 == NULL)
					continue;

					const Point &p1 = h->vertex()->point();
					const Point &p2 = h->next()->vertex()->point();
					const Point &p3 = h->next()->next()->vertex()->point();

					kernel k;
					Point d1 = k.construct_circumcenter_3_object()(p1,p2,p3);
					::glVertex3d(to_double(d1[0]),to_double(d1[1]),to_double(d1[2]));

					const Point &pp1 = h->opposite()->vertex()->point();
					const Point &pp2 = h->opposite()->next()->vertex()->point();
					const Point &pp3 = h->opposite()->next()->next()->vertex()->point();
					Point d2 = k.construct_circumcenter_3_object()(pp1,pp2,pp3);
					::glVertex3d(to_double(d2[0]),to_double(d2[1]),to_double(d2[2]));
				}
				else
				{
					// assembly and draw line segment
					const Point& p1 = h->prev()->vertex()->point();
					const Point& p2 = h->vertex()->point();
					::glVertex3d(to_double(p1[0]),to_double(p1[1]),to_double(p1[2]));
					::glVertex3d(to_double(p2[0]),to_double(p2[1]),to_double(p2[2]));
				}
			}
			::glEnd();
		}

		// superimpose vertices
		virtual void superimpose_vertices()
		{
			::glBegin(GL_POINTS);
				for (Point_iterator pPoint = this->points_begin(); pPoint != this->points_end(); pPoint++)
				::glVertex3d(to_double(pPoint->x()),to_double(pPoint->y()),to_double(pPoint->z()));
			::glEnd(); // // end point assembly
		}

		void drawCube()
		{
			::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
			glBegin(GL_QUAD_STRIP);
			  //Quads 1 2 3 4
				glVertex3f( 0.5f, 0.5f, 0.5f);   //V2
				glVertex3f( 0.5f,-0.5f, 0.5f);   //V1
				glVertex3f( 0.5f, 0.5f,-0.5f);   //V4
				glVertex3f( 0.5f,-0.5f,-0.5f);   //V3
				glVertex3f(-0.5f, 0.5f,-0.5f);   //V6
				glVertex3f(-0.5f,-0.5f,-0.5f);   //V5
				glVertex3f(-0.5f, 0.5f, 0.5f);   //V8
				glVertex3f(-0.5f,-0.5f, 0.5f);   //V7
				glVertex3f( 0.5f, 0.5f, 0.5f);   //V2
				glVertex3f( 0.5f,-0.5f, 0.5f);   //V1
			glEnd();
			glBegin(GL_QUADS);
			  //Quad 5
				glVertex3f(-0.5f, 0.5f,-0.5f);   //V6
				glVertex3f(-0.5f, 0.5f, 0.5f);   //V8
				glVertex3f( 0.5f, 0.5f, 0.5f);   //V2
				glVertex3f( 0.5f, 0.5f,-0.5f);   //V4
			  //Quad 6
				glVertex3f(-0.5f,-0.5f, 0.5f);   //V7
				glVertex3f(-0.5f,-0.5f,-0.5f);   //V5
				glVertex3f( 0.5f,-0.5f,-0.5f);   //V3
				glVertex3f( 0.5f,-0.5f, 0.5f);   //V1
			glEnd();
		}

		void gen_glListCube()
		{
			if (id_cube)
				return;
			id_cube = glGenLists(1);
			if (!id_cube)
				return;

			glNewList(id_cube, GL_COMPILE);
				drawCube();
			glEndList();
		}

		// superimpose vertices
		virtual void superimpose_spheres(bool glList, double scale)
		{
			//GLUquadricObj* pQuadric = gluNewQuadric();
			::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

			for (Vertex_iterator pVertex = this->vertices_begin(); pVertex !=  this->vertices_end(); pVertex++)
			{
				::glPushMatrix();
					double radius = to_double(average_edge_length_around(pVertex)*scale);
					::glTranslated(to_double(pVertex->point().x()), to_double(pVertex->point().y()), to_double(pVertex->point().z()));

					if (glList)
					{
						glScaled(radius, radius, radius);
						glCallList(id_cube);
					}
					else
					{
						glScaled(radius, radius, radius);
						drawCube();

						// or

						//::gluSphere(pQuadric,radius,24,24);
					}
				::glPopMatrix();
			}

			//gluDeleteQuadric(pQuadric);
		}

		// Render normals
		virtual void draw_normals()
		{
			//glColor3f(1.f, 0.f, 0.f);
			glBegin(GL_LINES);
			for (Vertex_iterator pVertex = this->vertices_begin(); pVertex !=  this->vertices_end(); pVertex++)
			{
					glVertex3d(to_double(pVertex->point().x()), to_double(pVertex->point().y()), to_double(pVertex->point().z()));
					glVertex3d(to_double(pVertex->point().x()) + to_double(pVertex->normal().x()), to_double(pVertex->point().y()) + to_double(pVertex->normal().y()),
								to_double(pVertex->point().z()) + to_double(pVertex->normal().z()));
			}
			glEnd(); // // end point assembly
		}

		// write in obj file format (OBJ).
		void write_obj(string output_name, int incr  = 1) // 1-based by default
		{
			std::ofstream stream(output_name.c_str());

			// output vertices
			for (Point_iterator pPoint = this->points_begin(); pPoint != this->points_end(); pPoint++)
			stream << 'v' << ' ' << pPoint->x() << ' ' <<
									pPoint->y() << ' ' <<
									pPoint->z() << std::endl;

			// precompute vertex indices
			this->set_index_vertices();

			// output facets
			for (Facet_iterator pFacet = this->facets_begin(); pFacet != this->facets_end(); pFacet++)
			{
				stream << 'f';
				Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
				do
					stream << ' ' << pHalfedge->vertex()->tag()+incr;
				while(++pHalfedge != pFacet->facet_begin());
				stream << std::endl;
			}

			stream.flush();
			stream.close();
		}

        // write in vrml file format (WRL).
		void write_vrml(string output_name)
		{
			/*std::ofstream stream(output_name.c_str());

			CGAL::VRML_2_ostream out(stream);
            out << (*this);

			stream.flush();
			stream.close();*/
		}

		void write_x3d(string output_name, bool write_color)
		{
			std::ofstream stream(output_name.c_str());
		}

		// write in off file format (OFF).
		void write_off(string output_name, bool write_color, bool write_normals)
		{
			std::ofstream file(output_name.c_str());

			if (write_normals)
			{
				file << "N";
			}

			if (write_color)
			{
				file << "C";
			}

			file << "OFF" << endl;

			file << "# MEPP output file" << endl;

			file << this->size_of_vertices() << " " << this->size_of_facets() << " " << (this->size_of_halfedges() >> 1) << endl; // MT

			// output vertices
			for (Vertex_iterator pVert = this->vertices_begin(); pVert != this->vertices_end(); pVert++)
			{
				file << pVert->point().x() << " " << pVert->point().y() << " " << pVert->point().z();

				if (write_normals)
				{
					file << "  ";
					file << (float)(to_double(pVert->normal().x())) << " " << (float)(to_double(pVert->normal().y())) << " " << (float)(to_double(pVert->normal().z()));
				}

				if (write_color)
				{
					file << "  ";
					file << pVert->color(0) << " " << pVert->color(1) << " " << pVert->color(2);
				}
				file << endl;
			}

			// precompute vertex indices
			this->set_index_vertices();

			// output facets
			for (Facet_iterator pFacet = this->facets_begin(); pFacet != this->facets_end(); pFacet++)
			{
				vector<int> indices;
				Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();

				do
				{
					indices.push_back(pHalfedge->vertex()->tag());
				}
				while(++pHalfedge != pFacet->facet_begin());

				file << indices.size();
				for (unsigned int i=0; i<indices.size(); i++)
				{
					file << " " << indices[i];
				}

				file << endl;
			}

			file.flush();
			file.close();
		}

		// draw bounding box
		virtual void gl_draw_bounding_box()
		{
			::glBegin(GL_LINES);

			// along x axis
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));

			// along y axis
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));

			// along z axis
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmin()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymin()),to_double(m_bbox.zmax()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmin()));
			::glVertex3d(to_double(m_bbox.xmax()),to_double(m_bbox.ymax()),to_double(m_bbox.zmax()));

			::glEnd();
		}

		// count #boundaries
		unsigned int calc_nb_boundaries()
		{
			unsigned int nb = 0;
			tag_halfedges(0);
			for (Halfedge_iterator he = this->halfedges_begin(); he != this->halfedges_end(); he++)
			{
				if (he->is_border() && he->tag() == 0)
				{
					nb++;
					Halfedge_handle curr = he;
					do
					{
						curr  = curr->next();
						curr->tag(1);
					}
					while(curr != he);
				}
			}
			m_nb_boundaries = nb;
			return nb;
		}

		// tag component
		void tag_component(Facet_handle pSeedFacet, const int tag_free, const int tag_done)
		{
			pSeedFacet->tag(tag_done);
			std::list<Facet_handle> facets;
			facets.push_front(pSeedFacet);
			while(!facets.empty())
			{
				Facet_handle pFacet = facets.front();
				facets.pop_front();
				pFacet->tag(tag_done);
				Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
				Halfedge_around_facet_circulator end = pHalfedge;
				CGAL_For_all(pHalfedge,end)
				{
					Facet_handle pNFacet = pHalfedge->opposite()->facet();
					if (pNFacet != NULL && pNFacet->tag() == tag_free)
					{
						facets.push_front(pNFacet);
						pNFacet->tag(tag_done);
					}
				}
			}
		}

		// count #components
		unsigned int calc_nb_components()
		{
			unsigned int nb = 0;
			tag_facets(0);
			for (Facet_iterator pFacet = this->facets_begin(); pFacet != this->facets_end(); pFacet++)
			{
				if (pFacet->tag() == 0)
				{
					nb++;
					tag_component(pFacet,0,1);
				}
			}
			m_nb_components = nb;
			return nb;
		}

		// compute the genus
		// V - E + F + B = 2 (C - G)
		// C -> #connected components
		// G : genus
		// B : #boundaries
		int genus()
		{
			int c = this->nb_components();
			int b = this->nb_boundaries();
			int v = this->size_of_vertices();
			int e = this->size_of_halfedges()/2;
			int f = this->size_of_facets();
			return genus(c,v,f,e,b);
		}

		int genus(int c, int v, int f, int e, int b)
		{
			return (2*c+e-b-f-v)/2;
		}

		void copy_from(MEPP_Common_Polyhedron<kernel, items> *new_mesh)
		{
			this->clear();

			CCopyPoly<MEPP_Common_Polyhedron<kernel, items>, kernel > copier;
			copier.copy(new_mesh, this);
			
			this->compute_bounding_box();

			this->compute_normals();
			this->compute_type();

			(void)this->calc_nb_components();
			(void)this->calc_nb_boundaries();

			copy_textures(new_mesh);
		}

		bool has_color() { return m_has_color; }
		void has_color(bool has_col) { m_has_color = has_col; }

		void fit_mesh()
		{
			// Transform mesh for standard size and position

			this->compute_bounding_box();

			// Translate
			FT xaverage = (this->xmax() + this->xmin())/2;
			FT yaverage = (this->ymax() + this->ymin())/2;
			FT zaverage = (this->zmax() + this->zmin())/2;

			Vector transl(-xaverage, -yaverage, -zaverage);

			// Scale
			FT scale = this->xmax()-this->xmin();
			if ( scale < (this->ymax()-this->ymin()) )
			{
				scale = this->ymax()-this->ymin();
			}
			if ( scale < (this->zmax()-this->zmin()) )
			{
				scale = this->zmax()-this->zmin();
			}
			scale = 10/scale;

			Affine_transformation t(CGAL::TRANSLATION, transl);
			std::transform( this->points_begin(), this->points_end(), this->points_begin(), t);

			Affine_transformation s(CGAL::SCALING, scale);
			std::transform( this->points_begin(), this->points_end(), this->points_begin(), s);

			// Compute new bounding box
			this->compute_bounding_box();
		}

#if (0)
		bool load_mesh(string path)
		{
			// Check extension
			size_t dot = path.find_last_of('.');
			if (dot == string::npos)
			{
				return false;
			}

			string extension = path.substr(dot+1);
			for (unsigned int i=0; i<extension.size(); i++)
			{
				extension[i] = toupper(extension[i]);
			}

			if ( (extension == "OFF") )
			{
				if ( !this->load_mesh_off(path) ) { return false; }
			}
			else if ( (extension == "OBJ") )
			{
				if ( !this->load_mesh_obj(path) ) { return false; }
			}
			else if ( (extension == "X3D") )
			{
				if ( !this->load_mesh_x3d(path) ) { return false; }
			}
#if (0)
			// Ajout Céline :
			else if ( (extension == "DAT") )
			{
				if ( !this->load_mesh_dat(path) ) { return false; }
			}
			else if ( (extension == "BW") )
			{
				if ( !this->load_mesh_bw(path) ) { return false; }
			}
#endif
			/*else if ( (extension == "SMF") )
			{
				if ( !this->load_mesh_smf(path) ) { return false; }
			}
			else if ( (extension == "PLY") )
			{
				if ( !this->load_mesh_ply(path) ) { return false; }
			}
			else
			{
				return false;
			}*/

			/*MEPP: if (extension != "DAT" && extension != "BW")  // ajout Céline
				this->fit_mesh();*/

			this->compute_normals();
			this->compute_type();

			return true;
		}
#endif

		int load_mesh_off(std::string filename)
		{
		    // read from stream
			std::ifstream stream(filename.c_str());
			if (!stream)
			{
				return -1;
			}
			stream.close();

			OFF_CGALImporter<HalfedgeDS, Point> poly_builder(filename);
			this->delegate(poly_builder);

			return 0;
		}

		int load_mesh_obj(std::string filename)
		{
			// read from stream
			std::ifstream stream(filename.c_str());
			if (!stream)
			{
				return -1;
			}
			stream.close();

			Builder_obj<HalfedgeDS> builder(filename);
			this->delegate(builder);
			m_has_texture_coordinates = builder.hasTextureCoordinates();
			m_use_halfedge_texture_coordinates = m_has_texture_coordinates;

			return 0;
		}

		int load_mesh_smf(std::string filename)
		{
			// read from stream
			std::ifstream stream(filename.c_str());
			if (!stream)
			{
				return -1;
			}
			stream.close();

			// build the polyhedron
			Builder_smf<HalfedgeDS> poly_builder(filename);
			this->delegate(poly_builder);

			return 0;
		}
		
		int load_mesh_ply(std::string filename)
		{
			// read from stream
			std::ifstream stream(filename.c_str());
			if (!stream)
			{
				return -1;
			}
			stream.close();

			// build the polyhedron
			Builder_ply<HalfedgeDS> poly_builder(filename);
			this->delegate(poly_builder);
			m_has_texture_coordinates = poly_builder.hasTextureCoordinates();

			return 0;
		}

		int load_mesh_x3d(string filename)
		{
			IFSData ifsdata;
			X3DMeshExtractor parser;
			parser.load(filename, ifsdata);

			X3D_CGALImporter<HalfedgeDS, Point> poly_builder(&(ifsdata.vertex), &(ifsdata.face));
			this->delegate(poly_builder);

			if (ifsdata.colorMode == COLOR_PER_FACE)
			{
				Facet_iterator f = this->facets_begin();

				for (unsigned int index=0; index<ifsdata.color.size(); index++)
				{
					Halfedge_around_facet_circulator pHalfedge = f->facet_begin();
					do
					{
						pHalfedge->vertex()->color(ifsdata.color[index][0], ifsdata.color[index][1], ifsdata.color[index][2]);
					}
					while(++pHalfedge != f->facet_begin());
					f++;
				}

				this->has_color(true);

			}
			else if (ifsdata.colorMode == COLOR_PER_VERTEX)
			{
				Vertex_iterator v = this->vertices_begin();
				for (unsigned int index=0; index<ifsdata.color.size(); index++)
				{
					v->color(ifsdata.color[index][0], ifsdata.color[index][1], ifsdata.color[index][2]);
					v++;
				}

				this->has_color(true);
			}
			else
			{
				if (ifsdata.color.size()) // color per vertex
				{
					Vertex_iterator v = this->vertices_begin();
					for (unsigned int index=0; index<ifsdata.color.size(); index++)
					{
						v->color(ifsdata.color[index][0], ifsdata.color[index][1], ifsdata.color[index][2]);
						v++;
					}

					this->has_color(true);
				}
			}

			return 0;
		}

#ifdef WITH_ASSIMP
		int load_mesh_assimp(aiMesh* mesh)
		{
			// build the polyhedron
			Builder_dae<HalfedgeDS> poly_builder(mesh);
			this->delegate(poly_builder);
			if (poly_builder.loadedSucess) {				
				m_has_texture_coordinates = poly_builder.hasTextureCoordinates();				
				return 0;
			}
			else
				return -1;
		}
#endif

        void triangulate()
        {
            Facet_iterator f = this->facets_begin();
            Facet_iterator f2 = this->facets_begin();
            do //for (; f != this->facets_end(); f++)
            {
                f = f2;
                if (f == this->facets_end())
                {
                    break;
                }
                f2++;

				if (!(f->is_triangle()))
				{					
					int num = (int)(f->facet_degree() - 3);
					Halfedge_handle h = f->halfedge();				

					h = this->make_hole(h);

					Halfedge_handle g = h->next();
					g = g->next();
					Halfedge_handle new_he = this->add_facet_to_border (h, g);
					new_he->texture_coordinates(h->texture_coordinates(0),h->texture_coordinates(1));
					new_he->opposite()->texture_coordinates(g->texture_coordinates(0),g->texture_coordinates(1));
					g=new_he;

					num--;
					while (num != 0)
					{
						g = g->opposite();
						g = g->next();				
						Halfedge_handle new_he = this->add_facet_to_border (h, g);		
						new_he->texture_coordinates(h->texture_coordinates(0),h->texture_coordinates(1));
						new_he->opposite()->texture_coordinates(g->texture_coordinates(0),g->texture_coordinates(1));
						g=new_he;
			
						num--;
					}

					this->fill_hole(h);
				}

            } while (true);

			this->compute_normals();
			this->compute_type();
        }

		void set_textures(QStringList files)
		{
			// MT : not required ???
			for (unsigned int i = 0; i < m_texture_array.size(); i++)
			{
				if ( m_texture_array[i].m_id > 0 )
					glDeleteTextures(1, &m_texture_array[i].m_id);
			}

			m_texture_array.clear();

			for (int i = 0; i < files.size(); i++)
			{
				Texture texture;

				texture.m_name = files[i];
				bool load_sucess = texture.m_data.load(texture.m_name);
				if (load_sucess)
				{
					texture.m_id = convertTextureToGLFormat(texture.m_data);

					m_texture_array.push_back(texture);
				}
			}
		}

		void set_textures(vector<Texture> textures)
		{
			// MT : not required ???
			for (unsigned int i = 0; i < m_texture_array.size(); i++)
			{
				if ( m_texture_array[i].m_id > 0 )
					glDeleteTextures(1, &m_texture_array[i].m_id);
			}

			m_texture_array.clear();

			for (unsigned int i = 0; i < textures.size(); i++) // MT add unsigned
			{
				Texture texture = textures[i];			
				texture.m_id = convertTextureToGLFormat(texture.m_data);
				m_texture_array.push_back(texture);				
			}
		}

		GLuint convertTextureToGLFormat( QImage& _texsrc )
		{
		  {
			// adjust texture size: 2^k * 2^l
			int tex_w, w( _texsrc.width()  );
			int tex_h, h( _texsrc.height() );

			for (tex_w=1; tex_w <= w; tex_w <<= 1) {};
			for (tex_h=1; tex_h <= h; tex_h <<= 1) {};
			tex_w >>= 1;
			tex_h >>= 1;
			_texsrc = _texsrc.scaled( tex_w, tex_h, Qt::IgnoreAspectRatio, Qt::SmoothTransformation );
		  }

		  QImage texture( QGLWidget::convertToGLFormat ( _texsrc ) );
		  
		  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		  glPixelStorei(GL_UNPACK_SKIP_ROWS,   0);
		  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
		  glPixelStorei(GL_UNPACK_ALIGNMENT,   1);
		  glPixelStorei(GL_PACK_ROW_LENGTH,    0);
		  glPixelStorei(GL_PACK_SKIP_ROWS,     0);
		  glPixelStorei(GL_PACK_SKIP_PIXELS,   0);
		  glPixelStorei(GL_PACK_ALIGNMENT,     1);    
		  
		  GLuint tex_id_ = 0;
		  glGenTextures(1, &tex_id_);
		  glBindTexture(GL_TEXTURE_2D, tex_id_);
		    
		  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);      
		  
		  glTexImage2D(GL_TEXTURE_2D,       // target
				   0,                   // level
				   GL_RGBA,             // internal format
				   texture.width(),     // width  (2^n)
				   texture.height(),    // height (2^m)
				   0,                   // border
				   GL_RGBA,             // format
				   GL_UNSIGNED_BYTE,    // type
				   texture.bits() );    // pointer to pixels

		  return tex_id_;
		}

		inline void use_halfedge_texture_coordinates(bool value)
		{
			 m_use_halfedge_texture_coordinates = value;
		}

		inline bool use_halfedge_texture_coordinates()
		{
			return ( m_use_halfedge_texture_coordinates );
		}

		inline bool has_texture()
		{
			return (m_texture_array.size() > 0);
		}

		void apply_texture_to_vertex_colors(unsigned int indexe = 0,
			float modifier_offset_texture_coordinate_s = 1.0f, float modifier_offset_texture_coordinate_t = 1.0f)
		{
			if (m_texture_array.size() > indexe && indexe)
				for (Vertex_iterator it = this->vertices_begin(); it != this->vertices_end(); ++it)
				{
					std::vector<int> position_coordinate_upper_left;
					std::vector<int> position_coordinate_upper_right;
					std::vector<int> position_coordinate_down_left;
					std::vector<int> position_coordinate_down_right;
					std::vector<int> position_coordinate_target;

					std::vector<int> color_upper_left;
					std::vector<int> color_upper_right;
					std::vector<int> color_down_left;
					std::vector<int> color_down_right;
					std::vector<int> color_target;

					std::vector<float> texture_coordinate_upper_left;
					std::vector<float> texture_coordinate_upper_right;
					std::vector<float> texture_coordinate_down_left;
					std::vector<float> texture_coordinate_down_right;
					std::vector<float> texture_coordinate_target;

					float texture_coordinate_s = it->texture_coordinates(0);
					float texture_coordinate_t = 1-it->texture_coordinates(1);

					float offset_texture_coordinate_s = (1.0f/m_texture_array[indexe].m_data.size().width())*modifier_offset_texture_coordinate_s;
					float offset_texture_coordinate_t = (1.0f/m_texture_array[indexe].m_data.size().height())*modifier_offset_texture_coordinate_t;

					texture_coordinate_upper_left.push_back(texture_coordinate_s-offset_texture_coordinate_s);
					texture_coordinate_upper_left.push_back(texture_coordinate_t+offset_texture_coordinate_t);

					texture_coordinate_upper_right.push_back(texture_coordinate_s+offset_texture_coordinate_s);
					texture_coordinate_upper_right.push_back(texture_coordinate_t+offset_texture_coordinate_t);

					texture_coordinate_down_left.push_back(texture_coordinate_s-offset_texture_coordinate_s);
					texture_coordinate_down_left.push_back(texture_coordinate_t-offset_texture_coordinate_t);

					texture_coordinate_down_right.push_back(texture_coordinate_s+offset_texture_coordinate_s);
					texture_coordinate_down_right.push_back(texture_coordinate_t-offset_texture_coordinate_t);

					texture_coordinate_target.push_back(texture_coordinate_s);
					texture_coordinate_target.push_back(texture_coordinate_t);

					GetPos(indexe, position_coordinate_upper_left, texture_coordinate_upper_left);
					GetPixel(indexe, position_coordinate_upper_left, color_upper_left);

					GetPos(indexe, position_coordinate_upper_right, texture_coordinate_upper_right);
					GetPixel(indexe, position_coordinate_upper_right, color_upper_right);

					GetPos(indexe, position_coordinate_down_left, texture_coordinate_down_left);
					GetPixel(indexe, position_coordinate_down_left, color_down_left);

					GetPos(indexe, position_coordinate_down_right, texture_coordinate_down_right);
					GetPixel(indexe, position_coordinate_down_right, color_down_right);

					GetPos(indexe, position_coordinate_target, texture_coordinate_target);
					GetPixel(indexe, position_coordinate_target, color_target);

					std::vector<int> color_p1;
					std::vector<int> color_p2;

					interpolate(color_p1,color_upper_left,color_down_left,texture_coordinate_target[1],texture_coordinate_upper_left[1],texture_coordinate_down_left[1]);

					interpolate(color_p2,color_upper_right,color_down_right,texture_coordinate_target[1],texture_coordinate_upper_right[1],texture_coordinate_down_right[1]);

					interpolate(color_target,color_p1,color_p2,texture_coordinate_target[0],texture_coordinate_upper_right[0],texture_coordinate_down_right[0]);

					it->color(color_target[0]/255.,color_target[1]/255.,color_target[2]/255.);
				}
		}

		void copy_textures(const MEPP_Common_Polyhedron* mesh)
		{
			m_has_texture_coordinates = mesh->m_has_texture_coordinates;
			m_use_halfedge_texture_coordinates = mesh->m_use_halfedge_texture_coordinates;

			m_texture_array = mesh->m_texture_array;
			for (unsigned int i = 0; i < m_texture_array.size(); i++)
				m_texture_array[i].m_id = convertTextureToGLFormat(m_texture_array[i].m_data);
		}

		Texture get_texture( int index )
		{			
			return m_texture_array[index];		
		}	
	protected:
		void GetPos(int indexe, std::vector<int> &position, const std::vector<float> &texture_coordinate_target) 
		{
			position.clear();
			position.push_back(texture_coordinate_target[0]*m_texture_array[indexe].m_data.size().width());
			position.push_back(texture_coordinate_target[1]*m_texture_array[indexe].m_data.size().height());

			if (position[0] >= m_texture_array[indexe].m_data.size().width())
				position[0] = m_texture_array[indexe].m_data.size().width() - 1;
			if (position[1] >= m_texture_array[indexe].m_data.size().height())
				position[1] = m_texture_array[indexe].m_data.size().height() - 1;

			if (position[0] < 0)
				position[0] = 0;
			if (position[1] < 0)
				position[1] = 0;
		}

		void GetPixel(int indexe, const std::vector<int> &position, std::vector<int> &color ) 
		{
			color.clear();
			color.resize(3);

			QRgb rgb = m_texture_array[indexe].m_data.pixel(position[0],position[1]);
			QColor pixel_color = QColor(rgb);
			pixel_color.getRgb(&color[0],&color[1],&color[2]);
		}

		void interpolate(std::vector<int> &color,const std::vector<int> &color1, const std::vector<int> &color2,
			const float &y, const float &y1, const float &y2)
		{
			color.clear();
			color.resize(3);
			for (unsigned int i = 0; i < color.size(); i++)
			{
				if ((y2 - y1)!=0)
					color[i] = color1[i] + (color2[i]-color1[i])*(y - y1)/(y2 - y1);
				else
					color[i] = color1[i];
			}
		}

#if (0)
		bool load_mesh_off(string filename)
		{
		    // Read from stream
			std::ifstream stream(filename.c_str());
			if (!stream)
			{
				return false;
			}

            // Check for color
			string s;
			stream >> s;
			while (s[0] == '#')
			{
			    stream >> s;
			}

			/*MEPP: if ((s == "COFF") || (s == "NCOFF")) // MT
			{
			    // Color: using custom loader
			    stream.close();
			    return this->load_mesh_off_color(filename);
			}*/

			stream.close();

			// No color: using CGAL loader
			stream.open(filename.c_str());
			if (!stream)
			{
				return false;
			}
			stream >> *this;

			return true;
		}

		/*MEPP: bool load_mesh_off_color(string filename)
		{
			this->has_color(true);

			// build the polyhedron
			OFF_CGALImporter<HalfedgeDS, Point> poly_builder(filename);
			this->delegate(poly_builder);

			return true;
		}*/
#endif

#if (0)
		// Ajout Céline :
		bool load_mesh_dat(string filename)
		{
			FILE *pFile = fopen(filename.c_str(),"rt");
			if (pFile == NULL)
				return false;

			Builder_dat<HalfedgeDS, MEPP_Common_Polyhedron<kernel, items>, kernel > builder(pFile);
			this->delegate(builder);

			// on récupère le nombre de niveaux de résolution du polyèdre semi-régulier :
			int nb_subd_levels = builder.nb_levels();
			this->tag(nb_subd_levels);	// on stocke dans m_tag pour l'instant !

			fclose(pFile);
			return true;
		}

		bool load_mesh_bw(string filename)
		{
			FILE *pFile = fopen(filename.c_str(),"rt");
			if (pFile == NULL)
				return false;

			Builder_bw<HalfedgeDS, MEPP_Common_Polyhedron<kernel, items>, kernel > builder(pFile);
			this->delegate(builder);

			// on récupère le nombre de niveaux de résolution du polyèdre semi-régulier :
			int nb_subd_levels = builder.nb_levels();
			this->tag(nb_subd_levels);	// on stocke dans m_tag pour l'instant !

			fclose(pFile);
			return true;
		}
#endif

};

#endif
