#ifndef CPOLYHEDRON_FROM_POLYGON_BUILDER_3
#define CPOLYHEDRON_FROM_POLYGON_BUILDER_3

/*!
 * \file CPolyhedron_from_polygon_builder_3.h
 * \brief An incremental builder to build a polyhedron
 * \author Cyril Leconte
 */

#include "../../../../mepp/Polyhedron/polyhedron.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>

/*!
 *
 */
template <class HDS>
class CPolyhedron_from_polygon_builder_3 : public CGAL::Modifier_base<HDS> {

private:
	typedef typename HDS::Traits::Point_3									Point_3;
	typedef typename std::vector<unsigned long>								Indices;
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS>			Builder;

	// Member variables
	std::vector<Point_3>													m_Sorted_vertices;
	std::vector<Indices>													m_Facets_indices;
	std::vector<bool>														m_Used_vertices;

public:
	// Constructors
	CPolyhedron_from_polygon_builder_3() {}
	void add_triangle(Facet_handle f, bool invert)
	{
		Indices	vi;
		Halfedge_handle he = f->facet_begin();

		if(he->vertex()->Label == 0xFFFFFFFF) add_vertex(he->vertex()->point(), he->vertex()->Label);
		vi.push_back(he->vertex()->Label);
		if(!invert)
		{
			if(he->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->vertex()->point(), he->next()->vertex()->Label);
			vi.push_back(he->next()->vertex()->Label);
			if(he->next()->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);
			vi.push_back(he->next()->next()->vertex()->Label);
		}
		else
		{
			if(he->next()->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);
			vi.push_back(he->next()->next()->vertex()->Label);
			if(he->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->vertex()->point(), he->next()->vertex()->Label);
			vi.push_back(he->next()->vertex()->Label);
		}
		m_Used_vertices[vi[0]] = true;
		m_Used_vertices[vi[1]] = true;
		m_Used_vertices[vi[2]] = true;
		m_Facets_indices.push_back(vi);
	}

	void add_triangle(vector<vector<unsigned long>> &T, Halfedge_handle &he)
	{
		for(unsigned int i = 0;i != T.size();++i)
		{
			for(unsigned int j = 0 ; j != 3 ; ++j)
			{
				switch (T[i][j])
				{
				case 0xFFFFFFFF:
					if(he->vertex()->Label != 0xFFFFFFFF)
					{
						T[i][j] = he->vertex()->Label;
					}
					else
					{
						T[i][j] = m_Sorted_vertices.size();
						he->vertex()->Label = T[i][j];
						m_Sorted_vertices.push_back(he->vertex()->point());
						m_Used_vertices.push_back(false);
					}
					break;
				case 0xFFFFFFFE:
					if(he->next()->vertex()->Label != 0xFFFFFFFF)
					{
						T[i][j] = he->next()->vertex()->Label;
					}
					else
					{
						T[i][j] = m_Sorted_vertices.size();
						he->next()->vertex()->Label = T[i][j];
						m_Sorted_vertices.push_back(he->next()->vertex()->point());
						m_Used_vertices.push_back(false);
					}
					break;
				case 0xFFFFFFFD:
					if(he->next()->next()->vertex()->Label != 0xFFFFFFFF)
					{
						T[i][j] = he->next()->next()->vertex()->Label;
					}
					else
					{
						T[i][j] = m_Sorted_vertices.size();
						he->next()->next()->vertex()->Label = T[i][j];
						m_Sorted_vertices.push_back(he->next()->next()->vertex()->point());
						m_Used_vertices.push_back(false);
					}
					break;
				}
			}
			m_Used_vertices[T[i][0]] = true;
			m_Used_vertices[T[i][1]] = true;
			m_Used_vertices[T[i][2]] = true;
			m_Facets_indices.push_back(T[i]);
		}
	}

	void add_vertex(Point3d &p, unsigned long &l)
	{
		l = m_Sorted_vertices.size();
		m_Sorted_vertices.push_back(p);
		m_Used_vertices.push_back(false);
	}

	void clear_unused_vertices()
	{
		unsigned long i = 0;
		unsigned long indice;

		while(!m_Used_vertices.empty() && !m_Used_vertices.back())
		{
			m_Used_vertices.pop_back();
			m_Sorted_vertices.pop_back();
		}

		while(i != m_Sorted_vertices.size())
		{
			if(!m_Used_vertices[i])
			{
				m_Used_vertices.pop_back();
				m_Sorted_vertices[i] = m_Sorted_vertices.back();
				m_Sorted_vertices.pop_back();

				indice = m_Sorted_vertices.size();
				for(unsigned int j = 0;j != m_Facets_indices.size();j++)
				{
					if(m_Facets_indices[j][0] == indice) m_Facets_indices[j][0] = i;
					else if(m_Facets_indices[j][1] == indice) m_Facets_indices[j][1] = i;
					else if(m_Facets_indices[j][2] == indice) m_Facets_indices[j][2] = i;
				}

				while(!m_Used_vertices.back())
				{
					m_Used_vertices.pop_back();
					m_Sorted_vertices.pop_back();
				}
			}
			++i;
		}
	}
	void operator()(HDS& hds)
	{
		Builder B(hds, true);
		B.begin_surface(3,1);
		add_vertices(B);	
		add_facets(B);
		B.end_surface();
	}
	
private:
	void add_vertices(Builder& B)
	{
		for (int i = 0; i != this->m_Sorted_vertices.size(); i++)
		{
			B.add_vertex(this->m_Sorted_vertices[i]);
		}
	}

	void add_facets(Builder &B)
	{
		for (int i = 0; i != this->m_Facets_indices.size(); i++)
		{
			B.begin_facet();
			for (int j = 0; j != this->m_Facets_indices[i].size(); j++)
			{
				B.add_vertex_to_facet(this->m_Facets_indices[i][j]);
			}
			B.end_facet();
		}
	}
};

#endif // CPOLYHEDRON_FROM_POLYGON_BUILDER_3
