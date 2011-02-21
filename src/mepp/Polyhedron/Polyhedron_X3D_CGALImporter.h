#ifndef X3D_CGALIMPORTER_H
#define X3D_CGALIMPORTER_H

#include <vector>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "Correct_CGAL_Structure.h"

// X3D Importer for CGAL
// Builds an polyhedron from a list of vertices and indices
template <class HDS, class Point>
class X3D_CGALImporter : public CGAL::Modifier_base<HDS> 
{
	public:

		X3D_CGALImporter()
		{
		
		}

		X3D_CGALImporter(std::vector< std::vector<float> >* vertices, std::vector< std::vector<int> >* facets)
		{
			m_vertices = vertices;
			m_facets = facets;
		}
		
		void operator()( HDS& hds)
		{
			CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true);

			builder.begin_surface((int)m_vertices->size(), (int)m_facets->size());
			
			for (unsigned int i=0; i<m_vertices->size(); i++)
			{
				builder.add_vertex(Point((*m_vertices)[i][0], (*m_vertices)[i][1], (*m_vertices)[i][2]));
			}

			// Correction of orientation
			//CORRECT_CGAL_STRUCTURE Structure_Corrector(m_vertices, NULL , m_facets);		
			//m_facets = Structure_Corrector.Correct_Facet_Orientation();
			
			for (unsigned int i=0; i<m_facets->size(); i++)
			{
				builder.begin_facet();

				for (unsigned int j=0; j<(*m_facets)[i].size(); j++)
				{
					builder.add_vertex_to_facet((*m_facets)[i][j]);
				}

				builder.end_facet();
			}
		
			builder.end_surface();

			if (builder.check_unconnected_vertices())
			{
				builder.remove_unconnected_vertices();
			}

		}
	private:
		std::vector< std::vector<float> >* m_vertices;
		std::vector< std::vector<int> >* m_facets;

};

#endif // X3D_CGALIMPORTER_H

