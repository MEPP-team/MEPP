#ifndef OFF_CGALIMPORTER_H
#define OFF_CGALIMPORTER_H

#include <vector>
#include <sstream>
#include <fstream>
#include <list>
#include <queue>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include "Correct_CGAL_Structure.h"

using namespace std;

// OFF Importer for CGAL
// Builds an polyhedron from a list of vertices and indices (and colors)
// Loads vertex color if available

// Structure:
// [C]OFF (for simple OFF or vertex color OFF)
// # comment lines (unlimited)
// num_vert num_faces num_edges(0)
// x y z r g b a
// ...
// n index0 ... indexn
// ...

template <class HDS, class Point>
class OFF_CGALImporter : public CGAL::Modifier_base<HDS>
{
	public:

		typedef typename HDS::Vertex_handle Vertex_handle;

		OFF_CGALImporter(std::string filename)
		{
			m_filename = filename;
		}

		void operator()( HDS& hds)
		{

			std::ifstream file(m_filename.c_str(), std::ios::in);
			if (!file.is_open())
			{
				return;
			}
			

			std::string s;
			bool color = false;
			bool normal = false; // MT
			unsigned int num_vert = 0;
			unsigned int num_faces = 0;
			
			getline(file,s);
			while (s[0] == '#')
			{
			    getline(file,s);
			}

			// name
			//file >> s;
			if (s.compare(0,4,"COFF")==0) 
				color = true;
			if (s.compare(0,4,"NOFF")==0)
				normal = true;
			if ( (s.compare(0,5,"NCOFF")==0) || (s.compare(0,5,"CNOFF")==0) )
			{
				normal = true; 
				color = true; 
			} // MT

			// comments
			file >> s;
			while (s[0] == '#')
			{
				getline(file, s);
				file >> s;
			}

			// header
			stringstream sstream(s);
			sstream >> num_vert;
			file >> num_faces;
			file >> s;

            CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

			builder.begin_surface((int)num_vert, (int)num_faces);

			float Temp_coord[3];
			float Temp_color[3];

			float Max_color_value = 0;
			
			bool Is_color_integer = false;
			// geom + color
			for (unsigned int i = 0; i < num_vert; i++)
			{
				vector<float> Vert_coord;
				vector<float> Vert_color;

				file >> Temp_coord[0];
				file >> Temp_coord[1];
				file >> Temp_coord[2];
				
				for (unsigned j = 0 ; j < 3; j++)
					Vert_coord.push_back(Temp_coord[j]);
				
				Vertices_position.push_back(Vert_coord);

				if (normal) // MT : on ignore les normales
				{
					file >> s;
					file >> s;
					file >> s;
				}

				if (color)
				{
					file >> Temp_color[0];
					file >> Temp_color[1];
					file >> Temp_color[2];
					for (unsigned j = 0 ; j < 3; j++)
						Vert_color.push_back(Temp_color[j]);
					
					Vertices_color.push_back(Vert_color);

					for (unsigned j = 0; j < 3; j++)
					{
						if (Temp_color[j] > Max_color_value)
							Max_color_value = Temp_color[j];
					}
				}
			}			

			// Color value can be integer between [0; 255]
			// or float between [0.0 ; 1.0]

			// if max of color value > 1.0, we consider that color values are integers;
			if (Max_color_value > 2.0)
				Is_color_integer = true;
			else
				Is_color_integer = false;

			// connectivity
			for (unsigned int i=0; i < num_faces; i++)
			{
				unsigned int face_size;
				unsigned int index;
				vector<int> vect_index;

				file >> face_size;				

				for (unsigned int j=0; j<face_size; j++)
				{
					file >> index;
					vect_index.push_back(index);					
				}

				Facets.push_back(vect_index);
			}

			file.close();


			// construction
			for (unsigned int i = 0; i < num_vert; i++)
			{
				Vertex_handle vertex = builder.add_vertex(Point(Vertices_position[i][0], Vertices_position[i][1], Vertices_position[i][2]));

				if (color)
				{
					if (Is_color_integer)
					{
						int RGB[3];
						RGB[0] = (int)floor(Vertices_color[i][0] + 0.5);
						RGB[1] = (int)floor(Vertices_color[i][1] + 0.5);
						RGB[2] = (int)floor(Vertices_color[i][2] + 0.5);

						vertex->color((float)RGB[0]/255.0, (float)RGB[1]/255.0, (float)RGB[2]/255.0);
					}
					else
						vertex->color(Vertices_color[i][0], Vertices_color[i][1], Vertices_color[i][2]);						
				}
			}

			// Correction of orientation
			//CORRECT_CGAL_STRUCTURE Structure_Corrector(&Vertices_position, &Vertices_color, &Facets);			
			//Facets = *(Structure_Corrector.Correct_Facet_Orientation());			
			
			// connectivity
			for (unsigned int i = 0; i < num_faces; i++)
			{
				builder.begin_facet();

				for (unsigned int j = 0; j < Facets[i].size(); j++)
				{					
					builder.add_vertex_to_facet(Facets[i][j]);
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

		std::string m_filename;		
		
                vector<vector<float> > Vertices_position;
                vector<vector<float> > Vertices_color;

                vector<vector<int> > Facets;
};

#endif // OFF_CGALIMPORTER_H

