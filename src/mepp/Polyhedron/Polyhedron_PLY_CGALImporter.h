#ifndef PLY_CGALIMPORTER_H
#define PLY_CGALIMPORTER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cstdio>
#include "ply.h"
#include "mesh.h"
#include "Correct_CGAL_Structure.h"

// Can only read triangular meshes 

using namespace std;



template <class HDS>
class Builder_ply : public CGAL::Modifier_base<HDS>
{

private:

	typedef typename HDS::Vertex_handle Vertex_handle;
    typedef typename HDS::Vertex::Point Point;
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    string m_filename;

    FILE *m_pFile;
	vector< vector<float> > Vertices_position;
	vector< vector<float> > Vertices_color;
	vector< vector<int> >   Facets;
	

	Mesh_ply mesh;

public:

    Builder_ply(string filename)
    {
        m_filename = filename;
        
    }

    bool is_ok()
    {
        return (m_pFile != NULL);
    }

    ~Builder_ply() {}

    void operator()(HDS& hds)
    {
        Builder builder(hds,true);
        builder.begin_surface(3,1,6);
            
		this->initialize();
		// Correction of orientation
		//CORRECT_CGAL_STRUCTURE Structure_Corrector(&Vertices_position, &Vertices_color, &Facets);		
		//Facets = *(Structure_Corrector.Correct_Facet_Orientation());

		this->construct(builder);
        builder.end_surface();
		
		if (builder.check_unconnected_vertices())
		{
			builder.remove_unconnected_vertices();
		}
    }


private:
	void initialize()
	{		
		Import_PLY(this->m_filename.c_str(), &mesh);
		
		int Number_vertices = (int)this->mesh.mVertices.size();

		
		int Number_faces = (int)this->mesh.mIndices.size()/3;
		
		for (int i = 0 ; i < Number_vertices; i++)
		{
			vector<float> Temp_position;
			Temp_position.push_back(this->mesh.mVertices[i].x);
			Temp_position.push_back(this->mesh.mVertices[i].y);
			Temp_position.push_back(this->mesh.mVertices[i].z);

			Vertices_position.push_back(Temp_position);

			//Vertex_handle vertex = builder.add_vertex(Point(coord[0],coord[1],coord[2]));
			if (this->mesh.mColors.size() != 0)
			{
				vector<float> Temp_color;
				Temp_color.push_back(mesh.mColors[i].x);
				Temp_color.push_back(mesh.mColors[i].y);
				Temp_color.push_back(mesh.mColors[i].z);

				Vertices_color.push_back(Temp_color);
				//vertex->color(rgb[0], rgb[1], rgb[2]);
			}
		}
		for (int i = 0 ; i < Number_faces; i++)
		{
			vector<int> Temp_facet;
			//builder.begin_facet();
			for (int j = 0; j < 3;j++)
			{				
				Temp_facet.push_back(this->mesh.mIndices[3*i+j]);
				//builder.add_vertex_to_facet(index);
			}
			Facets.push_back(Temp_facet);
			//builder.end_facet();
		}		


	}    

	void construct(Builder &builder)
	{
		int Number_vertices = (int)this->Vertices_position.size();
		int Number_faces = (int)this->Facets.size();

		for (int i = 0 ; i < Number_vertices; i++)
		{
			float coord[3];
			coord[0] = this->Vertices_position[i][0];
			coord[1] = this->Vertices_position[i][1];
			coord[2] = this->Vertices_position[i][2];

			Vertex_handle vertex = builder.add_vertex(Point(coord[0],coord[1],coord[2]));
			if (this->Vertices_color.size() != 0)
			{
				float rgb[3];
				rgb[0] = this->Vertices_color[i][0];
				rgb[1] = this->Vertices_color[i][1];
				rgb[2] = this->Vertices_color[i][2];
				vertex->color(rgb[0], rgb[1], rgb[2]);
			}
		}
		for (int i = 0 ; i < Number_faces; i++)
		{
			builder.begin_facet();
			for (int j = 0; j < (int)Facets[i].size();j++)
			{
				int index = this->Facets[i][j];
				builder.add_vertex_to_facet(index);
			}
			builder.end_facet();
		}
	}
};

#endif // PLY_CGALIMPORTER_H
