#ifndef SMF_CGALIMPORTER_H
#define SMF_CGALIMPORTER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cstdio>

#include "Correct_CGAL_Structure.h"

using namespace std;

template <class HDS>
class Builder_smf : public CGAL::Modifier_base<HDS>
{

private:

	typedef typename HDS::Vertex_handle Vertex_handle;
    typedef typename HDS::Vertex::Point Point;
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    string m_filename;
    FILE *m_pFile;
	vector< vector <float> > Vertices_coordinates;
	vector< vector <float> > Vertices_color;
	vector< vector <int> >   Facets;
	

public:

    Builder_smf(string filename)
    {
        m_filename = filename;
        m_pFile = fopen(filename.c_str(),"rt");
    }

    bool is_ok()
    {
        return (m_pFile != NULL);
    }

    ~Builder_smf() {}

    void operator()(HDS& hds)
    {
        Builder builder(hds,true);
        builder.begin_surface(3,1,6);         
		this->read_vertices();
		this->read_colors();
        this->read_facets();
		
		// Correction of orientation
		//CORRECT_CGAL_STRUCTURE Structure_Corrector(&Vertices_coordinates, &Vertices_color, &Facets);		
		//Facets = *(Structure_Corrector.Correct_Facet_Orientation());

		this->construct(builder);
        builder.end_surface();
		
		if (builder.check_unconnected_vertices())
		{
			builder.remove_unconnected_vertices();
		}

        fclose(m_pFile);
    }

private:

  // read vertex coordinates
    void read_vertices()
    {
        fseek(m_pFile,0,SEEK_SET);

        char pLine[512];
        while(fgets(pLine,512,m_pFile))
		{
			if (pLine[0] == 'v')
			{
				float x,y,z;
				if (sscanf(pLine,"v %f %f %f",&x,&y,&z) == 3)
				{
					vector<float> Temp_position;
					Temp_position.push_back(x);
					Temp_position.push_back(y);
					Temp_position.push_back(z);
					
					Vertices_coordinates.push_back(Temp_position);
				}
			}
		}
	}


	void read_colors()
	{
		fseek(m_pFile,0,SEEK_SET);

        char pLine[512];
        while(fgets(pLine,512,m_pFile))
		{
			if (pLine[0] == 'c')
			{
				float c0,c1,c2;
				if (sscanf(pLine,"c %f %f %f",&c0,&c1,&c2) == 3)
				{
					vector<float> Temp_color;
					Temp_color.push_back(c0);
					Temp_color.push_back(c1);
					Temp_color.push_back(c2);

					Vertices_color.push_back(Temp_color);
				}
			}
		}
	}
    // read facets and uv coordinates per halfedge
    void read_facets()
    {
        fseek(m_pFile,0,SEEK_SET);

        char pLine[512];
        while(fgets(pLine,512,m_pFile))
        {
			if (pLine[0] == 'f')
			{
				int i0,i1,i2;
				if (sscanf(pLine,"f %d %d %d",&i0,&i1,&i2) == 3)
				{
					vector<int> Temp_facet;
					Temp_facet.push_back(i0-1);
					Temp_facet.push_back(i1-1);
					Temp_facet.push_back(i2-1);

					Facets.push_back(Temp_facet);
				}
			}			
        }        
    }

	void construct(Builder &builder)
	{
		int Number_vertices = this->Vertices_coordinates.size();
		int Number_faces = this->Facets.size();

		for (int i = 0 ; i < Number_vertices; i++)
		{
			float coord[3];
			coord[0] = this->Vertices_coordinates[i][0];
			coord[1] = this->Vertices_coordinates[i][1];
			coord[2] = this->Vertices_coordinates[i][2];

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
			for (int j = 0; j < (int)Facets[i].size();j++)	// MT
			{
				int index = this->Facets[i][j];
				builder.add_vertex_to_facet(index);
			}
			builder.end_facet();
		}			
	}
};

#endif // SMF_CGALIMPORTER_H
