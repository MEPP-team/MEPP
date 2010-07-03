
#ifndef OBJ_CGALIMPORTER_H
#define OBJ_CGALIMPORTER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cstdio>

#include "Correct_CGAL_Structure.h"

using namespace std;

template <class HDS>
class Builder_obj : public CGAL::Modifier_base<HDS>
{

private:

	typedef typename HDS::Vertex_handle Vertex_handle;
    typedef typename HDS::Vertex::Point Point;
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    string m_filename;
    FILE *m_pFile;

	vector< vector<float> > Vertices_position;
	vector< vector<float> > Vertices_color;
	vector< vector<int> > Facets;

	bool Is_integer;


public:

    Builder_obj(string filename)
    {
        m_filename = filename;
        m_pFile = fopen(filename.c_str(),"rt");
    }

    bool is_ok()
    {
        return (m_pFile != NULL);
    }

    ~Builder_obj() {}

    void operator()(HDS& hds)
    {
        Builder builder(hds,true);
        builder.begin_surface(3,1,6);
            //this->read_vertices(builder);
            //this->read_facets(builder);
		this->read_vertices();
        this->read_facets();

		// Correction of orientation
		//CORRECT_CGAL_STRUCTURE Structure_Corrector(&Vertices_position, &Vertices_color, &Facets);		
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

		float Max_color_value = -5000; // to decide if integer between 0 and 255 or float between 0.0 and 1.0
		

        while(fgets(pLine,512,m_pFile))
		{

        /*char * line = NULL;
        size_t len = 0;
        ssize_t read;

        while ((read = getline(&line, &len, m_pFile)) != -1)*/
		
			if ((pLine[0] == 'v') && (pLine[1] != 'n'))
			{
				

				float x,y,z,rgb1,rgb2,rgb3;

				if (sscanf(pLine,"v %f %f %f %f %f %f",&x,&y,&z, &rgb1, &rgb2, &rgb3) == 3)
				{
					vector<float> Temp_position;
					Temp_position.push_back(x);
					Temp_position.push_back(y);
					Temp_position.push_back(z);

					Vertices_position.push_back(Temp_position);
				}
				  //builder.add_vertex(Point(x,y,z));
				

				if (sscanf(pLine,"v %f %f %f %f %f %f",&x,&y,&z, &rgb1, &rgb2, &rgb3) == 6)
				{
					vector<float> Temp_position;
					Temp_position.push_back(x);
					Temp_position.push_back(y);
					Temp_position.push_back(z);

					Vertices_position.push_back(Temp_position);

					vector<float> Temp_color;
					Temp_color.push_back(rgb1);
					Temp_color.push_back(rgb2);
					Temp_color.push_back(rgb3);

					Vertices_color.push_back(Temp_color);
					
					if (rgb1 > Max_color_value)  Max_color_value = rgb1;
					if (rgb2 > Max_color_value)  Max_color_value = rgb2;
					if (rgb3 > Max_color_value)  Max_color_value = rgb3;
				}
			}
		}

		if (Max_color_value > 2.0)
			this->Is_integer = true;
		else
			this->Is_integer = false;

        /*if (line)
           free(line);*/
    }

    // read facets and uv coordinates per halfedge
    void read_facets()
    {
        fseek(m_pFile,0,SEEK_SET);

        char pLine[512];
        while(fgets(pLine,512,m_pFile))

        {
            //char *pTmp = line;
            char *pTmp = pLine;
            if (pTmp[0] == 'f')
            {
                int index,n;
                char index_ascii[512],n_ascii[512];

                // create facet
                //builder.begin_facet();

                pTmp += 2; // jump after 'f '
				if (strstr(pTmp,"//")) // MT
                {
					vector<int> Temp_facet;

                    while(sscanf(pTmp,"%d//%d",&index,&n)) // MT
                    {
                        //itoa(index,index_ascii,10);
                        //itoa(n,n_ascii,10);						

                        sprintf(index_ascii,"%d",index);
                        sprintf(n_ascii,"%d",n);
                        
						Temp_facet.push_back(index-1);
						
						//builder.add_vertex_to_facet(index-1);
                        pTmp += (2 + strlen(index_ascii) + strlen(n_ascii)); // MT
                        if (strlen(pTmp) < 3)
                            break;
                        else
                            pTmp += 1;
                    }
					Facets.push_back(Temp_facet);
                }
                else if (strstr(pTmp,"/")) // MT
                {
					vector<int> Temp_facet;

                    while(sscanf(pTmp,"%d/%d",&index,&n)) // MT
                    {
                        //itoa(index,index_ascii,10);
                        //itoa(n,n_ascii,10);						

                        sprintf(index_ascii,"%d",index);
                        sprintf(n_ascii,"%d",n);
					
						Temp_facet.push_back(index-1);

                        //builder.add_vertex_to_facet(index-1);
                        pTmp += (1 + strlen(index_ascii) + strlen(n_ascii)); // MT
                        if (strlen(pTmp) < 3)
                            break;
                        else
                            pTmp += 1;
                    }
					Facets.push_back(Temp_facet);
					
                }
                else
                {
					vector<int> Temp_facet;

                    while(sscanf(pTmp,"%d",&index))
                    {
                        //itoa(index,index_ascii,10);
                        sprintf(index_ascii,"%d",index);
                        pTmp += strlen(index_ascii);

						Temp_facet.push_back(index-1);
                        //builder.add_vertex_to_facet(index-1);
                        if (strlen(pTmp) < 3)
                            break;
                        else
                            pTmp += 1;
                    }
					Facets.push_back(Temp_facet);
                }
                //builder.end_facet();
				
            }
        }

        /*if (line)
           free(line);*/
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
				if (Is_integer)
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
				/*float rgb[3];
				rgb[0] = this->Vertices_color[i][0];
				rgb[1] = this->Vertices_color[i][1];
				rgb[2] = this->Vertices_color[i][2];
				vertex->color(rgb[0], rgb[1], rgb[2]);
			}*/
		}
		for (int i = 0 ; i < Number_faces; i++)
		{
			builder.begin_facet();
			for (unsigned int j = 0; j < Facets[i].size();j++)
			{
				int index = this->Facets[i][j];
				builder.add_vertex_to_facet(index);
			}
			builder.end_facet();
		}

		
	}
};

#endif // OBJ_CGALIMPORTER_H
