#ifndef OBJ_CGALIMPORTER_H
#define OBJ_CGALIMPORTER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cstdio>

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/formatter.hpp>
#endif

#include "Correct_CGAL_Structure.h"

#include <QFileInfo>
#include <QDir>

using namespace std;

template <class HDS>
class Builder_obj : public CGAL::Modifier_base<HDS>
{

private:

    typedef typename HDS::Vertex_handle Vertex_handle;
    typedef typename HDS::Halfedge_handle Halfedge_handle;
    typedef typename HDS::Vertex::Point Point;
    typedef typename HDS::Vertex::Point Normal;
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    string m_filename;
    string m_mtlfilename;
    FILE *m_pFile;

	vector< vector<float> > Vertices_position;
	vector< vector<float> > Vertices_color;
	vector< vector<int> > Facets;
	
	vector< vector<float> > Vertices_uvs;
	vector< vector<float> > Vertices_normals;	
	vector< vector<int> > Facets_Uvs;
	vector< vector<int> > Facets_Normals;	

	bool Is_integer;
	bool m_has_texture_coordinates;

public:

    Builder_obj(string filename)
    {
		m_has_texture_coordinates = false;
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

		this->read_file();
			
		// Correction of orientation
		CORRECT_CGAL_STRUCTURE Structure_Corrector(&Vertices_position, &Vertices_color, &Facets);		
		Facets = *(Structure_Corrector.Correct_Facet_Orientation());

		this->construct(builder);

        builder.end_surface();

		if (builder.check_unconnected_vertices())
		{
			builder.remove_unconnected_vertices();
		}

        fclose(m_pFile);		
    }
	
	bool hasTextureCoordinates()
	{
		return m_has_texture_coordinates;
	}

	QStringList read_mtl_file()
	{
		QStringList list;
		FILE * mtlFile;
		vector< string > splitLine;			
		QFileInfo info( m_filename.c_str() );		
		
		QString path= QDir(info.dir()).filePath(m_mtlfilename.c_str());
		if ( (mtlFile = fopen(path.toStdString().c_str(),"rt")) )
		{
			char pLine[512];			
			while(fgets(pLine,512,m_pFile))
			{	
				string line(pLine);
				if ( line.find("map_Kd")!=std::string::npos)
				{
					boost::trim_right(line);
					boost::split( splitLine, line, boost::is_space(), boost::token_compress_on );
					list.push_back(QDir(info.dir()).filePath(splitLine[1].c_str()));
				}				
			}			
			fclose(mtlFile);
		}

		return list;
	}

private:
	// MT
	float s2f(string str)
	{
		 istringstream buffer(str);
		 float temp;
		 buffer >> temp;
		 return temp;
	}

	// MT
	int s2i(string str)
	{
		 istringstream buffer(str);
		 int temp;
		 buffer >> temp;
		 return temp;
	}

	void read_file()
	{
		fseek(m_pFile,0,SEEK_SET);

        char pLine[512];

		float Max_color_value = -5000; 		
		this->Is_integer = false;
		vector< string > splitLine,splitFacet;
		float rgb1,rgb2,rgb3;
		while(fgets(pLine,512,m_pFile))
		{			  	
			string line(pLine);
			boost::trim_right(line);
			boost::split( splitLine, line, boost::is_space(), boost::token_compress_on );
			if (splitLine.size() >0)
			{
				if (splitLine[0] == "v")	{
					vector<float> Temp_position;
					Temp_position.push_back( s2f( splitLine[1] ) );
					Temp_position.push_back( s2f( splitLine[2] ));
					Temp_position.push_back( s2f( splitLine[3] ));
					
					if ( splitLine.size() >= 7 )
					{
						vector<float> Temp_color;
						rgb1 = s2f( splitLine[4] );
						rgb2 = s2f( splitLine[5] );
						rgb3 = s2f( splitLine[6] );

						Temp_color.push_back( rgb1 );
						Temp_color.push_back( rgb2 );
						Temp_color.push_back( rgb3 );

						if ( rgb1 > 2 && rgb2 > 2 && rgb3 > 2)
							this->Is_integer = true;

						Vertices_color.push_back( Temp_color );
					}

					Vertices_position.push_back( Temp_position );
				}
				else if (splitLine[0] == "vn") {
					vector<float> Temp_normals;
					Temp_normals.push_back( s2f( splitLine[1] ) );
					Temp_normals.push_back( s2f( splitLine[2] ) );
					Temp_normals.push_back( s2f( splitLine[3] ) );

					Vertices_normals.push_back( Temp_normals );
				}
				else if (splitLine[0] == "vt") {
					vector<float> Temp_uvs;
					Temp_uvs.push_back( s2f( splitLine[1] ) );
					Temp_uvs.push_back( s2f( splitLine[2] ) );					

					Vertices_uvs.push_back( Temp_uvs );
				}
				else if (splitLine[0] == "f") {
					vector<int> Temp_vertex;
					vector<int> Temp_uv;
					vector<int> Temp_normal;
					
					for (unsigned int i = 1; i < splitLine.size(); i++) // MT add unsigned
					{
						
						boost::split( splitFacet, splitLine[i], boost::is_any_of("/"), boost::token_compress_off );
						int facetSize = splitFacet.size();
						Temp_vertex.push_back( s2i( splitFacet[0] )-1 );

						if ( facetSize >1 && splitFacet[1] != "")
							Temp_uv.push_back( s2i( splitFacet[1] ) -1);
						
						if ( facetSize > 2 )
							Temp_normal.push_back( s2i( splitFacet[2] )-1 );
					}
					
					Facets.push_back( Temp_vertex );
					if ( Temp_uv.size() > 0)
						Facets_Uvs.push_back( Temp_uv );
					if ( Temp_normal.size() > 0)
						Facets_Normals.push_back( Temp_normal );

				}
				else if (splitLine[0] == "mtllib")
				{
					m_mtlfilename = splitLine[1];
				}
				else{
					//comment or material...
				}
			}
		}

		
	}

	void construct(Builder &builder)
	{
		int Number_vertices = (int)this->Vertices_position.size();
		int Number_faces = (int)this->Facets.size();

		vector<Vertex_handle> vertexList;
		
		for (int i = 0 ; i < Number_vertices; i++)
		{
			float coord[3];
			coord[0] = this->Vertices_position[i][0];
			coord[1] = this->Vertices_position[i][1];
			coord[2] = this->Vertices_position[i][2];

			Vertex_handle vertex = builder.add_vertex(Point(coord[0],coord[1],coord[2]));
			vertexList.push_back( vertex );

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
		}
		
		for (int i = 0 ; i < Number_faces; i++)
		{
			if ( builder.test_facet(Facets[i].begin(),Facets[i].end()) )
			{
				builder.begin_facet();

				for (unsigned int j = 0; j < Facets[i].size();j++)
				{
					int index = this->Facets[i][j];
					builder.add_vertex_to_facet(index);

					if ((unsigned int)i < Facets_Normals.size()) // MT add unsigned int
					{
						int vertexIndex =  this->Facets[i][j];
						int normalIndex =  this->Facets_Normals[i][j];
						// vertexList[vertexIndex ]->normal().x = this->Vertices_normals[normalIndex][0] ;
						// vertexList[vertexIndex ]->normal().y = this->Vertices_normals[normalIndex][1] ;
						// vertexList[vertexIndex ]->normal().z = this->Vertices_normals[normalIndex][2] ;
					}	
				}			

				Halfedge_handle he = builder.end_facet();				

				if ((unsigned int)i < Facets_Uvs.size()) // MT add unsigned int
					for (unsigned int j = 0; j < Facets[i].size();j++) {
						int vertexIndex =  this->Facets[i][j];						
						int uvIndex =  this->Facets_Uvs[i][j];						

						he->texture_coordinates( this->Vertices_uvs[uvIndex][0], this->Vertices_uvs[uvIndex][1]);
						he = he->next();

						m_has_texture_coordinates = true;
				}			
			}
			else 
				qDebug("Facet %d can't be added.",i);
		}
	}
};

#endif // OBJ_CGALIMPORTER_H
