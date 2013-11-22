#ifndef ASSIMP_CGALIMPORTER_H
#define ASSIMP_CGALIMPORTER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cstdio>
#include "mesh.h"
#include "Correct_CGAL_Structure.h"

// Can only read triangular meshes 

// assimp include files
#include <assimp/Importer.hpp>	// OO version Header!
#include <assimp/postprocess.h>
#include <assimp/scene.h>

using namespace std;

template <class HDS>
class Builder_dae : public CGAL::Modifier_base<HDS>
{

private:

	typedef typename HDS::Vertex_handle Vertex_handle;
	typedef typename HDS::Vertex::Point Point;
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
	string m_filename;

	FILE *m_pFile;
	vector< vector<float> > Vertices_position;
	vector< vector<float> > Vertices_color;
	vector< vector<float> > Vertices_texture_coordinates;
	vector< vector<int> >   Facets;
	bool m_has_texture_coordinates;

	aiMesh* m_Mesh;	

public:
	bool loadedSucess;

	Builder_dae(aiMesh* mesh)
	{
		loadedSucess = false;
		m_Mesh = mesh;
		m_has_texture_coordinates = false;		
	}

	bool is_ok()
	{
		return (m_pFile != NULL);
	}

	~Builder_dae() {}

	void operator()(HDS& hds)
	{
		Builder builder(hds, true);
				
		builder.begin_surface(3, 1, 6);
		this->construct(builder);
		
		if (!this->loadedSucess)
			return;

		builder.end_surface();

		if (builder.check_unconnected_vertices())
		{
			builder.remove_unconnected_vertices();
		}
	}

	bool hasTextureCoordinates()
	{
		return m_has_texture_coordinates;
	}

private:
	/*bool initialize()
	{
		scene = (aiScene*) importer.ReadFile(m_filename, aiProcessPreset_TargetRealtime_Fast);

		// If the import failed, report it
		if (!scene) {
			const char * errorStr = importer.GetErrorString();
			printf("%s\n",errorStr);
			return false;
		}

		return true;
	}*/

	void construct(Builder &builder)
	{	
		for (unsigned int i = 0; i < m_Mesh->mNumVertices; i++)
		{
			Vertex_handle vertex = builder.add_vertex(Point(m_Mesh->mVertices[i].x, m_Mesh->mVertices[i].y, m_Mesh->mVertices[i].z));

			if (m_Mesh->mColors[0] != NULL)				
				vertex->color(m_Mesh->mColors[0][i].r, m_Mesh->mColors[0][i].g, m_Mesh->mColors[0][i].b);

			if (m_Mesh->mTextureCoords[0] != NULL)
			{
				m_has_texture_coordinates = true;
				vertex->texture_coordinates(m_Mesh->mTextureCoords[0][i].x, m_Mesh->mTextureCoords[0][i].y);
			}
		}

		for (unsigned int i = 0; i < m_Mesh->mNumFaces; i++)
		{
			builder.begin_facet();			
			for (unsigned int j = 0; j < m_Mesh->mFaces[i].mNumIndices; j++)
			{
				int index = m_Mesh->mFaces[i].mIndices[j];
				builder.add_vertex_to_facet(index);
				if (builder.error())
					return;
			}
			builder.end_facet();
		}

		loadedSucess = true;
	}
};

#endif // ASSIMP_CGALIMPORTER_H
