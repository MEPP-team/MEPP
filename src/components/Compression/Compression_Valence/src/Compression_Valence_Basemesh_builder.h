#ifndef Compression_Valence_BASEMESH_BUILDER_H
#define Compression_Valence_BASEMESH_BUILDER_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "../../../../mepp/Polyhedron/Polyhedron_Builder.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/circulator.h>
#include <CGAL/basic.h>
#include <iostream>

/**
 \property	template <class HDS,class Polyhedron,
 			class kernel> class CModifyBasemeshBuilder : public CGAL::Modifier_base<HDS>

 \brief	Gets the modifier base.

 \value	.
 */

template <class HDS,class Polyhedron,class kernel>

/**
 \class	CModifyBasemeshBuilder

 \brief	Modify basemesh builder. 

 */

class CModifyBasemeshBuilder : public CGAL::Modifier_base<HDS>
{
	// Access the Color_index_list

	/**
	 \fn	vector<int>* :::GetColor_index_list(void) const
	
	 \brief	Gets the color index list.
	
	 \return	null if it fails, else the color index list.
	 */

	vector<int>* GetColor_index_list(void) const				{ return(color_index_list);				};

	/**
	 \fn	void :::SetColor_index_list(vector<int>* _color_index_list)
	
	 \brief	Sets a color index list.
	
	 \param [in,out]	_color_index_list	If non-null, list of color indexes.
	 */

	void SetColor_index_list(vector<int>* _color_index_list)	{ color_index_list = _color_index_list;	};


private:

	/**
	 \typedef	typename HDS::Vertex Vertex
	
	 \brief	Defines an alias representing the vertex.
	 		///////////////////////////////////////////////////////////////////////////////////////////
	 		//////.
	 */

	typedef typename HDS::Vertex          Vertex;

	/**
	 \typedef	typename Vertex::Point Point
	
	 \brief	Defines an alias representing the point.
	 */

	typedef typename Vertex::Point        Point;

	/**
	 \typedef	typename CGAL::Enriched_polyhedron_incremental_builder_3<HDS> builder
	
	 \brief	Defines an alias representing the builder.
	 */

	typedef typename CGAL::Enriched_polyhedron_incremental_builder_3<HDS> builder;
	
	//private fileds :
	vector<Point> * Vertex_List;
	vector<int>   * Facet_List;
	vector<float> * Color_List;
	vector<int>   * Color_index_list;
	 

public:

	// life cycle

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	:::CModifyBasemeshBuilder(vector<Point> *_vertex_list, vector<int> *_facet_list,
	/// 	vector<float> *_color_list, vector<int> *_color_index_list)
	///
	/// \brief	Constructor.
	///
	/// \param [in,out]	_vertex_list	 	If non-null, list of vertices.
	/// \param [in,out]	_facet_list		 	If non-null, list of facets.
	/// \param [in,out]	_color_list		 	If non-null, list of colors.
	/// \param [in,out]	_color_index_list	If non-null, list of color indexes.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	CModifyBasemeshBuilder(vector<Point> *_vertex_list, vector<int> *_facet_list, vector<float> *_color_list, vector<int> *_color_index_list)
	{
		CGAL_assertion(_vertex_list->size() != 0);
		CGAL_assertion(_facet_list->size() != 0);
		Vertex_List = _vertex_list;
		Facet_List = _facet_list;
		Color_List = _color_list;
		Color_index_list = _color_index_list;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	:::~CModifyBasemeshBuilder()
	///
	/// \brief	Finaliser.
	///
	////////////////////////////////////////////////////////////////////////////////////////////////////

	~CModifyBasemeshBuilder() {}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	void :::operator()( HDS& hds)
	///
	/// \brief	Default constructor.
	///
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void operator()( HDS& hds)
	{
		builder B(hds,true);        //builder init
		B.begin_surface(3,1,6);     //initialisation of the new polyhedron
			add_vertices(B);        //create all vertices for the new polyhedron
			add_facets(B);
		B.end_surface();
	}

	// add vertices

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	void :::add_vertices(builder &B)
	///
	/// \brief	Adds a vertices. 
	///
	/// \param [in,out]	B	The.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void add_vertices(builder &B)
	{
		///< Number of vertices
		size_t Number_vertex = Vertex_List->size();
		for (int i = 0; i < (int)Number_vertex ; i++)
		{
			Point Pt = (*Vertex_List)[i];

			Vertex_handle v = B.add_vertex(Pt);
			
			if (!Color_List->empty())
			{
				float Temp_color[3];
				Temp_color[0] = (*Color_List)[3*i + 0];
				Temp_color[1] = (*Color_List)[3*i + 1];
				Temp_color[2] = (*Color_List)[3*i + 2];
			
				v->color(Temp_color[0], Temp_color[1], Temp_color[2]);
			}
			if (!Color_index_list->empty())
			{
				v->Vertex_Color_Index = (*Color_index_list)[i];
			}
			
		}		
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	void :::add_facets(builder &B)
	///
	/// \brief	Adds a facets. 
	///
	///
	/// \param [in,out]	B	The.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void add_facets(builder &B)
	{
		int Number_facet = (int)Facet_List->size() / 3;

		for (int i = 0; i < Number_facet ; i++)
		{
			B.begin_facet();
			int index_vertex1 = (*Facet_List)[3*i+0];
			int index_vertex2 = (*Facet_List)[3*i+1];
			int index_vertex3 = (*Facet_List)[3*i+2];
			B.add_vertex_to_facet(index_vertex1);
			B.add_vertex_to_facet(index_vertex2);
			B.add_vertex_to_facet(index_vertex3);
			B.end_facet();
		}		
		
		CGAL_assertion(!B.check_unconnected_vertices());
	}
};

#endif

#endif 
