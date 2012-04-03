#ifndef Compression_Valence_Web_ITEMS_H
#define Compression_Valence_Web_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence_Web

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class P, class Norm, class Plane>

/**
 \class	Compression_Valence_Web_Facet

 \brief	Compression valence facet. 

 */

class Compression_Valence_Web_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Compression_Valence_Web_Facet::Compression_Valence_Web_Facet()
		///
		/// \brief	Default constructor.
		///
		////////////////////////////////////////////////////////////////////////////////////////////////////

		Compression_Valence_Web_Facet() {}
	int Facet_Flag_web;
	int Component_Number_web;


};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>

/**
 \class	Compression_Valence_Web_Halfedge

 \brief	Compression valence web halfedge. 

 */

class Compression_Valence_Web_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Compression_Valence_Web_Halfedge::Compression_Valence_Web_Halfedge()
		///
		/// \brief	Default constructor.
		///
		////////////////////////////////////////////////////////////////////////////////////////////////////

		Compression_Valence_Web_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>

/**
 \class	Compression_Valence_Web_Vertex

 \brief	Compression valence web vertex. 

*/

class Compression_Valence_Web_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{	

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Compression_Valence_Web_Vertex::Compression_Valence_Web_Vertex()
		///
		/// \brief	Default constructor.
		///
		////////////////////////////////////////////////////////////////////////////////////////////////////

		Compression_Valence_Web_Vertex()
		{
			Seed_Edge_web = -1;
		}

		int m_color_int_web[3];

		int Vertex_Flag_web;
		int Vertex_Number_web;
		int Vertex_Color_Index_web;
		
		//For valence_driven
		int Vertex_Sign_web;
		int Q_Index_web;
		int Seed_Edge_web;
		int Component_Number_web;


		int color_int_web(int index) { return m_color_int_web[index]; };
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Compression_Valence_Vertex::color_int(int color0, int color1, int color2)
		///
		/// \brief	Color int.
		///
		/// \param	color0	The color 0.
		/// \param	color1	The first color.
		/// \param	color2	The second color.
		////////////////////////////////////////////////////////////////////////////////////////////////////

		void color_int_web(int color0, int color1, int color2) { m_color_int_web[0] = color0; m_color_int_web[1] = color1; m_color_int_web[2] = color2; };


};


template <class kernel, class items>

/**
 \class	Compression_Valence_Web_Polyhedron

 \brief	Compression valence web polyhedron. 

  */

class Compression_Valence_Web_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Compression_Valence_Web_Polyhedron::Compression_Valence_Web_Polyhedron()
		///
		/// \brief	Default constructor.
		///
		////////////////////////////////////////////////////////////////////////////////////////////////////

		Compression_Valence_Web_Polyhedron() {}
};

#endif

#endif // Compression_Valence_Web_ITEMS_H
