#ifndef Compression_Valence_ITEMS_H
#define Compression_Valence_ITEMS_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "../../../../mepp/Polyhedron/polyhedron_shared_items.h"

template <class Refs, class T, class Norm>
class Compression_Valence_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

		Compression_Valence_Facet() {}
		int Facet_Flag;
		//int Patch_Index;
		

};


template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Compression_Valence_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

		Compression_Valence_Halfedge() {}
};


template <class Refs, class T, class P, class Norm>
class Compression_Valence_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{	
	protected:
		

	public:

		Compression_Valence_Vertex() {}

		//float color_float(int index) { return m_color_float[index]; };
		//void color_float(float color0, float color1, float color2) {m_color_float[0] = color0; m_color_float[1] = color1; m_color_float[2] = color2;};
		
		int color_int(int index) { return m_color_int[index]; };
		void color_int(int color0, int color1, int color2) { m_color_int[0] = color0; m_color_int[1] = color1; m_color_int[2] = color2; };
		
		
		//float m_color_float[3];
		int m_color_int[3];
		
		int Vertex_Flag;
		int Vertex_Number;
		int Vertex_Color_Index;
		
		//For valence_driven
		int Vertex_Sign;
		int Q_Index;
		int Seed_Edge;

		int Component_Number;

		
};


template <class kernel, class items>
class Compression_Valence_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

		Compression_Valence_Polyhedron() {}
};

#endif

#endif // Compression_Valence_ITEMS_H
