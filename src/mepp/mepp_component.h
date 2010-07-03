///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef MEPP_COMPONENT_H
#define MEPP_COMPONENT_H

#include <QString>

#include "./Polyhedron/polyhedron.h"

class Viewer;

class mepp_component
{
	public:
		mepp_component(Viewer* v, PolyhedronPtr p)
		{
			init = 0;
			componentName = "nothing";

			m_viewer = v;
			m_polyhedron_ptr = p;
		}
		
		virtual ~mepp_component() {}

		Viewer* get_viewer() { return m_viewer; }
		PolyhedronPtr get_polyhedron_ptr() { return m_polyhedron_ptr; }

		void set_remove(bool r) { remove = r; }
		bool get_remove() { return remove; }

		void set_init(int i) { init = i; }
		int get_init() { return init; }

		void set_componentName(QString n) { componentName = n; }
		QString get_componentName() { return componentName; }

	protected:
		Viewer* m_viewer;
		PolyhedronPtr m_polyhedron_ptr;

		QString	componentName;
		int init;
		bool remove;
};

#endif // MEPP_COMPONENT_H
