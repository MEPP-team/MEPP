/*!
 * \file mepp_component.h
 * \brief mepp_component file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */
#ifndef MEPP_COMPONENT_H
#define MEPP_COMPONENT_H

#include <QString>

#include "./Polyhedron/polyhedron.h"

class Viewer;

/*! 
 * \class mepp_component
 * \brief mepp_component class.
 */
class mepp_component
{
	public:
		/*!
		 * \brief Constructor.
		 *
		 * \param v : Viewer pointer.
		 * \param p : Polyhedron pointer.
		 */
		mepp_component(Viewer* v, PolyhedronPtr p)
		{
			init = 0;
			componentName = "nothing";

			m_viewer = v;
			m_polyhedron_ptr = p;
		}
		
		/*!
		 * \brief Destructor.
		 */
		virtual ~mepp_component() {}

		/*!
		 * \fn Viewer* get_viewer()
		 * \brief Get the Viewer pointer.
		 *
		 * \return the Viewer pointer.
		 */
		Viewer* get_viewer() { return m_viewer; }
		/*!
		 * \fn PolyhedronPtr get_polyhedron_ptr()
		 * \brief Get the Polyhedron pointer.
		 *
		 * \return the Polyhedron pointer.
		 */
		PolyhedronPtr get_polyhedron_ptr() { return m_polyhedron_ptr; }

		/*!
		 * \fn set_remove(bool r)
		 * \brief Set the 'remove' flag for this component.
		 *
		 * \param r the 'remove' flag for this component (true or false).
		 */
		void set_remove(bool r) { remove = r; }
		/*!
		 * \fn bool get_remove()
		 * \brief Get the 'remove' flag for this component.
		 *
		 * \return the 'remove' flag for this component.
		 */
		bool get_remove() { return remove; }

		/*!
		 * \fn void set_init(int i)
		 * \brief Set the current (init) state.
		 *
		 * \param i the current (init) state.
		 */
		void set_init(int i) { init = i; }
		/*!
		 * \fn int get_init()
		 * \brief Get the current (init) state.
		 *
		 * \return the current (init) state.
		 */
		int get_init() { return init; }

		/*!
		 * \fn void set_componentName(QString n)
		 * \brief Set the component name.
		 *
		 * \param n the component name.
		 */
		void set_componentName(QString n) { componentName = n; }
		/*!
		 * \fn QString get_componentName()
		 * \brief Get the component name.
		 *
		 * \return the component name.
		 */
		QString get_componentName() { return componentName; }

	protected:
		Viewer* m_viewer;								//!< the Viewer pointer
		PolyhedronPtr m_polyhedron_ptr;					//!< the Polyhedron pointer

		QString	componentName;							//!< the component name
		int init;										//!< the current (init) state
		bool remove;									//!< the 'remove' flag for this component (true or false)
};

#endif // MEPP_COMPONENT_H
