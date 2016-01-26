/*!
 * \file mepp_component_plugin_interface.h
 * \brief mepp_component_plugin_interface file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */
#ifndef HEADER_MEPP_COMPONENT_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_PLUGIN_INTERFACE_H

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829

#include <QList>
#include <QtPlugin>

class QAction;
class QMdiSubWindow;

#include "mainwindow.hxx"

#include "mepp_component.h"

#endif

typedef boost::shared_ptr<mepp_component> mepp_componentPtr;

using namespace qglviewer; // for draw_link

/*! 
 * \class SleeperThread
 * \brief SleeperThread class.
 */
class SleeperThread : public QThread
{
	public:
		/*!
		 * \fn static void msleep(unsigned long msecs)
		 * \brief Sleep function.
		 *
		 * \param msecs sleeping time in ms.
		 */
		static void msleep(unsigned long msecs)
		{
			QThread::msleep(msecs);
		}
};

/*! 
 * \class mepp_component_plugin_interface
 * \brief mepp_component_plugin_interface class.
 */
class mepp_component_plugin_interface 
{
	public:
		/*!
		 * \brief Constructor.
		 */
		mepp_component_plugin_interface()
		{
			mMenu = NULL;
			mParentMenu = NULL;
			mPluginName = "Empty";
		}
		/*!
		 * \brief Destructor.
		 */
		virtual ~mepp_component_plugin_interface() {}

		/*!
		 * \fn QString getPluginName()
		 * \brief Get the plugin name.
		 *
		 * \return the plugin name.
		 */
		QString getPluginName() { return mPluginName; }

		QToolBar* mToolBar;		//!< the QToolBar pointer
		QMenu* mMenu;			//!< the QMenu pointer

		QMenu* mParentMenu;		//!< the QMenu parent pointer
								//!< choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
								//!< menuCompression, menuWatermaking, menuExamples

		/*!
		 * \fn QAction* getActionFromMainWindow(mainwindow* mw, QString action_name)
		 * \brief Get action object from its name.
		 *
		 * \param mw the mainwindow pointer.
 		 * \param action_name the action name.
		 *
		 * \return the action.
		 */
		QAction* getActionFromMainWindow(mainwindow* mw, QString action_name)
		{
			return mw->findChild<QAction*>(action_name);
		}

		/*!
		 * \fn void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		 * \brief Init the plugin.
		 *
		 * \param mainWindow the mainwindow pointer.
 		 * \param lw the subwindows list.
		 */
        virtual void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
        {
            this->mw = mainWindow;
            this->lwindow = lw;
        }

		/*!
		 * \fn QList<QAction*> actions()
		 * \brief Get action object from its name.
		 *
		 * \return the actions list.
		 */
		virtual QList<QAction*> actions() const = 0;

		/*!
		 * \fn void setListWindow(QList<QMdiSubWindow *> lw)
		 * \brief Set the subwindows list.
		 *
		 * \param lw the subwindows list.
		 */
		void setListWindow(QList<QMdiSubWindow *> lw)
		{
			this->lwindow = lw;

			cleanComponentForViewer();
		}

		/*!
		 * \fn bool doesExistComponentForViewer(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		 * \brief Return true if a component exists for this viewer. If not, return false.
		 *
		 * \param viewer : Viewer pointer.
		 * \param polyhedron_ptr : Polyhedron pointer.
		 *
		 * \return true or false.
		 */
		template <class COMPONENT_PTR, class COMPONENT>
		bool doesExistComponentForViewer(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		{
			COMPONENT_PTR component_ptr;

			bool find = false;
			for (int c=0; c<lcomponent.size(); c++)
			{
				Viewer* vc = lcomponent[c]->get_viewer();
				if (viewer == vc)
				{
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					if (component_ptr)
					{
						find = true;
						break;
					}
				}
			}

			return find;
		}

		/*!
		 * \fn COMPONENT_PTR findOrCreateComponentForViewer(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		 * \brief Return a component pointer for this viewer. Create one if no one exists.
		 *
		 * \param viewer : Viewer pointer.
		 * \param polyhedron_ptr : Polyhedron pointer.
		 *
		 * \return the component pointer.
		 */
		template <class COMPONENT_PTR, class COMPONENT>
		COMPONENT_PTR findOrCreateComponentForViewer(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		{
			COMPONENT_PTR component_ptr;

			bool find = false;
			for (int c=0; c<lcomponent.size(); c++)
			{
				Viewer* vc = lcomponent[c]->get_viewer();
				if (viewer == vc)
				{
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					if (component_ptr)
					{
						find = true;
						break;
					}
				}
			}

			if (!find)
			{
				component_ptr = COMPONENT_PTR(new COMPONENT(viewer, polyhedron_ptr));
				lcomponent << component_ptr;
			}

			return component_ptr;
		}

		/*!
		 * \fn bool doesExistComponentForPolyhedron(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		 * \brief Return true if a component exists for this polyhedron. If not, return false.
		 *
		 * \param viewer : Viewer pointer.
		 * \param polyhedron_ptr : Polyhedron pointer.
		 *
		 * \return true or false.
		 */
		template <class COMPONENT_PTR, class COMPONENT>
		bool doesExistComponentForPolyhedron(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		{
			COMPONENT_PTR component_ptr;

			bool find = false;
			for (int c=0; c<lcomponent.size(); c++)
			{
				PolyhedronPtr p = lcomponent[c]->get_polyhedron_ptr();
				if (polyhedron_ptr == p)
				{
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					if (component_ptr)
					{
						find = true;
						break;
					}
				}
			}

			return find;
		}

		/*!
		 * \fn COMPONENT_PTR findOrCreateComponentForPolyhedron(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		 * \brief Return a component pointer for this polyhedron. Create one if no one exists.
		 *
		 * \param viewer : Viewer pointer.
		 * \param polyhedron_ptr : Polyhedron pointer.
		 *
		 * \return the component pointer.
		 */
		template <class COMPONENT_PTR, class COMPONENT>
		COMPONENT_PTR findOrCreateComponentForPolyhedron(Viewer* viewer, PolyhedronPtr polyhedron_ptr)
		{
			COMPONENT_PTR component_ptr;

			bool find = false;
			for (int c=0; c<lcomponent.size(); c++)
			{
				PolyhedronPtr p = lcomponent[c]->get_polyhedron_ptr();
				if (polyhedron_ptr == p)
				{
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					if (component_ptr)
					{
						find = true;
						break;
					}
				}
			}

			if (!find)
			{
				component_ptr = COMPONENT_PTR(new COMPONENT(viewer, polyhedron_ptr));
				lcomponent << component_ptr;
			}

			return component_ptr;
		}

		/*!
		 * \fn void pre_draw()
		 * \brief Function called just before the draw.
		 */
		virtual void pre_draw() {}
		/*!
		 * \fn void post_draw()
		 * \brief Function called just after the draw.
		 */
		virtual void post_draw() {}

		/*!
		 * \fn void pre_draw_all_scene()
		 * \brief Function called just before all drawings (only in Space mode).
		 */
		virtual void pre_draw_all_scene() {}
		/*!
		 * \fn void post_draw_all_scene()
		 * \brief Function called just after all drawings (only in Space mode).
		 */
		virtual void post_draw_all_scene() {}

		/*!
		 * \fn void draw_link(Viewer* viewer, int frame_i, int frame_j, Vec p0, Vec p1)
		 * \brief Draw a link between 2 points of 2 frames in the viewer.
		 *
		 * \param viewer the Viewer pointer.
		 * \param frame_i first frame.
 		 * \param frame_j second frame.
 		 * \param p0 first point.
 		 * \param p1 second point.
		 */
		void draw_link(Viewer* viewer, int frame_i, int frame_j, Vec p0, Vec p1)
		{
			Vec vi = viewer->frame(frame_i)->position();
			Quaternion qi = viewer->frame(frame_i)->orientation();
			
			Vec vj = viewer->frame(frame_j)->position();
			Quaternion qj = viewer->frame(frame_j)->orientation();
			
			Vec vti, vtj;

			glBegin(GL_LINES);
				vti = qi*p0+vi;
				glVertex3f(vti.x, vti.y, vti.z);
				
				vtj = qj*p1+vj;
				glVertex3f(vtj.x, vtj.y, vtj.z);
			glEnd();
		}

		/*!
		 * \fn static int load_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer)
		 * \brief Function called when a mesh is loading from a specific component.
		 *
		 * \param polyhedron_ptr : the Polyhedron pointer.
		 * \param filename : the filename of the mesh.
		 * \param viewer the Viewer pointer.
		 */
		static int load_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer);
		/*!
		 * \fn static int save_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer)
		 * \brief Function called when a mesh is saving from a specific component.
		 *
		 * \param polyhedron_ptr : the Polyhedron pointer.
		 * \param filename : the filename of the mesh.
		 * \param viewer the Viewer pointer.
		 */
		static int save_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer);

		/*!
		 * \fn void OnMouseLeftDown(QMouseEvent *event)
		 * \brief Function called on mouse left down clic.
		 *
		 * \param event the QMouseEvent event.
		 */
		virtual void OnMouseLeftDown(QMouseEvent *event) {}
		/*!
		 * \fn void OnMouseLeftUp(QMouseEvent *event)
		 * \brief Function called on mouse left up clic.
		 *
		 * \param event the QMouseEvent event.
		 */
		virtual void OnMouseLeftUp(QMouseEvent *event) {}

		/*!
		 * \fn void OnMouseRightDown(QMouseEvent *event)
		 * \brief Function called on mouse right down clic.
		 *
		 * \param event the QMouseEvent event.
		 */
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		/*!
		 * \fn void OnMouseRightUp(QMouseEvent *event)
		 * \brief Function called on mouse right up clic.
		 *
		 * \param event the QMouseEvent event.
		 */
		virtual void OnMouseRightUp(QMouseEvent *event) {}

		/*!
		 * \fn void OnMouseMotion(QMouseEvent *event)
		 * \brief Function called on mouse motion.
		 *
		 * \param event the QMouseEvent event.
		 */
		virtual void OnMouseMotion(QMouseEvent *event) {}

		/*!
		 * \fn void OnMouseWheel(QWheelEvent *event)
		 * \brief Function called when wheel is moving.
		 *
		 * \param event the QWheelEvent event.
		 */
		virtual void OnMouseWheel(QWheelEvent *event) {}

		/*!
		 * \fn void OnKeyPress(QKeyEvent *event)
		 * \brief Function called on key press down.
		 *
		 * \param event the QKeyEvent event.
		 */
		virtual void OnKeyPress(QKeyEvent *event) {}
		/*!
		 * \fn void OnKeyRelease(QKeyEvent *event)
		 * \brief Function called on key press release.
		 *
		 * \param event the QKeyEvent event.
		 */
		virtual void OnKeyRelease(QKeyEvent *event) {}

	private:
		void cleanComponentForViewer()
		{
			QList<mepp_componentPtr>::iterator c;
			for (c = lcomponent.begin(); c < lcomponent.end(); c++)
				(*c)->set_remove(true);

			//for (c = lcomponent.begin(); c < lcomponent.end(); c++)
			{
				for (int i=0; i<lwindow.size(); i++)
				{
					Viewer* viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());

					for (c = lcomponent.begin(); c < lcomponent.end(); c++)
					{
						if ((*c)->get_viewer() == viewer)
						{
							(*c)->set_remove(false);
							//break;
						}
					}
				}
			}

			for (c = lcomponent.begin(); c < lcomponent.end(); c++)
			{
				if ((*c)->get_remove())
					lcomponent.erase(c);
			}
		}

	protected:
		mainwindow* mw;							//!< the mainwindow pointer

		QList<QMdiSubWindow *> lwindow;			//!< the subwindows list
		QList<mepp_componentPtr> lcomponent;	//!< the mepp_component pointer list

		QString mPluginName;					//!< the plugin name
};

Q_DECLARE_INTERFACE(mepp_component_plugin_interface,
                    "fr.liris.MEPP.PluginInterface/1.0")

#endif
