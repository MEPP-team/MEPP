///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_MEPP_COMPONENT_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_PLUGIN_INTERFACE_H

#include <QList>
#include <QtPlugin>

class QAction;
class QMdiSubWindow;

#include "mainwindow.hxx"

#include "mepp_component.h"
typedef boost::shared_ptr<mepp_component> mepp_componentPtr;

using namespace qglviewer; // for draw_link

class SleeperThread : public QThread
{
	public:
		static void msleep(unsigned long msecs)
		{
			QThread::msleep(msecs);
		}
};

class mepp_component_plugin_interface 
{
	public:
		mepp_component_plugin_interface()
		{
			mMenu = NULL;
			mParentMenu = NULL;
			mPluginName = "Empty";
		}
		virtual ~mepp_component_plugin_interface() {}

		QString getPluginName() { return mPluginName; }
		QToolBar* mToolBar;
		QMenu* mMenu;

		QMenu* mParentMenu; // choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
							// menuCompression, menuWatermaking, menuExamples

		// get action object from its name
		QAction* getActionFromMainWindow(mainwindow* mw, QString action_name)
		{
			return mw->findChild<QAction*>(action_name);
		}

        virtual void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
        {
            this->mw = mainWindow;
            this->lwindow = lw;
        }

		virtual QList<QAction*> actions() const = 0;

		void setListWindow(QList<QMdiSubWindow *> lw)
		{
			this->lwindow = lw;

			cleanComponentForViewer();
		}

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
					find = true;
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					break;
				}
			}

			return find;
		}

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
					find = true;
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					break;
				}
			}

			if (!find)
			{
				component_ptr = COMPONENT_PTR(new COMPONENT(viewer, polyhedron_ptr));
				lcomponent << component_ptr;
			}

			return component_ptr;
		}

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
					find = true;
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					break;
				}
			}

			return find;
		}

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
					find = true;
					component_ptr = boost::dynamic_pointer_cast<COMPONENT>(lcomponent[c]);
					break;
				}
			}

			if (!find)
			{
				component_ptr = COMPONENT_PTR(new COMPONENT(viewer, polyhedron_ptr));
				lcomponent << component_ptr;
			}

			return component_ptr;
		}

		virtual void pre_draw() {}
		virtual void post_draw() {}
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene() {}

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

		static int load_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer);
		static int save_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer);

		virtual void OnMouseLeftDown(QMouseEvent *event) {}
		virtual void OnMouseLeftUp(QMouseEvent *event) {}
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event) {}
		virtual void OnMouseMotion(QMouseEvent *event) {}
		virtual void OnMouseWheel(QWheelEvent *event) {}
		virtual void OnKeyPress(QKeyEvent *event) {}
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
		mainwindow* mw;

		QList<QMdiSubWindow *> lwindow;
		QList<mepp_componentPtr> lcomponent;

		QString mPluginName;
};

Q_DECLARE_INTERFACE(mepp_component_plugin_interface,
                    "fr.liris.MEPP.PluginInterface/1.0")

#endif
