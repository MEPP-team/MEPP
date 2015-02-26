///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////
#ifndef HEADER_MEPP_COMPONENT_CGAL_EXAMPLE_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_CGAL_EXAMPLE_PLUGIN_INTERFACE_H

#include <QtGlobal> // important, for QT_VERSION

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_CGAL_Example_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);
#if QT_VERSION >= 0x050000
	Q_PLUGIN_METADATA(IID "mepp_component_CGAL_Example_plugin")
#endif

	public:
		mepp_component_CGAL_Example_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_CGAL_Example_plugin()
		{
			delete actionStep_1; delete actionStep_2; delete actionStep_3;
			delete actionStep_4; delete actionStep_5; delete actionStep_6;
			delete actionStep_7; delete actionStep_8; delete actionStep_9;
			delete actionStep_10;
			delete actionStep_11;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			
			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuExamples;

			// début --- actions ---
			/*actionExample = new QAction(tr("Action example"), this);
			if (actionExample)
				connect(actionExample, SIGNAL(triggered()), this, SLOT(example()));*/

			actionStep_1 = new QAction(tr("Triangulate And Random Color Facets"), this);
			if (actionStep_1)
				connect(actionStep_1, SIGNAL(triggered()), this, SLOT(step1()));

			actionStep_2 = new QAction(tr("Create Center Vertex"), this);
			if (actionStep_2)
				connect(actionStep_2, SIGNAL(triggered()), this, SLOT(step2()));

			actionStep_3 = new QAction(tr("Show Black And White Facets"), this);
			if (actionStep_3)
				connect(actionStep_3, SIGNAL(triggered()), this, SLOT(step3()));

			actionStep_4 = new QAction(tr("Draw Connections"), this);
			if (actionStep_4)
				connect(actionStep_4, SIGNAL(triggered()), this, SLOT(step4()));

			actionStep_5 = new QAction(tr("Set Position And Orientation"), this);
			if (actionStep_5)
				connect(actionStep_5, SIGNAL(triggered()), this, SLOT(step5()));

			actionStep_6 = new QAction(tr("New/Add Polyhedron"), this);
			if (actionStep_6)
				connect(actionStep_6, SIGNAL(triggered()), this, SLOT(step6()));

			actionStep_7 = new QAction(tr("Load File From Component"), this);
			if (actionStep_7)
				connect(actionStep_7, SIGNAL(triggered()), this, SLOT(step7()));

			actionStep_8 = new QAction(tr("Save File From Component"), this);
			if (actionStep_8)
				connect(actionStep_8, SIGNAL(triggered()), this, SLOT(step8()));

			actionStep_9 = new QAction(tr("Sample to use Curvature component from this component"), this);
			if (actionStep_9)
				connect(actionStep_9, SIGNAL(triggered()), this, SLOT(step9()));

			actionStep_10 = new QAction(tr("Sample to load texture file from this component"), this);
			if (actionStep_10)
				connect(actionStep_10, SIGNAL(triggered()), this, SLOT(step10()));

			actionStep_11 = new QAction(tr("Init for picking vertex"), this);
			if (actionStep_11)
				connect(actionStep_11, SIGNAL(triggered()), this, SLOT(step11()));
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	//<< actionExample
										<< actionStep_1
										<< actionStep_2
										<< actionStep_3
										<< NULL			// menu separator
										<< actionStep_4
										<< actionStep_5
										<< NULL			// menu separator
										<< actionStep_6
										<< NULL			// menu separator
										<< actionStep_7
										<< actionStep_8
										<< NULL			// menu separator
										<< actionStep_9
										<< NULL			// menu separator
										<< actionStep_10
										<< NULL			// menu separator
										<< actionStep_11;
		}	
		
		virtual void pre_draw();
		virtual void post_draw();
		virtual void pre_draw_all_scene();
		virtual void post_draw_all_scene();

		virtual void OnMouseLeftDown(QMouseEvent *event);
		virtual void OnMouseLeftUp(QMouseEvent *event);
		virtual void OnMouseRightDown(QMouseEvent *event);
		virtual void OnMouseRightUp(QMouseEvent *event);
		virtual void OnMouseMotion(QMouseEvent *event);
		virtual void OnMouseWheel(QWheelEvent *event);
		virtual void OnKeyPress(QKeyEvent *event);
		virtual void OnKeyRelease(QKeyEvent *event);

	private:

		void draw_connections(Viewer* viewer, int frame_i, int frame_j);

	public slots:
		void example();

		void step1();
		void step2();
		void step3();
		void step4();
		void step5();
		void step6();
		void step7();
		void step8();
		void step9();
		void step10();
		void step11();

	private:
		QAction *actionExample;
		QAction *actionStep_1, *actionStep_2, *actionStep_3;
		QAction *actionStep_4, *actionStep_5;
		QAction *actionStep_6;
		QAction *actionStep_7, *actionStep_8;
		QAction *actionStep_9;
		QAction *actionStep_10;
		QAction *actionStep_11;
};

#endif

#endif
