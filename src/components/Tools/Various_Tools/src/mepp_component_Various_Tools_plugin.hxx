#ifndef HEADER_MEPP_COMPONENT_VARIOUS_TOOLS_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_VARIOUS_TOOLS_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Various_Tools

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_Various_Tools_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);
#if QT_VERSION >= 0x050000
	Q_PLUGIN_METADATA(IID "mepp_component_Various_Tools_plugin")
#endif

	public:
		mepp_component_Various_Tools_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_Various_Tools_plugin()
		{
			delete actionQuadTriangle; delete actionSqrt3; delete actionSqrt3Twice;
			delete actionDooSabin; delete actionLoop; delete actionCatmullClark;
			delete actionTriangulate;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			
			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuTools;

			// début --- actions ---
			actionQuadTriangle = new QAction(tr("Quad / triangle"), this);
			if (actionQuadTriangle)
				connect(actionQuadTriangle, SIGNAL(triggered()), this, SLOT(QuadTriangle()));

			actionSqrt3 = new QAction(tr("Sqrt3"), this);
			if (actionSqrt3)
				connect(actionSqrt3, SIGNAL(triggered()), this, SLOT(Sqrt3()));

			actionSqrt3Twice = new QAction(tr("Sqrt3 (apply twice)"), this);
			if (actionSqrt3Twice)
				connect(actionSqrt3Twice, SIGNAL(triggered()), this, SLOT(Sqrt3Twice()));

			actionDooSabin = new QAction(tr("Doo-Sabin"), this);
			if (actionDooSabin)
				connect(actionDooSabin, SIGNAL(triggered()), this, SLOT(DooSabin()));

			actionLoop = new QAction(tr("Loop"), this);
			if (actionLoop)
				connect(actionLoop, SIGNAL(triggered()), this, SLOT(Loop()));

			actionCatmullClark = new QAction(tr("Catmull-Clark"), this);
			if (actionCatmullClark)
				connect(actionCatmullClark, SIGNAL(triggered()), this, SLOT(CatmullClark()));

			actionTriangulate = new QAction(tr("Triangulate"), this);
			if (actionTriangulate)
				connect(actionTriangulate, SIGNAL(triggered()), this, SLOT(Triangulate()));
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionQuadTriangle
										<< actionSqrt3
										<< actionSqrt3Twice
										<< actionDooSabin
										<< actionLoop
										<< actionCatmullClark
										<< NULL			// menu separator
										<< actionTriangulate;
		}	
		
		virtual void pre_draw() {}
		virtual void post_draw() {}
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene() {}

		virtual void OnMouseLeftDown(QMouseEvent *event) {}
		virtual void OnMouseLeftUp(QMouseEvent *event) {}
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event) {}
		virtual void OnMouseMotion(QMouseEvent *event) {}
		virtual void OnMouseWheel(QWheelEvent *event) {}
		virtual void OnKeyPress(QKeyEvent *event) {}
		virtual void OnKeyRelease(QKeyEvent *event) {}

	public slots:
		void QuadTriangle();
		void Sqrt3();
		void Sqrt3Twice();
		void DooSabin();
		void Loop();
		void CatmullClark();
		void Triangulate();

	private:
		QAction *actionQuadTriangle, *actionSqrt3, *actionSqrt3Twice, *actionDooSabin, *actionLoop, *actionCatmullClark;
		QAction *actionTriangulate;
};

#endif

#endif
