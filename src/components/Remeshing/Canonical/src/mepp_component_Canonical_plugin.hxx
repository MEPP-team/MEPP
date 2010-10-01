#ifndef HEADER_MEPP_COMPONENT_CANONICAL_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_CANONICAL_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Canonical

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_Canonical_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:
		mepp_component_Canonical_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_Canonical_plugin()
		{
			delete actionSimplification;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			
			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuRemeshing_Subdivision;

			// début --- actions ---
			actionSimplification = new QAction(tr("Valence Driven Simplification"), this);
			if (actionSimplification)
				connect(actionSimplification, SIGNAL(triggered()), this, SLOT(ValenceDrivenSimplification()));
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionSimplification;
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
		void ValenceDrivenSimplification();

	private:
		QAction *actionSimplification;
};

#endif

#endif
