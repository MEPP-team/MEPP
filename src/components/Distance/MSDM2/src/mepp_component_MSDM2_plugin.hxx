/*!
	\file MSDM2_Component.h
	\brief MSDM2 Perceptual distance calculation
	\brief According to: Multiscale Metric for 3D Mesh Visual Quality Assessment
			G. Lavoue
			In Symposium on Geometry Processing (SGP) 2011
	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2011
 */

#ifndef HEADER_MEPP_COMPONENT_MSDM2_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_MSDM2_PLUGIN_INTERFACE_H

#include <QtGlobal> // important, for QT_VERSION

#include <QObject>

#include <mepp_config.h>
//#ifdef BUILD_component_MSDM2

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_MSDM2_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);
#if QT_VERSION >= 0x050000
	Q_PLUGIN_METADATA(IID "mepp_component_MSDM2_plugin")
#endif

	public:
		mepp_component_MSDM2_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_MSDM2_plugin()
		{
			delete action_distance; delete action_display;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			
			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuDistance_Quality_measure;

			// début --- actions ---
			/*actionExample = new QAction(tr("Action example"), this);
			if (actionExample)
				connect(actionExample, SIGNAL(triggered()), this, SLOT(example()));*/

			action_distance = new QAction(tr("Compute MSDM2 distance"), this);
			if (action_distance)
				connect(action_distance, SIGNAL(triggered()), this, SLOT(MSDM2_computation()));

			action_display = new QAction(tr("Local MSDM2 to Color Map "), this);
			if (action_display)
				connect(action_display, SIGNAL(triggered()), this, SLOT(DistanceToColorMap()));

			
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	//<< actionExample
										<< action_distance
										<< action_display
									
										;
		}	
		
		virtual void pre_draw() {}
		virtual void post_draw() {}
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene()  {}

		virtual void OnMouseLeftDown(QMouseEvent *event) {}
		virtual void OnMouseLeftUp(QMouseEvent *event) {}
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event) {}
		virtual void OnMouseMotion(QMouseEvent *event) {}
		virtual void OnMouseWheel(QWheelEvent *event) {}
		virtual void OnKeyPress(QKeyEvent *event) {}
		virtual void OnKeyRelease(QKeyEvent *event) {}

	private:

		void draw_connections(Viewer* viewer, int frame_i, int frame_j);

	public slots:
		void MSDM2_computation();
		void DistanceToColorMap();
	

	private:
		QAction *action_distance;
		QAction *action_display;

};

#endif

//#endif
