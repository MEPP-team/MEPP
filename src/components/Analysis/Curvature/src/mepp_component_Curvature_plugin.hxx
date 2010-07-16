#ifndef HEADER_MEPP_COMPONENT_CURVATURE_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_CURVATURE_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_Curvature_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:
		mepp_component_Curvature_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_Curvature_plugin()
		{
			delete actionCalculate_Normal_Cycle_Algorithm; delete actionMinimum_Curvature_to_Color_Map; delete actionMaximum_Curvature_to_Color_Map;
			delete actionDisplay_Minimum_Curvature_direction; delete actionDisplay_Maximum_Curvature_direction;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();

			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuAnalysis_Filtering;

			// début --- actions ---
			actionCalculate_Normal_Cycle_Algorithm = new QAction(tr("Calculate (Normal Cycle Algorithm)"), this);
			if (actionCalculate_Normal_Cycle_Algorithm)
				connect(actionCalculate_Normal_Cycle_Algorithm, SIGNAL(triggered()), this, SLOT(OnCurvature()));

			actionMinimum_Curvature_to_Color_Map= new QAction(tr("Minimum Curvature to Color Map"), this);
			if (actionMinimum_Curvature_to_Color_Map)
				connect(actionMinimum_Curvature_to_Color_Map, SIGNAL(triggered()), this, SLOT(OnDisplayMin()));

			actionMaximum_Curvature_to_Color_Map= new QAction(tr("Maximum Curvature to Color Map"), this);
			if (actionMaximum_Curvature_to_Color_Map)
				connect(actionMaximum_Curvature_to_Color_Map, SIGNAL(triggered()), this, SLOT(OnDisplayMax()));

			actionDisplay_Minimum_Curvature_direction= new QAction(tr("Display Minimum Curvature direction"), this);
			if (actionDisplay_Minimum_Curvature_direction)
				connect(actionDisplay_Minimum_Curvature_direction, SIGNAL(triggered()), this, SLOT(OnDisplayMinDir()));

			actionDisplay_Maximum_Curvature_direction= new QAction(tr("Display Maximum Curvature direction"), this);
			if (actionDisplay_Maximum_Curvature_direction)
				connect(actionDisplay_Maximum_Curvature_direction, SIGNAL(triggered()), this, SLOT(OnDisplayMaxDir()));
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionCalculate_Normal_Cycle_Algorithm
										<< NULL
										<< actionMinimum_Curvature_to_Color_Map
										<< actionMaximum_Curvature_to_Color_Map
										<< NULL
										<< actionDisplay_Minimum_Curvature_direction
										<< actionDisplay_Maximum_Curvature_direction;
		}	
		
		virtual void pre_draw() {}
		virtual void post_draw();
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene() {}

		virtual void OnMouseLeftDown(QMouseEvent *event);
		virtual void OnMouseLeftUp(QMouseEvent *event) {}
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event) {}
		virtual void OnMouseMotion(QMouseEvent *event) {}
		virtual void OnMouseWheel(QWheelEvent *event) {}
		virtual void OnKeyPress(QKeyEvent *event) {}
		virtual void OnKeyRelease(QKeyEvent *event) {}

	public slots:
		void OnCurvature();
		void OnDisplayMin();
		void OnDisplayMax();
		void OnDisplayMinDir();
		void OnDisplayMaxDir();

	private:
		QAction *actionCalculate_Normal_Cycle_Algorithm;
		QAction *actionMinimum_Curvature_to_Color_Map, *actionMaximum_Curvature_to_Color_Map;
		QAction *actionDisplay_Minimum_Curvature_direction, *actionDisplay_Maximum_Curvature_direction;
};

#endif

#endif
