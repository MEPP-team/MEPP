#ifndef HEADER_MEPP_COMPONENT_VSA_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_VSA_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

/**
 \class	mepp_component_VSA_plugin

 \brief	Mepp component curvature plugin. 

 */
class mepp_component_VSA_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:
		mepp_component_VSA_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_VSA_plugin()
		{
			delete actionSegmentation; delete actionColorMap;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();
			
			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuSegmentation;

			// début --- actions ---
			actionSegmentation = new QAction(tr("Variational Segmentation"), this);
			if (actionSegmentation)
				connect(actionSegmentation, SIGNAL(triggered()), this, SLOT(VariationalSegmentation()));

			actionColorMap = new QAction(tr("Face Labels to Color Map"), this);
			if (actionColorMap)
				connect(actionColorMap, SIGNAL(triggered()), this, SLOT(FaceLabelsToColorMap()));
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionSegmentation
										<< actionColorMap;
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

		/**
		 \fn	void mepp_component_VSA_plugin::VariationalSegmentation();
		
		 \brief	Launch the Variational segmentation.
		
		
		 */
		void VariationalSegmentation();

		/**
		 \fn	void FaceLabelsToColorMap();
		
		 \brief	Display the segmentation.
		 */
		void FaceLabelsToColorMap();

	private:
		QAction *actionSegmentation, *actionColorMap;
};

#endif

#endif
