///////////////////////////////////////////////////////////////////////////
// Author: Guillaume LAVOUE
// Year: 2010
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
//According to: Perceptually driven 3d distance metrics with application to watermarking
//			G. Lavoué, E. Drelie Gelasca, F. Dupont, A. Baskurt, and T. Ebrahimi 
//			In SPIE Applications of Digital Image Processing XXIX, vol. 6312, 2006.
//
///////////////////////////////////////////////////////////////////////////

#ifndef HEADER_MEPP_COMPONENT_MSDM_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_MSDM_PLUGIN_INTERFACE_H


#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

/**
 \class	mepp_component_MSDM_plugin

 \brief	Mepp component msdm plugin. 

 
 */

class mepp_component_MSDM_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:
		mepp_component_MSDM_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_MSDM_plugin()
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

			action_distance = new QAction(tr("Compute MSDM distance"), this);
			if (action_distance)
				connect(action_distance, SIGNAL(triggered()), this, SLOT(MSDM_computation()));

			action_display = new QAction(tr("Local MSDM to Color Map "), this);
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
		virtual void post_draw_all_scene() {}

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

		/**
		 \fn	void mepp_component_MSDM_plugin::MSDM_computation();
		
		 \brief MSDM Calculation.	
		
		 */
		void MSDM_computation();

		/**
		 \fn	void DistanceToColorMap();
		
		 \brief	 Display the MSDM color map.
		
		
		 */
		void DistanceToColorMap();

	private:
		QAction *action_distance;
		QAction *action_display;
};

#endif

#endif
