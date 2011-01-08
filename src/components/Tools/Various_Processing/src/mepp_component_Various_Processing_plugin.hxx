#ifndef HEADER_MEPP_COMPONENT_VARIOUS_PROCESSING_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_VARIOUS_PROCESSING_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_Various_Processing_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:
		mepp_component_Various_Processing_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_Various_Processing_plugin()
		{
			delete actionNoiseAddition; delete actionLaplacianSmoothing; delete actionCoordinateQuantization;
			delete actionTranslation; delete actionRotation; delete actionUniformScaling;
			delete actionSurfaceSubdivision; delete actionSurfaceSimplification;
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
			actionNoiseAddition = new QAction(tr("Noise addition"), this);
			if (actionNoiseAddition)
				connect(actionNoiseAddition, SIGNAL(triggered()), this, SLOT(NoiseAddition()));

			actionLaplacianSmoothing = new QAction(tr("Laplacian smoothing"), this);
			if (actionLaplacianSmoothing)
				connect(actionLaplacianSmoothing, SIGNAL(triggered()), this, SLOT(LaplacianSmoothing()));

			actionCoordinateQuantization = new QAction(tr("Coordinate quantization"), this);
			if (actionCoordinateQuantization)
				connect(actionCoordinateQuantization, SIGNAL(triggered()), this, SLOT(CoordinateQuantization()));

			actionTranslation = new QAction(tr("Translation"), this);
			if (actionTranslation)
				connect(actionTranslation, SIGNAL(triggered()), this, SLOT(Translation()));

			actionRotation = new QAction(tr("Rotation"), this);
			if (actionRotation)
				connect(actionRotation, SIGNAL(triggered()), this, SLOT(Rotation()));

			actionUniformScaling = new QAction(tr("Uniform scaling"), this);
			if (actionUniformScaling)
				connect(actionUniformScaling, SIGNAL(triggered()), this, SLOT(UniformScaling()));

			actionSurfaceSubdivision = new QAction(tr("Surface subdivision"), this);
			if (actionSurfaceSubdivision)
				connect(actionSurfaceSubdivision, SIGNAL(triggered()), this, SLOT(SurfaceSubdivision()));

			actionSurfaceSimplification = new QAction(tr("Surface simplification"), this);
			if (actionSurfaceSimplification)
				connect(actionSurfaceSimplification, SIGNAL(triggered()), this, SLOT(SurfaceSimplification()));
			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionNoiseAddition
										<< actionLaplacianSmoothing
										<< actionCoordinateQuantization
										<< NULL			// menu separator
										<< actionTranslation
										<< actionRotation
										<< actionUniformScaling
										<< NULL			// menu separator
										<< actionSurfaceSubdivision
#ifndef __linux__
                                                                                << actionSurfaceSimplification
#endif
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

	public slots:
		void NoiseAddition();
		void LaplacianSmoothing();
		void CoordinateQuantization();
		void Translation();
		void Rotation();
		void UniformScaling();
		void SurfaceSubdivision();
		void SurfaceSimplification();

	private:
		QAction *actionNoiseAddition, *actionLaplacianSmoothing, *actionCoordinateQuantization;
		QAction *actionTranslation, *actionRotation, *actionUniformScaling;
		QAction *actionSurfaceSubdivision, *actionSurfaceSimplification;
};

#endif

#endif
