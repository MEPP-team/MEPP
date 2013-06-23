///////////////////////////////////////////////////////////////////////////
// Author: Guillaume LAVOUE
// Year: 2010
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
//According to: Perceptually driven 3d distance metrics with application to watermarking
//			G. Lavoué, E. Drelie Gelasca, F. Dupont, A. Baskurt, and T. Ebrahimi 
//			In SPIE Applications of Digital Image Processing XXIX, vol. 6312, 2006.
//
///////////////////////////////////////////////////////////////////////////
#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#include "mepp_component_MSDM_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>
#include <CGAL/Timer.h>

#include "MSDM_Component.h"

#include "../../../Analysis/Curvature/src/Curvature_Component.h"

typedef boost::shared_ptr<MSDM_Component> MSDM_ComponentPtr;
typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;


void mepp_component_MSDM_plugin::MSDM_computation()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		MSDM_ComponentPtr component_ptr = findOrCreateComponentForViewer<MSDM_ComponentPtr, MSDM_Component>(viewer, polyhedron_ptr);

		float RadiusVal;
			

		SettingsDialog dial;
		if(viewer->getScenePtr()->get_nb_polyhedrons() == 2)
		{
			if (dial.exec() == QDialog::Accepted)
			{

				Timer timer;
				timer.start();	


				QApplication::setOverrideCursor(Qt::WaitCursor);

				RadiusVal=dial.Radius->value();

				PolyhedronPtr polyhedron_ptr_in1; 
				PolyhedronPtr polyhedron_ptr_in2;

				if (dial.radioGeo->isChecked())
				{
					polyhedron_ptr_in1= viewer->getScenePtr()->get_polyhedron(0);
					polyhedron_ptr_in2= viewer->getScenePtr()->get_polyhedron(1);
				}
				if (dial.radio1Ring->isChecked())
				{
					polyhedron_ptr_in1= viewer->getScenePtr()->get_polyhedron(1);
					polyhedron_ptr_in2= viewer->getScenePtr()->get_polyhedron(0);
				}
				
				mw->statusBar()->showMessage(tr("MSDM computation..."));

				//////////////processing here/////////////////////////////////////////////
				double maxdim=component_ptr->getMaxDim(polyhedron_ptr_in1);
		
				Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
			
				component_ptr_curvature->principal_curvature(polyhedron_ptr_in1,true,RadiusVal*maxdim);
				component_ptr_curvature->principal_curvature(polyhedron_ptr_in2,true,RadiusVal*maxdim);
				
				component_ptr->KmaxKmean(polyhedron_ptr_in1,maxdim);
				component_ptr->KmaxKmean(polyhedron_ptr_in2,maxdim);
		
				component_ptr->ComputeLocalCurvatureStatistics(polyhedron_ptr_in1,polyhedron_ptr_in2,0.015*maxdim,maxdim);

				double L;
				component_ptr->ComputeMSDM_FromStatistics(polyhedron_ptr_in1,polyhedron_ptr_in2,L);


				timer.stop();
				
				QApplication::restoreOverrideCursor();

				//////////////end processing here/////////////////////////////////////////////
				//mw->statusBar()->showMessage(tr("MSDM computation done"));
				QString time = QString("Processing time : %1 seconds \n").arg(timer.time(), 4, 'f', 3);
				QString value = QString("MSDM value : %1 \n").arg((float)L, 6, 'f', 5);
				QMessageBox::information(mw, APPLICATION, value+time);

				viewer->recreateListsAndUpdateGL();

				
			}
		}
		else
			QMessageBox::information(mw, APPLICATION, tr("MSDM computation needs 2 meshes opened in time or in space"));
	}

	
}



void mepp_component_MSDM_plugin::DistanceToColorMap()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);
	

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		MSDM_ComponentPtr component_ptr = findOrCreateComponentForViewer<MSDM_ComponentPtr, MSDM_Component>(viewer, polyhedron_ptr);
		component_ptr->ComputeMaxMin(polyhedron_ptr);	
		component_ptr->ConstructColorMap(polyhedron_ptr);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();

}


#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
Q_EXPORT_PLUGIN2(mepp_component_MSDM_plugin, mepp_component_MSDM_plugin);
#endif

#endif
