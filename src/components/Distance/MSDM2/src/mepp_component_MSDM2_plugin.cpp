/*!
	\file MSDM2_Component.h
	\brief MSDM2 Perceptual distance calculation
	\brief According to: Multiscale Metric for 3D Mesh Visual Quality Assessment
			G. Lavoue
			In Symposium on Geometry Processing (SGP) 2011
	\author Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
	\date	2011
 */
 
#include <mepp_config.h>
//#ifdef BUILD_component_MSDM2

#include "mepp_component_MSDM2_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>
#include <CGAL/Timer.h>

#include "MSDM2_Component.h"

#include "../../../Analysis/Curvature/src/Curvature_Component.h"

typedef boost::shared_ptr<MSDM2_Component> MSDM2_ComponentPtr;
typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;


void mepp_component_MSDM2_plugin::MSDM2_computation()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		
		MSDM2_ComponentPtr component_ptr = findOrCreateComponentForViewer<MSDM2_ComponentPtr, MSDM2_Component>(viewer, polyhedron_ptr);

		
		SettingsDialog dial;
		if(viewer->getScenePtr()->get_nb_polyhedrons() == 2)
		{
			if (dial.exec() == QDialog::Accepted)
			{

				Timer timer;
				timer.start();	

				viewer->getScenePtr()->get_polyhedron(0)->IsDistanceComputed=false;
				viewer->getScenePtr()->get_polyhedron(1)->IsDistanceComputed=false;

				QApplication::setOverrideCursor(Qt::WaitCursor);

				double ScaleNb=dial.Scale->value();

				Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);


				PolyhedronPtr polyhedron_ptr_in1; 
				PolyhedronPtr polyhedron_ptr_in2;
				if(!dial.radioSym->isChecked())
				{
					if (dial.radio12->isChecked())
					{
						polyhedron_ptr_in1= viewer->getScenePtr()->get_polyhedron(0);
						polyhedron_ptr_in2= viewer->getScenePtr()->get_polyhedron(1);
					}
					if (dial.radio21->isChecked())
					{
						polyhedron_ptr_in1= viewer->getScenePtr()->get_polyhedron(1);
						polyhedron_ptr_in2= viewer->getScenePtr()->get_polyhedron(0);
					}
				
					mw->statusBar()->showMessage(tr("MSDM2 computation..."));

					//////////////processing here/////////////////////////////////////////////
					double maxdim=component_ptr->getMaxDim(polyhedron_ptr_in2);
		
			
					double MSDM2Value;
					component_ptr->ProcessMSDM2_Multires(polyhedron_ptr_in2,polyhedron_ptr_in1,ScaleNb,maxdim,MSDM2Value,component_ptr_curvature);
	
					timer.stop();
				
					QApplication::restoreOverrideCursor();

					//////////////end processing here/////////////////////////////////////////////
					mw->statusBar()->showMessage(tr("MSDM2 computation done"));
					QString time = QString("Processing time : %1 seconds \n").arg(timer.time(), 4, 'f', 3);
					QString value = QString("MSDM2 value : %1 \n").arg((double)MSDM2Value, 6, 'f', 5);
					QMessageBox::information(mw, APPLICATION, value+time);

					viewer->recreateListsAndUpdateGL();
				}
				else //symmetric
				{
					mw->statusBar()->showMessage(tr("MSDM2 computation..."));

					//////////////processing here/////////////////////////////////////////////
					double maxdim0=component_ptr->getMaxDim(viewer->getScenePtr()->get_polyhedron(0));
					double maxdim1=component_ptr->getMaxDim(viewer->getScenePtr()->get_polyhedron(1));
					
					double MSDM2Value0;
					double MSDM2Value1;

					component_ptr->ProcessMSDM2_Multires(viewer->getScenePtr()->get_polyhedron(1),viewer->getScenePtr()->get_polyhedron(0),ScaleNb,maxdim1,MSDM2Value1,component_ptr_curvature);
	
					component_ptr->ProcessMSDM2_Multires(viewer->getScenePtr()->get_polyhedron(0),viewer->getScenePtr()->get_polyhedron(1),ScaleNb,maxdim0,MSDM2Value0,component_ptr_curvature);
	
					double MSDM2Value=(MSDM2Value0+MSDM2Value1)/2;

					timer.stop();
				
					QApplication::restoreOverrideCursor();

					//////////////end processing here/////////////////////////////////////////////
					mw->statusBar()->showMessage(tr("MSDM2 computation done"));
					QString time = QString("Processing time : %1 seconds \n").arg(timer.time(), 4, 'f', 3);
					QString value = QString("MSDM2 value : %1 \n").arg((double)MSDM2Value, 6, 'f', 5);
					QMessageBox::information(mw, APPLICATION, value+time);

					viewer->recreateListsAndUpdateGL();
				}
				
			}
		}
		else
			QMessageBox::information(mw, APPLICATION, tr("MSDM2 computation needs 2 meshes opened in time or in space"));
	}

	
}





void mepp_component_MSDM2_plugin::DistanceToColorMap()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);
	

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		MSDM2_ComponentPtr component_ptr = findOrCreateComponentForViewer<MSDM2_ComponentPtr, MSDM2_Component>(viewer, polyhedron_ptr);
		component_ptr->ComputeMaxMin(polyhedron_ptr,1);	
		component_ptr->ConstructColorMap(polyhedron_ptr,1);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();

}








Q_EXPORT_PLUGIN2(mepp_component_MSDM2_plugin, mepp_component_MSDM2_plugin);

//#endif
