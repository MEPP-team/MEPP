#include <mepp_config.h>
#ifdef BUILD_component_VSA

#include "mepp_component_VSA_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "VSA_Component.h"
typedef boost::shared_ptr<VSA_Component> VSA_ComponentPtr;

void mepp_component_VSA_plugin::VariationalSegmentation()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr);

		bool IsIncre;
		int RegNb,IterNb;
		//char RegNbChar[256];
		//char IterNbChar[256];

		if (!polyhedron_ptr->is_pure_triangle())
		{
			/*(void)wxMessageBox(_T("Segmentation not possible\n\n")
				_T("The mesh owns non-triangular facets")			  
					   
				   );*/
			QMessageBox::information(mw, APPLICATION, tr("Segmentation not possible: the mesh owns non-triangular facets."));
			return;
		}
		if (polyhedron_ptr->nb_components()!=1)
		{
			/*(void)wxMessageBox(_T("Segmentation not possible\n\n")
				_T("The mesh owns more than one component")			  
					   
				   );*/
			QMessageBox::information(mw, APPLICATION, tr("Segmentation not possible: the mesh owns more than one component."));
			return;
		}
		
		SettingsDialog dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			//strcpy(RegNbChar,dial.m_textCtrlReg->GetValue().ToAscii());
			//strcpy(IterNbChar,dial.m_textCtrlIter->GetValue().ToAscii());
			RegNb=dial.CtrlReg->value();//atoi(RegNbChar);
			IterNb=dial.CtrlIter->value();//atoi(IterNbChar);

			if (dial.radioInc->isChecked())
				IsIncre=true;
			else
				IsIncre=false;

			//wxBusyInfo busy(_T("Variational Segmentation..."));

			//m_frame->set_status_message(_T("Variational Segmentation..."));
			mw->statusBar()->showMessage(tr("Variational Segmentation..."));

			if(IsIncre==true)
				component_ptr->Variational_SegmentationIncr(polyhedron_ptr,RegNb,IterNb);
			else
				component_ptr->Variational_Segmentation(polyhedron_ptr,RegNb,IterNb);

			//m_frame->set_status_message(_T("Segmentation...done"));
			mw->statusBar()->showMessage(tr("Segmentation is done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_VSA_plugin::FaceLabelsToColorMap()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		VSA_ComponentPtr component_ptr = findOrCreateComponentForViewer<VSA_ComponentPtr, VSA_Component>(viewer, polyhedron_ptr);
		component_ptr->ConstructFaceColorMap(polyhedron_ptr);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
Q_EXPORT_PLUGIN2(mepp_component_VSA_plugin, mepp_component_VSA_plugin);
#endif

#endif
