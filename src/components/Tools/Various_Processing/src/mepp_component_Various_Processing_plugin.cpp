#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "mepp_component_Various_Processing_plugin.hxx"

#include "dialSettings_Various_Processing_Noise.hxx"
#include "dialSettings_Various_Processing_Smoothing.hxx"
#include "dialSettings_Various_Processing_Quantization.hxx"
#include "dialSettings_Various_Processing_Translation.hxx"
#include "dialSettings_Various_Processing_Rotation.hxx"
#include "dialSettings_Various_Processing_Uniform_Scaling.hxx"
#include "dialSettings_Various_Processing_Subdivision.hxx"
#include "dialSettings_Various_Processing_Simplification.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "Various_Processing_Component.h"
typedef boost::shared_ptr<Various_Processing_Component> Various_Processing_ComponentPtr;

void mepp_component_Various_Processing_plugin::NoiseAddition()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Noise dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char itensityChar[256];
			strcpy(itensityChar, dial.m_textCtrlIntensity->GetValue().ToAscii());*/
			double itensity = dial.Intensity->value();

			Noise_type noiseType = UNIFORM;
			if (dial.radioUniform->isChecked())
				noiseType = UNIFORM;
			else
				noiseType = GAUSSIAN;

			bool preserveBoundaries = dial.PreserveBoundaries->isChecked();

			//wxBusyInfo busy(_T("Add random noise..."));
			//m_frame->set_status_message(_T("Add random noise..."));
			mw->statusBar()->showMessage(tr("Add random noise..."));
			//double start = clock();

			component_ptr->NoiseAddition(polyhedron_ptr,noiseType,itensity,preserveBoundaries);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Add random noise...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Add random noise...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Processing_plugin::LaplacianSmoothing()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Smoothing dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char deformFactorChar[256];
			strcpy(deformFactorChar, dial.m_textCtrlDeformFactor->GetValue().ToAscii());*/
			double deformFactor = dial.DeformFactor->value();

			/*char IteraNumChar[256];
			strcpy(IteraNumChar, dial.m_textCtrlIteraNum->GetValue().ToAscii());*/
			int iteraNum = dial.IteraNum->value();

			bool preserveBoundaries = dial.PreserveBoundaries->isChecked();

			//wxBusyInfo busy(_T("Laplacian smoothing..."));
			//m_frame->set_status_message(_T("Laplacian smoothing..."));
			mw->statusBar()->showMessage(tr("Laplacian smoothing..."));
			//double start = clock();

			component_ptr->LaplacianSmoothing(polyhedron_ptr,deformFactor,iteraNum,preserveBoundaries);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Laplacian smoothing...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Laplacian smoothing...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Processing_plugin::CoordinateQuantization()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Quantization dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char bitDepthChar[256];
			strcpy(bitDepthChar, dial.m_textCtrlBitDepth->GetValue().ToAscii());*/
			int bitDepth = dial.CtrlBitDepth->value();

			//wxBusyInfo busy(_T("Coordinate quantization..."));
			//m_frame->set_status_message(_T("Coordinate quantization..."));
			mw->statusBar()->showMessage(tr("Coordinate quantization..."));
			//double start = clock();

			component_ptr->CoordinateQuantization(polyhedron_ptr,bitDepth);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Coordinate quantization...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Coordinate quantization...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Various_Processing_plugin::Translation()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Translation dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char xTranslationChar[256];
			strcpy(xTranslationChar, dial.m_textCtrlTranslationX->GetValue().ToAscii());*/
			double xTranslation = dial.CtrlTranslationX->value();

			/*char yTranslationChar[256];
			strcpy(yTranslationChar, dial.m_textCtrlTranslationY->GetValue().ToAscii());*/
			double yTranslation = dial.CtrlTranslationY->value();

			/*char zTranslationChar[256];
			strcpy(zTranslationChar, dial.m_textCtrlTranslationZ->GetValue().ToAscii());*/
			double zTranslation = dial.CtrlTranslationZ->value();

			//wxBusyInfo busy(_T("Object translation..."));
			//m_frame->set_status_message(_T("Object translation..."));
			mw->statusBar()->showMessage(tr("Object translation..."));
			//double start = clock();

			component_ptr->Translation(polyhedron_ptr,xTranslation,yTranslation,zTranslation);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Object translation...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Object translation...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Processing_plugin::Rotation()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Rotation dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char xAxisChar[256];
			strcpy(xAxisChar, dial.m_textCtrlXForAxis->GetValue().ToAscii());*/
			double xAxis = dial.CtrlXForAxis->value();

			/*char yAxisChar[256];
			strcpy(yAxisChar, dial.m_textCtrlYForAxis->GetValue().ToAscii());*/
			double yAxis = dial.CtrlYForAxis->value();

			/*char zAxisChar[256];
			strcpy(zAxisChar, dial.m_textCtrlZForAxis->GetValue().ToAscii());*/
			double zAxis = dial.CtrlZForAxis->value();

			/*char angleChar[256];
			strcpy(angleChar, dial.m_textCtrlAngle->GetValue().ToAscii());*/
			double angle = dial.CtrlAngle->value();

			//wxBusyInfo busy(_T("Object rotation..."));
			//m_frame->set_status_message(_T("Object rotation..."));
			mw->statusBar()->showMessage(tr("Object rotation..."));
			//double start = clock();

			component_ptr->Rotation(polyhedron_ptr, xAxis, yAxis, zAxis, angle);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Object rotation...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Object rotation...done"));

			viewer->recreateListsAndUpdateGL();
		}
}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Processing_plugin::UniformScaling()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Uniform_Scaling dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char scalingFactorChar[256];
			strcpy(scalingFactorChar, dial.m_textCtrlScalingFactor->GetValue().ToAscii());*/
			double scalingFactor = dial.ScalingFactor->value();

			//wxBusyInfo busy(_T("Object uniform scaling..."));
			//m_frame->set_status_message(_T("Object uniform scaling..."));
			mw->statusBar()->showMessage(tr("Object uniform scaling..."));
			//double start = clock();

			component_ptr->UniformScaling(polyhedron_ptr,scalingFactor);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Object uniform scaling...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Object uniform scaling...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Various_Processing_plugin::SurfaceSubdivision()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Subdivision dial;
		if (dial.exec() == QDialog::Accepted)
		{
			/*char depthChar[256];
			strcpy(depthChar, dial.m_textCtrlDepth->GetValue().ToAscii());*/
			int depth = dial.Depth->value();

			Subdivision_type subdivisionType = CATMULLCLARK;
			if (dial.radioCatmullClark->isChecked())
				subdivisionType = CATMULLCLARK;
			else if (dial.radioLoop->isChecked())
				subdivisionType = LOOP;
			else if (dial.radioDooSabin->isChecked())
				subdivisionType = DOOSABIN;
			else if (dial.radioSqrt3->isChecked())
				subdivisionType = SQRT3;
			else if (dial.radioMidpoint->isChecked())
				subdivisionType = MIDPOINT;

			if (subdivisionType==LOOP)
			{
				if (!polyhedron_ptr->is_pure_triangle())
				{
					/*(void)wxMessageBox(_T("Loop subdivision is not possible.\n\n")
						_T("The current mesh has non-triangular facets.")
						   );*/
					QMessageBox::information(mw, APPLICATION, tr("Loop subdivision is not possible: the current mesh has non-triangular facets."));
					return;
				}
			}

			if (subdivisionType==SQRT3)
			{
				if (!polyhedron_ptr->is_pure_triangle())
				{
					/*(void)wxMessageBox(_T("Sqrt3 subdivision is not possible.\n\n")
						_T("The current mesh has non-triangular facets.")
						   );*/
					QMessageBox::information(mw, APPLICATION, tr("Sqrt3 subdivision is not possible: the current mesh has non-triangular facets."));
					return;
				}
			}

			QApplication::setOverrideCursor(Qt::WaitCursor);

			//wxBusyInfo busy(_T("Surface subdivision..."));
			//m_frame->set_status_message(_T("Surface subdivision..."));
			mw->statusBar()->showMessage(tr("Surface subdivision..."));
			//double start = clock();

			component_ptr->Subdivision(polyhedron_ptr,subdivisionType,depth);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Surface subdivision...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Surface subdivision...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Processing_plugin::SurfaceSimplification()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		SettingsDialog_Various_Processing_Simplification dial;
		if (dial.exec() == QDialog::Accepted)
		{
			if (!polyhedron_ptr->is_pure_triangle())
			{
				/*(void)wxMessageBox(_T("Simplification is not possible.\n\n")
					_T("The current mesh has non-triangular facets.")
					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Simplification is not possible: the current mesh has non-triangular facets."));
				return;
			}

			QApplication::setOverrideCursor(Qt::WaitCursor);

			/*char targetEdgeNumChar[256];
			strcpy(targetEdgeNumChar, dial.m_textCtrlEdgeNum->GetValue().ToAscii());*/
			int targetEdgeNum = dial.EdgeNum->value();

			//wxBusyInfo busy(_T("Mesh simplification..."));
			//m_frame->set_status_message(_T("Mesh simplification..."));
			mw->statusBar()->showMessage(tr("Mesh simplification..."));
			//double start = clock();

			component_ptr->Simplification(polyhedron_ptr.get(),targetEdgeNum);

			/*float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
			m_frame->update_mesh_properties();
			m_frame->Refresh();
			wxString msg;
			msg.Printf(_T("Mesh simplification...done (%g s)"), duration);
			m_frame->set_status_message(msg);*/
			mw->statusBar()->showMessage(tr("Mesh simplification...done"));

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(mepp_component_Various_Processing_plugin, mepp_component_Various_Processing_plugin);

#endif
