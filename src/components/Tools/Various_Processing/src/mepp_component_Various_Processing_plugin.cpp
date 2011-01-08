#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "mepp_component_Various_Processing_plugin.hxx"

//#include "dialSettings_Various_Processing.hxx"

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
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		/*NoiseDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char itensityChar[256];
			strcpy(itensityChar, dial.m_textCtrlIntensity->GetValue().ToAscii());*/
			double itensity = 0.0005;//atof(itensityChar);

			Noise_type noiseType = UNIFORM;
			/*if (dial.m_radioBoxNoiseType->GetSelection()==0)
				noiseType = UNIFORM;
			else if (dial.m_radioBoxNoiseType->GetSelection()==1)
				noiseType = GAUSSIAN;*/

			bool preserveBoundaries = false;//dial.m_checkBoxPreserveBoundaries->GetValue();

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
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		/*SmoothingDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char deformFactorChar[256];
			strcpy(deformFactorChar, dial.m_textCtrlDeformFactor->GetValue().ToAscii());*/
			double deformFactor = 0.03;//atof(deformFactorChar);

			/*char IteraNumChar[256];
			strcpy(IteraNumChar, dial.m_textCtrlIteraNum->GetValue().ToAscii());*/
			int iteraNum = 10;//atoi(IteraNumChar);

			bool preserveBoundaries = true;//dial.m_checkBoxPreserveBoundaries->GetValue();

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
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		/*QuantizationDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char bitDepthChar[256];
			strcpy(bitDepthChar, dial.m_textCtrlBitDepth->GetValue().ToAscii());*/
			int bitDepth = 10;//atoi(bitDepthChar);

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
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		/*TranslationDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char xTranslationChar[256];
			strcpy(xTranslationChar, dial.m_textCtrlTranslationX->GetValue().ToAscii());*/
			double xTranslation = 0.2;//atof(xTranslationChar);

			/*char yTranslationChar[256];
			strcpy(yTranslationChar, dial.m_textCtrlTranslationY->GetValue().ToAscii());*/
			double yTranslation = 0.2;//atof(yTranslationChar);

			/*char zTranslationChar[256];
			strcpy(zTranslationChar, dial.m_textCtrlTranslationZ->GetValue().ToAscii());*/
			double zTranslation = 0.2;//atof(zTranslationChar);

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
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		/*RotationDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char xAxisChar[256];
			strcpy(xAxisChar, dial.m_textCtrlXForAxis->GetValue().ToAscii());*/
			double xAxis = 1.0;//atof(xAxisChar);

			/*char yAxisChar[256];
			strcpy(yAxisChar, dial.m_textCtrlYForAxis->GetValue().ToAscii());*/
			double yAxis = 0.0;//atof(yAxisChar);

			/*char zAxisChar[256];
			strcpy(zAxisChar, dial.m_textCtrlZForAxis->GetValue().ToAscii());*/
			double zAxis = 0.0;//atof(zAxisChar);

			/*char angleChar[256];
			strcpy(angleChar, dial.m_textCtrlAngle->GetValue().ToAscii());*/
			double angle = 90.0;//atof(angleChar);

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
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Processing_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Processing_ComponentPtr, Various_Processing_Component>(viewer, polyhedron_ptr);

		/*UniformScalingDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char scalingFactorChar[256];
			strcpy(scalingFactorChar, dial.m_textCtrlScalingFactor->GetValue().ToAscii());*/
			double scalingFactor = 1.5;//atof(scalingFactorChar);

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

		/*SubdivisionDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
		{
			/*char depthChar[256];
			strcpy(depthChar, dial.m_textCtrlDepth->GetValue().ToAscii());*/
			int depth = 1;//atoi(depthChar);

			Subdivision_type subdivisionType = CATMULLCLARK;
			/*if (dial.m_radioBoxSubdivisionType->GetSelection()==0)
				subdivisionType = CATMULLCLARK;
			else if (dial.m_radioBoxSubdivisionType->GetSelection()==1)
				subdivisionType = LOOP;
			else if (dial.m_radioBoxSubdivisionType->GetSelection()==2)
				subdivisionType = DOOSABIN;
			else if (dial.m_radioBoxSubdivisionType->GetSelection()==3)
				subdivisionType = SQRT3;
			else if (dial.m_radioBoxSubdivisionType->GetSelection()==4)
				subdivisionType = MIDPOINT;*/

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

		/*SimplificationDialogue dial(m_frame);
		if (dial.ShowModal() == wxID_OK)*/
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
			int targetEdgeNum = 80;//atoi(targetEdgeNumChar);

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
