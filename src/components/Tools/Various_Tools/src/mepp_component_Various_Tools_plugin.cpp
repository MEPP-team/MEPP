#include <mepp_config.h>
#ifdef BUILD_component_Various_Tools

#include "mepp_component_Various_Tools_plugin.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "Various_Tools_Component.h"
typedef boost::shared_ptr<Various_Tools_Component> Various_Tools_ComponentPtr;

void mepp_component_Various_Tools_plugin::QuadTriangle()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		//m_frame->set_status_message(_T("Quad/triangle subdivision..."));
		mw->statusBar()->showMessage(tr("Quad/triangle subdivision..."));

		//wxStopWatch timer;

		component_ptr->subdivide_quad(polyhedron_ptr);

		//m_frame->update_mesh_properties();

		/*float duration = ((float)timer.Time())/1000.;
		wxString msg;
		msg.Printf(_T("Quad/triangle subdivision...done (%g s)"),duration);
		m_frame->set_status_message(msg);*/
		mw->statusBar()->showMessage(tr("Quad/triangle subdivision...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Tools_plugin::Sqrt3()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		//m_frame->set_status_message(_T("Sqrt3 subdivision..."));
		mw->statusBar()->showMessage(tr("Sqrt3 subdivision..."));

		component_ptr->subdivide_sqrt3(polyhedron_ptr);

		//m_frame->update_mesh_properties();
		//m_frame->set_status_message(_T("Subdivide...done"));
		mw->statusBar()->showMessage(tr("Subdivide...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Tools_plugin::Sqrt3Twice()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		//m_frame->set_status_message(_T("Sqrt3 Twice subdivision..."));
		mw->statusBar()->showMessage(tr("Sqrt3 Twice subdivision..."));

		component_ptr->subdivide_sqrt3Twice(polyhedron_ptr);

		//m_frame->update_mesh_properties();
		//m_frame->set_status_message(_T("Subdivide...done"));
		mw->statusBar()->showMessage(tr("Subdivide...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Tools_plugin::DooSabin()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		//m_frame->set_status_message(_T("Doo-Sabin subdivision..."));
		mw->statusBar()->showMessage(tr("Doo-Sabin subdivision..."));

		//wxStopWatch timer;

		component_ptr->subdivide_doo(polyhedron_ptr);

		//m_frame->update_mesh_properties();

		/*float duration = ((float)timer.Time())/1000.;
		wxString msg;
		msg.Printf(_T("Doo-Sabin subdivision...done (%g s)"),duration);
		m_frame->set_status_message(msg);*/
		mw->statusBar()->showMessage(tr("Doo-Sabin subdivision...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Tools_plugin::Loop()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		//m_frame->set_status_message(_T("Loop subdivision..."));
		mw->statusBar()->showMessage(tr("Loop subdivision..."));

		//wxStopWatch timer;

		component_ptr->subdivide_loop(polyhedron_ptr);

		//m_frame->update_mesh_properties();

		/*float duration = ((float)timer.Time())/1000.;
		wxString msg;
		msg.Printf(_T("Loop subdivision...done (%g s)"),duration);
		m_frame->set_status_message(msg);*/
		mw->statusBar()->showMessage(tr("Loop subdivision...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}
void mepp_component_Various_Tools_plugin::CatmullClark()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		//m_frame->set_status_message(_T("Catmull-Clark subdivision..."));
		mw->statusBar()->showMessage(tr("Catmull-Clark subdivision..."));

		//wxStopWatch timer;

		component_ptr->subdivide_catmull(polyhedron_ptr);

		//m_frame->update_mesh_properties();

		/*float duration = ((float)timer.Time())/1000.;
		wxString msg;
		msg.Printf(_T("Catmull-Clark subdivision...done (%g s)"),duration);
		m_frame->set_status_message(msg);*/
		mw->statusBar()->showMessage(tr("Catmull-Clark subdivision...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Various_Tools_plugin::Triangulate()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Various_Tools_ComponentPtr component_ptr = findOrCreateComponentForViewer<Various_Tools_ComponentPtr, Various_Tools_Component>(viewer, polyhedron_ptr);

		if (polyhedron_ptr->is_pure_triangle())
		{
			//wxMessageBox(_T("Mesh is already triangular."), _T("Triangulation"));
			QMessageBox::information(mw, APPLICATION, tr("Mesh is already triangular."));
			return;
		}

		QApplication::setOverrideCursor(Qt::WaitCursor);

		//m_frame->set_status_message(_T("Triangulating..."));
		mw->statusBar()->showMessage(tr("Triangulating..."));

		component_ptr->triangulate(polyhedron_ptr);

		//m_frame->update_mesh_properties(false, false);

		//m_frame->Refresh();

		mw->statusBar()->showMessage(tr("Triangulating...done"));

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Various_Tools_plugin, mepp_component_Various_Tools_plugin);
#endif

#endif
