#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "mepp_component_Curvature_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "Curvature_Component.h"
typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;

void mepp_component_Curvature_plugin::post_draw()
{
	// active viewer
	//if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			Curvature_ComponentPtr component_ptr = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
			if (component_ptr->get_init() == 2 && polyhedron_ptr->curvature_is_calculated)
			{	
				glPushMatrix();
					// here your code
					if (component_ptr->get_displayMinDirections()==true)
					{
						glColor3f(1.0,0.0,0.0);
						::glBegin(GL_LINES);
						for (Vertex_iterator pVertex = polyhedron_ptr->vertices_begin();	pVertex !=	polyhedron_ptr->vertices_end();	pVertex++)
						{
							{

								float RayonMoyen=polyhedron_ptr->average_edge_length_around(pVertex);
								const Point3d& p1 = pVertex->point()-0.4*RayonMoyen*(pVertex->VKminCurv);
								const Point3d& p2 = pVertex->point()+0.4*RayonMoyen*(pVertex->VKminCurv);
								::glVertex3f(p1[0],p1[1],p1[2]);
								::glVertex3f(p2[0],p2[1],p2[2]);
							}


						}
						::glEnd();
					}

					if (component_ptr->get_displayMaxDirections()==true)
					{
						glColor3f(0.0,0.0,1.0);
						::glBegin(GL_LINES);
						for (Vertex_iterator pVertex = polyhedron_ptr->vertices_begin();	pVertex !=	polyhedron_ptr->vertices_end();	pVertex++)
						{
							{

								float RayonMoyen=polyhedron_ptr->average_edge_length_around(pVertex);
								const Point3d& p1 = pVertex->point()-0.4*RayonMoyen*(pVertex->VKmaxCurv);
								const Point3d& p2 = pVertex->point()+0.4*RayonMoyen*(pVertex->VKmaxCurv);
								::glVertex3f(p1[0],p1[1],p1[2]);
								::glVertex3f(p2[0],p2[1],p2[2]);
							}


						}
						::glEnd();
					}
				glPopMatrix();
			}
		}
	}
}

void mepp_component_Curvature_plugin::OnMouseLeftDown(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			//mw->statusBar()->showMessage(tr("mepp_component_Curvature_plugin: OnMouseLeftDown"), 1000);
		}
	}
}

void mepp_component_Curvature_plugin::OnCurvature()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Curvature_ComponentPtr component_ptr = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
		{
			bool IsGeo;
			double radius;
			//char radiusChar[256];

			SettingsDialog dial;
			if (dial.exec() == QDialog::Accepted)
			{
				QApplication::setOverrideCursor(Qt::WaitCursor);

				//strcpy(radiusChar,dial.m_textRadius->GetValue().ToAscii());
				//radius=atof(radiusChar);
				radius = dial.Radius->value();

				if (dial.radioGeo->isChecked())
					IsGeo=true;
				else
					IsGeo=false;

				//wxBusyInfo busy(_T("Calculating curvature..."));

				//m_frame->set_status_message(_T("Curvature..."));
				mw->statusBar()->showMessage(tr("Curvature..."));
				
				component_ptr->principal_curvature(polyhedron_ptr,IsGeo,radius);

				//m_frame->set_status_message(_T("Curvature...done"));
				mw->statusBar()->showMessage(tr("Curvature is done"));

				polyhedron_ptr->curvature_is_calculated=true;

				component_ptr->set_init(2);
				viewer->updateGL();
			}
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Curvature_plugin::OnDisplayMin()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Curvature_ComponentPtr component_ptr = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
		component_ptr->ConstructColorMap(polyhedron_ptr,1);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Curvature_plugin::OnDisplayMax()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Curvature_ComponentPtr component_ptr = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
		component_ptr->ConstructColorMap(polyhedron_ptr,2);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Curvature_plugin::OnDisplayMinDir()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Curvature_ComponentPtr component_ptr = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
		{
			if (component_ptr->get_displayMinDirections()==false)
				component_ptr->set_displayMinDirections(true);
			else
				component_ptr->set_displayMinDirections(false);
		}

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Curvature_plugin::OnDisplayMaxDir()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Curvature_ComponentPtr component_ptr = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
		{
			if (component_ptr->get_displayMaxDirections()==false)
				component_ptr->set_displayMaxDirections(true);
			else
				component_ptr->set_displayMaxDirections(false);
		}

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Curvature_plugin, mepp_component_Curvature_plugin);
#endif

#endif
