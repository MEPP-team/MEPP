///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

#include "mepp_component_CGAL_Example_plugin.hxx"

#include "dialSettings_CGAL_Example.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "CGAL_Example_Component.h"
typedef boost::shared_ptr<CGAL_Example_Component> CGAL_Example_ComponentPtr;

#if (0)
// we want to use Curvature component
#include "../../../Analysis/Curvature/src/Curvature_Component.h"
typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;
// we want to use Curvature component
#endif

void mepp_component_CGAL_Example_plugin::pre_draw()
{
	// active viewer
	//if (mw->activeMdiChild() != 0)
	/*{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
			if (component_ptr->get_init() == 2)
			{
				glPushMatrix();
					// here your code
					glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glPopMatrix();
			}
		}
	}*/
}

void mepp_component_CGAL_Example_plugin::post_draw()
{
	// active viewer
	//if (mw->activeMdiChild() != 0)
	/*{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
			if (component_ptr->get_init() == 2)
			{
				glPushMatrix();
					// here your code
				glPopMatrix();
			}
		}
	}*/
}

void mepp_component_CGAL_Example_plugin::pre_draw_all_scene()
{
	// active viewer
	//if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
			if (component_ptr->get_init() == 3)
			{
				glPushMatrix();
					// here your code
				glPopMatrix();
			}
		}
	}
}

void mepp_component_CGAL_Example_plugin::post_draw_all_scene()
{
	// active viewer
	//if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
			if (component_ptr->get_init() == 3)
			{
				int nbMesh = qMin(viewer->getScenePtr()->get_nb_polyhedrons(), viewer->get_nb_frames());
				if (nbMesh == 2)
				{
					glPushMatrix();
					glDisable(GL_LIGHTING);

						// here your code
						draw_connections(viewer, 0, 1); // link between first and second mesh (first and second frame)

					glEnable(GL_LIGHTING);
					glPopMatrix();
				}
			}
		}
	}
}

void mepp_component_CGAL_Example_plugin::draw_connections(Viewer* viewer, int frame_i, int frame_j)
{
#if (1)
	PolyhedronPtr pMesh_i = viewer->getScenePtr()->get_polyhedron(frame_i);
	PolyhedronPtr pMesh_j = viewer->getScenePtr()->get_polyhedron(frame_j);

	int nbp_i = pMesh_i->size_of_vertices();
	int nbp_j = pMesh_j->size_of_vertices();

	if (nbp_i == nbp_j)
	{
		glLineWidth(2);
		glColor3f(1., 0., 0.);

		Vertex_iterator pVertex_i = NULL;
		Vertex_iterator pVertex_j = NULL;
		for (pVertex_i = pMesh_i->vertices_begin(), pVertex_j = pMesh_j->vertices_begin(); pVertex_i != pMesh_i->vertices_end(); pVertex_i++, pVertex_j++)
		{
			Vec pi(pVertex_i->point().x(), pVertex_i->point().y(), pVertex_i->point().z());
			Vec pj(pVertex_j->point().x(), pVertex_j->point().y(), pVertex_j->point().z());

			draw_link(viewer, frame_i, frame_j, pi, pj);
		}
	}
#endif
}

void mepp_component_CGAL_Example_plugin::OnMouseLeftDown(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseLeftDown"), 1000);

			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);

			if (viewer->getScenePtr()->get_loadType() == Time)
			{
				PolyhedronPtr new_polyhedron_ptr(new Polyhedron(*(viewer->getScenePtr()->get_polyhedron())));

				component_ptr->GetClickedVertices(viewer->getScenePtr()->get_polyhedron(), event->x(), event->y(), 10);
				viewer->recreateListsAndUpdateGL();
				SleeperThread::msleep(300);

				viewer->getScenePtr()->add_polyhedron(new_polyhedron_ptr);
				viewer->getScenePtr()->set_current_polyhedron(viewer->getScenePtr()->get_nb_polyhedrons()-1);

				viewer->setDynTitle();

				viewer->recreateListsAndUpdateGL();
			}
			else if (viewer->getScenePtr()->get_loadType() == Normal)
			{
				component_ptr->GetClickedVertices(viewer->getScenePtr()->get_polyhedron(), event->x(), event->y(), 10);
				viewer->recreateListsAndUpdateGL();
			}
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnMouseLeftUp(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseLeftUp"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnMouseRightDown(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseRightDown"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnMouseRightUp(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseRightUp"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnMouseMotion(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseMotion"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnMouseWheel(QWheelEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			int rot = event->delta();

			if (rot<0)
				mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseWheel Up"), 1000);
			else
				mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnMouseWheel Down"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnKeyPress(QKeyEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnKeyPress"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::OnKeyRelease(QKeyEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			mw->statusBar()->showMessage(tr("mepp_component_CGAL_Example_plugin: OnKeyRelease"), 1000);
		}
	}
}

void mepp_component_CGAL_Example_plugin::step1()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		component_ptr->TriangulateAndRandomColorFacets(polyhedron_ptr);

		component_ptr->set_init(2);
		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_CGAL_Example_plugin::step2()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		
		//wxColour current_color(m_component->round(m_component->color(0)*255.), m_component->round(m_component->color(1)*255.), m_component->round(m_component->color(2)*255.));
		//wxColour new_color = ::wxGetColourFromUser(m_frame, current_color);
		QColor current_color(int(component_ptr->color(0)*255.), int(component_ptr->color(1)*255.), int(component_ptr->color(2)*255.));
		QColor new_color = QColorDialog::getColor(current_color, viewer);
		
		if (new_color.isValid())
		{
			component_ptr->color(float(new_color.red())/255., float(new_color.green())/255., float(new_color.blue())/255.);

			SettingsDialog_CGAL_Example dial;
			if (dial.exec() == QDialog::Accepted)
			{
				QApplication::setOverrideCursor(Qt::WaitCursor);

				//char iterationChar[256];
				//strcpy(iterationChar, dial.m_textIteration->GetValue().ToAscii());
				int iteration = dial.Iteration->value();//atoi(iterationChar);

				//wxBusyInfo busy(_T("Create center vertex..."));

				//m_frame->set_status_message(_T("Create center vertex..."));
				mw->statusBar()->showMessage(tr("Create center vertex..."));

				if (viewer->getScenePtr()->get_nb_polyhedrons()==1)
					component_ptr->CreateCenterVertex(polyhedron_ptr, true);
				else
				{
					for (int p=0; p<viewer->getScenePtr()->get_nb_polyhedrons(); p++)
						for (int i=0; i<iteration; i++)
							component_ptr->CreateCenterVertex(viewer->getScenePtr()->get_polyhedron(p), false);
				}

				/*m_frame->update_mesh_properties();
				m_frame->Refresh();
				m_frame->set_status_message(_T("Create center vertex...done"));*/
				mw->statusBar()->showMessage(tr("Create center vertex...done"));

				viewer->recreateListsAndUpdateGL();

				QApplication::restoreOverrideCursor();
				return;
			}
		}

		//m_frame->set_status_message(_T("Create center vertex...canceled"));
		mw->statusBar()->showMessage(tr("Create center vertex...canceled"));
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_CGAL_Example_plugin::step3()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		component_ptr->ShowBlackAndWhiteFacets(polyhedron_ptr);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_CGAL_Example_plugin::step4()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		
		component_ptr->set_init(3);
		viewer->updateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_CGAL_Example_plugin::step5()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		
		Vec pw(2, 0, 0); viewer->frame(viewer->getScenePtr()->get_current_polyhedron())->setPosition(pw); // position in world coordinate system
		//Vec pl(2, 0, 0); viewer->frame(viewer->getScenePtr()->get_current_polyhedron())->setTranslation(pl); // local frame translation

		Quaternion qw(0, 0, 0, 1); // identity quaternion (i.e., no rotation)
		viewer->frame(viewer->getScenePtr()->get_current_polyhedron())->setOrientation(qw); // rotation in world coordinate system
		/*Quaternion ql(0, 0, 0, 1); // identity quaternion (i.e., no rotation)
		viewer->frame(viewer->getScenePtr()->get_current_polyhedron())->setRotation(ql); // local frame rotation*/

		viewer->updateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_CGAL_Example_plugin::step6()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		// step 6a : begin
		emit(mw->get_actionAddEmpty()->trigger());

		int nb_polyhedrons = viewer->getScenePtr()->get_nb_polyhedrons();
		polyhedron_ptr = viewer->getScenePtr()->get_polyhedron(nb_polyhedrons-1);
		component_ptr->CreateTetrahedron(polyhedron_ptr);

		viewer->getScenePtr()->setcurrentFile(tr("internal mesh sample from empty"));
		viewer->setDynTitle();

		viewer->recreateListsAndUpdateGL();
		// step 6a : end
	}
	else
	{
		// step 6b : begin
		emit(mw->get_actionNewEmpty()->trigger());

		for (int i=0; i<lwindow.size(); i++) // all viewers
		{
			Viewer* viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
			if (viewer->getScenePtr()->get_polyhedron()->empty())
			{
				PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

				CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
				component_ptr->CreateTetrahedron(polyhedron_ptr);

				viewer->showAllScene();

				viewer->getScenePtr()->setcurrentFile(tr("internal mesh sample from empty"));
				viewer->setDynTitle();
			}
		}
		// step 6b : end
	}

	QApplication::restoreOverrideCursor();
}

int mepp_component_plugin_interface::load_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer)
{
	// here your code
	mepp_component_CGAL_Example_plugin *mepp_component_plugin = NULL;
	for (int i=0; i<viewer->lplugin.size(); ++i) {
		if (dynamic_cast<mepp_component_CGAL_Example_plugin*>(viewer->lplugin[i]) != 0) {
			mepp_component_plugin = dynamic_cast<mepp_component_CGAL_Example_plugin*>(viewer->lplugin[i]);
			//cout << "mepp_component_plugin found" << endl;
		}
	}

	int res;

	if (mepp_component_plugin)
	{
		CGAL_Example_ComponentPtr component_ptr = mepp_component_plugin->findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);

		res = polyhedron_ptr->load_mesh_off(filename.toStdString());
	}
	else
		res=-1;
	
	return res;
}
void mepp_component_CGAL_Example_plugin::step7()
{
	emit(mw->get_mainwindowActionOpen()->doSendParamsOpen(tr("Open Mesh File(s) - from CGAL_Example"), tr("OFF files (*.off)"), Normal, mepp_component_plugin_interface::load_file_from_component));

	//emit(mw->get_mainwindowActionOpen_space()->doSendParams(tr("Open Mesh File(s) (space) - from CGAL_Example"), tr("OFF files (*.off)"), mepp_component_plugin_interface::load_file_from_component));
	//emit(mw->get_mainwindowActionOpen_time()->doSendParams(tr("Open Mesh File(s) (time) - from CGAL_Example"), tr("OFF files (*.off)"), mepp_component_plugin_interface::load_file_from_component));
	//emit(mw->get_mainwindowActionOpen_and_Add_space()->doSendParams(tr("Open and Add Mesh File(s) (space) - from CGAL_Example"), tr("OFF files (*.off)"), mepp_component_plugin_interface::load_file_from_component));
	//emit(mw->get_mainwindowActionOpen_and_Add_time()->doSendParams(tr("Open and Add Mesh File(s) (time) - from CGAL_Example"), tr("OFF files (*.off)"), mepp_component_plugin_interface::load_file_from_component));
}

int mepp_component_plugin_interface::save_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer)
{
	// here your code
	return 0;
}
void mepp_component_CGAL_Example_plugin::step8()
{
	emit(mw->get_mainwindowActionSave_As()->doSendParams(tr("Save Mesh File - from CGAL_Example"), tr("OFF Files (*.off)"), mepp_component_plugin_interface::save_file_from_component));
}

void mepp_component_CGAL_Example_plugin::step9()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		QApplication::setOverrideCursor(Qt::WaitCursor);

		// we use CGAL_Example component here (as usual)
			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
			component_ptr->TriangulateAndRandomColorFacets(polyhedron_ptr);

			component_ptr->set_init(2);
		// we use CGAL_Example component here (as usual)

#if (0)
		// we use Curvature component here
			bool IsGeo;
			double radius;
			Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
				
			// params
			radius = 0.001;
			IsGeo=true;

			mw->statusBar()->showMessage(tr("Curvature..."));				
			component_ptr_curvature->principal_curvature(polyhedron_ptr,IsGeo,radius);
			mw->statusBar()->showMessage(tr("Curvature is done"));

			component_ptr_curvature->set_init(2);

			component_ptr_curvature->ConstructColorMap(polyhedron_ptr,1);
			viewer->recreateListsAndUpdateGL();
		// we use Curvature component here
#endif
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_CGAL_Example_plugin::example()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// all viewers
	for (int i=0; i<lwindow.size(); i++)
	{
		Viewer* viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		/*QMessageBox::information(mw, APPLICATION,
								   tr("Window: \"%1\".").
								   arg(viewer->userFriendlyCurrentFile()));*/

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		component_ptr->TriangulateAndRandomColorFacets(polyhedron_ptr);

		viewer->recreateListsAndUpdateGL();
	}	

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
		component_ptr->CreateCenterVertex(polyhedron_ptr, false);

		viewer->recreateListsAndUpdateGL();
	}

	QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(mepp_component_CGAL_Example_plugin, mepp_component_CGAL_Example_plugin);

#endif
