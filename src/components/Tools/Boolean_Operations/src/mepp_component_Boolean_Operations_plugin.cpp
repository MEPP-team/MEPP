#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

//#include "../../../../mepp/scene.h"

#include "mepp_component_Boolean_Operations_plugin.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "Boolean_Operations_Component.h"

typedef boost::shared_ptr<Boolean_Operations_Component> Boolean_Operations_ComponentPtr;

void mepp_component_Boolean_Operations_plugin::New_Position()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		float x, y, z;
		double a, b, c, w;
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		
		ScenePtr S = viewer->getScenePtr();

		for(int i = 0;i<viewer->getScenePtr()->get_nb_polyhedrons();i++)
		{
			Vertex_iterator pVertex = NULL;

			PolyhedronPtr P = S->get_polyhedron(i);

			viewer->frame(i)->getPosition(x,y,z);
			viewer->frame(i)->getOrientation(a,b,c,w);

			Vec T(x, y, z);
			Quaternion Q(a, b, c, w);

			for (pVertex = P->vertices_begin(); pVertex != P->vertices_end(); pVertex++)
			{
				Vec V = Q * Vec(pVertex->point().x(), pVertex->point().y(), pVertex->point().z()) + T;
				pVertex->point() = CGAL::ORIGIN + Vector(V[0], V[1], V[2]);
			}

			viewer->frame(i)->setPosition(0,0,0);
			viewer->frame(i)->setOrientation(0,0,0,1);
		}
		viewer->show();
		viewer->updateGL();
	}
	QApplication::restoreOverrideCursor();
}

void mepp_component_Boolean_Operations_plugin::Subdiviser()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Boolean_Operations_ComponentPtr component_ptr = findOrCreateComponentForViewer<Boolean_Operations_ComponentPtr, Boolean_Operations_Component>(viewer, polyhedron_ptr);
		component_ptr->SubdiviserPolyedre(polyhedron_ptr);

		viewer->updateGL();
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Boolean_Operations_plugin::Union()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		
		if(viewer->getScenePtr()->get_nb_polyhedrons() == 2)
		{
			ScenePtr S = viewer->getScenePtr();
			PolyhedronPtr polyhedron_ptr_in1 = viewer->getScenePtr()->get_polyhedron(0);
			PolyhedronPtr polyhedron_ptr_in2 = viewer->getScenePtr()->get_polyhedron(1);
			PolyhedronPtr polyhedron_ptr_out;

			emit(mw->get_actionNewEmpty()->trigger());
			
			for (int i=0; i<lwindow.size(); i++) // all viewers
			{
				viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
				if (viewer->getScenePtr()->get_polyhedron()->empty())
				{
					polyhedron_ptr_out = viewer->getScenePtr()->get_polyhedron();
				}
			}

			Boolean_Operations_ComponentPtr component_ptr = findOrCreateComponentForViewer<Boolean_Operations_ComponentPtr, Boolean_Operations_Component>(viewer, polyhedron_ptr_in1);
			component_ptr->Boolean_Union(polyhedron_ptr_in1, polyhedron_ptr_in2, polyhedron_ptr_out);

			polyhedron_ptr_out->compute_bounding_box();
			polyhedron_ptr_out->compute_normals();
			polyhedron_ptr_out->compute_type();
			(void)polyhedron_ptr_out->calc_nb_components();
			(void)polyhedron_ptr_out->calc_nb_boundaries();

			viewer->show();

			viewer->updateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Boolean_Operations_plugin::Inter()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		
		if(viewer->getScenePtr()->get_nb_polyhedrons() == 2)
		{
			ScenePtr S = viewer->getScenePtr();
			PolyhedronPtr polyhedron_ptr_in1 = viewer->getScenePtr()->get_polyhedron(0);
			PolyhedronPtr polyhedron_ptr_in2 = viewer->getScenePtr()->get_polyhedron(1);
			PolyhedronPtr polyhedron_ptr_out;

			emit(mw->get_actionNewEmpty()->trigger());
			
			for (int i=0; i<lwindow.size(); i++) // all viewers
			{
				viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
				if (viewer->getScenePtr()->get_polyhedron()->empty())
				{
					polyhedron_ptr_out = viewer->getScenePtr()->get_polyhedron();
				}
			}

			Boolean_Operations_ComponentPtr component_ptr = findOrCreateComponentForViewer<Boolean_Operations_ComponentPtr, Boolean_Operations_Component>(viewer, polyhedron_ptr_in1);
			component_ptr->Boolean_Inter(polyhedron_ptr_in1, polyhedron_ptr_in2, polyhedron_ptr_out);

			polyhedron_ptr_out->compute_bounding_box();
			polyhedron_ptr_out->compute_normals();
			polyhedron_ptr_out->compute_type();
			(void)polyhedron_ptr_out->calc_nb_components();
			(void)polyhedron_ptr_out->calc_nb_boundaries();

			viewer->show();

			viewer->updateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Boolean_Operations_plugin::Minus()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		
		if(viewer->getScenePtr()->get_nb_polyhedrons() == 2)
		{
			ScenePtr S = viewer->getScenePtr();
			PolyhedronPtr polyhedron_ptr_in1 = viewer->getScenePtr()->get_polyhedron(0);
			PolyhedronPtr polyhedron_ptr_in2 = viewer->getScenePtr()->get_polyhedron(1);
			PolyhedronPtr polyhedron_ptr_out;

			emit(mw->get_actionNewEmpty()->trigger());
			
			for (int i=0; i<lwindow.size(); i++) // all viewers
			{
				viewer = (Viewer *)qobject_cast<QWidget *>(lwindow[i]->widget());
				if (viewer->getScenePtr()->get_polyhedron()->empty())
				{
					polyhedron_ptr_out = viewer->getScenePtr()->get_polyhedron();
				}
			}

			Boolean_Operations_ComponentPtr component_ptr = findOrCreateComponentForViewer<Boolean_Operations_ComponentPtr, Boolean_Operations_Component>(viewer, polyhedron_ptr_in1);
			component_ptr->Boolean_Minus(polyhedron_ptr_in1, polyhedron_ptr_in2, polyhedron_ptr_out);

			polyhedron_ptr_out->compute_bounding_box();
			polyhedron_ptr_out->compute_normals();
			polyhedron_ptr_out->compute_type();
			(void)polyhedron_ptr_out->calc_nb_components();
			(void)polyhedron_ptr_out->calc_nb_boundaries();

			viewer->show();

			viewer->updateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}



Q_EXPORT_PLUGIN2(mepp_component_Boolean_Operations_plugin, mepp_component_Boolean_Operations_plugin);

#endif // BUILD_component_Boolean_Operations
