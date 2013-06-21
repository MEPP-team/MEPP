#include <mepp_config.h>
#ifdef BUILD_component_Canonical

#include "mepp_component_Canonical_plugin.hxx"

#include "dialSettings.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "Canonical_Component.h"
typedef boost::shared_ptr<Canonical_Component> Canonical_ComponentPtr;

void mepp_component_Canonical_plugin::ValenceDrivenSimplification()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Canonical_ComponentPtr component_ptr = findOrCreateComponentForViewer<Canonical_ComponentPtr, Canonical_Component>(viewer, polyhedron_ptr);

		polyhedron_ptr->normalize_border();
		if(polyhedron_ptr->size_of_border_edges() != 0)
		{
			/*(void)wxMessageBox(_T("Simplification not possible\n\n")
				_T("The mesh owns border edges")

				   );*/
			QMessageBox::information(mw, APPLICATION, tr("Simplification not possible: the mesh owns border edges."));
			return;
		}
		
		if(polyhedron_ptr->is_pure_triangle() == false)
		{
			/*(void)wxMessageBox(_T("Simplification not possible\n\n")
				_T("The mesh is not triangular")

				   );*/
			QMessageBox::information(mw, APPLICATION, tr("Simplification not possible: the mesh is not triangular."));
			return;
		}

		//MyDialog3 dial(m_frame);

		bool Normal_Flipping = false;
		bool Use_Metric = false;
		bool Use_forget_metric = false;

		unsigned int Number_Vertices = 0;
		int Forget_Value = 0;
		float Metric_threshold = 0.0;

		/*char NVerticesChar[256];
		char MValue[256];
		char FValue[256];*/

		SettingsDialog dial;
		if (dial.exec() == QDialog::Accepted)
		{
			QApplication::setOverrideCursor(Qt::WaitCursor);

			Normal_Flipping = dial.normal_flipping->isChecked();//dial.m_normal_flipping->GetValue();

			Use_Metric = dial.useMetric->isChecked();//dial.m_use_metric->GetValue();

			//strcpy(NVerticesChar,dial.m_number_vertices->GetValue().ToAscii());
			Number_Vertices = dial.number_vertices->value();//atoi(NVerticesChar);

			if (Use_Metric == true)
			{
				//strcpy(MValue,dial.m_metric_value->GetValue().ToAscii());
				Metric_threshold = dial.metric_value->value();//atof(MValue);

				Use_forget_metric = dial.forgetMetric->isChecked();//dial.m_forget_metric->GetValue();
				if (Use_forget_metric == true)
				{
					//strcpy(FValue,dial.m_forget_metric_value->GetValue().ToAscii());
					Forget_Value = dial.forget_metric_value->value();//atoi(FValue);
				}
			}
			int Number_Of_Vertices1 = polyhedron_ptr->size_of_vertices();
			int Number_Of_Vertices2 = polyhedron_ptr->size_of_vertices();

			while(polyhedron_ptr->size_of_vertices() > Number_Vertices)
			{
				Number_Of_Vertices1 = polyhedron_ptr->size_of_vertices();

				component_ptr->Decimation_Conquest(polyhedron_ptr,Normal_Flipping,Use_Metric,Metric_threshold,
					Use_forget_metric,Forget_Value);
				component_ptr->Regulation(polyhedron_ptr,Normal_Flipping,Use_Metric,Metric_threshold,
					Use_forget_metric,Forget_Value);

				Number_Of_Vertices2 = polyhedron_ptr->size_of_vertices();
				if (Number_Of_Vertices1 == Number_Of_Vertices2)
					break;
			}

			/*m_frame->update_mesh_properties();
			m_frame->Refresh();
			m_frame->set_status_message(_T("Simplification is done"));*/
			mw->statusBar()->showMessage(tr("Simplification is done"));

			viewer->recreateListsAndUpdateGL();

		}
	}

	QApplication::restoreOverrideCursor();
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Canonical_plugin, mepp_component_Canonical_plugin);
#endif

#endif
