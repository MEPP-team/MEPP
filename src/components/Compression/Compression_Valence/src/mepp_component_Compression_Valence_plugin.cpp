///////////////////////////////////////////////////////////////////////////
// Author: Ho LEE
// Year: 2011
// Month: MAY
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <time.h>

#include "mepp_component_Compression_Valence_plugin.hxx"

#include "dialSettingsComp.hxx"
#include "dialSettingsDecomp.hxx"
#include "dialSettingsJCW.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>
#include <CGAL/Timer.h>
#include "Compression_Valence_Component.h"
typedef boost::shared_ptr<Compression_Valence_Component> Compression_Valence_ComponentPtr;

void mepp_component_Compression_Valence_plugin::OnCompress()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	//QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (component_ptr->IsCompressed)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Compression not possible, the mesh is already compressed."));
				return;
			}
			if (component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Compression not possible, the mesh is being decompressed."));
				return;
			}
			
			if (!viewer->getScenePtr()->get_polyhedron()->is_pure_triangle())
			{				
				QMessageBox::information(mw, APPLICATION, tr("Compression not possible, the mesh owns non-triangular facets."));
				return;
			}				

			// read parameters
			SettingsDialogComp dial(mw, mw->saveLocation);
			if (dial.exec() == QDialog::Accepted)
			{	
				QString fileName = dial.file_name->text();
				if (fileName.isEmpty())
					return;

				QApplication::setOverrideCursor(Qt::WaitCursor);
				unsigned Qbit = dial.quanti->value();

				bool Is_normal_flipping_selected = dial.normal_flipping->isChecked();
				bool Is_adaptive_quantization_selected;

				if (dial.radioWithout->isChecked())
					Is_adaptive_quantization_selected = false;
				else
					Is_adaptive_quantization_selected = true;
				
				bool Is_compression_selected = false;
				if (dial.radioCompression->isChecked())
				{
					component_ptr->IsCompressed = true;
					Is_compression_selected = true;
				}
				else
					Is_compression_selected = false;

				bool Is_use_metric_selected = dial.useMetric->isChecked();				
				int Number_vertices = dial.number_vertices->value();

				bool Metric_threshold, Is_use_forget_metric_selected, Forget_metric_value;
				Metric_threshold=Is_use_forget_metric_selected=Forget_metric_value=false;
				if (Is_use_metric_selected == true) 
				{					
					Metric_threshold = dial.metric_value->value();

					Is_use_forget_metric_selected = dial.forgetMetric->isChecked();

					if (Is_use_forget_metric_selected == true)
					{
						Forget_metric_value = dial.forget_metric_value->value();
					}
				}

				bool Is_bijection_selected = dial.Use_Bijection->isChecked();
				char File_Name[256];
				strcpy(File_Name, fileName.toStdString().c_str());
				
				// compress
				QString Result = component_ptr->Main_Function(*viewer->getScenePtr()->get_polyhedron(),
															  viewer->getScenePtr()->currentFile().toStdString().c_str(),
														      File_Name,
															  Qbit,
															  Number_vertices,
															  Is_normal_flipping_selected,
			  												  Is_use_metric_selected,
															  Metric_threshold,
															  Is_use_forget_metric_selected,
															  Forget_metric_value,
															  Is_compression_selected,
															  Is_adaptive_quantization_selected,
															  Is_bijection_selected);															 
				
				QApplication::restoreOverrideCursor();

				if (Is_compression_selected)
				{				
					QMessageBox::information(mw, APPLICATION, Result);
				}				
				
				mw->statusBar()->showMessage(tr("Encoding is done"));
				viewer->recreateListsAndUpdateGL();	
			}
		}
	}	
}

int mepp_component_plugin_interface::load_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// here your code
	mepp_component_Compression_Valence_plugin *mepp_component_plugin = NULL;
	for (int i=0; i<viewer->lplugin.size(); ++i) {
		if (dynamic_cast<mepp_component_Compression_Valence_plugin*>(viewer->lplugin[i]) != 0) {
			mepp_component_plugin = dynamic_cast<mepp_component_Compression_Valence_plugin*>(viewer->lplugin[i]);			
		}
	}

	if (mepp_component_plugin)
	{
		Compression_Valence_ComponentPtr component_ptr = mepp_component_plugin->findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		
		component_ptr->File_name = filename.toStdString().c_str();				
							
		QString string = component_ptr->Decompress_Init(*polyhedron_ptr);
		viewer->getScenePtr()->set_loadType(Time);
		viewer->setFps(5);
		((mainwindow *)viewer->getParent())->statusBar()->showMessage(string);		
	}	
	
	return 0;	//todo: code de retour
}
void mepp_component_Compression_Valence_plugin::load_P3D_file()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	emit(mw->get_mainwindowActionOpen()->doSendParamsOpen(tr("Open P3D File - from Valence"), tr("P3D files (*.p3d)"), Specific, mepp_component_plugin_interface::load_file_from_component));
}

void mepp_component_Compression_Valence_plugin::OnMouseLeftUp(QMouseEvent *event)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);

			/*if (component_ptr->IsCompressed)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Compression not possible, the mesh is already compressed."));
				return;
			}*/
			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: please load .p3d file first."));
				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;

			// read from file
			if (!component_ptr->Sequence)
			{
				component_ptr->Decompression_From_File(*viewer->getScenePtr()->get_polyhedron(0));				
			}

			// from sequence
			else
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();

				if (component_ptr->Visu_level >= component_ptr->Total_layer)
					return;
				
				// if the next mesh already exists in the sequence				
				if (component_ptr->Visu_level < component_ptr->Process_level)
				{
					component_ptr->Visu_level++;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);
				}

				else
				{
					PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron; // MT

					component_ptr->Decompression_From_Sequence(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);					
					vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
					viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level-1, New_mesh);
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);
				}
			}
			ShowText();
			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnMouseRightUp(QMouseEvent *event)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
			
			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: please load .p3d file first."));
				return;
			}

			if (component_ptr->Possible_change_sequence == true)
			component_ptr->Possible_change_sequence = false;

			// read from file
			if (component_ptr->Sequence == false)
			{
				component_ptr->Decompression_Coarser_From_File(*viewer->getScenePtr()->get_polyhedron(0));				
			}
			// read from sequence
			else
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();

				if (component_ptr->Visu_level == 0)
					return;
				else
				{
					component_ptr->Visu_level--;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);				
				}
			}
			ShowText();
			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnMouseWheel(QWheelEvent *event)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			int rot = event->delta();

			if (rot<0)
				OnMouseLeftUp(NULL);	// todo: NULL
			else
				OnMouseRightUp(NULL);	// todo: NULL
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnDecompress_all()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	Timer timer;
	timer.start();

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{

			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: please load .p3d file first."));
				return;
			}			

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;

			// read from file
			if (!component_ptr->Sequence)
			{
				component_ptr->Decompression_All_From_File(*viewer->getScenePtr()->get_polyhedron(0));				
			}

			// read from sequence
			else
			{
				while(component_ptr->Process_level != component_ptr->Total_layer)
				{
					PolyhedronPtr New_mesh(new Polyhedron());
					component_ptr->Decompression_From_Sequence(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);					
					vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
					viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level-1, New_mesh);				
				}
				viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);				
			}
			ShowText();
			viewer->recreateListsAndUpdateGL();
		}
	}
	timer.stop();
	double Time = timer.time();
	QString Total_time = QString("Processing time : %1 s \n\n").arg(double(Time), 4, 'f', 3);			
	QMessageBox::information(mw, APPLICATION, Total_time);
	//QApplication::restoreOverrideCursor();
}

void mepp_component_Compression_Valence_plugin::OnDecompress_one_level()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;

	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			OnMouseLeftUp(NULL);	// todo: NULL				
		}
	}	
}

void mepp_component_Compression_Valence_plugin::OnDecompress_precedent_level()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;

	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			OnMouseRightUp(NULL);	
		}
	}	
}

void mepp_component_Compression_Valence_plugin::OnDecompress_go_to_specific_level()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: please load .p3d file first."));
				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;
			
			int Wanted_level = 0;
			int CLevel = 0;
			if (component_ptr->Sequence)
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();
				CLevel = component_ptr->Visu_level;
			}
			else
				CLevel = component_ptr->Current_level;			

			SettingsDialogDecomp dial;
			dial.currentLevel->setText(QString("%1").arg(int(CLevel)));
			dial.maxLevel->setText(QString("%1").arg(int(component_ptr->Total_layer)));

			if (dial.exec() == QDialog::Accepted)
			{				
				Wanted_level = dial.wantedLevel->value();
			}
			else
				return;	
			if(Wanted_level < 0)
				Wanted_level = 0;
			if(Wanted_level > component_ptr->Total_layer)
				Wanted_level = component_ptr->Total_layer;

			// read from file
			if (component_ptr->Sequence == false)
			{
				component_ptr->Decompression_Specific_Level_From_File(*viewer->getScenePtr()->get_polyhedron(0), Wanted_level);				
			}
			//read from sequence
			else
			{
				if (Wanted_level > component_ptr->Process_level)
				{
					while(component_ptr->Process_level != Wanted_level)
					{
						PolyhedronPtr New_mesh(new Polyhedron());
						component_ptr->Decompression_From_Sequence(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);					
						vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
						viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level-1, New_mesh);				
					}
					component_ptr->Visu_level = Wanted_level;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);						
				}
				else
				{
					component_ptr->Visu_level = Wanted_level;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);					
				}
			}
		}
		ShowText();
		viewer->recreateListsAndUpdateGL();		
	}
}

void mepp_component_Compression_Valence_plugin::OnDecompress_mesh_sequence_on_off()
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();		

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Switch of decompression mode not possible: please load .p3d file first."));
				return;
			}

			if ((component_ptr->Sequence) && (component_ptr->Possible_change_sequence))
			{
				component_ptr->Sequence = false;				
				mw->statusBar()->showMessage(tr("Generation of mesh sequence : OFF"));
			}
			else if ((!component_ptr->Sequence) && (component_ptr->Possible_change_sequence))
			{
				component_ptr->Sequence = true;
				mw->statusBar()->showMessage(tr("Generation of mesh sequence : ON"));
			}
			if (!component_ptr->Possible_change_sequence)
			{
				mw->statusBar()->showMessage(tr("Mode change can not be modified."));
			}

		}
	}
}
void mepp_component_Compression_Valence_plugin::ShowText(void)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			QString string = component_ptr->Show_Text();
			mw->statusBar()->showMessage(string);
		}
	}
}


void mepp_component_Compression_Valence_plugin::OnJCW(void)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// active viewer
	if (mw->activeMdiChild() != 0)
	{		
		Timer timer;
		timer.start();

		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		if (!viewer->getScenePtr()->get_polyhedron()->is_pure_triangle())
		{				
			QMessageBox::information(mw, APPLICATION, tr("JCW not possible: the mesh owns non-triangular facets."));
			return;
		}
		if (component_ptr->IsCompressed)
		{				
			QMessageBox::information(mw, APPLICATION, tr("JCW not possible, the mesh is already compressed."));
			return;
		}
		if (component_ptr->IsDecompress)
		{				
			QMessageBox::information(mw, APPLICATION, tr("JCW not possible, the mesh is being decompressed."));
			return;
		}
		
		//QApplication::setOverrideCursor(Qt::WaitCursor);

		SettingsDialogJCW dial(mw, mw->saveLocation);
		if (dial.exec() == QDialog::Accepted)
		{
			QString fileName = dial.file_name->text();

			QApplication::setOverrideCursor(Qt::WaitCursor);

			int Qbit = dial.quanti->value();
			int NVertices = dial.number_vertices->value();

			int Number_bins = dial.number_bins->value();
			int Number_regions = dial.number_regions->value();
			int Embedding_strength = dial.embedding_strength->value();
			QString Embedding_message = dial.embedding_message->text();

			bool Is_complete_reversibility_selected = dial.complete_reversibility->isChecked();				
			bool Is_divide_regions_selected = dial.divide_regions->isChecked();
			int Thres_divide_regions = 0;
			if (Is_divide_regions_selected)
				Thres_divide_regions = dial.thres_divide_regions->value();			

			bool Normal_flipping = false;
			bool Use_metric = false;
			float Metric_thread = 0.25;
			bool Use_forget_metric = false;
			int Forget_value = 3000;
			
			QString Resulting_message = component_ptr->Joint_Compression_Watermarking(*viewer->getScenePtr()->get_polyhedron(),
																					   viewer->getScenePtr()->currentFile().toStdString().c_str(),
																					   fileName.toStdString().c_str(),
																					   Number_bins,
																					   Number_regions,
 																					   Embedding_strength,
																					   Embedding_message.toStdString().c_str(),
																					   Is_complete_reversibility_selected,
																					   Is_divide_regions_selected,
																					   Thres_divide_regions,
																					   Qbit,
																					   NVertices,
												   									   Normal_flipping,
																					   Use_metric,
																					   Metric_thread,
																					   Use_forget_metric,
																					   Forget_value);														 
			//component_ptr->Compression(*viewer->getScenePtr()->get_polyhedron(), "output.p3d", Qbit, Connectivity_size, Color_size, Total_size, Initial_file_size);			
			
			QApplication::restoreOverrideCursor();		
			QMessageBox::information(mw, APPLICATION, Resulting_message);
			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnJCWdecompress(void)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: please load .p3d file first."));
				return;
			}		

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;


			// read from file
			if (component_ptr->Sequence == false)
			{
				component_ptr->JCW_Decompression_From_File(*viewer->getScenePtr()->get_polyhedron(0));				
			}

			// from sequence
			else
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();

				if (component_ptr->Visu_level >= component_ptr->Total_layer)
					return;

				// if the next mesh already exists in the sequence
				if (component_ptr->Visu_level < component_ptr->Process_level)
				{
					component_ptr->Visu_level++;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);					
				}

				// Insert a mesh into the sequence
				else
				{					

					PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron;	// MT

					component_ptr->JCW_Decompression_From_Sequence(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);					
					vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
					viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level-1, New_mesh);
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);					
				}
			}
			ShowText();
			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnJCWdecompress_without_extraction(void)
{
	std::cout << "entering " << __PRETTY_FUNCTION__ << std::endl;
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{				
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: please load .p3d file first."));
				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;


			// read from file
			if (component_ptr->Sequence == false)
			{
				component_ptr->JCW_Decompression_Without_Extraction_From_File(*viewer->getScenePtr()->get_polyhedron());//,component_ptr->File_name.c_str());				
			}

			// from sequence
			else
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();

				if (component_ptr->Visu_level >= component_ptr->Total_layer)
					return;

				// if the next mesh already exists in the sequence
				if (component_ptr->Visu_level < component_ptr->Process_level)
				{
					component_ptr->Visu_level++;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);					
				}

				// Insert a mesh into the sequence
				else
				{
					PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron;	// MT

					component_ptr->JCW_Decompression_Without_Extraction_From_Sequence(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);					
					vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
					viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level-1, New_mesh);
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);					
				}
			}
			ShowText();
			viewer->recreateListsAndUpdateGL();
		}
	}
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mepp_component_Compression_Valence_plugin, mepp_component_Compression_Valence_plugin);
#endif

#endif
