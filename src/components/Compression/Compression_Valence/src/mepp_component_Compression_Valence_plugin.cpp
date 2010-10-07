#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "mepp_component_Compression_Valence_plugin.hxx"

#include "dialSettingsComp.hxx"
#include "dialSettingsDecomp.hxx"

#include <QObject>
#include <QAction>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QMdiSubWindow>

#include "Compression_Valence_Component.h"
typedef boost::shared_ptr<Compression_Valence_Component> Compression_Valence_ComponentPtr;

void mepp_component_Compression_Valence_plugin::OnCompress()
{
	//QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (component_ptr->IsCompressed == true)
			{
				/*(void)wxMessageBox(_T("Compression not possible\n\n")
					_T("The mesh is already compressed")

					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Compression not possible: the mesh is already compressed."));
				return;
			}
			if (!viewer->getScenePtr()->get_polyhedron()->is_pure_triangle())
			{
				/*(void)wxMessageBox(_T("Compression not possible\n\n")
					_T("The mesh owns non-triangular facets")

					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Compression not possible: the mesh owns non-triangular facets."));
				return;
			}

			// todo IMPORTANT: taille du mesh d'origine à déterminer et stocker ici pour le calcul du ratio !!!
			component_ptr->Initial_file_size = 1;
			if (FILE *file = fopen(viewer->getScenePtr()->currentFile().toStdString().c_str(), "r"))
			{
				fseek(file,0,SEEK_END);
				component_ptr->Initial_file_size = ftell(file);
				fclose(file);
			}

			////////////////////////////////////////////////////
			unsigned Number_vertices = 0;
			unsigned Qbit = 0;		
			unsigned Forget_metric_value = -1;	
			unsigned Connectivity_size = 0;
			unsigned Color_size = 0;
			unsigned Total_size = 0;
			unsigned Init_number_vertices = 0;
			unsigned Final_number_vertices = 0;
			unsigned Number_layers = 0;
			
			float Metric_threshold = -1.0;
			double Processing_time;	
			double Connectivity_rate = 0;
			double Color_rate = 0;
			double Geometry_rate = 0;
			double Total_rate = 0;	
			
			bool Is_normal_flipping_selected = false;	
			bool Is_use_metric_selected = false;
			bool Is_use_forget_metric_selected = false;
			bool Is_adaptive_quantization_selected;
			bool Is_compression_selected;

			/*char String_number_vertices[256];
			char String_MValue[256];
			char String_FValue[256];*/

			/*wxFileName Temp_fileName;

			wxString Path;
			Path = Temp_fileName.GetCwd();	

			MyDialog3 dial(m_frame);
			dial.m_directory->SetPath( Path );
			size_t slash = this->File_name.find_last_of("/\\");
			this->File_name = this->File_name.substr(slash+1);
			dial.m_file_name->SetValue( wxString::FromAscii( this->File_name.c_str()));


			if (dial.ShowModal() == wxID_OK)*/
			SettingsDialogComp dial;
			if (dial.exec() == QDialog::Accepted)
			{	
				QString fileName = dial.file_name->text();
				if (fileName.isEmpty())
					return;

				QApplication::setOverrideCursor(Qt::WaitCursor);

				/*wxString N_Path = dial.m_directory->GetPath();
				wxString N_File_Name = dial.m_file_name->GetValue();


				wxString File = N_Path + _T("\\") + N_File_Name; // File name is chosenS*/
				//wxString File = N_File_Name; // File name is chosenS
				Qbit = dial.quanti->value();

				Is_normal_flipping_selected = dial.normal_flipping->isChecked();

				if (dial.radioWithout->isChecked())
					Is_adaptive_quantization_selected = false;
				else
					Is_adaptive_quantization_selected = true;
				
				if (dial.radioCompression->isChecked())
				{
					component_ptr->IsCompressed = true;
					Is_compression_selected = true;
				}
				else
					Is_compression_selected = false;

				Is_use_metric_selected = dial.useMetric->isChecked();
				//strcpy(String_number_vertices, dial.m_number_vertices->GetValue().ToAscii());
				Number_vertices = dial.number_vertices->value();

				if (Is_use_metric_selected == true) // Decided to use metric to select vertices to be removed
				{
					//strcpy(String_MValue, dial.m_metric_value->GetValue().ToAscii());
					Metric_threshold = dial.metric_value->value();

					Is_use_forget_metric_selected = dial.forgetMetric->isChecked();

					if (Is_use_forget_metric_selected == true)
					{
						//strcpy(String_FValue, dial.m_forget_metric_value->GetValue().ToAscii());
						Forget_metric_value = dial.forget_metric_value->value();
					}
				}

				char File_Name[256];
				strcpy(File_Name, fileName.toStdString().c_str()); //strcpy(File_Name,File.ToAscii());

				//wxBusyInfo busy(_T("Proceeding Encoding"));

				Processing_time = component_ptr->Main_Function(*viewer->getScenePtr()->get_polyhedron(),
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
															  Number_layers, 
															  Init_number_vertices,
															  Final_number_vertices,
															  Connectivity_size,
															  Color_size,
															  Total_size, 
															  component_ptr->Initial_file_size);

				// To show result
				Connectivity_rate = (double)Connectivity_size / Init_number_vertices;
				Color_rate = (double)Color_size / Init_number_vertices;
				Total_rate = (double)Total_size * 8 / Init_number_vertices;
				Geometry_rate = Total_rate - Connectivity_rate - Color_rate;


				QString string1 = QString("Base mesh : %1 vertices \n").arg(Final_number_vertices, 3);
				//string1.Printf(_T("Base mesh : %3d vertices \n"), Final_number_vertices);

				QString string2 = QString("Connectivity : %1 b/v \n").arg(float(Connectivity_rate), 4, 'f', 3);
				//string2.Printf(_T("Connectivity : %4.3f b/v \n"),float(Connectivity_rate));

				QString string3("Geometry : ");
				//string3 += wxString::Format(_T("%4.3f"),float(Geometry_rate));
				string3 += QString("%1").arg(float(Geometry_rate), 4, 'f', 3);
				string3 += " b/v\n";

				
				QString string4("Color : ");
				//string4 += wxString::Format(_T("%4.3f"),float(Color_rate));
				string4 += QString("%1").arg(float(Color_rate), 4, 'f', 3);
				string4 += " b/v\n";
				
				QString string5("Total size : ");
				//string5 += wxString::Format(_T("%4.3f"),float(Total_rate));
				string5 += QString("%1").arg(float(Total_rate), 4, 'f', 3);
				string5 += " b/v\n";

				QString string6("Ratio : ");
				//string6 += wxString::Format(_T("%3.3f %% \n\n"),(float)Total_size / this->Initial_file_size * 100);
				string6 += QString("%1 % \n\n").arg((float)Total_size / component_ptr->Initial_file_size * 100, 3, 'f', 3);


				QString string7("Number of layers : ");
				//string7 += wxString::Format(_T("%d"),Number_layers);
				string7 += QString("%1").arg(Number_layers);
				string7 += "\n";

				QString string8("Calculation time : ");
				//string8 += wxString::Format(_T("%3.2f seconds \n"),float(Processing_time));
				string8 += QString("%1 seconds \n").arg(float(Processing_time), 3, 'f', 2);
				
				QApplication::restoreOverrideCursor();

				if (Is_compression_selected)
				{
					QString t = string1 + string2 + string3 + string4 + string5 + string6 + string7 + string8;
					/*(void)wxMessageBox((t)
					);*/
					QMessageBox::information(mw, APPLICATION, t);
				}

				/*m_frame->update_mesh_properties();
				m_frame->Refresh();
				m_frame->set_status_message(_T("Encoding is done"));*/
				mw->statusBar()->showMessage(tr("Encoding is done"));

				viewer->recreateListsAndUpdateGL();
			}
		}
	}

	//QApplication::restoreOverrideCursor();
}

int mepp_component_plugin_interface::load_file_from_component(PolyhedronPtr polyhedron_ptr, QString filename, Viewer* viewer)
{
	// here your code
	mepp_component_Compression_Valence_plugin *mepp_component_plugin = NULL;
	for (int i=0; i<viewer->lplugin.size(); ++i) {
		if (dynamic_cast<mepp_component_Compression_Valence_plugin*>(viewer->lplugin[i]) != 0) {
			mepp_component_plugin = dynamic_cast<mepp_component_Compression_Valence_plugin*>(viewer->lplugin[i]);
			//cout << "mepp_component_plugin found" << endl;
		}
	}

	if (mepp_component_plugin)
	{
		Compression_Valence_ComponentPtr component_ptr = mepp_component_plugin->findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);

		bool Is_decompression = false;
		unsigned Init_size = 1;

		/*fseek(file,0,SEEK_END);
		component_ptr->Initial_file_size = ftell(file);
		size_t slash = spath.find_last_of("/\\");
		component_ptr->File_name = spath.substr(slash+1);
		size_t point = comp->File_name.find_last_of('.');
		component_ptr->File_name.replace(point+1,3,"p3d");*/

		// MT
		component_ptr->File_name = filename.toStdString().c_str();

		//if (FILE *file2 = fopen(spath.c_str(), "r"))
		if (FILE *file2 = fopen(filename.toStdString().c_str(), "r"))
		{
			fseek(file2,0,SEEK_END);
			component_ptr->Compressed_file_size = ftell(file2);
			fclose(file2);
		}

		Is_decompression = true;
		//if (component_ptr)
		{					
			component_ptr->Total_layer=component_ptr->Decompress_Init(*polyhedron_ptr,Init_size,filename.toStdString().c_str());
			component_ptr->IsDecompress = true;
			component_ptr->Current_level = 0;					
			component_ptr->Initial_file_size = Init_size;
			//res_load=true;
		}
		/*else 
		  assert(false);	//todo: temp*/

		if (Is_decompression) // Ajout Ho.
		{
			float prog = (float)component_ptr->Calculate_Current_File_Size() / component_ptr->Compressed_file_size * 100;
			float ratio = 1/ ((float)component_ptr->Calculate_Current_File_Size() / Init_size);

			QString string("Number of all levels : ");
			//string += wxString::Format(_T("%d"),int(comp->Total_layer));
			string += QString("%1").arg(int(component_ptr->Total_layer));
			string += "   |   ";

			//string += wxString::Format(_T("Prog : %3.3f %%"), prog);
			string += QString("Prog : %1 %").arg(prog, 3, 'f', 3);
			string += "   |   ";

			//string += wxString::Format(_T("Ratio : %.3f \n"), ratio);
			string += QString("Ratio : %1 \n").arg(ratio, 0, 'f', 3);

			component_ptr->Prog.push_back(prog);
			component_ptr->Ratio.push_back(ratio);

			/*child->update_mesh_properties();
			child->Refresh();					
			child->set_status_message(string);*/
			((mainwindow *)viewer->getParent())->statusBar()->showMessage(string);
		}
	}
	else
	{
		//todo: temp
	}
	
	return 0;	//todo: code de retour
}
void mepp_component_Compression_Valence_plugin::load_P3D_file()
{
	emit(mw->get_mainwindowActionOpen()->doSendParamsOpen(tr("Open P3D File - from Valence"), tr("P3D files (*.p3d)"), Specific, mepp_component_plugin_interface::load_file_from_component));
}

void mepp_component_Compression_Valence_plugin::OnMouseLeftUp(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);

			if (!component_ptr->IsDecompress)
			{
				/*(void)wxMessageBox(_T("Decompression not possible\n\n")
					_T("You must execute Initialization before")
					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: you must execute Initialization before."));
				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;

			// read from file
			if (component_ptr->Sequence == false)
			{
				if (component_ptr->Current_level >= component_ptr->Total_layer)
					return;

				if (component_ptr->Process_level == 0)
					WriteInfo();

				component_ptr->Current_level = component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(0), component_ptr->File_name.c_str());
				if (component_ptr->Current_level > component_ptr->Process_level)
				{
					component_ptr->Process_level++;
					WriteInfo();
				}

				ShowText();
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

					ShowText();
				}

				else
				{
					if (component_ptr->Process_level == 0)
						WriteInfo();

					PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron; // MT
					component_ptr->Copy_Polyhedron.copy(viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level).get(), New_mesh.get());	// MT: get() 2 fois !!!

					component_ptr->Attibute_Seed_Gate_Flag(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);

					New_mesh->compute_normals();			
					
					vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
					viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level, New_mesh);
					
					component_ptr->Process_level++;
					component_ptr->Visu_level++;
					
					component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level),component_ptr->File_name.c_str());

					WriteInfo();

					float prog = (float)component_ptr->Calculate_Current_File_Size() / component_ptr->Compressed_file_size * 100;
					float ratio = 1/((float)component_ptr->Calculate_Current_File_Size() / component_ptr->Initial_file_size);

					component_ptr->Prog.push_back(prog);
					component_ptr->Ratio.push_back(ratio);

					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);

					ShowText();
				}

			}

			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnMouseRightUp(QMouseEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);

			if (component_ptr->Possible_change_sequence == true)
			component_ptr->Possible_change_sequence = false;

			// read from file
			if (component_ptr->Sequence == false)
			{

				if (component_ptr->Current_level <= 0)
					return;

				component_ptr->Stop_Decoder();

				component_ptr->Total_layer = component_ptr->Decompress_Init(*viewer->getScenePtr()->get_polyhedron(0), component_ptr->Initial_file_size, component_ptr->File_name.c_str());
				int Temp_level = component_ptr->Current_level - 1;

				component_ptr->Current_level = 0;
				while(component_ptr->Current_level != Temp_level)
					component_ptr->Current_level = component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(0), component_ptr->File_name.c_str());

				ShowText();
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

					ShowText();
				}
			}

			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnMouseWheel(QWheelEvent *event)
{
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		if (doesExistComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr)) // important !!!
		{
			int rot = event->delta();

			if (rot<0)
			{
				OnMouseLeftUp(NULL);	// todo: NULL
			}
			else
			{
				OnMouseRightUp(NULL);	// todo: NULL
			}
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnDecompress_all()
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{
				/*(void)wxMessageBox(_T("Decompression not possible\n\n")
					_T("You must execute Initialization before"));*/
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: you must execute Initialization before."));

				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;

			// read from file
			if (component_ptr->Sequence == false)
			{
				if (component_ptr->Process_level == 0)
					WriteInfo();
				while(component_ptr->Current_level != component_ptr->Total_layer)
				{
					component_ptr->Current_level = component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(0), component_ptr->File_name.c_str());
					/*if (this->Current_level > this->Process_level)
					{
						this->Process_level++;
						this->WriteInfo();
					}*/
				}

				ShowText();
			}
			// read from sequence
			else
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();

				while(component_ptr->Visu_level != component_ptr->Total_layer)
				{
					// if the next mesh already exists in the sequence
					if (component_ptr->Visu_level < component_ptr->Process_level)
						component_ptr->Visu_level++;

					else
					{
						if (component_ptr->Process_level == 0)
							WriteInfo();

						PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron;	// MT
						component_ptr->Copy_Polyhedron.copy(viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level).get(), New_mesh.get());	// MT: get() 2 fois !!!
						
						component_ptr->Attibute_Seed_Gate_Flag(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);
						
						New_mesh->compute_normals();

						vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
						viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level, New_mesh);

						component_ptr->Process_level++;
						component_ptr->Visu_level++;

						component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), component_ptr->File_name.c_str());

						WriteInfo();

						float prog = (float)component_ptr->Calculate_Current_File_Size() / component_ptr->Compressed_file_size * 100;
						float ratio = 1/((float)component_ptr->Calculate_Current_File_Size() / component_ptr->Initial_file_size);

						component_ptr->Prog.push_back(prog);
						component_ptr->Ratio.push_back(ratio);
					}

				}

				viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);
				ShowText();
			}

			viewer->recreateListsAndUpdateGL();
		}
	}

	QApplication::restoreOverrideCursor();
}

void mepp_component_Compression_Valence_plugin::OnDecompress_one_level()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{
				/*(void)wxMessageBox(_T("Decompression not possible\n\n")
					_T("You must execute Initialization before")
					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: you must execute Initialization before."));

				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;


			// read from file
			if (component_ptr->Sequence == false)
			{
				if (component_ptr->Current_level >= component_ptr->Total_layer)
					return;

				if (component_ptr->Process_level == 0)
					WriteInfo();

				component_ptr->Current_level = component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(),component_ptr->File_name.c_str());
				if (component_ptr->Current_level > component_ptr->Process_level)
				{
					component_ptr->Process_level++;
					WriteInfo();
				}

				ShowText();
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

					ShowText();
				}

				// Insert a mesh into the sequence
				else
				{
					if (component_ptr->Process_level == 0)
						WriteInfo();

					PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron;	// MT
					component_ptr->Copy_Polyhedron.copy(viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level).get(), New_mesh.get());	// MT: get() 2 fois !!!

					component_ptr->Attibute_Seed_Gate_Flag(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);
					New_mesh->compute_normals();

					vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT
					viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level, New_mesh);

					component_ptr->Process_level++;
					component_ptr->Visu_level++;

					component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level),component_ptr->File_name.c_str());

					WriteInfo();

					float prog = (float)component_ptr->Calculate_Current_File_Size() / component_ptr->Compressed_file_size * 100;
					float ratio = 1/((float)component_ptr->Calculate_Current_File_Size() / component_ptr->Initial_file_size);

					component_ptr->Prog.push_back(prog);
					component_ptr->Ratio.push_back(ratio);

					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);

					ShowText();
				}
			}

			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnDecompress_precedent_level()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{
				/*(void)wxMessageBox(_T("Decompression not possible\n\n")
					_T("You must execute Initialization before")
					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: you must execute Initialization before."));

				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;

			// read from file
			if (component_ptr->Sequence == false)
			{

				if (component_ptr->Current_level <= 0)
					return;

				component_ptr->Stop_Decoder();

				component_ptr->Total_layer = component_ptr->Decompress_Init(*viewer->getScenePtr()->get_polyhedron(), component_ptr->Initial_file_size, component_ptr->File_name.c_str());
				int Temp_level = component_ptr->Current_level - 1;

				component_ptr->Current_level = 0;
		
				while(component_ptr->Current_level != Temp_level)
					component_ptr->Current_level = component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(0), component_ptr->File_name.c_str());

				ShowText();
			}


			// read from sequence
			else
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();

				if (component_ptr->Visu_level <= 0)
					return;

				component_ptr->Visu_level--;
				viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);

				ShowText();

			}

			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnDecompress_go_to_specific_level()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (!component_ptr->IsDecompress)
			{
				/*(void)wxMessageBox(_T("Decompression not possible\n\n")
					_T("You must execute Initialization before")
					   );*/
				QMessageBox::information(mw, APPLICATION, tr("Decompression not possible: you must execute Initialization before."));

				return;
			}

			if (component_ptr->Possible_change_sequence == true)
				component_ptr->Possible_change_sequence = false;


			//char WL[256];
			int Wanted_level = 0;


			int CLevel = 0;
			if (component_ptr->Sequence)
			{
				component_ptr->Visu_level = viewer->getScenePtr()->get_current_polyhedron();
				CLevel = component_ptr->Visu_level;
			}
			else
				CLevel = component_ptr->Current_level;

			/*Dial_decomp dial(m_frame);
			wxString Current;
			Current.Printf(_T("%d"), CLevel);

			wxString Total;
			Total.Printf(_T("%d"), component_ptr->Total_layer);

			dial.m_staticText4->SetLabel(Current);
			dial.m_staticText5->SetLabel(Total);*/

			SettingsDialogDecomp dial;
			dial.currentLevel->setText(QString("%1").arg(int(CLevel)));
			dial.maxLevel->setText(QString("%1").arg(int(component_ptr->Total_layer)));
			if (dial.exec() == QDialog::Accepted)
			{
				//strcpy(WL,dial.m_textCtrl3->GetValue().ToAscii());
				Wanted_level = dial.wantedLevel->value();
			}
			else
				return;


			if (Wanted_level < 0)
				Wanted_level = 0;

			if (Wanted_level > component_ptr->Total_layer)
				Wanted_level = component_ptr->Total_layer;

			if (Wanted_level == CLevel)
				return;

			// read from file
			if (component_ptr->Sequence == false)
			{
				if (Wanted_level > CLevel)
				{

					while(component_ptr->Current_level != Wanted_level)
					{
						component_ptr->Current_level = component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(0),component_ptr->File_name.c_str());
					}
				}

				else if (Wanted_level < CLevel)
				{
					component_ptr->Stop_Decoder();

					component_ptr->Total_layer = component_ptr->Decompress_Init(*viewer->getScenePtr()->get_polyhedron(0), component_ptr->Initial_file_size, component_ptr->File_name.c_str());

					component_ptr->Decompress_To_Level(*viewer->getScenePtr()->get_polyhedron(0), Wanted_level);
					component_ptr->Current_level = Wanted_level;
				}
				ShowText();

			}
			//read from sequence
			else
			{
				if (Wanted_level > CLevel)
				{
					while(component_ptr->Visu_level != Wanted_level)
					{
						// if the next mesh already exists in the sequence
						if (component_ptr->Visu_level < component_ptr->Process_level)
							component_ptr->Visu_level++;

						else
						{
							if (component_ptr->Process_level == 0)
								WriteInfo();

							PolyhedronPtr New_mesh(new Polyhedron());	//Polyhedron * New_mesh = new Polyhedron;	// MT
							component_ptr->Copy_Polyhedron.copy(viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level).get(), New_mesh.get());	// MT: get() 2 fois !!!

							component_ptr->Attibute_Seed_Gate_Flag(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level), *New_mesh);
							New_mesh->compute_normals();
							//this->m_polyhedron_list->push_back(New_mesh);

							vector<PolyhedronPtr/*Polyhedron**/>::iterator it = viewer->getScenePtr()->get_begin_polyhedrons();	// MT			
							viewer->getScenePtr()->insert_polyhedron(it + component_ptr->Process_level, New_mesh);

							component_ptr->Process_level++;
							component_ptr->Visu_level++;

							component_ptr->Decompress_Each_Step(*viewer->getScenePtr()->get_polyhedron(component_ptr->Process_level),component_ptr->File_name.c_str());

							WriteInfo();

							float prog = (float)component_ptr->Calculate_Current_File_Size() / component_ptr->Compressed_file_size * 100;
							float ratio = 1/((float)component_ptr->Calculate_Current_File_Size() / component_ptr->Initial_file_size);

							component_ptr->Prog.push_back(prog);
							component_ptr->Ratio.push_back(ratio);
						}

					}
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Process_level);
					ShowText();
				}
				else
				{
					component_ptr->Visu_level = Wanted_level;
					viewer->getScenePtr()->set_current_polyhedron(component_ptr->Visu_level);

					ShowText();
				}

			}

			viewer->recreateListsAndUpdateGL();
		}
	}
}

void mepp_component_Compression_Valence_plugin::OnDecompress_mesh_sequence_on_off()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if ((component_ptr->Sequence) && (component_ptr->Possible_change_sequence))
			{
				component_ptr->Sequence = false;
				//m_frame->set_status_message(_T("Sequence : OFF"));
				mw->statusBar()->showMessage(tr("Sequence : OFF"));
			}
			else if ((!component_ptr->Sequence) && (component_ptr->Possible_change_sequence))
			{
				component_ptr->Sequence = true;
				//m_frame->set_status_message(_T("Sequence : ON"));
				mw->statusBar()->showMessage(tr("Sequence : ON"));
			}
		}
	}
}

void mepp_component_Compression_Valence_plugin::WriteInfo(void)
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			if (component_ptr->Process_level == 0)
			{
				component_ptr->Dec_File_Info = component_ptr->File_name;
				size_t point = component_ptr->Dec_File_Info.find('.');
				component_ptr->Dec_File_Info.replace(point+1,3,"txt");

				component_ptr->Dec_Info = fopen(component_ptr->Dec_File_Info.c_str(),"w");
			}
			else
				component_ptr->Dec_Info = fopen(component_ptr->Dec_File_Info.c_str(),"a");

			int CLevel = 0, Number_vertices = 0;

			if (component_ptr->Sequence)
			{
				CLevel = component_ptr->Visu_level;
				Number_vertices = (int)viewer->getScenePtr()->get_polyhedron(CLevel)->size_of_vertices();
			}
			else
			{
				CLevel = component_ptr->Current_level;
				Number_vertices = (int)viewer->getScenePtr()->get_polyhedron(0)->size_of_vertices();
			}

			unsigned Current_file_size = component_ptr->Calculate_Current_File_Size();

			float prog = (float)Current_file_size / component_ptr->Compressed_file_size * 100;
			
			
			if (component_ptr->Process_level == component_ptr->Total_layer)
				prog = 100.0;

			float ratio = 1 / ((float)Current_file_size / component_ptr->Initial_file_size);

			fprintf(component_ptr->Dec_Info,"Level %2d   #v : %8d      %6u bytes     Prog : %7.3f %%    Ratio : %9.3f\n", CLevel, Number_vertices, Current_file_size, prog, ratio);
			fclose(component_ptr->Dec_Info);
		}
	}
}

void mepp_component_Compression_Valence_plugin::ShowText(void)
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		Compression_Valence_ComponentPtr component_ptr = findOrCreateComponentForViewer<Compression_Valence_ComponentPtr, Compression_Valence_Component>(viewer, polyhedron_ptr);
		{
			int CLevel = -1;
			if (component_ptr->Sequence)
				CLevel = component_ptr->Visu_level;
			else
				CLevel = component_ptr->Current_level;


			QString string = QString("Current level : %1/%2").arg(CLevel, 2).arg(component_ptr->Total_layer, 2);
			//string += wxString::Format(_T("%2d"),CLevel);
			//string += _T("/");
			//string += wxString::Format(_T("%2d"),component_ptr->Total_layer);
			string += "   |   ";

			float prog =0, ratio = 0;

			if (component_ptr->Sequence)
			{
				prog = component_ptr->Prog[CLevel];
				ratio = component_ptr->Ratio[CLevel];
			}
			else
			{

				prog = (float)component_ptr->Calculate_Current_File_Size() / component_ptr->Compressed_file_size * 100;
				ratio = 1/((float)component_ptr->Calculate_Current_File_Size() / component_ptr->Initial_file_size);
			}

			if (CLevel != component_ptr->Total_layer)
				//string += wxString::Format(_T("Prog : %3.3f %% "), prog);
				string += QString("Prog : %1 % ").arg(prog, 3, 'f', 3);
			else
				string += "Prog : 100.000 %";

			string += "   |   ";
			//string += wxString::Format(_T("Ratio : %3.3f\n"), ratio);
			string += QString("Ratio : %1\n").arg(ratio, 3, 'f', 3);

			/*m_frame->update_mesh_properties();
			m_frame->Refresh();

			m_frame->set_status_message(string);*/
			mw->statusBar()->showMessage(string);
		}
	}
}

Q_EXPORT_PLUGIN2(mepp_component_Compression_Valence_plugin, mepp_component_Compression_Valence_plugin);

#endif
