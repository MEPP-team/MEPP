///////////////////////////////////////////////////////////////////////////
// Author: Ho LEE
// Year: 2011
// Month: MAY
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "Compression_Valence_Component.h"
#include "Compression_Valence_Polyhedron.h"
#include "Matrix3X3.h"
#include "Compression_Valence_Common.h"
#include <CGAL/Timer.h>
#include "Processing_Kai.h"
#include <map>
#include <set>
#include <bitset>

#define COLOR_NUMBER 10000
#define USE_COLOR_METRIC
#define AC_BUFFER 1024 * 10000

const int MINIMUM_PREDICTION_NUMBER = 3;
const int LIMIT_NUMBER = 50;

QString Compression_Valence_Component::Main_Function(Polyhedron     & _pMesh,
													const char *     _Input_File_Name,
													const char*      _File_Name,
													const int      & _Qbit,
													const int      & _NVertices,
													const bool       _Normal_flipping,
													const bool       _Use_metric,
													const float    & _Metric_thread,
													const bool       _Use_forget_metric,
													const int      & _Forget_value, 
													const bool       _Compression_selected,
													const bool       _Adaptive_quantization,
													const bool       _Is_bijection_selected)												
{	
	Timer timer;
	timer.start();	

	// read size of input mesh (file size)
	if (FILE *file = fopen(_Input_File_Name, "r"))
	{
		fseek(file,0,SEEK_END);
		this->Initial_file_size = ftell(file);
		fclose(file);
	}

	// Use of bijection or not for the geometry encoding
	this->Is_Bijection_Enabled = _Is_bijection_selected;
	
	unsigned Init_number_vertices = (unsigned)_pMesh.size_of_vertices();
	
	// Initialization - Quantization, Color, Multiple components
	this->Global_Initialization(_pMesh, _Qbit, _File_Name);
		
	// When use of adaptive quantization is selected
	if (_Adaptive_quantization)	
		this->Adaptive_Quantization(_pMesh, _NVertices, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value, _Qbit);
	else
		this->Simplification(_pMesh, _NVertices, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value);		
	
	unsigned Connectivity_size=0, Color_size=0, Total_size=0;

	// Compression
	if (_Compression_selected)
		this->Compression(_pMesh, _File_Name, _Qbit, Connectivity_size, Color_size, Total_size);//, this->Initial_file_size);		

	unsigned Number_layers = this->GlobalCountOperation;
	unsigned Final_number_vertices = (unsigned)_pMesh.size_of_vertices();

	timer.stop();

	// To show result
	double Connectivity_rate = (double)Connectivity_size / Init_number_vertices;
	double Color_rate = (double)Color_size / Init_number_vertices;
	double Total_rate = (double)Total_size * 8 / Init_number_vertices;
	double Geometry_rate = Total_rate - Connectivity_rate - Color_rate;

	QString Res = QString("Base mesh : %1 vertices \n").arg(Final_number_vertices, 3);
	Res += QString("Connectivity : %1 b/v \n").arg(float(Connectivity_rate), 4, 'f', 3);	
	Res += QString("Geometry : ");
	Res += QString("%1").arg(float(Geometry_rate), 4, 'f', 3);
	Res += " b/v\n";	
	Res += QString("Color : ");				
	Res += QString("%1").arg(float(Color_rate), 4, 'f', 3);
	Res += " b/v\n";	
	Res += QString("Total size : ");				
	Res += QString("%1").arg(float(Total_rate), 4, 'f', 3);
	Res += " b/v\n";
	Res += QString("Ratio : ");
	Res += QString("%1 % \n\n").arg((float)Total_size / this->Initial_file_size * 100, 3, 'f', 3);
	Res += QString("Number of layers : ");
	Res += QString("%1").arg(Number_layers);
	Res += "\n";
	Res += QString("Calculation time : ");
	Res += QString("%1 seconds \n").arg(float(timer.time()), 3, 'f', 2);
	return Res;
}

// Description : To select the input gate.
void Compression_Valence_Component::Global_Initialization(Polyhedron & _pMesh, 
														  const int  & _Qbit,
														  const char * _File_Name)
{
	// (1) Determination if the mesh is colored.
	// (2) Conversion of color space.
	// (3) Quantization of converted colors into "vertex->color_int()"
	// (4) Re-paint mesh with re_calculated rgb colors from quantized converted color.
	// (5) Establish color palette.
	// (6) Color quantization - Descrease of possible color numbers 
	this->Color_Initialization(_pMesh);	
	
	// Initialization of multiple components (quantization is performed separately for each component)
	this->Multiple_Components_Initialization(_pMesh, _Qbit);		

	// Quantization of each component
	this->Quantization(_pMesh);
}	

void Compression_Valence_Component::Multiple_Components_Initialization(Polyhedron & _pMesh, 
																	   int const  & _Qbit)
{
	// Initialize vertex flags
	for (Vertex_iterator pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
	{
		pVertex->Seed_Edge = OTHER_COORDINATE;	
		pVertex->Component_Number = -1;
	}	
	
	// Several components
	// (1) To get the number of components;
	// (2) To get the volume, area and the number of vertices of each component;
	// (3) To know if each component is closed or not.
	// (4) Tag vertices and facets to correspoding component number

	_pMesh.tag_facets(-1);
	int Component_index = 0;

	for (Facet_iterator pFacet = _pMesh.facets_begin(); pFacet != _pMesh.facets_end(); pFacet++)
	{
		if (pFacet->tag() == -1)
		{	
			bool Is_closed = true;
			bool Is_seed_edge_found = false;

			float xmin =  50000., ymin =  50000., zmin =  50000.;
			float xmax = -50000., ymax = -50000., zmax = -50000.;
			double area = 0;
			int Number_vertices = 0; 

			pFacet->tag(Component_index);

			std::list<Facet_handle> facets;
			facets.push_front(pFacet);

			while(!facets.empty())
			{
				Facet_handle F = facets.front();
				facets.pop_front();
				
				F->tag(Component_index);
				
				// tag component number to facet
				F->Component_Number = Component_index;
				
				area += Area_Facet_Triangle(F->halfedge());			

				Halfedge_around_facet_circulator pHalfedge = F->facet_begin();
				Halfedge_around_facet_circulator end = pHalfedge;
				
				CGAL_For_all(pHalfedge, end)
				{
					// tag the vertex to its corresponding component number
					if (pHalfedge->vertex()->Component_Number == -1)
					{
						pHalfedge->vertex()->Component_Number = Component_index;
						Number_vertices++;
					}			
					
					// To calculate the bounding box of each component
					if (pHalfedge->vertex()->point().x() > xmax)
						xmax = pHalfedge->vertex()->point().x();

					if (pHalfedge->vertex()->point().y() > ymax)
						ymax = pHalfedge->vertex()->point().y();

					if (pHalfedge->vertex()->point().z() > zmax)
						zmax = pHalfedge->vertex()->point().z();

					if (pHalfedge->vertex()->point().x() < xmin)
						xmin = pHalfedge->vertex()->point().x();

					if (pHalfedge->vertex()->point().y() < ymin)
						ymin = pHalfedge->vertex()->point().y();

					if (pHalfedge->vertex()->point().z() < zmin)
						zmin = pHalfedge->vertex()->point().z();
					
					// To know if the component is closed or not
					if (pHalfedge->is_border_edge())
						Is_closed = false;						
					
					if (!Is_seed_edge_found)
					{
						if ( ( (!pHalfedge->is_border_edge()) && (!Is_Border_Vertex(pHalfedge)) && (!Is_Border_Vertex(pHalfedge->opposite())) ) || (pHalfedge->next()->vertex_degree() != 6) )
						{
							Is_seed_edge_found = true;
							
							//Seed edge of each component
							pHalfedge->vertex()->Seed_Edge = 2 * Component_index;
							pHalfedge->opposite()->vertex()->Seed_Edge = 2 * Component_index + 1;							
						}
					}

					Facet_handle pNFacet = pHalfedge->opposite()->facet();
					if (pNFacet != NULL && pNFacet->tag() == -1)
					{
						facets.push_front(pNFacet);
						pNFacet->tag(Component_index);
						pNFacet->Component_Number = Component_index;
					}
				}
			}

			this->xmin.push_back(xmin);
			this->ymin.push_back(ymin);
			this->zmin.push_back(zmin);

			this->xmax.push_back(xmax);
			this->ymax.push_back(ymax);
			this->zmax.push_back(zmax);
			
			double HighestBB = -5000;
			if(xmax - xmin > HighestBB)
				HighestBB = xmax - xmin;
			if(ymax - ymin > HighestBB)
				HighestBB = ymax - ymin;
			if(zmax - zmin > HighestBB)
				HighestBB = zmax - zmin;
			
			this->HighestLengthBB.push_back(HighestBB);

			// volume, area, number of vertices of each component
			Vector e1(xmax-xmin, 0., 0.);
			Vector e2(0., ymax-ymin, 0.);

			Vector normal = CGAL::cross_product(e1,e2);
			double base_area = sqrt(normal * normal);
			double volume = base_area * (zmax - zmin);	
			
			// For quasi-optimal determination of quantization precision. 
			// Originally, the highest dimension of bounding box was considered as 10.
			// We shoud take this point into account.
			volume *= pow(((double)10.0/(double)HighestBB), 3.0);
			area *= pow(((double)10.0/(double)HighestBB), 2.0);
			
			// Stock information for each component
			this->ComponentVolume.push_back(volume);
			this->ComponentArea.push_back(area);
			this->ComponentNumberVertices.push_back(Number_vertices);

			// To know if each component is open or closed.
			this->IsClosed.push_back(Is_closed);

			Component_index++;
		}
	}

	this->NumberComponents = Component_index;	

	list<int> li;
	list<Point_Int> pi;
	list<Color_Unit> cu;
	
	// Initilization of containers for each component
	for (int i = 0; i < this->NumberComponents; i++)
	{		
		//operation (decimation or quantization)
		this->ListOperation.push_back(li);
		
		//connectivity
		this->Connectivity.push_back(li);

		//geometry
		this->Geometry.push_back(pi);

		//vertex color
		this->VertexColor.push_back(cu);
		
		// Number of connectivity and geometry symbols
		this->NumberSymbol.push_back(li);
		this->NumberVertices.push_back(li);
		
		// Range of alpha and gamma
		this->AlphaRange.push_back(li);
		this->AlphaOffset.push_back(li);
		this->GammaRange.push_back(li);
		this->GammaOffset.push_back(li);
	
		// Number of operations for each component
		this->ComponentOperations.push_back(0);

		// Number of decimations for each component
		this->NumberDecimation.push_back(0);
		
		// Number of under quantization
		this->NumberChangeQuantization.push_back(0);		
		this->NumberColorQuantization.push_back(0);

		// Qbit for each component
		this->Qbit.push_back(_Qbit);
		
		// Displacement vector for under quantization for each component
		this->QuantizationCorrectVector.push_back(li);

		// The number of vertices for each under_quantization
		this->NumberQuantizationLayer.push_back(li);
		
		// For color quantization change
		this->NumberProcessedVertices.push_back(li);
		this->ColorChildcellIndex.push_back(li);
		this->ColorEncoderIndex.push_back(li);
	}
}

/*
	Description : Quantize all vertices so that the new positions 
	are reguliraly spaced in the 3D space. */
void Compression_Valence_Component::Quantization(Polyhedron & _pMesh) 
{
	// Quantization step for each component
	for (int i = 0; i < this->NumberComponents; i++)
	{
		float max = this->xmax[i] - this->xmin[i];
		
		if (this->ymax[i] - this->ymin[i] > max)
			max = this->ymax[i] - this->ymin[i];
		
		if (this->zmax[i] - this->zmin[i] > max)
			max = this->zmax[i] - this->zmin[i];
		
		int NbInteraval = pow(2., (int)this->Qbit[i]);

		float Q_Step = max / (float)NbInteraval;
		this->Quantization_Step.push_back(Q_Step);				
	}	

	// Vertex quantization
	for (Vertex_iterator pVert = _pMesh.vertices_begin();pVert != _pMesh.vertices_end(); pVert++)
	{
		double x = pVert->point().x();
		double y = pVert->point().y();
		double z = pVert->point().z();

		int Component_ID = pVert->Component_Number;

		int Qx = (int)(ceil((x - (double)this->xmin[Component_ID]) / (double)this->Quantization_Step[Component_ID])) - 1;
		if (Qx == -1)
			Qx = 0;
		int Qy = (int)(ceil((y - (double)this->ymin[Component_ID]) / (double)this->Quantization_Step[Component_ID])) - 1;
		if (Qy == -1)
			Qy = 0;
		int Qz = (int)(ceil((z - (double)this->zmin[Component_ID]) / (double)this->Quantization_Step[Component_ID])) - 1;
		if (Qz == -1)
			Qz = 0;
		
		pVert->point() = Point3d(this->xmin[Component_ID] + (Qx + 0.5) * this->Quantization_Step[Component_ID],
							     this->ymin[Component_ID] + (Qy + 0.5) * this->Quantization_Step[Component_ID],
							     this->zmin[Component_ID] + (Qz + 0.5) * this->Quantization_Step[Component_ID]);		
	}
}



// this->ColorArray -> contains all initial colors present in the input mesh.
void Compression_Valence_Component::Color_Initialization(Polyhedron &_pMesh)
{	
	
	// (1) To determine if the mesh is colored.
	//     We consider the mesh is colored if the number of colors are > 2.	
	Vertex_iterator pVertex = _pMesh.vertices_begin();	
	
	this->OnlyColor[0] = pVertex->color(0);
	this->OnlyColor[1] = pVertex->color(1);
	this->OnlyColor[2] = pVertex->color(2);
	
	for (; pVertex != _pMesh.vertices_end(); pVertex++)
	{	
		if ((pVertex->color(0) != this->OnlyColor[0]) || (pVertex->color(1) != this->OnlyColor[1]) || (pVertex->color(2) != this->OnlyColor[2]))
		{
			this->IsColored = true;
			break;
		}		
	}
	
	// Even there is not > 2 colors, one unique color can exist.
	if (!this->IsColored)
	{
		if ((this->OnlyColor[0] != 0.5)  || (this->OnlyColor[1] != 0.5) || (this->OnlyColor[2] != 0.5))
		{
			this->IsColored = true;
			this->IsOneColor = true;
		}
	}

	// (2) If the mesh is colored -> Color initialization.
	//     Conversion to Lab colors and determine max and min for each components.
	if ((this->IsColored) && (!this->IsOneColor))
	{
		float Temp_color[3];
		
		float C0_min = 5000, C1_min = 5000, C2_min = 5000;
		float C0_max = -5000, C1_max = -5000, C2_max = -5000;
		
			// calculate new color and find min values.
		for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
		{
			Temp_color[0] = pVertex->color(0);
			Temp_color[1] = pVertex->color(1);
			Temp_color[2] = pVertex->color(2);			
				
			// Color space conversion.
			float New_color[3];							
			RGB_To_LAB(Temp_color[0], Temp_color[1], Temp_color[2], New_color);				
			
			// assignment of color in new color space in order to quantize later.
			pVertex->color(New_color[0], New_color[1], New_color[2]);
			
			// calculate min and max to identify quantization step.
			if (New_color[0] > C0_max)
				C0_max = New_color[0];		
			if (New_color[0] < C0_min)
				C0_min = New_color[0];
					
			if (New_color[1] > C1_max)
				C1_max = New_color[1];
			if (New_color[1] < C1_min)
				C1_min = New_color[1];
					
			if (New_color[2] > C2_max)
				C2_max = New_color[2];
			if (New_color[2] < C2_min)
				C2_min = New_color[2];			
		}
		
		const int Nb_interaval = (int)pow(2.0, C0_QUANTIZATION) - 1;
		
		float Color_max = C0_max - C0_min;
		if ( C1_max - C1_min > Color_max)
			Color_max = C1_max - C1_min;
		if ( C2_max - C2_min > Color_max)
			Color_max = C2_max - C2_min;

		//Information needed to quantize converted color.
		this->C0_Min = C0_min;
		this->C1_Min = C1_min;
		this->C2_Min = C2_min;				

		// Step size of color quantization
		this->Color_Quantization_Step = (float)((Color_max) / Nb_interaval);
		
		//Enter quantized color valued into "vertex->color_int" to use lated.
		//Also, recalculate vertex color using these quantized color values.
		Color_Unit Resulting_color;		
		float New_vertex_color[3];
		float Reconstructed_color[3];

		//int Color_index = 0;
		for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
		{
			Temp_color[0] = pVertex->color(0);
			Temp_color[1] = pVertex->color(1);
			Temp_color[2] = pVertex->color(2);

			int Qc0 = (int)(floor((Temp_color[0] - C0_min) / this->Color_Quantization_Step + 0.5));		
			int Qc1 = (int)(floor((Temp_color[1] - C1_min) / this->Color_Quantization_Step + 0.5));
			int Qc2 = (int)(floor((Temp_color[2] - C2_min) / this->Color_Quantization_Step + 0.5));				
			
			Resulting_color.c0 = Qc0;
			Resulting_color.c1 = Qc1;
			Resulting_color.c2 = Qc2;
			pVertex->color_int(Resulting_color.c0, Resulting_color.c1, Resulting_color.c2);			

			New_vertex_color[0] = this->C0_Min + Resulting_color.c0 * this->Color_Quantization_Step;
			New_vertex_color[1] = this->C1_Min + Resulting_color.c1 * this->Color_Quantization_Step;
			New_vertex_color[2] = this->C2_Min + Resulting_color.c2 * this->Color_Quantization_Step;

			LAB_To_RGB(New_vertex_color[0], New_vertex_color[1], New_vertex_color[2], Reconstructed_color);			
			for (int i = 0 ; i < 3; i++)
			{
				if (Reconstructed_color[i] < 0.)
					Reconstructed_color[i] = 0.;
				if (Reconstructed_color[i] > 1.)
					Reconstructed_color[i] = 1.;
			}
			
			// re-paint the input mesh with reconstructed colors from Lab to RGB transformation.
			pVertex->color(Reconstructed_color[0], Reconstructed_color[1], Reconstructed_color[2]);									

		}	
	}
}


void Compression_Valence_Component::Adaptive_Quantization(Polyhedron  & _pMesh, 
													      const int   & _NVertices, 
														  const bool    _Normal_flipping,
														  const bool    _Use_metric,
														  const float & _Metric_thread,
														  const bool    _Use_forget_metric,
														  const int   & _Forget_value,
														  const int   & _Qbit)
{


	if(!this->IsColored)
	{
		//bool Is_any_vertex_removed = true;
		unsigned Last_Number = 0;
		unsigned Current_Number = _pMesh.size_of_vertices();			
		int Operation_choice = -1;

		//int Old_Q;
		//double mrms, mrmswrtBB, hausdorff, hausdorffwrtBB;
		bool Continue;
		do
		{	
			Last_Number = Current_Number;		
			Continue = false;

			for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
			{
				// if the ith component did not remove any vertex in last loop, it is not necessary to simplifiy it.
				if (this->ComponentOperations[Component_ID] == this->GlobalCountOperation)
				{
					int temp = 0;

					// Calculate area of each component.
					this->Recalculate_Component_Area(_pMesh, Component_ID, temp);
					
					// Estimated quantization precision of geometry
					int QG = Estimate_Geometry_Quantization(_pMesh, this->ComponentVolume[Component_ID], this->ComponentArea[Component_ID], this->ComponentNumberVertices[Component_ID]);

					if (QG < 4)
						QG = 4;		
					
					// if the current precision is > QG then we apply decrease of quantization precision
					if (QG < (int)this->Qbit[Component_ID])	// MT
					{
						Continue = true;
						Operation_choice = 1;
						
						// Reducing of quantization precision of 1 bit.
						this->Diminush_Geometry_Quantization_Precision(_pMesh, Component_ID);

						this->Qbit[Component_ID]--;
						this->NumberChangeQuantization[Component_ID]++;			

						this->ComponentOperations[Component_ID]++;
						this->ListOperation[Component_ID].push_front(Operation_choice);
					}		
					// else -> Decimation
					else
					{					
						Operation_choice = 0;

						unsigned Initial_number_vertices = _pMesh.size_of_vertices();
						
						// Decimation and regulation conquests.
						this->Decimation_Conquest(_pMesh, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value, Component_ID);
						this->Regulation(_pMesh, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value, Component_ID);
						
						this->NumberDecimation[Component_ID] += 1;						
						this->ComponentOperations[Component_ID] += 1;
						this->ListOperation[Component_ID].push_front(Operation_choice);	
						
						unsigned Diff_number_vertices = _pMesh.size_of_vertices() - Initial_number_vertices;
						this->ComponentNumberVertices[Component_ID] += Diff_number_vertices;						
						
						if (Diff_number_vertices == 0)
							this->Remove_Last_Phase_Elements(Component_ID);

						else
							Continue = true;	
					}				
				}				
			}

			Current_Number = _pMesh.size_of_vertices();	
			
			if (Continue)
				this->GlobalCountOperation++;			

			
			if (Current_Number < (unsigned)_NVertices)	// MT
				break;
			
		}while((Current_Number != Last_Number) || (Continue));
		
		_pMesh.compute_normals();
	}

	// if mesh is colored
	else
	{
		FILE * Operation_order = fopen("Operation_order.txt","w");
		fclose(Operation_order);

	    //bool Is_any_vertex_removed = true;
		unsigned Last_Number = 0;
		unsigned Current_Number = _pMesh.size_of_vertices();			
		int Operation_choice = -1;
		
		vector<int> QC_Initials;
		vector<int> QC_Finals;
		
		//To get QC_Initial and QC_Final
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{
			double Max_color, Mean_color;
			int Number_vertices = 0;
			
			// Calculate max and mean color error.
			this->Calculate_Edge_Color_Difference(_pMesh, Component_ID, Max_color, Mean_color, Number_vertices);
			double Mean_Max = Mean_color / Max_color;
	
			//int QC_Initial = floor(-56.334*Mean_Max*Mean_Max +4.6838*Mean_Max +7.8632 + 0.5);
			//int QC_Final = floor(-0.893 * log(Mean_Max*Area) + 8.0957 + 0.5);

			int QC_init = floor(-57.966*Mean_Max*Mean_Max + 5.311*Mean_Max +7.8062 + 0.5);			
			int QC_fin = floor(84.548*Mean_Max*Mean_Max - 33.723*Mean_Max + 7.8222 + 0.5);
			
			QC_Initials.push_back(QC_init);
			QC_Finals.push_back(QC_fin);
		}


		bool Continue;
		do
		{		  
			Last_Number = Current_Number;		
			Continue = false;

			for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
			{
				// if the i-th component did not remove any vertex in last loop, it is not necessary to simplifiy it.
				if (this->ComponentOperations[Component_ID] == this->GlobalCountOperation)
				{
					int temp = 0;
					this->Recalculate_Component_Area(_pMesh, Component_ID, temp);
					int Number_of_vertices = 0;
					double Max_color = 0., Mean_color = 0.;
					this->Calculate_Edge_Color_Difference(_pMesh, Component_ID, Max_color, Mean_color, Number_of_vertices);	
					double Mean_Max = Mean_color / Max_color;

					// Estimated color quantization(QC) and geometry quantization (QG)
					int QC = floor(-1.181 * log(Mean_Max*this->ComponentArea[Component_ID] / (double)Number_of_vertices) + 0.3281 + 0.5);
					int QG = Estimate_Geometry_Quantization(_pMesh, this->ComponentVolume[Component_ID], this->ComponentArea[Component_ID], this->ComponentNumberVertices[Component_ID]);
					
					// If the current color quantization precision > QC_INIT -> Decrease one bit from color quantization resolution.					
					if ((8 - this->NumberColorQuantization[Component_ID] > QC_Initials[Component_ID]) && (!this->IsOneColor))// && (QC_Final < 8 - this->NumberColorQuantization))
					{						
						Operation_choice = 2;
						Continue = true;
						
						// Reducing color quantization precision of 1 bit.
						this->Diminush_Color_Quantization_Precision(_pMesh, Component_ID);

						this->NumberColorQuantization[Component_ID]++;
						this->ComponentOperations[Component_ID]++;
						this->ListOperation[Component_ID].push_front(Operation_choice);
					}
					// If the current color quantization precision > Estimated precision -> Decrease
					else if (((8 - this->NumberColorQuantization[Component_ID] > QC ) && (QC >= QC_Finals[Component_ID])) && (!this->IsOneColor))
					{		
						Operation_choice = 2;
						Continue = true;
						
						// Reducing color quantization precision of 1 bit.
						this->Diminush_Color_Quantization_Precision(_pMesh, Component_ID);

						this->NumberColorQuantization[Component_ID]++;
						this->ComponentOperations[Component_ID]++;
						this->ListOperation[Component_ID].push_front(Operation_choice);
					}
					
					// Geometry quantization precision decrease
					else if (((int)this->Qbit[Component_ID] > QG) && (QG >= 5))
					{
						Continue = true;
						Operation_choice = 1;
						
						// Reducuing geometry quantization precision of 1 bit.
						this->Diminush_Geometry_Quantization_Precision(_pMesh, Component_ID);

						this->Qbit[Component_ID]--;
						this->NumberChangeQuantization[Component_ID]++;			

						this->ComponentOperations[Component_ID]++;
						this->ListOperation[Component_ID].push_front(Operation_choice);		
					}
					
					// Else we perform decimation
					else
					{					
						Operation_choice = 0;

						unsigned Initial_number_vertices = _pMesh.size_of_vertices();
						
						// Decimation and regulation conquests.
						this->Decimation_Conquest(_pMesh, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value, Component_ID);
						this->Regulation(_pMesh, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value, Component_ID);

						this->NumberDecimation[Component_ID] += 1;
						unsigned Diff_number_vertices = _pMesh.size_of_vertices() - Initial_number_vertices;
						
						this->ComponentOperations[Component_ID] += 1;
						this->ListOperation[Component_ID].push_front(Operation_choice);	
						
						this->ComponentNumberVertices[Component_ID] += Diff_number_vertices;

						if (Diff_number_vertices == 0)
							this->Remove_Last_Phase_Elements(Component_ID);

						else
							Continue = true;	
					}
					FILE * Operation_order = fopen("Operation_order.txt", "a");
					fprintf(Operation_order, "Component_ID = %d    Operation = %d \n", Component_ID, Operation_choice);
					fclose(Operation_order);
				}
			}

			Current_Number = _pMesh.size_of_vertices();	
			
			if (Continue)
				this->GlobalCountOperation++;	
			
			
			if (Current_Number < (unsigned)_NVertices)	// MT
				break;
			
		}while((Current_Number != Last_Number) || (Continue));
		
		_pMesh.compute_normals();

		fclose(Operation_order);
	}


}

// Description : This function select a set of independent vertices to be removed
int Compression_Valence_Component::Decimation_Conquest(Polyhedron  & _pMesh,
													   const bool    Normal_flipping,
													   const bool    Use_metric,
													   const float & Metric_thread,
													   const bool    Use_forget_metric,
													   const int   & Forget_value,
													   const int   & Component_ID)
{	
	
	// Calculate mean color and meah area for color metric
	double Max_color, Mean_color;
	int Temp_NV = 0;
	int Number_facets;
	this->Calculate_Edge_Color_Difference(_pMesh, Component_ID, Max_color, Mean_color, Temp_NV);
	this->Recalculate_Component_Area(_pMesh, Component_ID, Number_facets);
	double Mean_area = (double)this->ComponentArea[Component_ID] / (double)Number_facets;


	// Initialize vertex and face flags.
	Init(_pMesh);
		
	// to count number of independent vertices and number of symbol of connectivity.
	int Number_vertices = 0; 
	int Number_symbol = 0;	

	// To find first edge.
	Halfedge_iterator hi = _pMesh.halfedges_begin();			
	
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID+1))
		hi++;

	Halfedge_handle First_halfedge = &(*(hi));
	
	// Two vertices of seed edge are flaged CONQUERED
	First_halfedge->vertex()->Vertex_Flag = CONQUERED;
	First_halfedge->opposite()->vertex()->Vertex_Flag = CONQUERED;

	// These vertices are also flaged with sign flags for retriangulation
	First_halfedge->vertex()->Vertex_Sign = PLUS;
	First_halfedge->opposite()->vertex()->Vertex_Sign = MINUS;	
	
	std::queue<Halfedge_handle> Halfedges; // halfedge queue. 	
	Halfedges.push(First_halfedge);	// push the first halfedge in the queue.

	Halfedge_handle h; // The current gate
	
	/// Main loop
	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();		
		
		unsigned type = 0; // define type of retriangulation		
		
		int valence = (int)h->next()->vertex_degree(); // valence		

		if (h->is_border() == true)
			continue;		
		
		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED, do nothing!!
		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;

		// if its front vertex is free and has a valence <= 6 and it is not a border vertex.
		else if ( (h->next()->vertex()->Vertex_Flag == FREE) && 
				  (valence >= 3) && 
				  (valence <= 6) && 
				  (!Is_Border_Vertex(h->next())))
		{

			type = Find_Type(h, valence);
			
			// Check if the manifold property is violated.
			bool Manifold_condition = Is_Manifold_Property_Violated(h, type, valence);
			
            bool Normal_flipping_condition = false;
            
			if (Normal_flipping)
                Normal_flipping_condition = Is_Normal_Flipping_Occured(h, valence);
			
			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Geometric_metric_condition = false;
			
			if (Use_metric)
			{
				if (Use_forget_metric)
				{
					if (_pMesh.size_of_vertices() > (unsigned)Forget_value)
						Geometric_metric_condition = false;
					else
						Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
				}
				else
					Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
			}					
			bool Is_Color_Too_Important = false;			
			
			#ifdef USE_COLOR_METRIC
			if(this->IsColored)
				Is_Color_Too_Important = this->Error_Projected_Surface(_pMesh, h, Component_ID, Mean_color, Mean_area);
			#endif

			// remove the front vertex if its removal does not viloate the manifold property and some metrics			
			bool Check_all_condition = false;

			if ((!Manifold_condition) && (!Geometric_metric_condition) && (!Normal_flipping_condition) && (!Is_Color_Too_Important))
				Check_all_condition = true;			

			if (Check_all_condition) // All conditions are good. -> Remove the center vertex.
			{

				//increase number of vertices and symbols
				Number_vertices++;
				Number_symbol++;				
				
				Halfedge_handle g = h;				

				Point3d Coordinates_removed_vertex = g->next()->vertex()->point();
				Point_Int CRV = Change_Real_Int(Coordinates_removed_vertex, Component_ID);				

				// encoding of vertex color
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Removed_vertex_color;
					Removed_vertex_color = Get_Vertex_Color(g->next());					

					Color_Unit Average_color;

					Average_color = Get_Average_Vertex_Color_Lee(g, valence);			


					// Color difference from average color of neighbors					
					Color_Unit Color_diff = Removed_vertex_color - Average_color;					
					this->InterVertexColor.push_front(Color_diff);	

				}			
								

				// Enter symbol 'VALENCE CODE' into the list of symbols
				this->InterConnectivity.push_front(valence - 3);
				
				// Calculate the position of barycenter.		
				Point3d Barycenter = Barycenter_Patch_Before_Removal(g);

				Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

				// remove the front vertex
				_pMesh.erase_center_vertex(g->next()); 

				g = h;
				g->facet()->Facet_Flag = TEMP_FLAG;				
				
				// Flags and fill queue
				for (int j = 0; j < (valence-1); j++)
				{
					g = g->next();
					g->vertex()->Vertex_Flag = CONQUERED;
					
					g->opposite()->vertex()->Vertex_Flag = CONQUERED;
					
					if (!g->is_border_edge())
						Halfedges.push(g->opposite());
				}
				g = h;							
												
				Retriangulation(_pMesh, g, valence, 0, Component_ID);
				
				Vector normal = Normal_Patch(g, valence);
				Vector T2 = CGAL::NULL_VECTOR;
				Vector T1 = Calculate_T1_T2(h, normal, T2);

				if (normal == CGAL::NULL_VECTOR)
				{					
					T1 =     Vector(1.,0.,0.);
					T2 =     Vector(0.,1.,0.);
					normal = Vector(0.,0.,1.);
				}
				else if (T1 == CGAL::NULL_VECTOR)
				{					
					T1 =     Vector(1.,0.,0.);
					T2 =     Vector(0.,1.,0.);
					normal = Vector(0.,0.,1.);
				}	
				else if (T2 == CGAL::NULL_VECTOR)
				{					
					T1 =     Vector(1.,0.,0.);
					T2 =     Vector(0.,1.,0.);
					normal = Vector(0.,0.,1.);
				}	
				
				Point_Int Dist = CRV - BC;
				Point_Int Frenet_Coordinates;
				//Bijection
				if (this->Is_Bijection_Enabled)
					Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2,normal);
				else
					Frenet_Coordinates = Dist;

				this->InterGeometry.push_front(Frenet_Coordinates);				
			}
			
			// Conditions are not satisfied -> NULL CODE
			else 
			{
				// Enter symbol 'NULL PATCH' into the list of symbols
				this->InterConnectivity.push_front(4);

				Number_symbol++;
				h->facet()->Facet_Flag = CONQUERED;
				h->next()->vertex()->Vertex_Flag = CONQUERED;
				
				if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());				
				
				if (!h->prev()->is_border_edge())					
					Halfedges.push(h->prev()->opposite());				
			}
		}
		// border edge.
		else if ( (h->next()->vertex()->Vertex_Flag == FREE) && 
				  (valence >= 3) && 
				  (valence <= 4) && 
				  (Is_Border_Vertex(h->next()) == true))
		{
			/*****    conditions of vertex removal (based on area) will be added    *****/
			/*****    conditions of vertex removal (based on area) will be added    *****/		

			type = Find_Type(h, valence);
			
			Halfedge_handle g = h;
			
			bool Check_border_structure = true;
			
			Halfedge_handle First_border_edge;
			Halfedge_handle Standard_edge;			
			
					
			
			// To find first border edge.
			Halfedge_around_vertex_circulator Hvc = g->next()->vertex()->vertex_begin();
			Halfedge_around_vertex_circulator Hvc_end = Hvc;
			
			int Number_neighboring_border_vertices = 0;
			CGAL_For_all(Hvc, Hvc_end)
			{				
				if (Is_Border_Vertex(Hvc->opposite()))
					Number_neighboring_border_vertices++;					
				
				if (Hvc->is_border())
					First_border_edge = Hvc;
			}					
			
			if (Number_neighboring_border_vertices > 2)
				Check_border_structure = false;	
				
			if (Is_Manifold_Property_Violated(g, type, valence))
				Check_border_structure = false;
			// Added by Adrien Maglo.
			// Test if we are not in the case of a mesh pimple.
			if (Check_border_structure)
			{
				// Get the two patch border vertices.
				Vertex_handle vh1 = First_border_edge->next()->vertex();
				Vertex_handle vh2 = First_border_edge->opposite()->vertex();
				Halfedge_around_vertex_circulator vh_it = vh1->vertex_begin();
				Halfedge_around_vertex_circulator vh_it_end = vh_it;

				// Test if the two patch border vertices are not connected
				// by an edge.
				CGAL_For_all(vh_it, vh_it_end)
				{
					if (vh_it->opposite()->vertex() == vh2)
					{
						Check_border_structure = false;
						break;
					}
				}
			}			


			// if the hole is a triangle(Count_border_edges == 3), we should not remove the vertex.
			int Count_border_edges = 1;
			Standard_edge = First_border_edge;
			First_border_edge = First_border_edge->next();
			
			for (; First_border_edge != Standard_edge; First_border_edge = First_border_edge->next())
			{
				Count_border_edges++;				
				if (Count_border_edges > 5)
					break;
			}
			
			if (Count_border_edges <= 3)
				Check_border_structure = false;

			// Exception which causes a manifold violation
			if (Count_border_edges == 4)
				if (Standard_edge->prev()->opposite()->vertex_degree() == 2)
					Check_border_structure = false;

			/*if (Check_border_structure)
			{
				if(valence == 3)
				{
					if(Is_Border_Vertex(g))
					{
						if(g->opposite()->vertex_degree() == 3)
							Check_border_structure = false;   
					}
					else
					{
						if(g->vertex_degree() == 3)
							Check_border_structure = false;   
					}
				}
			}*/

		

			// Structure of triangles is good -> decimation.
			if (Check_border_structure) 
			{	
				Number_vertices++;
				Number_symbol++;
				
				// Connectivity Code = 5 for border valence 3 vertex
				//         "         = 6 for border valence 4 vertex
				this->InterConnectivity.push_front(valence + 2); 
			
				Point3d Real_vertex_position = h->next()->vertex()->point();
				Point_Int Vertex_position = Change_Real_Int(Real_vertex_position, Component_ID);				
				
				Color_Unit Removed_vertex_color;
				//int Vertex_color_index = -1;
				if ((this->IsColored) && (!this->IsOneColor))
				{	
					Removed_vertex_color = Get_Vertex_Color(h->next());
				}				
				

				vector<Halfedge_handle> Border_edges;
				int Number_jump = 0;

				while(g->next()->is_border_edge() == false)
				{
					g = g->next()->opposite()->next();
					Number_jump++;
				}
				
				Border_edges.push_back(g->opposite());
				for (int i = 0; i < (valence - 2); i++)
				{
					g = g->prev()->opposite()->prev();	
					Border_edges.push_back(g->opposite());													
				}		

				for (int i = 0; i < (valence - 1); i++)
				{
					Halfedge_handle Temp = Border_edges[i];
					Temp = Temp->opposite();
					
					CGAL_assertion(!Temp->is_border());
					_pMesh.erase_facet(Temp);
				}
								
				if (valence == 3)
				{
					Halfedge_handle Retriangulation_edge = Border_edges[valence - 2]->opposite()->prev();
					
					// One triangle has to be created to smooth the mesh boundary					
					_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());					
										
					Halfedge_handle Input_gate = Border_edges[Number_jump]->opposite();					

					// its front face is tagged CONQUERED
					Input_gate->facet()->Facet_Flag = CONQUERED;					
					Input_gate->next()->vertex()->Vertex_Flag = CONQUERED;

					if ((type == 1) || (type == 2) || (type == 4))
					{
						if (Input_gate->next()->vertex()->Vertex_Sign == NOSIGN)
							Input_gate->next()->vertex()->Vertex_Sign = PLUS;
					}
					else if (type == 3)
					{
						if (Input_gate->next()->vertex()->Vertex_Sign == NOSIGN)
							Input_gate->next()->vertex()->Vertex_Sign = MINUS;
					}
					

					Point3d Barycenter = Barycenter_Patch_After_Removal(Input_gate, valence);
					Point_Int BC = Change_Real_Int(Barycenter, Component_ID);					
					
					if ((this->IsColored) && (!this->IsOneColor))
					{
						Color_Unit Average_color;

						Average_color = Get_Average_Vertex_Color_After_Removal(g, valence); // g is the most left placed edge
						Color_Unit Color_diff = Removed_vertex_color - Average_color;						

						this->InterVertexColor.push_front(Color_diff);
					}

					Halfedge_handle g = Input_gate;
					Vector normal = Normal_Patch(g, valence);
					Vector T2 = CGAL::NULL_VECTOR;					
					
					Vector T1 = Calculate_T1_T2(g, normal, T2);

					if (normal == CGAL::NULL_VECTOR)
					{
						normal = Vector(0,0,1);	
						T1 = Vector(1,0,0);
						T2 = Vector(0,1,0);					
					}					
					
					Point_Int Dist = Vertex_position - BC;

					//bijection

					Point_Int Frenet_Coordinates;
					if (this->Is_Bijection_Enabled)
						Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2,normal);
					else
						Frenet_Coordinates = Dist;


					//Point_Int Frenet_Coordinates = Frenet_Rotation(Dist,T1,T2,normal);
					
					this->InterGeometry.push_front(Frenet_Coordinates);

					if (Number_jump == 0)
						Halfedges.push(Border_edges[1]);
					else
						Halfedges.push(Border_edges[0]);
				}		
				
				else if (valence == 4)
				{		
					Halfedge_handle Retriangulation_edge = Border_edges[valence - 2]->opposite()->prev();					
					
					if (   ( (Number_jump == 0) && ((type == 5) || (type == 8)) ) ||
						   ( (Number_jump == 1) && ((type == 6) || (type == 7)) ) ||
						   ( (Number_jump == 2) && ((type == 5) || (type == 8)) )  )

					{						
						_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[1]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag = CONQUERED;						

						_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());
						Border_edges[0]->opposite()->facet()->Facet_Flag = CONQUERED;						
					}

					else if (   ( (Number_jump == 0) && ((type == 6) || (type == 7)) ) ||
								( (Number_jump == 1) && ((type == 5) || (type == 8)) ) ||
								( (Number_jump == 2) && ((type == 6) || (type == 7)) )  )
					{												
						_pMesh.add_facet_to_border(Border_edges[2]->opposite(), Border_edges[0]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag = CONQUERED;						
						Halfedge_handle Temp_border = Border_edges[2]->opposite()->next();
						
						_pMesh.add_facet_to_border(Retriangulation_edge, Temp_border);
						Border_edges[2]->opposite()->facet()->Facet_Flag = CONQUERED;						
					}					

					Halfedge_handle Input_gate = Border_edges[Number_jump]->opposite();

					
					// Vertex Signs
					if ((type == 5) || (type == 8))
					{
						Halfedge_handle g = Input_gate;
						g = g->prev()->opposite()->next();
						if (g->vertex()->Vertex_Sign == NOSIGN)
							g->vertex()->Vertex_Sign = MINUS;
						if (g->opposite()->vertex()->Vertex_Sign == NOSIGN)
							g->opposite()->vertex()->Vertex_Sign = PLUS;
					}
					else if ((type == 6) || (type == 7))
					{
						Halfedge_handle g = Input_gate;
						g = g->next()->opposite()->prev();
						if (g->vertex()->Vertex_Sign == NOSIGN)
							g->vertex()->Vertex_Sign = PLUS;
						if (g->opposite()->vertex()->Vertex_Sign == NOSIGN)
							g->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					
					//Vertex Flags + Fill Halfedges queue.
					for (int i = 0; i < valence - 1; i++)
					{
						Border_edges[i]->vertex()->Vertex_Flag = CONQUERED;
						Border_edges[i]->opposite()->vertex()->Vertex_Flag = CONQUERED;

						if (i != (unsigned)Number_jump)
							Halfedges.push(Border_edges[i]);
					}
					
					if ((this->IsColored) && (!this->IsOneColor))
					{
						Color_Unit Average_color;

						Color_Unit Color_0 = Get_Vertex_Color(Border_edges[0]->opposite());
						Color_Unit Color_1 = Get_Vertex_Color(Border_edges[0]);
						Color_Unit Color_2 = Get_Vertex_Color(Border_edges[1]);
						Color_Unit Color_3 = Get_Vertex_Color(Border_edges[2]);
						Average_color.c0 = (int)(Color_0.c0 + Color_1.c0 + Color_2.c0 + Color_3.c0)/4.0;
						Average_color.c1 = (int)(Color_0.c1 + Color_1.c1 + Color_2.c1 + Color_3.c1)/4.0;
						Average_color.c2 = (int)(Color_0.c2 + Color_1.c2 + Color_2.c2 + Color_3.c2)/4.0;
						
						Color_Unit Color_diff = Removed_vertex_color - Average_color;						
						
						this->InterVertexColor.push_front(Color_diff);

					}				
					

					//encoding of geometry
					Point3d P0,P1,P2,P3;
					
					P0 = Border_edges[0]->opposite()->vertex()->point();				
					P1 = Border_edges[0]->vertex()->point();
					P2 = Border_edges[1]->vertex()->point();				
					P3 = Border_edges[2]->vertex()->point();					

					Point3d Barycenter((P0.x()+P1.x()+P2.x()+P3.x())/4,
									   (P0.y()+P1.y()+P2.y()+P3.y())/4,
									   (P0.z()+P1.z()+P2.z()+P3.z())/4);				
					
					Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

					Halfedge_handle g = Input_gate;
					Vector normal = Normal_Patch(g, valence);

					Vector T2 = CGAL::NULL_VECTOR;
					Vector T1 = Calculate_T1_T2(g, normal, T2);

					if (normal == CGAL::NULL_VECTOR)
					{
						T1 = Vector(1,0,0);
						T2 = Vector(0,1,0);	
						normal = Vector(0,0,1);										
					}					
					
					Point_Int Dist = Vertex_position - BC;
					Point_Int Frenet_Coordinates;
					if (this->Is_Bijection_Enabled)
						Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2,normal);
					else
						Frenet_Coordinates = Dist;
					
					this->InterGeometry.push_front(Frenet_Coordinates);					
				}					
			}

			else // Border vertex can not be removed
			{
				Number_symbol++;
				this->InterConnectivity.push_front(4);
				
				h->facet()->Facet_Flag = CONQUERED;				
				h->next()->vertex()->Vertex_Flag = CONQUERED;

				if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}				
				
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());

				if (!h->prev()->is_border_edge())
					Halfedges.push(h->prev()->opposite());			
			}	
		}
		

		

		// if its front face sholud be labelled as NULL PATCH		
		else
		{
			// Enter symbol 'NULL PATCH' into the list of symbols
			this->InterConnectivity.push_front(4);
			Number_symbol++;

			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;///////////////////////////

			if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}

			// two other output gates are pused into the fifo queue.
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)					
				Halfedges.push(h->prev()->opposite());
		}
	}	
	if ((this->IsColored) && (!this->IsOneColor))
	{		
		while(!this->InterVertexColor.empty())
		{
			Color_Unit Col = this->InterVertexColor.front();
			this->InterVertexColor.pop_front();
			this->VertexColor[Component_ID].push_front(Col);
		}
	}
	// InterConnectivity is the intermediate connectivity symbol container.
	// We put all symbols in the main container which is this->Connectivity
	while(!InterConnectivity.empty())
	{
		int Symbol = InterConnectivity.front();
		InterConnectivity.pop_front();
		this->Connectivity[Component_ID].push_front(Symbol);
	}
	
	// same operation than connectivity.
	while(!InterGeometry.empty())
	{
		Point_Int Geo = InterGeometry.front();
		InterGeometry.pop_front();
		this->Geometry[Component_ID].push_front(Geo);
	}
	
	this->NumberVertices[Component_ID].push_front(Number_vertices);
	this->NumberSymbol[Component_ID].push_front(Number_symbol);
	
	// if this decimation didn't remove any vertex, we should remove all connectivity symbols. for that we store number of symbols.
	this->DumpSymbolDecimation = Number_symbol;	


	return Number_vertices;
}


// Description : Regulation conquest
int Compression_Valence_Component::Regulation(Polyhedron  & _pMesh,
								              const bool    Normal_flipping,
								              const bool    Use_metric,
											  const float & Metric_thread,
											  const bool    Use_forget_metric,
											  const int   & Forget_value,
											  const int   & Component_ID)
{    
	double Max_color, Mean_color;
	int Temp_NV = 0;
	int Number_facets;
	this->Calculate_Edge_Color_Difference(_pMesh, Component_ID, Max_color, Mean_color, Temp_NV);
	this->Recalculate_Component_Area(_pMesh, Component_ID,Number_facets);
	double Mean_area = (double)this->ComponentArea[Component_ID] / (double)Number_facets;

	// Initialization
	Init(_pMesh);
	
	
	
	// Number of removed vertices and number of connectivity symbols
	int Number_vertices = 0;
	int Number_symbol = 0;

	Halfedge_iterator hi = _pMesh.halfedges_begin();
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))
		hi++;	

	Halfedge_handle First_halfedge = &(*(hi));	

	First_halfedge->vertex()->Vertex_Flag = CONQUERED;
	First_halfedge->opposite()->vertex()->Vertex_Flag = CONQUERED;
		
	std::queue<Halfedge_handle> Halfedges;
	Halfedges.push(First_halfedge);

	Halfedge_handle h;

	while(!Halfedges.empty())
	{		
		h = Halfedges.front();
		Halfedges.pop();
		
		size_t valence = h->next()->vertex_degree();
		//int Component_ID = h->next()->vertex()->Component_Number;

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;

		else if ((h->next()->vertex()->Vertex_Flag == FREE) && (valence == 3) && (Is_Border_Vertex(h->next()) == false)) // if valence is 3, remove the front vertex.
		{
			Halfedge_handle g = h;
			int type = 1; // ant type of valence 3

			// Check if the manifold property is violated.
			//bool Manifold_condition = Is_Manifold_Property_Violated(h, type, valence);			
           
			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Geometric_metric_condition = false;
			if (Use_metric == true)
			{
				if (Use_forget_metric == true)
				{
					if (_pMesh.size_of_vertices() > (unsigned)Forget_value)
						Geometric_metric_condition = false;
					else
						Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
				}

				else
					Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
			}	
			bool Is_Color_Too_Important = false;	

			#ifdef USE_COLOR_METRIC
			if(this->IsColored)
			{
				Is_Color_Too_Important = this->Error_Projected_Surface(_pMesh, h, Component_ID, Mean_color, Mean_area);
			}
			#endif

			// remove the front vertex if its removal does not viloate the manifold property and some metrics			
			bool Check_all_condition = false;

			if ((!Geometric_metric_condition) && (!Is_Color_Too_Important))
				Check_all_condition = true;			
			
			// All conditions are goods. Remove the front vertex
			if (Check_all_condition)
			{
				Number_vertices++;
				Number_symbol++;
				
				// Insertion of a symbol 3 to the list
				this->InterConnectivity.push_front(0);
				
				Point3d Geo_info = h->next()->vertex()->point();
				Point_Int Geo = Change_Real_Int(Geo_info, Component_ID);				
				
				Point3d Barycenter = Barycenter_Patch_Before_Removal(g);								
				Point_Int BC = Change_Real_Int(Barycenter, Component_ID);		
				
				
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Removed_vertex_color;
					
					Removed_vertex_color = Get_Vertex_Color(h->next());

					Color_Unit Average_color;
					Average_color = Get_Average_Vertex_Color_Lee(g, valence);
					
					Color_Unit Color_diff = Removed_vertex_color - Average_color;					

					this->InterVertexColor.push_front(Color_diff);
				}

				g = h;

				// P, Q, R are used to compute the normal of triangle after removal
				Point3d P = g->vertex()->point();
				g = g->next()->opposite()->next();
				Point3d Q = g->vertex()->point();
				g = g->next()->opposite()->next();
				Point3d R = g->vertex()->point();				
				g = h;

				Vector normal = Triangle_Normal(P,Q,R);
				Vector T2 = CGAL::NULL_VECTOR;
				Vector T1 = Calculate_T1_T2(h,normal,T2);
				
				if (T1 == CGAL::NULL_VECTOR)
				{				
					T1 = Vector(1,0,0);
					T2 = Vector(0,1,0);					
					normal = Vector(0,0,1);	
				}	

				else if (normal == CGAL::NULL_VECTOR)
				{				
					T1 = Vector(1,0,0);
					T2 = Vector(0,1,0);					
					normal = Vector(0,0,1);	
				}
				else if (T2 == CGAL::NULL_VECTOR)
				{				
					T1 = Vector(1,0,0);
					T2 = Vector(0,1,0);					
					normal = Vector(0,0,1);	
				}
				
				Point_Int Dist = Geo - BC;					
				Point_Int Frenet_Coordinates;
				
				if(this->Is_Bijection_Enabled)
					Frenet_Coordinates = Frenet_Rotation(Dist,T1,T2,normal);
				else
					Frenet_Coordinates = Dist;
				
				this->InterGeometry.push_front(Frenet_Coordinates);
								
				g = h->next();
				g->vertex()->Vertex_Flag = TO_BE_REMOVED;
				g->facet()->Facet_Flag = TO_BE_REMOVED;

				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				
				if (!g->prev()->is_border_edge())
				{
					Halfedge_handle h1 = g->prev()->opposite();
					h1->facet()->Facet_Flag = CONQUERED;
					h1->next()->vertex()->Vertex_Flag = CONQUERED;
					if (!h1->next()->is_border_edge())
						Halfedges.push(h1->next()->opposite());
					if (!h1->prev()->is_border_edge())
						Halfedges.push(h1->prev()->opposite());					
				}
				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				if (!g->prev()->is_border_edge())
				{
					Halfedge_handle h2 = g->prev()->opposite();
					h2->facet()->Facet_Flag = CONQUERED;
					h2->next()->vertex()->Vertex_Flag = CONQUERED;
					if (!h2->next()->is_border_edge())
						Halfedges.push(h2->next()->opposite());
					if (!h2->prev()->is_border_edge())
						Halfedges.push(h2->prev()->opposite());
				}
			}

			else
			{
				this->InterConnectivity.push_front(1);
				Number_symbol++;
				h->facet()->Facet_Flag = CONQUERED;
				h->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());
				if (!h->prev()->is_border_edge())					
					Halfedges.push(h->prev()->opposite());				
			}
		}
		else  // NULL triangle
		{
			this->InterConnectivity.push_front(1);
			Number_symbol++;
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())					
				Halfedges.push(h->prev()->opposite());			
		}
	}
	
	
	// Removal of all vertices with TO_BE_REMOVED flag
	Vertex_iterator pVertex = NULL;
	for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end();)
	{
		Vertex_handle vh = pVertex;
		pVertex++;

		if (vh->Vertex_Flag == TO_BE_REMOVED)
		{
			Halfedge_handle temp = vh->halfedge();
			temp = _pMesh.erase_center_vertex(temp);
			temp->facet()->Component_Number = Component_ID;
		}
	}	

	if ((this->IsColored) && (!this->IsOneColor))
	{		
		while(!this->InterVertexColor.empty())
		{
			Color_Unit Col = this->InterVertexColor.front();
			this->InterVertexColor.pop_front();
			this->VertexColor[Component_ID].push_front(Col);
		}	
	}
	while(!InterConnectivity.empty())
	{
		int Symbol = this->InterConnectivity.front();
		this->InterConnectivity.pop_front();
		this->Connectivity[Component_ID].push_front(Symbol);
	}
	
	while(!InterGeometry.empty())
	{
		Point_Int Geo = InterGeometry.front();
		InterGeometry.pop_front();
		this->Geometry[Component_ID].push_front(Geo);
	}
	
	this->NumberSymbol[Component_ID].push_front(Number_symbol);
	this->NumberVertices[Component_ID].push_front(Number_vertices);

	this->DumpSymbolRegulation = Number_symbol;
	
	

	return Number_vertices;
}



// Description : Decoding of the regulation conquest
void Compression_Valence_Component::Un_Regulation(Polyhedron &_pMesh, Arithmetic_Codec & Decoder, const int & Component_ID)
{	
	Init(_pMesh);

	Adaptive_Data_Model Connectivity(2);
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;

	unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];
	
	int Alpha_range = Decoder.get_bits(Qbit+1);
	int Alpha_offset = Decoder.get_bits(Qbit+1);
	int Gamma_range = Decoder.get_bits(Qbit+1);
	int Gamma_offset = Decoder.get_bits(Qbit+1);

	if (this->Smallest_Alpha < 0)
		Alpha_offset = Alpha_offset + this->Smallest_Alpha;

	if (this->Smallest_Gamma < 0)
		Gamma_offset = Gamma_offset + this->Smallest_Gamma;

	bool check_alpha = false;
	bool check_gamma = false;

	if ((Alpha_range == 0) || (Alpha_range == 1))
	{
		check_alpha = true;
		Alpha_range = 2;
	}

	if ((Gamma_range == 0) || (Gamma_range == 1))
	{
		check_gamma = true;
		Gamma_range = 2;
	}

	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);

	Adaptive_Data_Model alpha(Alpha_range);
	Adaptive_Data_Model gamma(Gamma_range);
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))	
		hi++;	
	
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))// already visited.
			continue;
		
		// read connectivity information
		int valence = Decoder.decode(Connectivity) + 3;
		
		// if valence is 3, remove the front vertex.
		if (valence == 3)
		{
			Halfedge_handle g = h;			

			// Insertion of a vertex
			Halfedge_handle pass = h;
			
			vector<Point3d> Vertices; //contains 1-ring and 2-rings vertices
			Point3d Barycenter = Barycenter_Patch_After_Removal(pass,3);
			Point_Int BC = Change_Real_Int(Barycenter, Component_ID);
			
			Vector normal = Triangle_Normal(pass);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(pass,normal,T2);

			if (T1 == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}	

			else if (normal == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}
			else if (T2 == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}
			
			Point_Int Frenet;
			
			if (check_alpha == false)
			{
				Frenet.x = Decoder.decode(alpha);
				Frenet.y = Decoder.decode(alpha);
			}
			else 
			{
				Frenet.x = 0;
				Frenet.y = 0;
			}

			if (check_gamma == false)
				Frenet.z = Decoder.decode(gamma);
			
			else 
				Frenet.z = 0;			

			Frenet.x -= Alpha_offset;
			Frenet.y -= Alpha_offset;
			Frenet.z -= Gamma_offset;

			Point_Int Diff;			
			if(this->Is_Bijection_Enabled)
				Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);
			else
				Diff = Frenet;
			Point_Int Center = BC + Diff;		
			
			Point3d Center_vertex = this->Change_Int_Real(Center,Component_ID);		
			
			
			// Assign the region number to inserted vertex
			Halfedge_handle reg = h;

			int Selected_region = 500000;
			//int Number_vertices = 500000;
			vector<int> T_Bin;
			vector<int> T_Number;

			for(int i = 0; i < valence; i++)
			{
				int N1 = reg->vertex()->Region_Number;
				bool Is_existed = false;
				for(unsigned int j = 0; j < T_Bin.size(); j++)
				{
					if(N1 == T_Bin[j])
					{
						T_Number[j]++;
						Is_existed = true;
					}
				}
				if(!Is_existed)
				{
					T_Bin.push_back(N1);
					T_Number.push_back(1);
				}
				reg = reg->next();
			}
			int Max = -5000;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] > Max)
					Max = T_Number[i];
			}
			vector<int> T_possible_bin;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] == Max)
					T_possible_bin.push_back(T_Bin[i]);
			}

			if(T_possible_bin.size() == 1)
			{
				Selected_region = T_possible_bin[0];
			}
			else
			{
				Selected_region = 5000;
				for(unsigned int i = 0; i < T_possible_bin.size(); i++)
				{
					if(T_possible_bin[i] < Selected_region)
						Selected_region = T_possible_bin[i];
				}
			}

			
			// Vertex insertion
			g = _pMesh.create_center_vertex(g);
			
			g->vertex()->point() = Center_vertex;
			g->vertex()->Seed_Edge = -1;		
			
			int RO = this->GlobalCountOperation - Decompress_count;

			g->vertex()->Removal_Order = RO;	

			g->vertex()->Region_Number = Selected_region;
			
			if(Selected_region != -1)
				this->m_Number_Vertices_Per_Regions[Selected_region]++;

			// Vertex flags
			g->vertex()->Vertex_Flag = CONQUERED;
			g = h;
			g->facet()->Facet_Flag = CONQUERED;
			g = g->next();
			g = g->prev_on_vertex();
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			g->facet()->Facet_Flag = CONQUERED;			
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h1 = g->prev()->opposite();
				h1->facet()->Facet_Flag = CONQUERED;
				h1->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h1->next()->is_border_edge())
					Halfedges.push(h1->next()->opposite());
				if (!h1->prev()->is_border_edge())
					Halfedges.push(h1->prev()->opposite());
			}
			g = g->prev_on_vertex();
			g->facet()->Facet_Flag = CONQUERED;
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h2 = g->prev()->opposite();
				h2->facet()->Facet_Flag = CONQUERED;
				h2->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h2->next()->is_border_edge())
					Halfedges.push(h2->next()->opposite());
				if (!h2->prev()->is_border_edge())
					Halfedges.push(h2->prev()->opposite());
			}
			
			
			if ((this->IsColored) && (!this->IsOneColor))
			{				
				g = h;				

				Color_Unit Predicted_color;

				Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
				
				Color_Unit Color_difference;
				Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
				Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
				Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;								
				
				Color_Unit CV = Predicted_color + Color_difference;				

				g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);

				float LAB[3];
				LAB[0] = this->C0_Min + CV.c0 * Color_step;
				LAB[1] = this->C1_Min + CV.c1 * Color_step;
				LAB[2] = this->C2_Min + CV.c2 * Color_step;

				float RGB[3];
				LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

				g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);
			}
			
			if ((this->IsColored) && (this->IsOneColor))
			{
				g = h->next();
				g->vertex()->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
			}
		}

		else
		{
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())
				Halfedges.push(h->prev()->opposite());
		}
	}
}



// Description : Decoding function of decimation conquest
void Compression_Valence_Component::Un_Decimation_Conquest(Polyhedron       & _pMesh, 
														   Arithmetic_Codec & Decoder,
														   const int        & Component_ID)
{
	Init(_pMesh);
	
	int Number_connectivity_symbols;
	if (this->IsClosed[Component_ID])
		Number_connectivity_symbols = 5;
	else
		Number_connectivity_symbols = 7;

	Adaptive_Data_Model Connectivity(Number_connectivity_symbols);	
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;	

	unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];

	int Alpha_range = Decoder.get_bits(Qbit+1);
	int Alpha_offset = Decoder.get_bits(Qbit+1);
	int Gamma_range = Decoder.get_bits(Qbit+1);
	int Gamma_offset = Decoder.get_bits(Qbit+1);

	if (this->Smallest_Alpha < 0)
		Alpha_offset = Alpha_offset + this->Smallest_Alpha;
	
	if (this->Smallest_Gamma < 0)
		Gamma_offset = Gamma_offset + this->Smallest_Gamma;

	bool check_alpha = false;
	bool check_gamma = false;

	if ((Alpha_range == 0) || (Alpha_range == 1))
	{
		check_alpha = true;
		Alpha_range = 2;
	}

	if ((Gamma_range == 0) || (Gamma_range == 1))
	{
		check_gamma = true;
		Gamma_range = 2;
	}
	

	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);


	Adaptive_Data_Model alpha(Alpha_range);
	Adaptive_Data_Model gamma(Gamma_range);
	
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID+1))
		hi++;	

	// Two vertices of seed edges are flaged CONQUERED
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;	
		
	// These vertices are also flages with sign flags for retriangulation
	hi->vertex()->Vertex_Sign = PLUS;
	hi->opposite()->vertex()->Vertex_Sign = MINUS;
	
	// Two vertices of seed edges are flaged CONQUERED
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;	

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		unsigned int valence = 0,type = 0; // define type of retriangulation

		if (h->is_border() == true)
			continue;	
		
		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED
		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;
			
		valence = Decoder.decode(Connectivity) + 3;		
				
		// if its front vertex is free
		if ((valence >= 3) && (valence <= 6))
		{
			type = Find_Type(h, valence);
			
			// remove the front vertex if its removal does not viloate the manifold property
			Halfedge_handle pass = h;
			Halfedge_handle g = h;
			
			Vector normal = Normal_Patch(pass, valence);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(h,normal, T2);			
			
			if (T1 == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}	

			else if (normal == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}
			else if (T2 == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}

			// remove edges to re_find polygon and (check. and attribute sign flag.)
			bool Check_Validity = Remove_Edges(_pMesh, g, type);
			Check_Validity = false;

			//Check_Validity = false;// * * * * * * * * * * *////
			if (Check_Validity == false)
			{
				g = h;
				Halfedge_handle pass = h;
				
				vector<Point3d> Vertices; //contains 1-ring and 2-ring vertices;
				Point3d Barycenter = Barycenter_Patch_After_Removal(pass,valence);				
				Point_Int BC = Change_Real_Int(Barycenter, Component_ID);
				
				Point_Int Frenet;
				if (check_alpha == false)
				{
					Frenet.x = Decoder.decode(alpha);
					Frenet.y = Decoder.decode(alpha);
				}
				else if (check_alpha == true)
				{
					Frenet.x = 0;
					Frenet.y = 0;
				}
				if (check_gamma == false)
					Frenet.z = Decoder.decode(gamma);
				else if (check_gamma == true)
					Frenet.z = 0;

				Frenet.x -= Alpha_offset;
				Frenet.y -= Alpha_offset;
				Frenet.z -= Gamma_offset;

				Point_Int Diff;			
				if(this->Is_Bijection_Enabled)
					Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);
				else
					Diff = Frenet;
				/*Point_Int Center = BC + Diff;
				Point_Int Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);*/

				Point_Int Center = BC + Diff;				
				
				Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);

				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				for(int i = 0; i < (int)valence; i++)
				{
					int N1 = reg->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}
					reg = reg->next();
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
				{
					Selected_region = T_possible_bin[0];
				}
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}
				
				g = _pMesh.create_center_vertex(g);
				g->vertex()->point() = Center_vertex;

				g->vertex()->Region_Number = Selected_region;				

				if(Selected_region != -1)
					this->m_Number_Vertices_Per_Regions[Selected_region]++;
				
				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;

				g->vertex()->Vertex_Flag = CONQUERED;

				g = h;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				g->vertex()->Vertex_Flag = CONQUERED;
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;

				for (unsigned int k = 0; k < (valence - 1); k++)
				{
					g = g->next()->opposite()->next();
					g->facet()->Facet_Flag = TO_BE_REMOVED;
					g->vertex()->Vertex_Flag = CONQUERED;
					g->opposite()->vertex()->Vertex_Flag = CONQUERED;
					if (g->is_border_edge() == false)
						Halfedges.push(g->opposite());
				}
				
				g->next()->vertex()->Seed_Edge = -1;

				if ((this->IsColored) && (!this->IsOneColor))
				{						
					g = h;				

					Color_Unit Predicted_color;
					Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);					
					
					Color_Unit Color_difference;
					Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
					Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
					Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;								
					
					Color_Unit CV = Predicted_color + Color_difference;					

					g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * Color_step;
					LAB[1] = this->C1_Min + CV.c1 * Color_step;
					LAB[2] = this->C2_Min + CV.c2 * Color_step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);
				
					//#endif
				}
				if ((this->IsColored) && (this->IsOneColor))
				{
					g = h->next();
					g->vertex()->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
				}

			}
		}
		// In case of border edge.
		else if ((valence == 8) || (valence == 9))
		{
			
			type = Find_Type(h, valence - 5);
		
			Halfedge_handle pass = h;			
			
			Vector normal = Normal_Patch(pass, valence - 5);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(h, normal, T2);			
			
			if (normal == CGAL::NULL_VECTOR)
			{							
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);	
				normal = Vector(0,0,1);				
			}						
			
			pass = h;				
			Point_Int Frenet;
			if (check_alpha == false)
			{
				Frenet.x = Decoder.decode(alpha);
				Frenet.y = Decoder.decode(alpha);
			}
			else if (check_alpha == true)
			{
				Frenet.x = 0;
				Frenet.y = 0;
			}

			if (check_gamma == false)
				Frenet.z = Decoder.decode(gamma);
			else if (check_gamma == true)
				Frenet.z = 0;

			Frenet.x -= Alpha_offset;
			Frenet.y -= Alpha_offset;
			Frenet.z -= Gamma_offset;
			
			Point_Int Diff;			
			if(this->Is_Bijection_Enabled)
				Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);
			else
				Diff = Frenet;

			//Point_Int Center = BC + Diff;

			//Point_Int Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);			
			Halfedge_handle g = h;			
						
			Color_Unit Predicted_color;
			if ((this->IsColored) && (!this->IsOneColor))
			{
				Predicted_color.c0 = Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
				Predicted_color.c1 = Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
				Predicted_color.c2 = Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
			}			

			// border edge with valence == 3
			if (valence == 8)
			{				
				Point3d Barycenter = Barycenter_Patch_After_Removal(pass, 3);				
				Point_Int BC = Change_Real_Int(Barycenter,Component_ID);

				Point_Int Center = BC + Diff;				
				Point3d Center_vertex = this->Change_Int_Real(Center,Component_ID);
				
				//#ifdef PREDICTION_METHOD
				Color_Unit Average_color;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Average_color = Get_Average_Vertex_Color_After_Removal(g, 3);
				}
				//#endif

				vector<Halfedge_handle> Border_edges;
				int Number_jump = 0;
				
				if (g->next()->is_border_edge())
					Number_jump = 0;
				if (g->prev()->is_border_edge())
				{
					Number_jump = 1;
					g = g->next();
				}
				
				Border_edges.push_back(g->opposite());
				
				g = g->prev();
				Border_edges.push_back(g->opposite());

				_pMesh.erase_facet(g);

				Halfedge_handle Prev_edge = Border_edges[1]->opposite()->prev();
					
				// g points the new vertex
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				g->vertex()->point() = Center_vertex;
				g->vertex()->Vertex_Flag = CONQUERED;

				//#ifdef PREDICTION_METHOD
				Color_Unit CV;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					CV = Average_color + Predicted_color;				

					g->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->vertex()->color(RGB[0], RGB[1], RGB[2]);
				}
				//#endif
				
				if ((this->IsColored) && (this->IsOneColor))
				{
					g->vertex()->color(this->OnlyColor[0],this->OnlyColor[1],this->OnlyColor[2]);
				}


				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());			
				
				Halfedge_handle Tag_handle;

				if (Number_jump == 0)
					Tag_handle = Border_edges[1];
				if (Number_jump == 1)
					Tag_handle = Border_edges[0]->opposite();				
				Tag_handle->vertex()->Vertex_Flag = CONQUERED;

				if ((type == 1) || (type == 2) || (type == 4))
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = PLUS;
				}
				else if (type == 3)
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = MINUS;
				}
				for (int i=0; i <2; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}				
			}
			// border edge with valence == 4
			else if (valence == 9)
			{				
				int Number_jump = -1;				
				vector<Halfedge_handle> Border_edges;

				if ((type == 5) || (type == 8))
				{
					// jump == 0;
					if (g->next()->is_border_edge())
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());
						Border_edges.push_back(g->prev()->opposite()->next()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					// jump == 1;
					else if (g->prev()->opposite()->next()->is_border_edge())
					{
						Number_jump = 1;
						
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[0]->opposite());
					}
					
					// jump == 2;
					else
					{
						Number_jump = 2;

						Border_edges.push_back(g->prev()->opposite()->next()->opposite());
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
				}
				else
				{
					if (g->prev()->is_border_edge())
					{
						Number_jump = 2;

						Border_edges.push_back(g->next()->opposite()->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					else if (g->next()->opposite()->prev()->is_border_edge())
					{
						Number_jump = 1;

						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
					else
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

				}
				
				g = h;
				
				Point3d p0 = Border_edges[0]->vertex()->point();
				Point3d p1 = Border_edges[0]->opposite()->vertex()->point();
				Point3d p2 = Border_edges[1]->vertex()->point();
				Point3d p3 = Border_edges[2]->vertex()->point();
				
				
				Point3d Barycenter((p0.x() + p1.x() + p2.x() + p3.x())/4.,
								   (p0.y() + p1.y() + p2.y() + p3.y())/4.,
								   (p0.z() + p1.z() + p2.z() + p3.z())/4.);
				
				Point_Int BC = Change_Real_Int(Barycenter,Component_ID);

				Point_Int Center = BC + Diff;				
				
				Point3d Center_vertex = this->Change_Int_Real(Center,Component_ID);				
				
				// to create the new facets
				Halfedge_handle Prev_edge = Border_edges[2]->opposite()->prev();
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[2]->opposite());
				g->vertex()->point() = Center_vertex;				
					
				//#ifdef PREDICTION_METHOD
				Color_Unit CV;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Average_color;

					Color_Unit Color_0 = Get_Vertex_Color(Border_edges[0]->opposite());
					Color_Unit Color_1 = Get_Vertex_Color(Border_edges[0]);
					Color_Unit Color_2 = Get_Vertex_Color(Border_edges[1]);
					Color_Unit Color_3 = Get_Vertex_Color(Border_edges[2]);

					Average_color.c0 = (int)(Color_0.c0 + Color_1.c0 + Color_2.c0 + Color_3.c0)/4.0;
					Average_color.c1 = (int)(Color_0.c1 + Color_1.c1 + Color_2.c1 + Color_3.c1)/4.0;
					Average_color.c2 = (int)(Color_0.c2 + Color_1.c2 + Color_2.c2 + Color_3.c2)/4.0;					
					
					CV = Average_color + Predicted_color;					
			
					g->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->vertex()->color(RGB[0], RGB[1], RGB[2]);
				}
				//#endif
				if ((this->IsColored) && (this->IsOneColor))
				{
					g->vertex()->color(this->OnlyColor[0],this->OnlyColor[1],this->OnlyColor[2]);
				}
				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());


				//vertex_tag
				if (Number_jump == 0)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}

				else if (Number_jump == 1)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}

				}
				else // jump == 2
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}
				
				for (int i=0; i < 3; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					
					Border_edges[i]->vertex()->Vertex_Flag = CONQUERED;
					Border_edges[i]->opposite()->vertex()->Vertex_Flag = CONQUERED;
					
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}					
			}
		}

			
		
		
		//  the symbol == N
		else if (valence == 7)
		{
			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;

			if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)
				Halfedges.push(h->prev()->opposite());	
		}
	}	
}


// When the last loop doen not remove any vertex, all information are deleted
void Compression_Valence_Component::Remove_Last_Phase_Elements(const int & Component_ID)
{
	this->NumberDecimation[Component_ID] -= 1; //this->NumberDecimation[Component_ID] - 1;
	this->ComponentOperations[Component_ID] -= 1;	
	this->ListOperation[Component_ID].pop_front();

	for (int i = 0; i < (this->DumpSymbolDecimation + this->DumpSymbolRegulation);i++)	
		this->Connectivity[Component_ID].pop_front();		

	for (int i = 0; i < 2; i++)
	{			
		this->NumberSymbol[Component_ID].pop_front();
		this->NumberVertices[Component_ID].pop_front();
	}	
}



// To store information needed for the reconstruction of the base mesh.
void Compression_Valence_Component::Write_Base_Mesh(Polyhedron & _pMesh, Arithmetic_Codec & enc, unsigned &Connectivity_size, unsigned & Color_size, const int & Num_color_base_mesh)
{

	unsigned int Max_Qbit = 0;
	for (int i =0; i< this->NumberComponents; i++)
	{
		enc.put_bits(this->ComponentOperations[i], 8);	 									   // number of operations < 256	
		enc.put_bits(this->Qbit[i]-4, 4);															// quantization bit < 16
		enc.put_bits(this->NumberChangeQuantization[i], 4);					// number of decrease of quantization resolution
		enc.put_bits(this->NumberColorQuantization[i], 3);

		if (this->Qbit[i] > Max_Qbit)
			Max_Qbit = this->Qbit[i];
	}
		
	enc.put_bits(_pMesh.size_of_vertices(), 15);									    // number of vertices of base mesh < 4096
	enc.put_bits(_pMesh.size_of_facets(), 16);										  // number of facets of base mesh < 8192	

	//int Base_color_index_bit = 0;

	int Basemesh_vertex_number = 0;
	vector<int> Seed_edges(2 * this->NumberComponents, -1);

		
	// Encoding of vertex information of base mesh //
	for (Vertex_iterator pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); Basemesh_vertex_number++, pVertex++)
	{			
		pVertex->Vertex_Number = Basemesh_vertex_number;
		
		//int SEDD = pVertex->Seed_Edge;

		if (pVertex->Seed_Edge != OTHER_COORDINATE)
			Seed_edges[pVertex->Seed_Edge] = Basemesh_vertex_number;		
		
		//int cid = pVertex->Component_Number;

		Point_Int Vertex = Change_Real_Int(pVertex->point(), pVertex->Component_Number);
		
		
		enc.put_bits(Vertex.x, Max_Qbit + 1);
		enc.put_bits(Vertex.y, Max_Qbit + 1);
		enc.put_bits(Vertex.z, Max_Qbit + 1);		
		
		if ((this->IsColored) && (!this->IsOneColor))
		{
			//#ifdef PREDICTION_METHOD
			int C0 = pVertex->color_int(0);
			int C1 = pVertex->color_int(1);
			int C2 = pVertex->color_int(2);
			
			enc.put_bits(C0, C0_QUANTIZATION);
			enc.put_bits(C1, C1_QUANTIZATION);
			enc.put_bits(C2, C2_QUANTIZATION);
					
			Color_size += 3 * C0_QUANTIZATION;
			//#endif
		}		
	}
	
	// Bits needed for each edge.
	int Facet_index_bit = (int)ceil(log((double)(_pMesh.size_of_vertices()+1))/log((double)2));
	
	int Count_facet_index = 0;	 
	for (Facet_iterator pFacet = _pMesh.facets_begin() ; pFacet != _pMesh.facets_end() ; pFacet++)
	{
		Count_facet_index++;

		Halfedge_handle pHalfedge = pFacet->halfedge();				
		do
		{
			enc.put_bits(pHalfedge->vertex()->Vertex_Number, Facet_index_bit);			
			pHalfedge = pHalfedge->next();
			

		} while (pHalfedge != pFacet->halfedge());		
	}
	
	// Store seed edge information.
	for (int i = 0 ; i < (int)Seed_edges.size(); i++)	// MT
		enc.put_bits(Seed_edges[i], Facet_index_bit);	
	
	Connectivity_size += Facet_index_bit * (Count_facet_index + Seed_edges.size());	
}




void Compression_Valence_Component::Calculate_Geometry_Color_Offset_Range()
{
	for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
	{
		list<int>::iterator Number_iterator = this->NumberVertices[Component_ID].begin();
		list<Point_Int>::iterator Vertex_iterator = this->Geometry[Component_ID].begin();
		
		unsigned Number_phases = this->NumberVertices[Component_ID].size();
		
		
		for (unsigned i = 0; i < Number_phases; i++)
		{
			int Number_vertices_layer = *Number_iterator;
			Number_iterator++;

			int alpha_max = -10000;
			int alpha_min = 10000;
			int gamma_max = -10000;
			int gamma_min = 10000;
			int alpha = 0,beta = 0,gamma = 0;	
		

			if (Number_vertices_layer != 0)
			{
				for (int j = 0; j < Number_vertices_layer; j++)
				{
					alpha = Vertex_iterator->x;
					if (alpha > alpha_max)
						alpha_max = alpha;
					if (alpha < alpha_min)
						alpha_min = alpha;

					beta = Vertex_iterator->y;
					if (beta > alpha_max)
						alpha_max = beta;
					if (beta < alpha_min)
						alpha_min = beta;

					gamma = Vertex_iterator->z;
					if (gamma > gamma_max)
						gamma_max = gamma;
					if (gamma < gamma_min)
						gamma_min = gamma;

					Vertex_iterator++;				
				}

				this->AlphaRange[Component_ID].push_back(alpha_max - alpha_min + 1);
				this->AlphaOffset[Component_ID].push_back(-alpha_min);

				this->GammaRange[Component_ID].push_back(gamma_max - gamma_min + 1);
				this->GammaOffset[Component_ID].push_back(-gamma_min);			
			}
			else
			{
				this->AlphaRange[Component_ID].push_back(0);
				this->AlphaOffset[Component_ID].push_back(0);

				this->GammaRange[Component_ID].push_back(0);
				this->GammaOffset[Component_ID].push_back(0);			
			}
		}
	}

	
	
	/* Calculate alpha_min and gamma_min of all coefficients 
	 * in order to prevent negative symbols.Compression is not possible. */
	list<int>::iterator it_gamma, it_alpha;
	for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
	{
		for (it_alpha = this->AlphaOffset[Component_ID].begin();it_alpha != this->AlphaOffset[Component_ID].end();it_alpha++)
		{
			if (*it_alpha < this->Smallest_Alpha)
				this->Smallest_Alpha = *it_alpha;
		}
		for (it_gamma = this->GammaOffset[Component_ID].begin();it_gamma != this->GammaOffset[Component_ID].end();it_gamma++)
		{
			if (*it_gamma < this->Smallest_Gamma)
				this->Smallest_Gamma = *it_gamma;
		}
	}
	


	if (this->IsColored)
	{

		int C0_min = 50000, C1_min = 50000, C2_min = 50000;
		int C0_max = -50000, C1_max = -50000, C2_max = -50000;

		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{
			//#ifdef PREDICTION_METHOD
			list<Color_Unit>::iterator Vertex_color_iterator;
			for (Vertex_color_iterator = this->VertexColor[Component_ID].begin(); Vertex_color_iterator != this->VertexColor[Component_ID].end(); Vertex_color_iterator++)
			//#endif
			{
				if (Vertex_color_iterator->c0 < C0_min)
					C0_min = Vertex_color_iterator->c0;
				if (Vertex_color_iterator->c0 > C0_max)
					C0_max = Vertex_color_iterator->c0;
				
				if (Vertex_color_iterator->c1 < C1_min)
					C1_min = Vertex_color_iterator->c1;
				if (Vertex_color_iterator->c1 > C1_max)
					C1_max = Vertex_color_iterator->c1;

				if (Vertex_color_iterator->c2 < C2_min)
					C2_min = Vertex_color_iterator->c2;
				if (Vertex_color_iterator->c2 > C2_max)
					C2_max = Vertex_color_iterator->c2;
			}

		}

		
		this->C0_Range = C0_max - C0_min + 1;
		this->C1_Range = C1_max - C1_min + 1;
		this->C2_Range = C2_max - C2_min + 1;

		if (this->C0_Range <= 1) this->C0_Range = 2;
		if (this->C1_Range <= 1) this->C1_Range = 2;
		if (this->C2_Range <= 1) this->C2_Range = 2;
			
		this->Smallest_C0 = C0_min;
		this->Smallest_C1 = C1_min;
		this->Smallest_C2 = C2_min;
	}
}


//void Compression_Valence_Component::Color_Metric_Roy(Polyhedron &_pMesh, double &min, double &max,double &mean, double &rms)
//{
//	Mesh_roy * Simplified = new Mesh_roy;
//	Mesh_roy * Temp = new Mesh_roy;	
//	
//	Vector3d Pos, Color;
//	Vector3i Face;
//
//	for (int i = 0; i < this->Original->VertexNumber(); i++)
//	{
//		Pos = this->Original->Vertex(i);
//		Temp->AddVertex(Pos);
//		Color = this->Original->Color(i);
//		Temp->AddColor(Color);
//	}
//	for (int i = 0; i < this->Original->FaceNumber(); i++)
//	{
//		Face = this->Original->Face(i);
//		Temp->AddFace(Face);
//	}
//
//
//	for (Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
//	{
//		Point3d pt = pVert->point();
//
//		Pos[0] = pt.x();
//		Pos[1] = pt.y();
//		Pos[2] = pt.z();
//
//		Color[0] = pVert->color(0);
//		Color[1] = pVert->color(1);
//		Color[2] = pVert->color(2);
//
//		Simplified->AddVertex(Pos);
//		Simplified->AddColor(Color);
//	}
//	
//	_pMesh.set_index_vertices();
//
//	for (Facet_iterator pFacet = _pMesh.facets_begin(); pFacet != _pMesh.facets_end(); pFacet++)
//	{
//		int count = 0;
//
//		Halfedge_around_facet_circulator pH = pFacet->facet_begin();
//		do
//		{
//			Face[count] = pH->vertex()->tag();
//			count++;
//		}while(++pH != pFacet->facet_begin());
//		
//		Simplified->AddFace(Face);
//	}
//
//	
//	Deviation * dev = new Deviation;
//
//	dev->Initialization(Temp, Simplified, 0, 0.5);
//	dev->SetDeviationColorBound(0);
//
//	bool IsOK = dev->Compute(COLOR_DEVIATION);
//
//	min = dev->Min();
//	max = dev->Max();
//	mean = dev->Mean();
//	rms = dev->Rms();	
//	
//	delete dev;
//	delete Simplified;
//	delete Temp;
//}


void Compression_Valence_Component::Simplification(Polyhedron  & _pMesh,
									   const int   & NVertices, 
									   const bool    Normal_flipping,
									   const bool    Use_metric,
									   const float & Metric_thread, 
									   const bool    Use_forget_metric,
									   const int   & Forget_value)
{		
	

	//bool Is_any_vertex_removed = true;
	
	unsigned Last_Number = 0;
	unsigned Current_Number = _pMesh.size_of_vertices();

	//int Operation_choice = -1;
	
	do
	{		
		Last_Number = Current_Number;
		
		// Simplify component by component
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{
			// if the ith component did not remove any vertex in last loop, it is not necessary to simplifiy it.
			if (this->ComponentOperations[Component_ID] == this->GlobalCountOperation)
			{
				unsigned Initial_number_vertices = _pMesh.size_of_vertices();

				this->Decimation_Conquest(_pMesh, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value, Component_ID);				
				this->Regulation(_pMesh, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value, Component_ID);				

				int Diff_number_vertices = _pMesh.size_of_vertices() - Initial_number_vertices;
				
				this->ComponentOperations[Component_ID] += 1;
				this->NumberDecimation[Component_ID] += 1;
				
				this->ListOperation[Component_ID].push_front(0);	
				
				if (Diff_number_vertices == 0)	
					this->Remove_Last_Phase_Elements(Component_ID);
			}
		}		

		Current_Number = _pMesh.size_of_vertices();
		if (Current_Number != Last_Number)
			this->GlobalCountOperation++;

		if (Current_Number < (unsigned)NVertices) // MT
			break;
		
	}while((Current_Number != Last_Number));
	
	_pMesh.compute_normals();
}

void Compression_Valence_Component::Compression(Polyhedron     & _pMesh, 
												const char     * File_Name, 
												const int      & _Qbit, 
												unsigned       & Connectivity_size, 
												unsigned       & Color_size, 
												unsigned       & Total_size)
												//const unsigned & Initial_file_size)
{
	// Calculate offset and range for the compression.
	this->Calculate_Geometry_Color_Offset_Range();	
	
	FILE * fp  = fopen(File_Name, "wb");													//Main FILE to save compression information.		
	
	int res;

	res=fwrite(&this->Smallest_Alpha, sizeof(int), 1, fp);				  			 // smallest value of alpha (to save the negative value)
	res=fwrite(&this->Smallest_Gamma, sizeof(int), 1, fp);						 	 // smallest value of gamma (to save the negative value)
	res=fwrite(&this->Initial_file_size, sizeof(unsigned), 1, fp);			// Intial size of the input file (To visualize during decompression)	
	res=fwrite(&this->NumberComponents, sizeof(int), 1, fp);
	
	for (int i =0; i<this->NumberComponents; i++)
	{
		res=fwrite(&this->Quantization_Step[i], sizeof(float), 1, fp);							    // Quantization_Step(step of quantization)
		res=fwrite(&this->xmin[i], sizeof(float), 1, fp);																	   // xmin value
		res=fwrite(&this->ymin[i], sizeof(float), 1, fp);																	   // ymin value
		res=fwrite(&this->zmin[i], sizeof(float), 1, fp);																	   // zmin value
	}
	
	
	// Type of mesh
	char Colored;
	char OneColor;	

	if (this->IsColored)
	{
		Colored = 1;
		if (this->IsOneColor)
			OneColor = 1;
		else
			OneColor = 0;
	}
	else			
		Colored = 0;	
	
	res=fwrite(&Colored, sizeof(char), 1, fp);
	if (this->IsColored)
		res=fwrite(&OneColor, sizeof(char), 1, fp);




	int Num_color_base_mesh = 0;
	if ((this->IsColored) && (!this->IsOneColor))
	{

		res=fwrite(&this->Color_Quantization_Step, sizeof(float), 1, fp);

		// En-tete pour la couleur
		res=fwrite(&this->C0_Min, sizeof(float), 1, fp); // smallest value of c0 
		res=fwrite(&this->C1_Min, sizeof(float), 1, fp); // smallest value of c1 
		res=fwrite(&this->C2_Min, sizeof(float), 1, fp); // smallest value of c2		
			
		res=fwrite(&this->Smallest_C0, sizeof(int), 1, fp);  
		res=fwrite(&this->Smallest_C1, sizeof(int), 1, fp); 
		res=fwrite(&this->Smallest_C2, sizeof(int), 1, fp);
		
		Color_size += sizeof(int) * 8 * 9;		
	
	}	
	
	if ((this->IsColored) && (this->IsOneColor))
	{

		res=fwrite(&this->OnlyColor[0], sizeof(float), 1, fp); // smallest value of c0 
		res=fwrite(&this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1 
		res=fwrite(&this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2		
	}

	// Declaration du codeur.
	Arithmetic_Codec enc(AC_BUFFER); 
	enc.start_encoder(); 	

	// To calculate connectivity rate
	Arithmetic_Codec Connectivity_encoder(AC_BUFFER); 
	Connectivity_encoder.start_encoder(); 	
		
	for (int i = 0; i < this->NumberComponents; i++)
	{
		if (this->IsClosed[i])
			enc.put_bits(0,1);			
					
		else
			enc.put_bits(1,1);			
	
	}

	if(this->Is_Bijection_Enabled)
		enc.put_bits(1,1);
	else
		enc.put_bits(0,1);

	/*	Write information of base mesh.
		geometry + connectivity + color information(if the mesh is colored). */
	this->Write_Base_Mesh(_pMesh, enc, Connectivity_size, Color_size, Num_color_base_mesh);	

	// To calculate color rate
	Arithmetic_Codec Color_enc(AC_BUFFER); 
	Color_enc.start_encoder();


	//#ifdef PREDICTION_METHOD

	Adaptive_Data_Model C0_model;
	Adaptive_Data_Model C1_model;
	Adaptive_Data_Model C2_model;
	
	//To calculate color rate	
	Adaptive_Data_Model PC_C0_model;
	Adaptive_Data_Model PC_C1_model;
	Adaptive_Data_Model PC_C2_model;	
	
	if ((this->IsColored) && (!this->IsOneColor))
	{
		enc.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
		enc.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
		enc.put_bits(this->C2_Range, C2_QUANTIZATION + 1);		

		C0_model.set_alphabet(this->C0_Range);
		C1_model.set_alphabet(this->C1_Range);
		C2_model.set_alphabet(this->C2_Range);		
		
		//To calculate color rate
		Color_enc.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
		Color_enc.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
		Color_enc.put_bits(this->C2_Range, C2_QUANTIZATION + 1);

		PC_C0_model.set_alphabet(this->C0_Range);
		PC_C1_model.set_alphabet(this->C1_Range);
		PC_C2_model.set_alphabet(this->C2_Range);		

	}
		
	this->DM_JCW_MOVE_ERROR.set_alphabet(3);	

	// Main loop of compression //
	for (int i = 0; i < this->GlobalCountOperation; i++)
	{		

		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (i < this->ComponentOperations[Component_ID])
			{
				int Number_connectivity_symbols;

				if (this->IsClosed[Component_ID])
					Number_connectivity_symbols = 5;
				else
					Number_connectivity_symbols = 7;

				int Type_operation = this->ListOperation[Component_ID].front();
				this->ListOperation[Component_ID].pop_front();				

				// decimation is chosen to be applied.
				if (Type_operation == 0) 
				{			
					enc.put_bits(0, 2);					
					
					for (int j = 0; j < 2 ; j++) //Decimation and regulation
					{
						Adaptive_Data_Model Connectivity;				
						
						if (j == 0)	Connectivity.set_alphabet(2);
						else		Connectivity.set_alphabet(Number_connectivity_symbols);
							
						Adaptive_Data_Model Temp_connectivity;
						if (j == 0)	Temp_connectivity.set_alphabet(2);
						else		Temp_connectivity.set_alphabet(Number_connectivity_symbols);						
						
						int Alpha_range = this->AlphaRange[Component_ID].front();
						this->AlphaRange[Component_ID].pop_front();
						int Alpha_offset = this->AlphaOffset[Component_ID].front();
						this->AlphaOffset[Component_ID].pop_front();

						int Gamma_range = this->GammaRange[Component_ID].front();
						this->GammaRange[Component_ID].pop_front();
						int Gamma_offset = this->GammaOffset[Component_ID].front();
						this->GammaOffset[Component_ID].pop_front();				

						enc.put_bits(Alpha_range, _Qbit+1);
						if (this->Smallest_Alpha < 0) 	enc.put_bits(Alpha_offset - this->Smallest_Alpha, _Qbit+1);
						else							enc.put_bits(Alpha_offset, _Qbit+1);

						enc.put_bits(Gamma_range, _Qbit+1);				
						if (this->Smallest_Gamma < 0)	enc.put_bits(Gamma_offset - this->Smallest_Gamma, _Qbit+1);
						else							enc.put_bits(Gamma_offset, _Qbit+1);				
						
						bool check_alpha = false;
						bool check_gamma = false;

						if ((Alpha_range == 0) || (Alpha_range == 1))
						{
							check_alpha = true;
							Alpha_range = 2;
						}
						if ((Gamma_range == 0) || (Gamma_range == 1))
						{
							check_gamma = true;
							Gamma_range = 2;
						}				
						
						Adaptive_Data_Model alpha(Alpha_range);
						Adaptive_Data_Model gamma(Gamma_range);						

						int Number_symbols = this->NumberSymbol[Component_ID].front();
						this->NumberSymbol[Component_ID].pop_front();

						for (unsigned k = 0; k < (unsigned)Number_symbols; k++)	// MT
						{
							unsigned symbol = this->Connectivity[Component_ID].front();
							this->Connectivity[Component_ID].pop_front();			

							enc.encode(symbol, Connectivity);

							// To calculare connectivity rate
							Connectivity_encoder.encode(symbol, Temp_connectivity);
														
							if (((j == 0) && (symbol != 1)) || ((j == 1) && (symbol != 4)))
							{
								Point_Int Coeff = this->Geometry[Component_ID].front();
								this->Geometry[Component_ID].pop_front();
								
								int x = Coeff.x + Alpha_offset;
								int y = Coeff.y + Alpha_offset;
								int z = Coeff.z + Gamma_offset;					
								
								if (check_alpha == false)
								{						
									enc.encode(x, alpha);
									enc.encode(y, alpha);						
								}
								if (check_gamma == false)
								{	
									enc.encode(z, gamma);						
								}								
							
								
								//#ifdef PREDICTION_METHOD
								if ((this->IsColored) && (!this->IsOneColor))
								{
									Color_Unit VC = this->VertexColor[Component_ID].front();
									this->VertexColor[Component_ID].pop_front();

									enc.encode(VC.c0 - this->Smallest_C0, C0_model);
									enc.encode(VC.c1 - this->Smallest_C1, C1_model);
									enc.encode(VC.c2 - this->Smallest_C2, C2_model);								

									// To calculate color rate
									Color_enc.encode(VC.c0 - this->Smallest_C0, PC_C0_model);
									Color_enc.encode(VC.c1 - this->Smallest_C1, PC_C1_model);
									Color_enc.encode(VC.c2 - this->Smallest_C2, PC_C2_model);
								}
								//#endif							
							}
						}
						alpha.reset();
						gamma.reset();										
					}
					
					if(!this->m_N_Errors.empty())
					{
						int tt_count = this->m_N_Errors.front();
						this->m_N_Errors.pop_front();
							
						//if(!this->m_JCW_Move_Error.empty())
						for(int tt_i = 0; tt_i < tt_count; tt_i++)
						{
							vector<int> Temp_error = this->m_JCW_Move_Error.front();
							this->m_JCW_Move_Error.pop_front();

							for(int m = 0; m < 3; m++)
							{
								int Te = Temp_error[m];

								if (Te == -1)
									Te = 2;

								enc.encode(Te, this->DM_JCW_MOVE_ERROR);
							}									
						}
					}
					
				}
		
				// Decrease of geometry quantization resolution.						
				else if(Type_operation == 1)
				{
					enc.put_bits(1, 2);				

					int Number_vertices = this->NumberQuantizationLayer[Component_ID].front();
					this->NumberQuantizationLayer[Component_ID].pop_front();
					
					Adaptive_Data_Model Under_quantization_model(8);
					
					for (int i = 0; i < Number_vertices; i++)
					{
						int Under_quantization_coeff = this->QuantizationCorrectVector[Component_ID].front();
						this->QuantizationCorrectVector[Component_ID].pop_front();

						enc.encode(Under_quantization_coeff, Under_quantization_model);						
					}		
				}
				
				// Decrease of color quantization resolution.
				else
				{
					enc.put_bits(2,2);
					Arithmetic_Codec Color_enc(AC_BUFFER); 
					Color_enc.start_encoder();

					int Number_vertices = this->NumberProcessedVertices[Component_ID].front();
					this->NumberProcessedVertices[Component_ID].pop_front();					
					
					Adaptive_Data_Model *Color_quantization_model = new Adaptive_Data_Model[COLOR_NUMBER];
					for(int i = 0 ; i < COLOR_NUMBER; i++)
						Color_quantization_model[i].set_alphabet(8);
					
					// To measure color rates
					Adaptive_Data_Model *Temp_quantization_model = new Adaptive_Data_Model[COLOR_NUMBER];
					for(int i = 0 ; i < COLOR_NUMBER; i++)
						Temp_quantization_model[i].set_alphabet(8);

					for(int i = 0; i < Number_vertices; i++)
					{
						int Color_index = this->ColorEncoderIndex[Component_ID].front();
						this->ColorEncoderIndex[Component_ID].pop_front();

						int Childcell_index = this->ColorChildcellIndex[Component_ID].front();
						this->ColorChildcellIndex[Component_ID].pop_front();				

						enc.encode(Childcell_index, Color_quantization_model[Color_index]);			

						Color_enc.encode(Childcell_index, Temp_quantization_model[Color_index]);
					}
					Color_size += Color_enc.stop_encoder() * 8;

					delete[] Color_quantization_model;
					delete[] Temp_quantization_model;
				}				
			}			
		}

	}			

	Connectivity_size += Connectivity_encoder.stop_encoder() * 8;
	Color_size += Color_enc.stop_encoder() * 8;

	enc.write_to_file(fp);	
	fclose(fp);

	FILE* f_size = fopen(File_Name,"rb");
	fseek(f_size, 0, SEEK_END);
	Total_size = ftell(f_size);
}



void Compression_Valence_Component::Calculate_Edge_Color_Difference(Polyhedron & _pMesh, 
																	const int & _Component_ID, 
																	double & _Max_color, 
																	double & _Mean_color,
																	int & Number_of_vertices)
{	

	float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;

	float Color_small_step = 0.0;

	if(this->NumberColorQuantization[_Component_ID] == 0)
	{
		Color_small_step = this->Color_Quantization_Step;
	}
	
	else
	{
		Color_small_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[_Component_ID]);	
	}

	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge != 2*_Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*_Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
	}	

	std::queue<Vertex*> vertices;		
	
	// push a input gate to begin loop 
	// in order to treat all vertices.	
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));	
	
	// To count number of vertices;;
	int Count_vertices = 0;
	double Max_color = -5000.;
	double Mean_color = 0.;
	double L_min =  5000, A_min =  5000, B_min =  5000;
	double L_max = -5000, A_max = -5000, B_max = -5000;

	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;		
		else
		{			
			v->Vertex_Flag = CONQUERED;	
			v->Vertex_Number = Count_vertices;			
						
			Halfedge_around_vertex_circulator hvc = v->vertex_begin();
			Halfedge_around_vertex_circulator phvc = hvc;			
			
			double Mean = 0.0;
			int Count = 0;
			CGAL_For_all(hvc, phvc)
			{
				if(hvc->opposite()->vertex()->Vertex_Flag == FREE)
				{
					Color_Unit Color_0, Color_1;
					Color_0.c0 = hvc->vertex()->color_int(0);
					Color_0.c1 = hvc->vertex()->color_int(1);
					Color_0.c2 = hvc->vertex()->color_int(2);

					Color_1.c0 = hvc->opposite()->vertex()->color_int(0);
					Color_1.c1 = hvc->opposite()->vertex()->color_int(1);
					Color_1.c2 = hvc->opposite()->vertex()->color_int(2);

					double LAB_0[3], LAB_1[3];

					LAB_0[0] = C0_min + Color_0.c0 * Color_small_step;
					LAB_0[1] = C1_min + Color_0.c1 * Color_small_step;
					LAB_0[2] = C2_min + Color_0.c2 * Color_small_step;

					LAB_1[0] = C0_min + Color_1.c0 * Color_small_step;
					LAB_1[1] = C1_min + Color_1.c1 * Color_small_step;
					LAB_1[2] = C2_min + Color_1.c2 * Color_small_step;

					if(LAB_0[0] > L_max)
						L_max = LAB_0[0];
					if(LAB_0[1] > A_max)
						A_max = LAB_0[1];
					if(LAB_0[2] > B_max)
						B_max = LAB_0[2];
					
					if(LAB_1[0] > L_max)
						L_max = LAB_1[0];
					if(LAB_1[1] > A_max)
						A_max = LAB_1[1];
					if(LAB_1[2] > B_max)
						B_max = LAB_1[2];
					
					if(LAB_0[0] < L_min)
						L_min = LAB_0[0];
					if(LAB_0[1] < A_min)
						A_min = LAB_0[1];
					if(LAB_0[2] < B_min)
						B_min = LAB_0[2];
					
					if(LAB_1[0] < L_min)
						L_min = LAB_1[0];
					if(LAB_1[1] < A_min)
						A_min = LAB_1[1];
					if(LAB_1[2] < B_min)
						B_min = LAB_1[2];
					double diff = 0.0;
					for(int i = 0; i < 3; i++)
					{
						diff += (LAB_0[i] - LAB_1[i]) * (LAB_0[i] - LAB_1[i]);
					}

					diff = sqrt(diff);
					
					Mean += diff;
					Count++;
				}
			}
			
			if (Count != 0)
			{
				Mean /= Count;
				Mean_color += Mean;
			}

			Halfedge_around_vertex_circulator h = v->vertex_begin();			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}

			//increment number of vertices;
			Count_vertices++;
		}
	}

	Max_color = (L_max - L_min) * (L_max - L_min) + 
		        (A_max - A_min) * (A_max - A_min) + 
		        (B_max - B_min) * (B_max - B_min);

	Max_color = sqrt(Max_color);
	_Max_color = Max_color;
	_Mean_color = 3 * Mean_color / Count_vertices;
	Number_of_vertices = Count_vertices;
}

// Differentes histogrammes pour chaque couleur.
void Compression_Valence_Component::Diminush_Color_Quantization_Precision(Polyhedron &_pMesh, const int _Component_ID)
{

	//int N = 16; // Number of neighboring colors to use;
	//float Color_diff_seuil = 30.0;

	/*float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;*/

	float Color_small_step = 0.0;

	if(this->NumberColorQuantization[_Component_ID] == 0)
		Color_small_step = this->Color_Quantization_Step;		
	
	else
		Color_small_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[_Component_ID]);		

	double Color_large_step = Color_small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();		
	
	while((hi->vertex()->Seed_Edge != 2*_Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*_Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
	}	

	//premiere passe pour sous quantifier et donne l'indice de symbol a chaque sommet
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		if (pVert->Component_Number == _Component_ID)
		{
			int LAB_init_C0 = pVert->color_int(0);
			int LAB_init_C1 = pVert->color_int(1);
			int LAB_init_C2 = pVert->color_int(2);				
			
			// Les coordonnees apres under-quantization
			int LAB_final_C0 = LAB_init_C0 / 2;
			int LAB_final_C1 = LAB_init_C1 / 2;
			int LAB_final_C2 = LAB_init_C2 / 2;
			
			int i = LAB_init_C0 % 2;
			int j = LAB_init_C1 % 2;
			int k = LAB_init_C2 % 2;
			
			int Q_index = Get_Correct_Vector(i, j, k);
			
			pVert->color_int(LAB_final_C0, LAB_final_C1, LAB_final_C2);
			pVert->Q_Index = Q_index;

			float RGB[3];
			
			// To re-colorify the vertex with decreased color precision
			vector<float> LAB;
			
			LAB.push_back(this->C0_Min + LAB_final_C0 * Color_large_step);
			LAB.push_back(this->C1_Min + LAB_final_C1 * Color_large_step);
			LAB.push_back(this->C2_Min + LAB_final_C2 * Color_large_step);

			LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

			pVert->color(RGB[0], RGB[1], RGB[2]);		
		}
	}	

	std::queue<Vertex*> vertices;	

	std::list<int> Inter_color_quantization;
	std::list<int> Inter_index_per_color;
	
	// push a input gate to begin loop 
	// in order to treat all vertices.	
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));

	int Vertex_index = 0;

	// Color table which stocks all present colors
        vector<vector<int> > Color_table;

	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;			
			v->Vertex_Number = Vertex_index;

			int Q_index = v->Q_Index;		
			vector<int> LAB;
		
			LAB.push_back(v->color_int(0));
			LAB.push_back(v->color_int(1));
			LAB.push_back(v->color_int(2));

			bool Is_existed_color = false;
			
			// C'est une couleur deja existante?
			for(unsigned i = 0 ; i < Color_table.size(); i++)
			{
				if((LAB[0] == Color_table[i][0]) && (LAB[1] == Color_table[i][1]) && (LAB[2] == Color_table[i][2]))
				{
					Is_existed_color = true;
					Inter_color_quantization.push_front(Q_index);
					Inter_index_per_color.push_front(i);
					break;
				}			
			}
			
			// Sinon, agrandir le tabeau.
			if(Is_existed_color == false)
			{
				Inter_color_quantization.push_front(Q_index);
				Inter_index_per_color.push_front(Color_table.size());

				Color_table.push_back(LAB);				
			}
			Halfedge_around_vertex_circulator h;

			// First_vertex -> To find the seed edge
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}

			// To find an deterministic order = a given vertex and a vertex 
			//with the highest value of vertex number
			else 
			{
				int Comp_number = -2;
				
				h = v->vertex_begin();

				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h, h2)
				{
					if(h->opposite()->vertex()->Vertex_Number > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if(h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}

			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}

			Vertex_index++;
		}
	}

	

	int ssss = Inter_color_quantization.size();
	this->NumberProcessedVertices[_Component_ID].push_front(ssss);

	while(!Inter_color_quantization.empty())
	{
		int index = Inter_color_quantization.front();
		Inter_color_quantization.pop_front();

		this->ColorChildcellIndex[_Component_ID].push_front(index);
	}
	while(!Inter_index_per_color.empty())
	{
		int Color = Inter_index_per_color.front();
		Inter_index_per_color.pop_front();
		
		this->ColorEncoderIndex[_Component_ID].push_front(Color);
	}
}

void Compression_Valence_Component::Recalculate_Component_Area(Polyhedron & _pMesh, const int & _Component_ID, int & Number_facets)
{
	double Area = 0.0;	
	Number_facets = 0;
	for(Facet_iterator pFacet = _pMesh.facets_begin(); pFacet != _pMesh.facets_end(); pFacet++)
	{
		if(pFacet->Component_Number == _Component_ID)
		{
			Number_facets += 1;
			Area += Area_Facet_Triangle(pFacet->halfedge());
		}
	}
	Area *= pow(((double)10.0/(double)this->HighestLengthBB[_Component_ID]), 2.0);
	this->ComponentArea[_Component_ID] = Area;
}


double Compression_Valence_Component::Calculate_Area(Polyhedron & _pMesh)
{	
	double mesh_area = 0;
	for(Facet_iterator pFacet = _pMesh.facets_begin(); pFacet != _pMesh.facets_end(); pFacet++)
	{
		Halfedge_handle h = pFacet->halfedge();
		mesh_area += Area_Facet_Triangle(h);
	}

	return mesh_area;
}

void Compression_Valence_Component::Augment_Geometry_Quantization_Precision(Polyhedron &_pMesh,Arithmetic_Codec & Decoder, const int & Component_ID)
{
	Adaptive_Data_Model Under_quantization_model(8);	

	// Premier_Pas == distance d'un Quantization_Step de grill de quantification (Q)
	double Small_step = 0.0;
			
	if (this->NumberChangeQuantization[Component_ID] == 0)	
		Small_step = this->Quantization_Step[Component_ID];		
	
	else
		Small_step = this->Quantization_Step[Component_ID] * pow(2.0, this->NumberChangeQuantization[Component_ID] - 1);		
	
	//double Large_step = Small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();			
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID +1))
		hi++;	

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin();pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;		
		pVert->Vertex_Number = -1;
		pVert->Component_Number = -1;
	}	
	//_pMesh.compute_normals();	
	
	std::queue<Vertex*> vertices;
	std::list<int> Layer_symbols;	
	
	// push first vertex to begin loop 
	// to treat all vertices.
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));
	
	int Vertex_index = 0;
	
	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;
			v->Vertex_Number = Vertex_index;
			v->Component_Number = Component_ID;
									
			unsigned valence = v->vertex_degree();
			
			Point3d *Neighbors = new Point3d[valence];			
			Halfedge_around_vertex_circulator hvc = v->vertex_begin();
			Halfedge_around_vertex_circulator hvc2 = hvc;		
			unsigned count_neighbor = 0;
			CGAL_For_all(hvc,hvc2)
			{
				Neighbors[count_neighbor] = hvc->opposite()->vertex()->point();
				count_neighbor++;
			}
			
			Point3d Center = v->point();
			vector<int> D(6,0);

			for (unsigned i = 0; i < valence; i++)
			{
				Vector Vect = Neighbors[i] - Center;				
				
				if (Vect.x() <= 0)					
				{
					int I_dist = Vect.x() / Small_step / 2;
					D[0] += abs(I_dist);
				}
				else					
				{
					int I_dist = Vect.x() / Small_step / 2;
					D[1] += abs(I_dist);
				}
					
				if (Vect.y() <= 0)					
				{
					int I_dist = Vect.y() / Small_step / 2;
					D[2] += abs(I_dist);
				}
				else					
				{
					int I_dist = Vect.y() / Small_step / 2;
					D[3] += abs(I_dist);
				}
					
				if (Vect.z() <= 0)					
				{
					int I_dist = Vect.z() / Small_step / 2;
					D[4] += abs(I_dist);
				}
				else							
				{
					int I_dist = Vect.z() / Small_step / 2;
					D[5] += abs(I_dist);
				}
			}

			vector<double> U(3,0);
			U[0] = abs((double)D[0] / (double)(D[0] + D[1]) - 0.5);
			U[1] = abs((double)D[2] / (double)(D[2] + D[3]) - 0.5);
			U[2] = abs((double)D[4] / (double)(D[4] + D[5]) - 0.5);		
			
			multimap<double, int> U_order;
			U_order.insert(pair<double, int>(U[0],0));
			U_order.insert(pair<double, int>(U[1],1));
			U_order.insert(pair<double, int>(U[2],2));			
			
			multimap<double,int>::iterator it;

			vector<int> Weight_axe(3,0);
			int Temp_weight = 1;
			for (it = U_order.begin(); it != U_order.end(); it++)
			{
				Weight_axe[it->second] = Temp_weight;
				Temp_weight++;
			}			
			
			vector<int> Priority(8,0);
			Priority[0] = Weight_axe[0]*D[0] + Weight_axe[1]*D[2] + Weight_axe[2]*D[4];
			Priority[1] = Weight_axe[0]*D[1] + Weight_axe[1]*D[2] + Weight_axe[2]*D[4];
			Priority[2] = Weight_axe[0]*D[0] + Weight_axe[1]*D[3] + Weight_axe[2]*D[4];
			Priority[3] = Weight_axe[0]*D[1] + Weight_axe[1]*D[3] + Weight_axe[2]*D[4];
			Priority[4] = Weight_axe[0]*D[0] + Weight_axe[1]*D[2] + Weight_axe[2]*D[5];
			Priority[5] = Weight_axe[0]*D[1] + Weight_axe[1]*D[2] + Weight_axe[2]*D[5];
			Priority[6] = Weight_axe[0]*D[0] + Weight_axe[1]*D[3] + Weight_axe[2]*D[5];
			Priority[7] = Weight_axe[0]*D[1] + Weight_axe[1]*D[3] + Weight_axe[2]*D[5];
			
			multimap<int, int> Priority_map;
			multimap<int, int>::iterator it_reorder;
			for (int i =0; i < 8; i++)
				Priority_map.insert(pair<int, int>(Priority[i],i));
			
			vector<int> Priority_reorder(8,0);
			
			int Temp_priority = 0;
			///////////////////
			bool Is_same_priority_value = false;

			it_reorder = Priority_map.begin(); 
			for(int i = 0; i < 7; i++)
			{
				int P0 = it_reorder->first;
				it_reorder++;
				int P1 = it_reorder->first;

				if (P0 == P1)
					Is_same_priority_value = true;
			}

			if(!Is_same_priority_value)
			{			
				for (it_reorder = Priority_map.begin(); it_reorder != Priority_map.end(); it_reorder++)
				{
					Priority_reorder[it_reorder->second] = 7 - Temp_priority;
					Temp_priority++;					
				}			
			}
			else
			{
				for(int i = 0; i < 8; i++)
				{
					Priority_reorder[i] = i;
				}
			}								
			
			int Reordered_number = Decoder.decode(Under_quantization_model);
			
			for (int i = 0 ; i < 8; i++)
			{
				if (Priority_reorder[i] == Reordered_number)
				{
					v->Q_Index = i;
					break;
				}
			}			
			
			Halfedge_around_vertex_circulator h;
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}
			else
			{
				int Comp_number = -2;
				h = v->vertex_begin();
				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number > Comp_number)
					{
						Comp_number = h->opposite()->vertex()->Vertex_Number;
					}
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}			
			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
				{
					vertices.push(&(*(h->opposite()->vertex())));
				}
			}
			Vertex_index++;

			delete []Neighbors;
		}
	}
	Vertex_iterator pVertex = NULL;
	for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
	{
		if (pVertex->Component_Number == Component_ID)
		{
			int Coeff[3];
			Get_Coefficient_Up_Quantization(pVertex->Q_Index,Coeff);
			Point3d New_position(pVertex->point().x() + Coeff[0] * Small_step / 2,
								 pVertex->point().y() + Coeff[1] * Small_step / 2,
								 pVertex->point().z() + Coeff[2] * Small_step / 2);
			pVertex->point() = New_position;
		}
		
	}	
	this->Qbit[Component_ID]++;
	this->NumberChangeQuantization[Component_ID]--;	
}



/*  Description : ADAPTIVE_QUANTIZATION
Decreasing of quantization resolution based on the prediction of PENG.
Opposite function is up_quantization. */
void Compression_Valence_Component::Diminush_Geometry_Quantization_Precision(Polyhedron &_pMesh, const int & Component_ID)
{		

	// stock three mins for the reconstruction.
	float _xmin = this->xmin[Component_ID];
	float _ymin = this->ymin[Component_ID];
	float _zmin = this->zmin[Component_ID];

	// Premier_Pas == distance d'un Quantization_Step de grill de quantification (Q)
	double Small_step = 0.0;
			
	if (this->NumberChangeQuantization[Component_ID] == 0)		
		Small_step = this->Quantization_Step[Component_ID];		
	
	else
		Small_step = this->Quantization_Step[Component_ID] * pow(2.0, this->NumberChangeQuantization[Component_ID]);		
		
	// Large_step == distance d'un Quantization_Step de grille de quantification(Q - 1)
	double Large_step = Small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin();pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
	}

	
	//premiere passe pour sous quantifie et donne l'indice de symbol a chaque sommet
	for (Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		if (pVert->Component_Number == Component_ID)
		{
			double x = pVert->point().x();
			double y = pVert->point().y();
			double z = pVert->point().z();
			
			// Qx, Qy, Qz sont des coordonnees avant under-Quantization
			int Qx = (int)(ceil((x - _xmin)/Small_step)) - 1;
			if (Qx == -1)
				Qx = 0;
			int Qy = (int)(ceil((y - _ymin)/Small_step)) - 1;
			if (Qy == -1)
				Qy = 0;
			int Qz = (int)(ceil((z - _zmin)/Small_step)) - 1;
			if (Qz == -1)
				Qz = 0;				

			// Les coordonnees apres under-quantization
			int Second_Qx = Qx / 2;
			int Second_Qy = Qy / 2;
			int Second_Qz = Qz / 2;
			
			int i = Qx % 2;
			int j = Qy % 2;
			int k = Qz % 2;
			
			int Q_index = Get_Correct_Vector(i, j, k);
			pVert->Q_Index = Q_index;

			// Les sommets sont deplaces vers le centre de cellule mere.
			pVert->point() = Point3d(_xmin+(Second_Qx + 0.5) * Large_step,
									 _ymin+(Second_Qy + 0.5) * Large_step,
									 _zmin+(Second_Qz + 0.5) * Large_step);			
		}
	}
	
	/////// appliquer calcul des normales?? ou Quantization_Step???	
	//_pMesh.compute_normals();	

	std::queue<Vertex*> vertices;	

	std::list<int> Layer_symbols;
	
	// push a input gate to begin loop 
	// in order to treat all vertices.
	
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));

	int Vertex_index = 0;


	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;			
			v->Vertex_Number = Vertex_index;			

			int Q_index = v->Q_Index;			
			unsigned Count_neighbor = 0;
			unsigned valence = v->vertex_degree();
			
			Point3d *Neighbors = new Point3d[valence];			
			Halfedge_around_vertex_circulator hvc = v->vertex_begin();
			Halfedge_around_vertex_circulator phvc = hvc;			
			
			CGAL_For_all(hvc,phvc)
			{
				Neighbors[Count_neighbor] = hvc->opposite()->vertex()->point();
				Count_neighbor++;
			}
			
			Point3d Center = v->point();			
			
			vector<int> D(6,0);

			for (unsigned i = 0; i < valence; i++)
			{
				Vector Vect = Neighbors[i] - Center;				
				
				if (Vect.x() <= 0)					
				{
					int I_dist = Vect.x() / Small_step / 2;
					D[0] += abs(I_dist);
				}
				else					
				{
					int I_dist = Vect.x() / Small_step / 2;
					D[1] += abs(I_dist);
				}
					
				if (Vect.y() <= 0)					
				{
					int I_dist = Vect.y() / Small_step / 2;
					D[2] += abs(I_dist);
				}
				else					
				{
					int I_dist = Vect.y() / Small_step / 2;
					D[3] += abs(I_dist);
				}
					
				if (Vect.z() <= 0)					
				{
					int I_dist = Vect.z() / Small_step / 2;
					D[4] += abs(I_dist);
				}
				else							
				{
					int I_dist = Vect.z() / Small_step / 2;
					D[5] += abs(I_dist);
				}
			}

			vector<double> U(3,0);
			U[0] = abs((double)D[0] / (double)(D[0] + D[1]) - 0.5);
			U[1] = abs((double)D[2] / (double)(D[2] + D[3]) - 0.5);
			U[2] = abs((double)D[4] / (double)(D[4] + D[5]) - 0.5);		
			
			multimap<double, int> U_order;
			U_order.insert(pair<double, int>(U[0],0));
			U_order.insert(pair<double, int>(U[1],1));
			U_order.insert(pair<double, int>(U[2],2));			
			
			multimap<double,int>::iterator it;

			vector<int> Weight_axe(3,0);
			int Temp_weight = 1;
			for (it = U_order.begin(); it != U_order.end(); it++)
			{
				Weight_axe[it->second] = Temp_weight;
				Temp_weight++;
			}			
			
			vector<int> Priority(8,0);
			Priority[0] = Weight_axe[0]*D[0] + Weight_axe[1]*D[2] + Weight_axe[2]*D[4];
			Priority[1] = Weight_axe[0]*D[1] + Weight_axe[1]*D[2] + Weight_axe[2]*D[4];
			Priority[2] = Weight_axe[0]*D[0] + Weight_axe[1]*D[3] + Weight_axe[2]*D[4];
			Priority[3] = Weight_axe[0]*D[1] + Weight_axe[1]*D[3] + Weight_axe[2]*D[4];
			Priority[4] = Weight_axe[0]*D[0] + Weight_axe[1]*D[2] + Weight_axe[2]*D[5];
			Priority[5] = Weight_axe[0]*D[1] + Weight_axe[1]*D[2] + Weight_axe[2]*D[5];
			Priority[6] = Weight_axe[0]*D[0] + Weight_axe[1]*D[3] + Weight_axe[2]*D[5];
			Priority[7] = Weight_axe[0]*D[1] + Weight_axe[1]*D[3] + Weight_axe[2]*D[5];
			
			multimap<int, int> Priority_map;
			multimap<int, int>::iterator it_reorder;
			for (int i =0; i < 8; i++)
				Priority_map.insert(pair<int, int>(Priority[i],i));
			
			vector<int> Priority_reorder(8,0);
			
			int Temp_priority = 0;
			bool Is_same_priority_value = false;

			it_reorder = Priority_map.begin(); 
			for(int i = 0; i < 7; i++)
			{
				int P0 = it_reorder->first;
				it_reorder++;
				int P1 = it_reorder->first;

				if (P0 == P1)
					Is_same_priority_value = true;
			}

			if(!Is_same_priority_value)
			{			
				for (it_reorder = Priority_map.begin(); it_reorder != Priority_map.end(); it_reorder++)
				{
					Priority_reorder[it_reorder->second] = 7 - Temp_priority;
					Temp_priority++;				
				}			
			}
			else
			{
				for(int i = 0; i < 8; i++)
				{
					Priority_reorder[i] = i;
				}
			}

			int Reordered_number = Priority_reorder[Q_index];			
			
			Halfedge_around_vertex_circulator h;

			// First_vertex -> To find the seed edge
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}

			else
			{
				int Comp_number = -2;
				h = v->vertex_begin();
				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}	
			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}
			Vertex_index++;

			delete []Neighbors;
			Layer_symbols.push_front(Reordered_number);		
		}
	}
	
	this->NumberQuantizationLayer[Component_ID].push_front(Layer_symbols.size());
	while(!Layer_symbols.empty())
	{
		int index = Layer_symbols.front();
		Layer_symbols.pop_front();

		this->QuantizationCorrectVector[Component_ID].push_front(index);
	}		
}


QString Compression_Valence_Component::Decompress_Init(Polyhedron &_pMesh)//, unsigned & _Initial_file_size, const char* File_Name)
{

	if (FILE *file2 = fopen(this->File_name.c_str(), "r"))
	{
		fseek(file2,0,SEEK_END);
		this->Compressed_file_size = ftell(file2);
		fclose(file2);
	}	

	_pMesh.clear();
	this->IsClosed.clear();
	this->Qbit.clear();
	this->xmin.clear();
	this->ymin.clear();
	this->zmin.clear();
	this->Quantization_Step.clear();
	this->NumberColorQuantization.clear();
	this->NumberChangeQuantization.clear();
	this->ComponentOperations.clear();
	
	this->Decompress_count = 0;	
	
	FILE* fp = fopen(this->File_name.c_str(),"rb");
	
	this->DM_JCW_MOVE_ERROR.set_alphabet(3);

	int res;		
			
	res=fread(&this->Smallest_Alpha, sizeof(int), 1, fp);
	res=fread(&this->Smallest_Gamma, sizeof(int), 1, fp);
	res=fread(&this->Initial_file_size, sizeof(unsigned), 1, fp);	 
	res=fread(&this->NumberComponents, sizeof(int), 1, fp);
	
	float Qpas;
	float t_xmin, t_ymin, t_zmin;
	
	// Read geometry information for each component;
	for (int i =0; i < this->NumberComponents; i++)
	{
		res=fread(&Qpas, sizeof(float), 1, fp); // quantization step
		res=fread(&t_xmin, sizeof(float), 1, fp); // x_min
		res=fread(&t_ymin, sizeof(float), 1, fp); // y_min
		res=fread(&t_zmin, sizeof(float), 1, fp); // z_min

		this->Quantization_Step.push_back(Qpas);
		this->xmin.push_back(t_xmin);
		this->ymin.push_back(t_ymin);
		this->zmin.push_back(t_zmin);
	}
	
	char Colored, One_color;
	res=fread(&Colored, sizeof(char), 1 , fp);		
	
	this->IsOneColor = false;
	if (Colored == 1)
	{
		this->IsColored = true;	
		res=fread(&One_color, sizeof(char), 1 , fp);
		if (One_color == 1)
			this->IsOneColor = true;		
	}
	else				
		this->IsColored = false;
	
	// read color information for each component
	if ((this->IsColored) && (!this->IsOneColor))
	{	
		res=fread(&this->Color_Quantization_Step, sizeof(float), 1, fp);	
		
		res=fread(&this->C0_Min, sizeof(float), 1, fp); // smallest absolute position of c0 
		res=fread(&this->C1_Min, sizeof(float), 1, fp); // smallest absolute position of c1
		res=fread(&this->C2_Min, sizeof(float), 1, fp); // smallest absolute position of c2
			
		res=fread(&this->Smallest_C0, sizeof(int), 1, fp);  //smallest quantized positions
		res=fread(&this->Smallest_C1, sizeof(int), 1, fp); 
		res=fread(&this->Smallest_C2, sizeof(int), 1, fp);			
	}
	
	if ((this->IsColored) && (this->IsOneColor))
	{
		res=fread(&this->OnlyColor[0], sizeof(float), 1, fp); // smallest absolute position of c0 
		res=fread(&this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1 
		res=fread(&this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2
	}
	this->Decoder.set_buffer(AC_BUFFER);
	this->Decoder.read_from_file(fp);
	
	// To know if each component is colored or not, and closed of not.
	for (int i = 0; i < this->NumberComponents; i++)
	{
		if (Decoder.get_bits(1) == 0)
			this->IsClosed.push_back(true);					
		else
			this->IsClosed.push_back(false);		
	}

	if (Decoder.get_bits(1) == 1)
		this->Is_Bijection_Enabled = true;
	else
		this->Is_Bijection_Enabled = false;
	
	this->GlobalCountOperation = -1;
	
	unsigned Max_Qbit = 0;
	for (int i = 0; i < this->NumberComponents; i++)
	{
		int Number_operation = Decoder.get_bits(8);  // Number of total operations.
		this->ComponentOperations.push_back(Number_operation);
		
		if (Number_operation > this->GlobalCountOperation)
			this->GlobalCountOperation = Number_operation;
			
		int Qbit = Decoder.get_bits(4); // Initial quantization bit of geometry
		Qbit += 4;
		this->Qbit.push_back(Qbit);
		
		int NCQ = Decoder.get_bits(4); // Number of change of geometry quantization
		this->NumberChangeQuantization.push_back(NCQ);
		
		int Number_color_quantization_change = Decoder.get_bits(3); // Number of change of color quantization
		this->NumberColorQuantization.push_back(Number_color_quantization_change);

		if (this->Qbit[i] > Max_Qbit)
			Max_Qbit = this->Qbit[i];
	}
	
	int Number_basemesh_vertex = Decoder.get_bits(15); // Number of vertices of base mesh
	int Number_basemesh_facet  = Decoder.get_bits(16); // Number of facets of base mesh
		
	// Vectors for generation of base mesh
	vector<Point3d> vlist;
	vector<int>     flist;
	vector<float>   clist;
	vector<int>     Color_index_list;

	for (int i = 0; i < Number_basemesh_vertex; i++)
	{
		Point_Int Pt_int;
		Pt_int.x = Decoder.get_bits(Max_Qbit + 1); // Read geometry info
		Pt_int.y = Decoder.get_bits(Max_Qbit + 1);
		Pt_int.z = Decoder.get_bits(Max_Qbit + 1);
		
		// All vertices have quantization precision of component 0
		// That'll be corrected below.
		Point3d Pt_real = Change_Int_Real(Pt_int, 0); 
		vlist.push_back(Pt_real);		

		if ((this->IsColored) && (!this->IsOneColor))
		{
			//#ifdef PREDICTION_METHOD
		
			Color_Unit TC;
			TC.c0 = Decoder.get_bits(C0_QUANTIZATION); // Read color info.
			TC.c1 = Decoder.get_bits(C1_QUANTIZATION);
			TC.c2 = Decoder.get_bits(C2_QUANTIZATION);

			float L = this->C0_Min + TC.c0 * this->Color_Quantization_Step;
			float a = this->C1_Min + TC.c1 * this->Color_Quantization_Step;
			float b = this->C2_Min + TC.c2 * this->Color_Quantization_Step;
			
			clist.push_back(L);
			clist.push_back(a);
			clist.push_back(b);							
		//	#endif			
		}		
	}	
	
	// Read connectivity information
	int Facet_index_bit = (int)ceil(log((double)(Number_basemesh_vertex + 1)) / log((double)2));
	
	for (int i = 0; i < (Number_basemesh_facet * 3); i++)
	{
		int v = Decoder.get_bits(Facet_index_bit);
		flist.push_back(v);
	}
	
	// Generation of base mesh using builder()
	CModifyBasemeshBuilder<HalfedgeDS, Polyhedron, Enriched_kernel> builder(&vlist, &flist, &clist, &Color_index_list);
	_pMesh.delegate(builder);
	
	_pMesh.compute_normals();

	// Seed Edges;
	map<int, int> Seed_Edges;

	for (int i = 0; i < 2*this->NumberComponents; i++)
	{
		int Vertex_number = Decoder.get_bits(Facet_index_bit);
		Seed_Edges.insert(pair<int,int>(Vertex_number, i));
	}	

	int Basemesh_vertex_number = 0;
	
	Vertex_iterator pVertex = NULL;
	
	map<int,int>::iterator Seed_edge_iterator = Seed_Edges.begin();
	
	int Count_detected_vertices = 0;
	for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); Basemesh_vertex_number++, pVertex++)
	{	
		if (Count_detected_vertices < this->NumberComponents * 2)
		{
			if (Basemesh_vertex_number == Seed_edge_iterator->first) 
			{
				pVertex->Seed_Edge = Seed_edge_iterator->second;
				Seed_edge_iterator++;
				Count_detected_vertices++;
			}	
		}
		else
			pVertex->Seed_Edge = OTHER_COORDINATE;

		pVertex->Component_Number = -1;
	}

	/*float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;*/

	float Color_small_step = 0.0;	

	if ((this->IsColored) && (!this->IsOneColor))
	{
		for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
		{
			/*float LAB[3];
			LAB[0] = pVertex->color(0);
			LAB[1] = pVertex->color(1);
			LAB[2] = pVertex->color(2);
			
			float RGB[3];
			LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
			pVertex->color(RGB[0], RGB[1], RGB[2]);		*/		
		}
	
		this->C0_Range = Decoder.get_bits(C0_QUANTIZATION + 1);
		this->Color_0_Model.set_alphabet(this->C0_Range);		
		this->C1_Range = Decoder.get_bits(C1_QUANTIZATION + 1);
		this->Color_1_Model.set_alphabet(this->C1_Range);		
		this->C2_Range = Decoder.get_bits(C2_QUANTIZATION + 1);
		this->Color_2_Model.set_alphabet(this->C2_Range);		
	}
	if ((this->IsColored) && (this->IsOneColor))
	{
		for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
		{
			pVertex->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
		}		
	}

	// To get number of vertices of each component and restore the real position of vertices
	//if (this->NumberComponents != 1)
	{
		_pMesh.tag_facets(-1);
		int Component_index = 0;
		
		for (int Component_number = 0; Component_number < this->NumberComponents; Component_number++)
		{
			Halfedge_iterator hi = _pMesh.halfedges_begin();			
	
			while((hi->vertex()->Seed_Edge != 2 * Component_number) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_number + 1))
				hi++;			

			Facet_handle fh = hi->facet();
			fh->tag(Component_index);

			std::list<Facet_handle> facets;
			facets.push_front(fh);

			while(!facets.empty())
			{
				Facet_handle F = facets.front();
				facets.pop_front();				

				F->tag(Component_index);

				Halfedge_around_facet_circulator pHalfedge = F->facet_begin();
				Halfedge_around_facet_circulator end = pHalfedge;
					
				CGAL_For_all(pHalfedge, end)
				{
					// tag the vertex to its corresponding component number
					if (pHalfedge->vertex()->Component_Number == -1)
					{
						pHalfedge->vertex()->Component_Number = Component_index;
						Point3d Wrong_position = pHalfedge->vertex()->point();
						
						// The correct position of vertex is restored.
						Point_Int Temp_pos = Change_Real_Int(Wrong_position, 0); 
						Point3d Real_position = Change_Int_Real(Temp_pos, Component_number);

						pHalfedge->vertex()->point() = Real_position;
						
						float Wrong_LAB[3];
						Wrong_LAB[0] = pHalfedge->vertex()->color(0);
						Wrong_LAB[1] = pHalfedge->vertex()->color(1);
						Wrong_LAB[2] = pHalfedge->vertex()->color(2);												

						int Original_LAB[3];
						Original_LAB[0] = (int)floor((Wrong_LAB[0] - this->C0_Min)/this->Color_Quantization_Step + 0.5);
						Original_LAB[1] = (int)floor((Wrong_LAB[1] - this->C1_Min)/this->Color_Quantization_Step + 0.5);
						Original_LAB[2] = (int)floor((Wrong_LAB[2] - this->C2_Min)/this->Color_Quantization_Step + 0.5);
												
						pHalfedge->vertex()->color_int(Original_LAB[0], Original_LAB[1], Original_LAB[2]);						

						if(this->NumberColorQuantization[Component_number] == 0)
							Color_small_step = this->Color_Quantization_Step;		
						else
							Color_small_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_number]);		

						float RGB[3];
						float LAB[3];
						LAB[0] = this->C0_Min + Original_LAB[0] * Color_small_step;
						LAB[1] = this->C1_Min + Original_LAB[1] * Color_small_step;
						LAB[2] = this->C2_Min + Original_LAB[2] * Color_small_step;

						LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
						pHalfedge->vertex()->color(RGB[0], RGB[1], RGB[2]);
					}

					Facet_handle pNFacet = pHalfedge->opposite()->facet();
					if (pNFacet != NULL && pNFacet->tag() == -1)
					{
						facets.push_front(pNFacet);
						pNFacet->tag(Component_index);
					}
				}
			}
			//this->ComponentNumberVertices.push_back(Number_vertices);
			Component_index++;
		}
	}

	this->IsDecompress = true;
	this->Current_level = 0;					
	//this->Initial_file_size = _Initial_file_size;
	this->Total_layer = this->GlobalCountOperation;

		
	float prog = (float)this->Calculate_Current_File_Size() / this->Compressed_file_size * 100;
	float ratio = 1/ ((float)this->Calculate_Current_File_Size() / this->Initial_file_size);

	QString string("Number of all levels : ");
	string += QString("%1").arg(int(this->Total_layer));
	string += "   |   ";			
	string += QString("Prog : %1 %").arg(prog, 3, 'f', 3);
	string += "   |   ";			
	string += QString("Ratio : %1 \n").arg(ratio, 0, 'f', 3);

	this->Prog.clear();
	this->Prog.push_back(prog);
	this->Ratio.clear();
	this->Ratio.push_back(ratio);
	return string;
}


// Description : To decode step by step - show intermediate meshes
int Compression_Valence_Component::Decompress_Each_Step(Polyhedron &_pMesh, const char* File_Name)
{		
	if (this->Decompress_count < this->GlobalCountOperation)
	{
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (this->Decompress_count < this->ComponentOperations[Component_ID])
			{	
				int Operation = Decoder.get_bits(2);
				if (Operation == 0)
				{
					this->Un_Regulation(_pMesh, Decoder, Component_ID);					
					this->Un_Decimation_Conquest(_pMesh, Decoder, Component_ID);
				}
				else if (Operation == 1)
					this->Augment_Geometry_Quantization_Precision(_pMesh, Decoder, Component_ID);		
				else if (Operation == 2)
					this->Augment_Color_Quantization_Precision(_pMesh, Decoder, Component_ID);				
			}
		}			
	}	
		
			
	this->Decompress_count++;	
	_pMesh.compute_normals();	

	return this->Decompress_count;
	
}


void Compression_Valence_Component::Augment_Color_Quantization_Precision(Polyhedron &_pMesh, Arithmetic_Codec & Decoder, const int & _Component_ID)
{

	Adaptive_Data_Model *Color_quantization_model = new Adaptive_Data_Model[COLOR_NUMBER];
	for(int i = 0 ; i < COLOR_NUMBER; i++)
		Color_quantization_model[i].set_alphabet(8);

	/*float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;*/

	float Color_small_step = 0.0;
	float Color_large_step = 0.0;	

	// define quantization step for (QC = large step) and (QC+1 = small_step)
	if(this->NumberColorQuantization[_Component_ID] == 0)		
		Color_large_step = this->Color_Quantization_Step;			
	else
		Color_large_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[_Component_ID]);
		
	Color_small_step = Color_large_step / (float)2.0;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();		
	
	while((hi->vertex()->Seed_Edge != 2*_Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*_Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
		pVert->Component_Number = -1;
	}	
	
	// Generation of color table
        vector<vector<int> > Color_table;

	std::queue<Vertex*> vertices;				
	// push a input gate to begin loop 
	// in order to treat all vertices.
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));

	int Vertex_index = 0;

	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;			
			v->Vertex_Number = Vertex_index;
			v->Component_Number = _Component_ID;

			vector<int> LAB;
		
			LAB.push_back(v->color_int(0));
			LAB.push_back(v->color_int(1));
			LAB.push_back(v->color_int(2));

			bool Is_existed_color = false;

			int Color_index = -1;
			
			//Color table creation
			for(int i = 0 ; i < (int)Color_table.size(); i++)
			{
				if((LAB[0] == Color_table[i][0]) && (LAB[1] == Color_table[i][1]) && (LAB[2] == Color_table[i][2]))
				{
					Is_existed_color = true;
					Color_index = Decoder.decode(Color_quantization_model[i]);

					break; 
				}			
			}

			if(Is_existed_color == false)
			{				
				Color_index = Decoder.decode(Color_quantization_model[Color_table.size()]);
				Color_table.push_back(LAB);
			}			
			
			v->Q_Index = Color_index;
			
			//unsigned Count_neighbor = 0;
			//unsigned valence = v->vertex_degree();

			Halfedge_around_vertex_circulator h;

			// First_vertex -> To find the seed edge
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}

			// Synchronization proble = find an deterministic order using the vertex with the highest vertex number
			else 
			{
				int Comp_number = -2;				
				h = v->vertex_begin();
				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h,h2)
				{
					if(h->opposite()->vertex()->Vertex_Number > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if(h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}

			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}
			Vertex_index++;			
		}

	}

	// Here, we increase a bit of quantization precision of corresponding component 
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		if (pVert->Component_Number == _Component_ID)
		{
			int Coeff[3];
			Get_Coefficient_Up_Quantization(pVert->Q_Index, Coeff);

			int LAB_init_C0 = pVert->color_int(0);
			int LAB_init_C1 = pVert->color_int(1);
			int LAB_init_C2 = pVert->color_int(2);				
			
			// Each coordinate is multiplied by 2		
			int LAB_final_C0 = LAB_init_C0 * 2;
			int LAB_final_C1 = LAB_init_C1 * 2;
			int LAB_final_C2 = LAB_init_C2 * 2;
			
			// We sum the childcell index to restore true color
			if(Coeff[0] == 1)	LAB_final_C0 += 1;
			if(Coeff[1] == 1)	LAB_final_C1 += 1;
			if(Coeff[2] == 1)	LAB_final_C2 += 1;			

			pVert->color_int(LAB_final_C0, LAB_final_C1, LAB_final_C2);		
			
			float Inter_color[3];
			float RGB[3];

			Inter_color[0] = this->C0_Min + LAB_final_C0 * Color_small_step;
			Inter_color[1] = this->C1_Min + LAB_final_C1 * Color_small_step;
			Inter_color[2] = this->C2_Min + LAB_final_C2 * Color_small_step;

			LAB_To_RGB(Inter_color[0], Inter_color[1], Inter_color[2], RGB);

			pVert->color(RGB[0], RGB[1], RGB[2]);			
		}
	}

	//this->Q_color++;
	this->NumberColorQuantization[_Component_ID]--;

	delete [] Color_quantization_model;
}





// Description : Change a point coordinates in real to integer coordinates
Point_Int Compression_Valence_Component::Change_Real_Int(const Point3d &pt, const int & Component_ID)
{
	Point_Int Point;

	double Quantization_step = 0.0;

	// If the quantization resolution is decreased, 
	// we increase the step of quantization by a power of two.
	if (this->NumberChangeQuantization[Component_ID] == 0)
		Quantization_step = this->Quantization_Step[Component_ID];
	else
		Quantization_step = this->Quantization_Step[Component_ID] * pow(2.0, (int)this->NumberChangeQuantization[Component_ID]);	
	
	float xmin = this->xmin[Component_ID];
	float ymin = this->ymin[Component_ID];
	float zmin = this->zmin[Component_ID];
	
	double x = pt.x();
	double y = pt.y();
	double z = pt.z();

	Point.x = (int)(ceil((x-xmin)/Quantization_step)) - 1;
	if (Point.x == -1)
		Point.x = 0;

	Point.y = (int)(ceil((y-ymin)/Quantization_step)) - 1;
	if (Point.y == -1)
		Point.y = 0;

	Point.z = (int)(ceil((z-zmin)/Quantization_step)) - 1;
	if (Point.z == -1)
		Point.z = 0;

	return Point;
}



// Description : Change a point coordinates in integer to real coordinates
Point3d Compression_Valence_Component::Change_Int_Real(const Point_Int &Point, const int & Component_ID)
{
	float Quantization_step = 0;
	
	// If the quantization resolution is decreased, 
	// we increase the step of quantization by a power of two.
	if (this->NumberChangeQuantization[Component_ID] == 0)
		Quantization_step = this->Quantization_Step[Component_ID];
	else
		Quantization_step = this->Quantization_Step[Component_ID] * pow(2.0,(int)this->NumberChangeQuantization[Component_ID]);		
		
	float xmin = this->xmin[Component_ID];
	float ymin = this->ymin[Component_ID];
	float zmin = this->zmin[Component_ID];	

	double x = xmin + (Point.x + 0.5) * Quantization_step;
	double y = ymin + (Point.y + 0.5) * Quantization_step;
	double z = zmin + (Point.z + 0.5) * Quantization_step;

	Point3d pt(x, y, z);

	return pt;
}




// To match first vertices flags between two meshes. Must be used when copying a mesh.
void Compression_Valence_Component::Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh)
{
	Vertex_iterator pVertex = NULL;
	Vertex_iterator pVertex2 = New_mesh.vertices_begin();

	for (pVertex = Original_mesh.vertices_begin(); pVertex != Original_mesh.vertices_end(); pVertex++, pVertex2++)
	{
		//pVertex2->Seed_Edge = pVertex->Seed_Edge;
		
		int LAB[3];
		LAB[0] = pVertex->color_int(0);
		LAB[1] = pVertex->color_int(1);
		LAB[2] = pVertex->color_int(2);
		pVertex2->color_int(LAB[0], LAB[1], LAB[2]);

		float RGB[3];
		RGB[0] = pVertex->color(0);
		RGB[1] = pVertex->color(1);
		RGB[2] = pVertex->color(2);

		pVertex2->color(RGB[0], RGB[1], RGB[2]);
	}	
}


//// To reordering the color index for mapping table method.
//int Compression_Valence_Component::Mapping_Table_Index_Reordering(Polyhedron &_pMesh)
//{
//	int			     New_index = 0;
//	int				 Num_color_base_mesh;
//	vector<int>	     Reordered_color_index;
//	Vertex_iterator  pVertex;
//	
//	this->ReorderingColorIndex.clear();
//
//	for (int i = 0; i < NUMBER_SEEDS; i++)
//		this->ReorderingColorIndex.push_back(-1);
//	
//	for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
//	{
//		int Current_vertex_index = pVertex->Vertex_Color_Index;
//		
//		if (this->ReorderingColorIndex[Current_vertex_index] == -1)
//		{
//			this->ReorderingColorIndex[Current_vertex_index] = New_index;
//			New_index++;
//		}
//	}	
//	Num_color_base_mesh = New_index;
//
//
//	list<int>::iterator Color_iter = this->ColorIndex.begin();
//
//	for (;Color_iter != this->ColorIndex.end(); Color_iter++)
//	{
//		int Current_vertex_index = *Color_iter;
//		
//		if (this->ReorderingColorIndex[Current_vertex_index] == -1)
//		{
//			this->ReorderingColorIndex[Current_vertex_index] = New_index;
//			New_index++;
//		}
//	}	
//
//	return Num_color_base_mesh;
//}
//


//void Compression_Valence_Component::Separate_Components(Polyhedron &_pMesh)
//{
//	CCopyPoly<Polyhedron, Enriched_kernel> Copy_Polyhedron;
//
//	int Number_components = _pMesh.nb_components();
//	
//	if (Number_components == 1)
//		return;
//
//	for (int i = 0; i < Number_components; i++)
//	{
//		Polyhedron * New_mesh = new Polyhedron;		
//		Copy_Polyhedron.copy(&(_pMesh), New_mesh);
//		New_mesh->tag_facets(-1);
//		int Component_index = 0;
//
//		for (Facet_iterator pFacet = New_mesh->facets_begin(); pFacet != New_mesh->facets_end(); pFacet++)
//		{
//			if (pFacet->tag() == -1)
//			{				
//				New_mesh->tag_component(pFacet, -1, Component_index);
//				Component_index++;
//			}
//		}
//		
//		for (int j = 0; j < Number_components - 1; j++)
//		{
//
//			for (Halfedge_iterator pHedge = New_mesh->halfedges_begin(); pHedge != New_mesh->halfedges_end(); )
//			{
//				Halfedge_handle h = pHedge;
//				pHedge++;							
//				
//				if (!h->is_border())
//				{
//					int Facet_tag = h->facet()->tag();
//
//					if (Facet_tag != i)
//					{
//						New_mesh->erase_connected_component(h);							
//						break;
//					}
//				}
//			}
//		}
//
//		QString Outputfile = QString("Separate_%1").arg(i);	//Outputfile += wxString::Format("%d",i);
//		Outputfile += ".off";
//		New_mesh->write_off(Outputfile.toStdString(), true, false);
//		delete New_mesh;
//	}	
//}

int Compression_Valence_Component::GetResolutionChange(Polyhedron *_pMesh, float Prec)
{
	int MinX,MinY,MaxX,MaxY;
	MinX = MinY = 1000;
	MaxX = MaxY = 0;
	GLdouble *model ;  GLdouble *proj ;  GLint *view;
	view  = new int[4096];
	model = new double[4096];
	proj  = new double[4096];
	glGetIntegerv (GL_VIEWPORT, view);
	glGetDoublev (GL_MODELVIEW_MATRIX, model);
	glGetDoublev (GL_PROJECTION_MATRIX, proj);

	GLdouble wx ; 
	GLdouble wy ; 
	GLdouble wz;

	for (Vertex_iterator pVertex = _pMesh->vertices_begin(); pVertex!= _pMesh->vertices_end(); pVertex++)
	{
		gluProject (pVertex->point().x(),pVertex->point().y(),pVertex->point().z(), model, proj, view, &wx, &wy, &wz);  // on simule la projection du sommet dans l'espace window
		if (wz>0. && wz<1)
		{
			if (wx < MinX)
				MinX = wx;
			if (wx > MaxX)
				MaxX = wx;
			if (wy < MinY)
				MinY = wy;
			if (wy > MaxY)
				MaxY = wy;
		}		
	}

	float aire = (MaxX - MinX)*(MaxY - MinY);
	delete[] view;
	delete[] model;
	delete[] proj;

	////~la moitie des triangle sont affich? on choisi n pixel/triangle
	
	return 2*floor(aire/Prec);
}

void Compression_Valence_Component::Convert_To_Spherical(const Point3d & Pt, double * Spheric)
{	
	double r = 0., theta = 0., phi = 0.;	

	r = ((double)Pt.x() - (double)this->m_VC[0]) * ((double)Pt.x() - (double)this->m_VC[0]) +
		((double)Pt.y() - (double)this->m_VC[1]) * ((double)Pt.y() - (double)this->m_VC[1]) +
		((double)Pt.z() - (double)this->m_VC[2]) * ((double)Pt.z() - (double)this->m_VC[2]);	
	
	r = sqrt(r);	
	
	theta = atan2((double)Pt.y() - (double)this->m_VC[1],(double)Pt.x() - (double)this->m_VC[0]);
	phi   = acos((double)(Pt.z() - (double)this->m_VC[2])/r);

	Spheric[0] = r;
	Spheric[1] = theta;
	Spheric[2] = phi;
}

void Compression_Valence_Component::Convert_To_Cartesian(const double * Spheric, double * Cartesian)
{
	double r = Spheric[0];
	double theta = Spheric[1];
	double phi = Spheric[2];

	Cartesian[0] = r * sin(phi) * cos(theta) + (double)m_VC[0];
	Cartesian[1] = r * sin(phi) * sin(theta) + (double)m_VC[1];	
	Cartesian[2] = r * cos(phi) + (double)m_VC[2];
}

QString Compression_Valence_Component::Joint_Compression_Watermarking(Polyhedron & _pMesh,
																	 const char  * _Input_File_Name, 
																	 const char  * _Output_File_Name, 
																	 const int   & _Number_bins, 
																	 const int   & _Number_regions,
																	 const int   & _Embedding_strength,
																	 const char  * _Embedding_message, 
																	 const bool    _Is_complete_reversibility_selected, 
																	 const bool    _Is_divide_regions_selected,
																	 const int   & _Thres_divide_regions, 
																	 const int   & _Qbit, 
																	 const int   & _NVertices, 
																	 const bool    _Normal_flipping, 
																	 const bool    _Use_metric, 
																	 const float & _Metric_thread,
																	 const bool    _Use_forget_metric, 
																	 const int   & _Forget_value)	
{

	Timer timer;
	timer.start();

	if (FILE *file = fopen(_Input_File_Name, "r"))
	{
		fseek(file,0,SEEK_END);
		this->Initial_file_size = ftell(file);
		fclose(file);
	}
	
	int Init_number_vertices = (int)_pMesh.size_of_vertices();

	// Current version deals only with the one-component mesh.
	int NB_components = _pMesh.calc_nb_components();
	if (NB_components != 1)
		return 0;

	// read bits to insert
	this->Read_Information_To_Hide(_Embedding_message);
	
	// Set number of histogram bis.
	this->Set_Number_Bin(_Number_bins);

	// Set shifting level of histogram.
	this->Set_Embedding_Level(_Embedding_strength);

	// Set number of regions for segmentation.
	this->Set_Number_Region(_Number_regions);
	
	this->Number_Save_Over_bins = _Embedding_strength * 2 + 3;
	
	// Calculate mesh center.
	this->JCW_Calculate_Mesh_Center(_pMesh);
	this->Multiple_Components_Initialization(_pMesh, _Qbit);

	Polyhedron * Temp_mesh = new Polyhedron;		
	Copy_Polyhedron.copy(&(_pMesh), Temp_mesh);
	
	// Calculate radius of each vertices
	// Rmin and Rmax are used to determine delta_r.
	this->JCW_Calculate_Radius(*Temp_mesh);

	// Expansion of mesh
	// After the displacement of vertex
	// its position can be outside of initial mesh bounding box.
	this->JCW_Expand_Mesh(*Temp_mesh);

	this->Is_Division_Big_Regions_Enabled = _Is_divide_regions_selected;
	this->Is_Complete_Reversibility_Enabled = _Is_complete_reversibility_selected;

	int division, reversibility;
	if(this->Is_Division_Big_Regions_Enabled)
		division = 1;
	else
		division = 0;
	if (this->Is_Complete_Reversibility_Enabled)
		reversibility = 1;
	else
		reversibility = 0;	

	FILE * f_Dist = fopen("Keys.txt", "ab");

	int res;

	res=fwrite(&this->m_Dist, sizeof(double), 1, f_Dist);
	res=fwrite(&this->m_NumberBin, sizeof(int), 1, f_Dist);
	res=fwrite(&this->m_NumberRegion, sizeof(int), 1, f_Dist);
	res=fwrite(&this->m_EmbeddingStrength, sizeof(int), 1, f_Dist);
	res=fwrite(&division, sizeof(int), 1, f_Dist);
	res=fwrite(&reversibility, sizeof(int), 1, f_Dist);
	res=fwrite(&_Thres_divide_regions, sizeof(int), 1, f_Dist);

	fclose(f_Dist);
	delete Temp_mesh;
	
	this->Color_Initialization(_pMesh);	

	// Quantization of the bounding box obtained by previous expansion.
	this->Quantization(_pMesh);	
	
	//_pMesh.write_off("Quantized.off",false,false);

	//bool Is_any_vertex_removed = true;	
	unsigned Last_Number = 0;
	unsigned Current_Number = _pMesh.size_of_vertices();
	//int Operation_choice = -1;
	
	do
	{		
		Last_Number = Current_Number;		
		
		// if the ith component did not remove any vertex in last loop, it is not necessary to simplifiy it.
		if (this->ComponentOperations[0] == this->GlobalCountOperation)
		{
			unsigned Initial_number_vertices = _pMesh.size_of_vertices();			
			
			// Simplification of input mesh for segmentation.
			// Compression is 'not' followed.
			this->JCW_Decimation_For_Segmentation(_pMesh, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric,_Forget_value, 0);
			this->JCW_Regulation_For_Segmentation(_pMesh, _Normal_flipping, _Use_metric, _Metric_thread, _Use_forget_metric, _Forget_value, 0);

			int Diff_number_vertices = _pMesh.size_of_vertices() - Initial_number_vertices;
			
			this->ComponentOperations[0] += 1;
			this->NumberDecimation[0] += 1;
			
			this->ListOperation[0].push_front(0);	
			
			if (Diff_number_vertices == 0)	
				this->Remove_Last_Phase_Elements(0);
		}

		Current_Number = _pMesh.size_of_vertices();
		if (Current_Number != Last_Number)
			this->GlobalCountOperation++;

		if (Current_Number < (unsigned)_NVertices) // MT
			break;
		
	}while((Current_Number != Last_Number));		
	
	//_pMesh.write_off("Base_mesh.off", true, false);
	
	

	// Generation of region on base mesh.
	this->JCW_Generate_Regions_On_Base_Mesh(_pMesh);
	
	// create one another version of the base mesh and copy also vertices flags. 
	Polyhedron * Base_Mesh = new Polyhedron;
	Copy_Polyhedron.copy(&(_pMesh), Base_Mesh);
	
	Vertex_iterator pV_BM = Base_Mesh->vertices_begin();
	for(Vertex_iterator pV = _pMesh.vertices_begin(); pV != _pMesh.vertices_end(); pV++)
	{
		//pV->color(0.5, 0.5, 0.5);

		pV_BM->Seed_Edge = pV->Seed_Edge;
		pV_BM->Component_Number = pV->Component_Number;		

		pV_BM->m_color_int[0] = pV->m_color_int[0];
		pV_BM->m_color_int[1] = pV->m_color_int[1];
		pV_BM->m_color_int[2] = pV->m_color_int[2];

		pV_BM++;
	}	

	while(this->Decompress_count < this->ComponentOperations[0])
	{
		// Mesh 'TM' is used to get region index of each inserted vertex.
		Polyhedron * TM = new Polyhedron;		
		Copy_Polyhedron.copy(&(_pMesh), TM);
		
		Vertex_iterator pV_temp = TM->vertices_begin();
		for(Vertex_iterator pV = _pMesh.vertices_begin(); pV != _pMesh.vertices_end(); pV++)
		{
			pV_temp->Seed_Edge = pV->Seed_Edge;
			pV_temp->Region_Number = pV->Region_Number;
			pV_temp->Component_Number = pV->Component_Number;
			pV_temp++;
		}

		list<int> FP_Connectivity;
		list<Point3d> FP_Geometry;
		list<int> FP_Region_Number;		
		
		// These refine operations are used to get region information of each inserted vertex.
		this->JCW_Un_Regulation_For_Region_Detection(*TM, 0, FP_Connectivity, FP_Geometry, FP_Region_Number);
		this->JCW_Un_Decimation_For_Region_Detection(*TM, 0, FP_Connectivity, FP_Geometry, FP_Region_Number);
		delete TM;

		list<Point3d> SP_Moved_Position;
		list<Point_Int> SP_Watermarked_Position;
		list<Point3d> SP_Original_Position;

        list<vector<int> > JCW_ERROR;		

		// Insert operation gives the watermarked position of each inserted vertex and its position when we extract the watermark.
		this->JCW_Region_Mass_Center_Insert_Watermark(_pMesh, FP_Geometry, FP_Region_Number, SP_Watermarked_Position, SP_Moved_Position, SP_Original_Position, JCW_ERROR);
		
		// Real refinement operations.
		this->JCW_Un_Regulation_For_Insertion(_pMesh, 0, FP_Connectivity, SP_Moved_Position, SP_Original_Position, SP_Watermarked_Position, JCW_ERROR);
		this->JCW_Un_Decimation_For_Insertion(_pMesh, 0, FP_Connectivity, SP_Moved_Position, SP_Original_Position, SP_Watermarked_Position, JCW_ERROR);		

		// We move inserted vertices to its position when extracting watermark.
		// Synchronization with the decoding process.
		for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
		{
			if(pVert->Removal_Order == this->GlobalCountOperation - this->Decompress_count)
			{
				Point3d WP(pVert->Watermarked_Position[0], pVert->Watermarked_Position[1], pVert->Watermarked_Position[2]);
				pVert->point() = WP;
			}
		}

		if (this->Is_Division_Big_Regions_Enabled)						
			this->JCW_Divide_Big_Regions(_pMesh, _Thres_divide_regions);
		
		if (this->Is_Complete_Reversibility_Enabled)
			this->JCW_Code_Difference_Histogram_Shifting(_pMesh,0);

		this->Decompress_count++;
	}

	_pMesh.compute_normals();

	while(!this->JCW_Connectivity.empty())
	{
		int symbol = this->JCW_Connectivity.front() - 3;
		this->JCW_Connectivity.pop_front();

		this->Connectivity[0].push_back(symbol);
	}
	while(!this->JCW_Geometry.empty())
	{
		Point_Int Geo = this->JCW_Geometry.front();
		this->JCW_Geometry.pop_front();

		this->Geometry[0].push_back(Geo);
	}
	
	this->IsCompressed = true;
		
	_pMesh.compute_normals();

	unsigned Connectivity_size =0, Color_size=0, Total_size=0;
	int Number_inserted_bits = this->N_Inserted_Watermarks;

	this->Compression(*Base_Mesh, _Output_File_Name, _Qbit, Connectivity_size, Color_size, Total_size);//, Initial_file_size);

	double Connectivity_rate = (double)Connectivity_size / Init_number_vertices;
	double Color_rate = (double)Color_size / Init_number_vertices;
	double Total_rate = (double)Total_size * 8 / Init_number_vertices;
	double Geometry_rate = Total_rate - Connectivity_rate - Color_rate;

	timer.stop();
	double Time = timer.time();

	QString Res = QString("JCW done!\n\n");			
	Res += QString("Processing time : %1 s \n\n").arg(double(Time), 4, 'f', 3);			
	Res += QString("%1 bits inserted \n\n").arg(int(Number_inserted_bits));
	Res += QString("Connectivity : %1 b/v \n").arg(float(Connectivity_rate), 4, 'f', 3);					
	Res += QString("Geometry : ");
	Res += QString("%1").arg(float(Geometry_rate), 4, 'f', 3);
	Res += " b/v\n";
	Res += QString("Color : ");
	Res += QString("%1").arg(float(Color_rate), 4, 'f', 3);
	Res += " b/v\n";		
	Res += QString("Total size : ");				
	Res += QString("%1").arg(float(Total_rate), 4, 'f', 3);
	Res += " b/v\n";

	//QString t = string1 + string2 + string3 + string4 + string5 + string6 + string7;	

	return Res;	
}

void Compression_Valence_Component::JCW_Calculate_Mesh_Center(Polyhedron &_pMesh)
{
	FILE * f_Center = fopen("Keys.txt", "wb");
	double xc = 0., yc = 0., zc = 0.;
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		xc += pVert->point().x();
		yc += pVert->point().y();
		zc += pVert->point().z();
	}
	
	xc /= _pMesh.size_of_vertices();
	yc /= _pMesh.size_of_vertices();
	zc /= _pMesh.size_of_vertices();

	this->m_VC[0] = xc;
	this->m_VC[1] = yc;
	this->m_VC[2] = zc;

	int res;
	
	res=fwrite(&xc, sizeof(double), 1, f_Center);
	res=fwrite(&yc, sizeof(double), 1, f_Center);
	res=fwrite(&zc, sizeof(double), 1, f_Center);

	fclose(f_Center);
}

void Compression_Valence_Component::JCW_Calculate_Radius(Polyhedron & _pMesh)
{
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		double Spheric[3];		
		Point3d Cartesian = pVert->point();
		
		this->Convert_To_Spherical(Cartesian, Spheric);
		pVert->Spherical_Coordinates(Spheric);		

		if (Spheric[0] < this->m_Rmin)
			this->m_Rmin = Spheric[0];
		if (Spheric[0] > this->m_Rmax)
			this->m_Rmax = Spheric[0];
	}

	this->m_Dist = (this->m_Rmax) / this->m_NumberBin;

	FILE * dist = fopen("dist.txt","w");
	fprintf(dist, "%lf", this->m_Rmax);
	fclose(dist);
}



void Compression_Valence_Component::JCW_Expand_Mesh(Polyhedron & _pMesh)
{
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		double Spheric[3];		
		
		Spheric[0] = pVert->Spherical_Coordinates(0);
		Spheric[1] = pVert->Spherical_Coordinates(1);
		Spheric[2] = pVert->Spherical_Coordinates(2);

		Spheric[0] += this->Number_Save_Over_bins * this->m_EmbeddingStrength * this->m_Dist;

		double Cart[3];
		
		this->Convert_To_Cartesian(Spheric, Cart);

		Point3d temp(Cart[0], Cart[1], Cart[2]);
		pVert->point() = temp;
	}
	
	_pMesh.compute_bounding_box();
	
	this->xmin[0] = _pMesh.xmin();
	this->ymin[0] = _pMesh.ymin();
	this->zmin[0] = _pMesh.zmin();
	this->xmax[0] = _pMesh.xmax();
	this->ymax[0] = _pMesh.ymax();
	this->zmax[0] = _pMesh.zmax();
}

void Compression_Valence_Component::JCW_Generate_Regions_On_Base_Mesh(Polyhedron & _pMesh)
{
	// Initialization of vertices
	int Vertex_count = 0;
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Number = Vertex_count;
		pVert->Region_Number = -1;
		Vertex_count++;
	}
	
	// Init of face flag
	for(Facet_iterator pFace = _pMesh.facets_begin(); pFace != _pMesh.facets_end(); pFace++)
		pFace->Facet_Flag = FREE;
		
	std::queue<Halfedge_handle> Halfedges;

	int NVertices = (int)_pMesh.size_of_vertices();
	
	// Pas is the step used to select the seed vertex.
	float Pas = (float)NVertices / (float)this->m_NumberRegion;
	
	set<int> List;
	set<int>::iterator it_List;

	for(int i = 0; i < this->m_NumberRegion; i++)
	{
		int Seed_vertex_index = floor(i * Pas + 0.5);
		List.insert(Seed_vertex_index);
	}
		

	// Find seed edges for each region.
	int Seed_count = 0;
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		int VN = pVert->Vertex_Number;
		it_List = List.find(VN);
		
		// If VN is in the seed list
		// We have to attribute a region index.
		if(it_List != List.end())
		{
			Halfedge_around_vertex_circulator hvc = pVert->vertex_begin();
			Halfedge_around_vertex_circulator hvc2 = hvc;

			int Minimal_index = 5000;
			
			CGAL_For_all(hvc,hvc2)
			{
				int VN2 = hvc->opposite()->vertex()->Vertex_Number;

				if(hvc->is_border_edge())
					continue;

				if(VN2 != -1)
				{
					it_List = List.find(VN2);
					if(it_List == List.end())
					{
						if (VN2 < Minimal_index)
							Minimal_index = VN2;
					}
				}
			}

			// if one seed vertex is surround by other seed vertices,
			//this vertex is not used as a seed vertex.
			if (Minimal_index != 5000)
			{
				pVert->Region_Number = Seed_count;

				hvc = pVert->vertex_begin();
				hvc2 = hvc;
				CGAL_For_all(hvc,hvc2)
				{
					if(hvc->opposite()->vertex()->Vertex_Number == Minimal_index)
					{
						hvc->opposite()->vertex()->Region_Number = Seed_count;

						Halfedges.push(hvc);
						Halfedges.push(hvc->opposite());						
					}
				}
				Seed_count++;
			}
		}
	}
	_pMesh.write_off("base.off",true,false);

	// vector to count number of vertices in each region
	for(int i = 0; i < this->m_NumberRegion; i++)
		this->m_Number_Vertices_Per_Regions.push_back(0);	

	while(!Halfedges.empty())
	{
		Halfedge_handle h = Halfedges.front();
		Halfedges.pop();

		int N1 = h->vertex()->Region_Number;
		int N2 = h->opposite()->vertex()->Region_Number;
		
		if (h->next()->vertex()->Region_Number != -1)
			continue;
		if (h->facet()==NULL)
			continue;
		if (h->facet()->Facet_Flag != FREE)
			continue;
		if ((N1 == -1) || (N2 == -1))
			continue;

		if (N1 == N2)
		{
			h->next()->vertex()->Region_Number = N1;
			this->m_Number_Vertices_Per_Regions[N1]++;

			h->facet()->Facet_Flag = CONQUERED;

			Halfedges.push(h->next()->opposite());
			Halfedges.push(h->prev()->opposite());			
		}
	}
	

	// If there are vertices without region index,
	// We run a loop.
	int Count;
	do
	{
		Count = 0;
		// To deal with the vertices which do not belong to any region.
		for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
		{
			if(pVert->Region_Number == -1)
			{
				int Selected_reg = 5000;
				int Number_vertices = 5000;

				Halfedge_around_vertex_circulator hvc = pVert->vertex_begin();
				Halfedge_around_vertex_circulator hvc2 = hvc;		
				
				CGAL_For_all(hvc,hvc2)
				{
					int N_reg = hvc->opposite()->vertex()->Region_Number;
					if (N_reg != -1)
					{
						if(this->m_Number_Vertices_Per_Regions[N_reg] < Number_vertices)
						{
							Number_vertices = this->m_Number_Vertices_Per_Regions[N_reg];
							Selected_reg = N_reg;								
						}
						if(this->m_Number_Vertices_Per_Regions[N_reg] == Number_vertices)
						{
							if(N_reg < Selected_reg)
								Selected_reg = N_reg;
						}
					}
				}
				if((Selected_reg >= 0) && (Selected_reg < this->m_NumberRegion))
				{
					pVert->Region_Number = Selected_reg;
					this->m_Number_Vertices_Per_Regions[Selected_reg]++;
				}
				else
					Count++;
			}				
		}	
	}while(Count != 0);
}

void Compression_Valence_Component::Initialize_Spherical_Coordinates(Polyhedron & _pMesh)
{
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		double Spheric[3];		
		Point3d Cartesian = pVert->point();
		
		this->Convert_To_Spherical(Cartesian, Spheric);
		pVert->Spherical_Coordinates(Spheric);
	}	
}

void Compression_Valence_Component::JCW_Region_Mass_Center_Insert_Watermark(Polyhedron & _pMesh, 
																			list<Point3d> & FP_Geometry, 
																			list<int> & FP_Region_Number, 
																			list<Point_Int> & SP_Watermarked_Position, 
																			list<Point3d> & SP_Moved_Position,
																			list<Point3d> & SP_Original_Position,
                                                                            list<vector<int> > & JCW_Error)
{		
    vector<vector<int> > Hist_inserted_vertices;
    vector<vector<int> > Hist_remained_vertices;

	vector<int> Temp(this->m_NumberBin + this->Number_Save_Over_bins, 0);
	
	for(int i = 0; i < this->m_NumberRegion; i++)
	{
		Hist_inserted_vertices.push_back(Temp);
		Hist_remained_vertices.push_back(Temp);
	}

	vector<Point3d> I_Geo;
	vector<int> I_RN;	
	
	while(!FP_Geometry.empty())
	{
		Point3d Geo = FP_Geometry.front();
		FP_Geometry.pop_front();
		
		I_Geo.push_back(Geo);
	}
	while(!FP_Region_Number.empty())
	{
		int RN = FP_Region_Number.front();
		FP_Region_Number.pop_front();

		I_RN.push_back(RN);
	}
	
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{		
		int RN = pVert->Region_Number;
		double Spheric[3];

		Point3d Pos = pVert->point();
		this->Convert_To_Spherical(Pos, Spheric);
		pVert->Spherical_Coordinates(Spheric);

		double R = Spheric[0];			

		int Hist_number = ceil(R / this->m_Dist) - 1;
		if(Hist_number == -1)
			Hist_number = 0;

		Hist_remained_vertices[RN][Hist_number]++;
	}	

	for(unsigned i = 0; i < I_Geo.size(); i++)
	{
		Point3d Pt = I_Geo[i];
		int RN = I_RN[i];

		double Spheric[3];
		this->Convert_To_Spherical(Pt, Spheric);

		double R = Spheric[0];

		int Hist_number = ceil(R / this->m_Dist) - 1;
		if(Hist_number == -1)
			Hist_number = 0;
			
		Hist_inserted_vertices[RN][Hist_number]++;
	}
	
	// random generaion of watermark bits
	vector<int> Watermarks_to_insert(this->m_NumberRegion, -1);
	
	//srand(time(NULL));

	//for(int i = 0; i < this->m_NumberRegion; i++)
	//{				
	//	int Bit_to_insert = rand() % 2;		
	//	//Bit_to_insert = 0;
	//	Watermarks_to_insert.push_back(Bit_to_insert);
	//}		

	vector<int> Cases;
	vector<int> Directions(this->m_NumberRegion, -1);	
	
	srand(time(NULL));

	// Calculate center of mass to determine the case of each region
	for(int i = 0; i < this->m_NumberRegion; i++)
	{
		double Center_mass_decimated = 0, Center_mass_remained = 0;
		int Count_decimated = 0, Count_remained = 0;

		for(int j = 0; j < this->m_NumberBin + this->Number_Save_Over_bins; j++)
		{
			Center_mass_decimated += j * Hist_inserted_vertices[i][j];
			Center_mass_remained += j * Hist_remained_vertices[i][j];

			Count_decimated += Hist_inserted_vertices[i][j];
			Count_remained += Hist_remained_vertices[i][j];		
		}
		
		// If the number of vertices < LIMIT_NUMBER, we do not embed watermark. 
		// It is because of the histogram perturbation caused by the quantization.
		if(Count_remained < LIMIT_NUMBER)
		{
			Cases.push_back(2);
			continue;
		}

		Center_mass_decimated /= (double)Count_decimated;
		Center_mass_remained /= (double)Count_remained;		
		
		double Diff_CM = Center_mass_decimated - Center_mass_remained;		
		
		// If the difference of center of mass < 2* embedding level, we can embed a bit
		if(abs(Diff_CM) < this->m_EmbeddingStrength)
		{
			Cases.push_back(0);
			//

			//for(int i = 0; i < this->m_NumberRegion; i++)
			//{				
			//	int Bit_to_insert = rand() % 2;		
			//	//Bit_to_insert = 0;
			//	Watermarks_to_insert.push_back(Bit_to_insert);
			//}	
			int Bit_to_insert;
			if(!this->m_Watermarks.empty())
			{
				Bit_to_insert = this->m_Watermarks.front();
				this->m_Watermarks.pop_front();

				this->N_Inserted_Watermarks++;				
			}
			else
			{
				this->N_Inserted_Watermarks++;				
				Bit_to_insert = rand() % 2;				
			}

			FILE * f_insert = fopen("Inserted_watermarks.txt", "a");				
			fprintf(f_insert, "%d\t%d\n", Bit_to_insert, this->Decompress_count);
			fclose(f_insert);

			Watermarks_to_insert[i] = Bit_to_insert;			
		}
		else
		{
			Cases.push_back(1);
			if(Center_mass_decimated - Center_mass_remained > 0)
				Directions[i] = 1;			
			else
				Directions[i] = 0;	
		}		
	}		
	
	// Move vertices according to bit to insert
	// We move only vertices which will be removed
	// To minimize deformation.
	for(unsigned i = 0; i < I_Geo.size(); i++)
	{	
		int Region_index = I_RN[i];			

		if((Cases[Region_index] == 0)  || (Cases[Region_index] == 1))
		{
			Point3d P_bar_R = I_Geo[i];
			Point_Int P_bar_I = Change_Real_Int(P_bar_R, 0);
			
			SP_Original_Position.push_back(P_bar_R);

			double P_bar_S[3];				
			this->Convert_To_Spherical(P_bar_R, P_bar_S);	

			int Old_hist_number = ceil(P_bar_S[0] / this->m_Dist) - 1;
			if(Old_hist_number == -1)
				Old_hist_number = 0;
		
			int Direction = -1;

			if (Cases[Region_index] == 0)
			{			
				if(Watermarks_to_insert[Region_index] == 0)
				{					
					P_bar_S[0] -= this->m_EmbeddingStrength * this->m_Dist;
					Direction = 0;										
				}
				else
				{					
					P_bar_S[0] += this->m_EmbeddingStrength * this->m_Dist;
					Direction = 1;
				}					
			}			
			else if(Cases[Region_index] == 1)
			{
				if(Directions[Region_index] == 0)
				{
					P_bar_S[0] -= this->m_EmbeddingStrength * this->m_Dist;
					Direction = 0;
				}		
				else
				{						
					P_bar_S[0] += this->m_EmbeddingStrength * this->m_Dist;
					Direction = 1;
				}
			}
			double Q_R[3];
			double Q_bar_R[3];
			this->Convert_To_Cartesian(P_bar_S, Q_R);
	
			Q_bar_R[0] = Q_R[0];
			Q_bar_R[1] = Q_R[1];
			Q_bar_R[2] = Q_R[2];		
		
			int Qx = ceil((Q_bar_R[0]-(double)this->xmin[0])/(double)this->Quantization_Step[0]) -1;		
			if(Qx == -1)
				Qx = 0;
			int Qy = ceil((Q_bar_R[1]-(double)this->ymin[0])/(double)this->Quantization_Step[0]) -1;		
			if(Qy == -1)
				Qy = 0;
			int Qz = ceil((Q_bar_R[2]-(double)this->zmin[0])/(double)this->Quantization_Step[0]) -1;		
			if(Qz == -1)
				Qz = 0;

			Point_Int T_Q;
			T_Q.x = Qx;
			T_Q.y = Qy;
			T_Q.z = Qz;
			SP_Watermarked_Position.push_back(T_Q);			

			Q_bar_R[0] = (double)this->xmin[0] + (Qx + 0.5)*(double)Quantization_Step[0];
			Q_bar_R[1] = (double)this->ymin[0] + (Qy + 0.5)*(double)Quantization_Step[0];
			Q_bar_R[2] = (double)this->zmin[0] + (Qz + 0.5)*(double)Quantization_Step[0];		

			Point3d Pt_Q_bar(Q_bar_R[0], Q_bar_R[1], Q_bar_R[2]);
		
			double Q_bar_S[3];
			Convert_To_Spherical(Pt_Q_bar, Q_bar_S);				

			// To verify if we can not obtain the original position
			
			if(Direction == 0)
				Q_bar_S[0] += this->m_EmbeddingStrength * this->m_Dist;
			else
				Q_bar_S[0] -= this->m_EmbeddingStrength * this->m_Dist;					

			double Pp_R[3];		

			Convert_To_Cartesian(Q_bar_S, Pp_R);		

			int Qx2 = ceil((Pp_R[0] - (double)xmin[0])/(double)this->Quantization_Step[0]) -1;
			if(Qx2 == -1)
				Qx2 = 0;
			int Qy2 = ceil((Pp_R[1] - (double)ymin[0])/(double)this->Quantization_Step[0]) -1;			
			if(Qy2 == -1)
				Qy2 = 0;
			int Qz2 = ceil((Pp_R[2] - (double)zmin[0])/(double)this->Quantization_Step[0]) -1;
			if(Qz2 == -1)
				Qz2 = 0;
			
			double MP_real[3];
			MP_real[0] = (double)this->xmin[0] + (Qx2 + 0.5)*(double)Quantization_Step[0];
			MP_real[1] = (double)this->ymin[0] + (Qy2 + 0.5)*(double)Quantization_Step[0];
			MP_real[2] = (double)this->zmin[0] + (Qz2 + 0.5)*(double)Quantization_Step[0];
			
			Point3d MP(MP_real[0], MP_real[1], MP_real[2]);
			SP_Moved_Position.push_back(MP);

			vector<int> Err;
			Err.push_back(Qx2 - P_bar_I.x);
			Err.push_back(Qy2 - P_bar_I.y);
			Err.push_back(Qz2 - P_bar_I.z);			
			
			JCW_Error.push_back(Err);
		}
		else
		{
			Point3d MP = I_Geo[i];
			Point_Int T_Q = this->Change_Real_Int(MP, 0);
			SP_Watermarked_Position.push_back(T_Q);
			SP_Moved_Position.push_back(MP);
			SP_Original_Position.push_back(MP);
			
			vector<int> Err(3, 0);
			JCW_Error.push_back(Err);
		}
	}
}

void Compression_Valence_Component::JCW_Region_Mass_Center_Extract_Watermark(Polyhedron & _pMesh)
{
	FILE * f_extract;
	if(this->Decompress_count == 0)
	{
		f_extract = fopen("Extracted_watermark.txt", "w");		
	}
	else
	{
		f_extract = fopen("Extracted_watermark.txt", "a");			
	}

        vector<vector<int> > Hist_inserted_vertices;
        vector<vector<int> > Hist_remained_vertices;

	vector<int> Temp(this->m_NumberBin + this->Number_Save_Over_bins, 0);
	for(int i = 0; i < this->m_NumberRegion; i++)
	{
		Hist_inserted_vertices.push_back(Temp);
		Hist_remained_vertices.push_back(Temp);
	}


	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		
		double R = pVert->Spherical_Coordinates(0);

		int Region_index = pVert->Region_Number;		

		int Hist_number = ceil(R / this->m_Dist) - 1;
		if(Hist_number == -1)
			Hist_number = 0;

		if(Hist_number >= this->m_NumberBin + this->Number_Save_Over_bins)
			Hist_number =  this->m_NumberBin + this->Number_Save_Over_bins - 1;
		
		int RO = this->GlobalCountOperation - this->Decompress_count;			
		
		if(pVert->Removal_Order == RO)
			Hist_inserted_vertices[Region_index][Hist_number]++;
		else
			Hist_remained_vertices[Region_index][Hist_number]++;
	
	}

	vector<int> Cases(this->m_NumberRegion, -1);
	vector<int> Watermarks(this->m_NumberRegion, -1);
	vector<int> Direction(this->m_NumberRegion, -1);		

	for(int i = 0; i < this->m_NumberRegion; i++)
	{
		float Center_mass_inserted = 0.;
		float Center_mass_remained = 0.;

		int Count_inserted = 0;
		int Count_remained = 0;

		for(int j = 0; j < this->m_NumberBin + this->Number_Save_Over_bins; j++)
		{
			Count_inserted += Hist_inserted_vertices[i][j];
			Count_remained += Hist_remained_vertices[i][j];

			Center_mass_inserted += j * Hist_inserted_vertices[i][j];
			Center_mass_remained += j * Hist_remained_vertices[i][j];
		}		

		
		if(Count_remained < LIMIT_NUMBER)
		{
			Cases[i] = 2;			
			continue;
		}

		Center_mass_inserted /= (float)Count_inserted;
		Center_mass_remained /= (float)Count_remained;
		
		double Diff_CM = Center_mass_inserted - Center_mass_remained;

		if (abs(Diff_CM) > 2 * this->m_EmbeddingStrength)
		{
			Cases[i] = 1;

			if(Diff_CM > 0)
				Direction[i] = 1;
			else
				Direction[i] = 0; 
		}
		else
		{
			Cases[i] = 0;			
			if(Diff_CM > 0)
			{
				Direction[i] = 1;
				Watermarks[i] = 1;
				fprintf(f_extract, "1\t%d\n", this->Decompress_count);
				this->m_Watermarks.push_back(1);
			}
			else
			{
				Direction[i] = 0;
				Watermarks[i] = 0;
				fprintf(f_extract, "0\t%d\n", this->Decompress_count);
				this->m_Watermarks.push_back(0);
			}
		}		
	}	
	fclose(f_extract);


	
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		
		double Spheric[3];
		Spheric[0] = pVert->Spherical_Coordinates(0);
		Spheric[1] = pVert->Spherical_Coordinates(1);
		Spheric[2] = pVert->Spherical_Coordinates(2);

		int Bin_index = pVert->Region_Number;
		
		if(Cases[Bin_index] != 2)
		{
			
			if(pVert->Removal_Order == this->GlobalCountOperation - this->Decompress_count)
			{
				if(Direction[Bin_index] == 0)
					Spheric[0] += this->m_EmbeddingStrength * this->m_Dist;
				else
					Spheric[0] -= this->m_EmbeddingStrength * this->m_Dist;
			}						
			double Cart[3];
			Convert_To_Cartesian(Spheric, Cart);

			int Qx = ceil((Cart[0] - this->xmin[0])/this->Quantization_Step[0]) - 1;
			if (Qx == -1)
				Qx = 0;
			int Qy = ceil((Cart[1] - this->ymin[0])/this->Quantization_Step[0]) - 1;
			if (Qy == -1)
				Qy = 0;
			int Qz = ceil((Cart[2] - this->zmin[0])/this->Quantization_Step[0]) - 1;
			if (Qz == -1)
				Qz = 0;

			Cart[0] = this->xmin[0] + (Qx + 0.5)*Quantization_Step[0];
			Cart[1] = this->ymin[0] + (Qy + 0.5)*Quantization_Step[0];
			Cart[2] = this->zmin[0] + (Qz + 0.5)*Quantization_Step[0];
			Point3d Temp(Cart[0], Cart[1], Cart[2]);

			pVert->point() = Temp;
		}
	}	
	//vector<int> res;
	//return res;
}

//void Compression_Valence_Component::JCW_Choose_Valid_Vertices(Polyhedron &_pMesh)
//{
	//// Initialize vertex flag.
	//// Initially all vertices are valid.
	//for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	//{
	//	pVert->Valid_Vertex = true;
	//}

	//vector<double> Comp_angles;

	//double Pas = (double)this->Quantization_Step[0] / (double)this->m_Dist;
	//double Ang = Pas / 2.0;

	//while(Ang <= 1.0)
	//{
	//	Comp_angles.push_back(Ang);
	//	Ang += Pas;
	//}

	//Ang = -Pas/2.0;
	//while(Ang >= -1.0)
	//{
	//	Comp_angles.push_back(Ang);
	//	Ang -= Pas;
	//}
	//
	//int Count = 0;

	//for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	//{
	//	double Spheric[3];
	//	Spheric[0] = pVert->Spherical_Coordinates(0);
	//	Spheric[1] = pVert->Spherical_Coordinates(1);
	//	Spheric[2] = pVert->Spherical_Coordinates(2);
	//		
	//	double Check_x = sin(Spheric[2]) * cos(Spheric[1]);
	//	double Check_y = sin(Spheric[2]) * sin(Spheric[1]);
	//	double Check_z = cos(Spheric[2]);

	//	double Min_x = 5000., Min_y = 5000., Min_z = 5000.;

	//	for(unsigned i = 0; i < Comp_angles.size(); i++)
	//	{
	//		if (abs(Check_x - Comp_angles[i]) < Min_x)
	//			Min_x = abs(Check_x - Comp_angles[i]);
	//		
	//		if (abs(Check_y - Comp_angles[i]) < Min_y)
	//			Min_y = abs(Check_y - Comp_angles[i]);
	//		
	//		if (abs(Check_z - Comp_angles[i]) < Min_z)
	//			Min_z = abs(Check_z - Comp_angles[i]);
	//	}

	//	if((Min_x < THRES_ANGLE) || (Min_y < THRES_ANGLE) || (Min_z < THRES_ANGLE))		
	//		pVert->Valid_Vertex = false;	
	//	if(pVert->Valid_Vertex == false)
	//		Count++;
	//}

//}

// Description : This function select a set of independent vertices to be removed
int Compression_Valence_Component::JCW_Decimation_For_Segmentation(Polyhedron  & _pMesh,
																   const bool    Normal_flipping,
																   const bool    Use_metric,
													               const float & Metric_thread,
													               const bool    Use_forget_metric,
													               const int   & Forget_value,
													               const int   & Component_ID)
{	
	
	// Calculate mean color and meah area for color metric
	double Max_color, Mean_color;
	int Temp_NV = 0;
	int Number_facets;
	double Mean_area;
	if(this->IsColored)
	{
		this->Calculate_Edge_Color_Difference(_pMesh, Component_ID, Max_color, Mean_color, Temp_NV);
		this->Recalculate_Component_Area(_pMesh, Component_ID, Number_facets);
		Mean_area= (double)this->ComponentArea[Component_ID] / (double)Number_facets;
	}


	// Initialize vertex and face flags.
	Init(_pMesh);
		
	// to count number of independent vertices and number of symbol of connectivity.
	int Number_vertices = 0; 
	int Number_symbol = 0;	

	// To find first edge.
	Halfedge_iterator hi = _pMesh.halfedges_begin();			
	
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID+1))
		hi++;

	Halfedge_handle First_halfedge = &(*(hi));
	
	// Two vertices of seed edge are flaged CONQUERED
	First_halfedge->vertex()->Vertex_Flag = CONQUERED;
	First_halfedge->opposite()->vertex()->Vertex_Flag = CONQUERED;

	// These vertices are also flaged with sign flags for retriangulation
	First_halfedge->vertex()->Vertex_Sign = PLUS;
	First_halfedge->opposite()->vertex()->Vertex_Sign = MINUS;	
	
	std::queue<Halfedge_handle> Halfedges; // halfedge queue. 	
	Halfedges.push(First_halfedge);	// push the first halfedge in the queue.

	Halfedge_handle h; // The current gate
	
	/// Main loop
	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		// define type for retriangulation
		unsigned type = 0;

		// valence
		int valence = (int)h->next()->vertex_degree(); 		

		if (h->is_border() == true)
			continue;		
		
		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED, do nothing!!
		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;

		// if its front vertex is free and has a valence <= 6 and it is not a border vertex.
		else if ( (h->next()->vertex()->Vertex_Flag == FREE) && 
				  (valence >= 3) && 
				  (valence <= 6) && 
				  (!Is_Border_Vertex(h->next())))
		{

			type = Find_Type(h, valence);
			
			// Check if the manifold property is violated.
			bool Manifold_condition = Is_Manifold_Property_Violated(h, type, valence);
			
            bool Normal_flipping_condition = false;
            
			if (Normal_flipping)
                Normal_flipping_condition = Is_Normal_Flipping_Occured(h, valence);
			
			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Geometric_metric_condition = false;
			
			if (Use_metric)
			{
				if (Use_forget_metric)
				{
					if (_pMesh.size_of_vertices() > (unsigned)Forget_value)
						Geometric_metric_condition = false;
					else
						Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
				}
				else
					Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
			}					
			bool Is_Color_Too_Important = false;			
			
			#ifdef USE_COLOR_METRIC
			if(this->IsColored)
				Is_Color_Too_Important = this->Error_Projected_Surface(_pMesh, h, Component_ID, Mean_color, Mean_area);
			#endif

			// remove the front vertex if its removal does not viloate the manifold property and some metrics			
			bool Check_all_condition = false;

			if ((!Manifold_condition) && (!Geometric_metric_condition) && (!Normal_flipping_condition) && (!Is_Color_Too_Important))
				Check_all_condition = true;			

			if (Check_all_condition) // All conditions are good. -> Remove the center vertex.
			{

				//increase number of vertices and symbols
				Number_vertices++;
				Number_symbol++;				
				
				Halfedge_handle g = h;				

				Point3d Coordinates_removed_vertex = g->next()->vertex()->point();
				Point_Int CRV = Change_Real_Int(Coordinates_removed_vertex, Component_ID);
				
				// encoding of vertex color
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Removed_vertex_color;
					Removed_vertex_color = Get_Vertex_Color(g->next());					

					Color_Unit Average_color;

					Average_color = Get_Average_Vertex_Color_Lee(g, valence);

					// Color difference from average color of neighbors					
					Color_Unit Color_diff = Removed_vertex_color - Average_color;					
					
					//#ifdef PREDICTION_METHOD
					this->InterVertexColor.push_front(Color_diff);
					//#endif
					
				}								

				// Enter symbol 'VALENCE CODE' into the list of symbols
				this->InterConnectivity.push_front(valence - 3);				

				// remove the front vertex
				_pMesh.erase_center_vertex(g->next()); 

				g = h;
				g->facet()->Facet_Flag = TEMP_FLAG;				
				
				// Flags and fill queue
				for (int j = 0; j < (valence-1); j++)
				{
					g = g->next();
					g->vertex()->Vertex_Flag = CONQUERED;
					
					g->opposite()->vertex()->Vertex_Flag = CONQUERED;
					
					if (!g->is_border_edge())
						Halfedges.push(g->opposite());
				}
				g = h;							
												
				Retriangulation(_pMesh, g, valence, 0, Component_ID);				
				
				this->InterGeometry.push_front(CRV);
			}
			
			// Conditions are not satisfied -> NULL CODE
			else 
			{
				// Enter symbol 'NULL PATCH' into the list of symbols
				this->InterConnectivity.push_front(4);

				Number_symbol++;
				h->facet()->Facet_Flag = CONQUERED;
				h->next()->vertex()->Vertex_Flag = CONQUERED;
				
				if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());				
				
				if (!h->prev()->is_border_edge())					
					Halfedges.push(h->prev()->opposite());				
			}
		}
		// border edge.
		else if ( (h->next()->vertex()->Vertex_Flag == FREE) && 
				  (valence >= 3) && 
				  (valence <= 4) && 
				  (Is_Border_Vertex(h->next())) )
		{
			/*****    conditions of vertex removal (based on area) will be added    *****/			

			type = Find_Type(h, valence);
			
			Halfedge_handle g = h;
			
			bool Check_border_structure = true;
			
			Halfedge_handle First_border_edge;
			Halfedge_handle Standard_edge;					
			
			// To find first border edge.
			Halfedge_around_vertex_circulator Hvc = g->next()->vertex()->vertex_begin();
			Halfedge_around_vertex_circulator Hvc_end = Hvc;
			
			int Number_neighboring_border_vertices = 0;
			CGAL_For_all(Hvc, Hvc_end)
			{				
				if (Is_Border_Vertex(Hvc->opposite()))
					Number_neighboring_border_vertices++;					
				
				if (Hvc->is_border())
					First_border_edge = Hvc;
			}					
			
			if (Number_neighboring_border_vertices > 2)
				Check_border_structure = false;	
				
			if (Is_Manifold_Property_Violated(g, type, valence))
				Check_border_structure = false;


			// Added by Adrien Maglo.
			// Test if we are not in the case of a mesh pimple.
			if (Check_border_structure)
			{
				// Get the two patch border vertices.
				Vertex_handle vh1 = First_border_edge->next()->vertex();
				Vertex_handle vh2 = First_border_edge->opposite()->vertex();
				Halfedge_around_vertex_circulator vh_it = vh1->vertex_begin();
				Halfedge_around_vertex_circulator vh_it_end = vh_it;

				// Test if the two patch border vertices are not connected by an edge.
				CGAL_For_all(vh_it, vh_it_end)
				{
					if (vh_it->opposite()->vertex() == vh2)
					{
						Check_border_structure = false;
						break;
					}
				}
			}			


			// if the hole is a triangle(Count_border_edges == 3), we should not remove the vertex.
			int Count_border_edges = 1;
			Standard_edge = First_border_edge;
			First_border_edge = First_border_edge->next();
			
			for (; First_border_edge != Standard_edge; First_border_edge = First_border_edge->next())
			{
				Count_border_edges++;				
				if (Count_border_edges > 5)
					break;
			}
			
			if (Count_border_edges <= 3)
				Check_border_structure = false;

			// Exception which causes a manifold violation
			if (Count_border_edges == 4)
				if (Standard_edge->prev()->opposite()->vertex_degree() == 2)
					Check_border_structure = false;		

			// Structure of triangles is good -> decimation.
			if (Check_border_structure) 
			{	
				Number_vertices++;
				Number_symbol++;
				
				// Connectivity Code = 5 for border valence 3 vertex
				//         "         = 6 for border valence 4 vertex
				this->InterConnectivity.push_front(valence + 2); 
			
				Point3d Real_vertex_position = g->next()->vertex()->point();
				Point_Int Vertex_position = Change_Real_Int(Real_vertex_position, Component_ID);				
				
				Color_Unit Removed_vertex_color;
				//int Vertex_color_index = -1;
				
				if ((this->IsColored) && (!this->IsOneColor))
				{	
					Removed_vertex_color = Get_Vertex_Color(h->next());
				}				
				

				vector<Halfedge_handle> Border_edges;
				int Number_jump = 0;

				while(g->next()->is_border_edge() == false)
				{
					g = g->next()->opposite()->next();
					Number_jump++;
				}
				
				Border_edges.push_back(g->opposite());
				for (int i = 0; i < (valence - 2); i++)
				{
					g = g->prev()->opposite()->prev();	
					Border_edges.push_back(g->opposite());													
				}		

				for (int i = 0; i < (valence - 1); i++)
				{
					Halfedge_handle Temp = Border_edges[i];
					Temp = Temp->opposite();
					
					CGAL_assertion(!Temp->is_border());
					_pMesh.erase_facet(Temp);
				}
								
				if (valence == 3)
				{
					Halfedge_handle Retriangulation_edge = Border_edges[valence - 2]->opposite()->prev();
					
					// One triangle has to be created to smooth the mesh boundary					
					_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());					
										
					Halfedge_handle Input_gate = Border_edges[Number_jump]->opposite();					

					// its front face is tagged CONQUERED
					Input_gate->facet()->Facet_Flag = CONQUERED;					
					Input_gate->next()->vertex()->Vertex_Flag = CONQUERED;

					if ((type == 1) || (type == 2) || (type == 4))
					{
						if (Input_gate->next()->vertex()->Vertex_Sign == NOSIGN)
							Input_gate->next()->vertex()->Vertex_Sign = PLUS;
					}
					else if (type == 3)
					{
						if (Input_gate->next()->vertex()->Vertex_Sign == NOSIGN)
							Input_gate->next()->vertex()->Vertex_Sign = MINUS;
					}			
					
					if ((this->IsColored) && (!this->IsOneColor))
					{
						Color_Unit Average_color;

						Average_color = Get_Average_Vertex_Color_After_Removal(g, valence); // g is the most left placed edge
						Color_Unit Color_diff = Removed_vertex_color - Average_color;						

						//#ifdef PREDICTION_METHOD
						this->InterVertexColor.push_front(Color_diff);
						//#endif
					}

					Halfedge_handle g = Input_gate;

					this->InterGeometry.push_front(Vertex_position);

					if (Number_jump == 0)
						Halfedges.push(Border_edges[1]);
					else
						Halfedges.push(Border_edges[0]);
				}		
				
				else if (valence == 4)
				{		
					Halfedge_handle Retriangulation_edge = Border_edges[valence - 2]->opposite()->prev();					
					
					if (   ( (Number_jump == 0) && ((type == 5) || (type == 8)) ) ||
						   ( (Number_jump == 1) && ((type == 6) || (type == 7)) ) ||
						   ( (Number_jump == 2) && ((type == 5) || (type == 8)) )  )

					{						
						_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[1]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag = CONQUERED;						

						_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());
						Border_edges[0]->opposite()->facet()->Facet_Flag = CONQUERED;						
					}

					else if (   ( (Number_jump == 0) && ((type == 6) || (type == 7)) ) ||
								( (Number_jump == 1) && ((type == 5) || (type == 8)) ) ||
								( (Number_jump == 2) && ((type == 6) || (type == 7)) )  )
					{												
						_pMesh.add_facet_to_border(Border_edges[2]->opposite(), Border_edges[0]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag = CONQUERED;						
						Halfedge_handle Temp_border = Border_edges[2]->opposite()->next();
						
						_pMesh.add_facet_to_border(Retriangulation_edge, Temp_border);
						Border_edges[2]->opposite()->facet()->Facet_Flag = CONQUERED;						
					}					

					Halfedge_handle Input_gate = Border_edges[Number_jump]->opposite();

					
					// Vertex Signs
					if ((type == 5) || (type == 8))
					{
						Halfedge_handle g = Input_gate;
						g = g->prev()->opposite()->next();
						if (g->vertex()->Vertex_Sign == NOSIGN)
							g->vertex()->Vertex_Sign = MINUS;
						if (g->opposite()->vertex()->Vertex_Sign == NOSIGN)
							g->opposite()->vertex()->Vertex_Sign = PLUS;
					}
					else if ((type == 6) || (type == 7))
					{
						Halfedge_handle g = Input_gate;
						g = g->next()->opposite()->prev();
						if (g->vertex()->Vertex_Sign == NOSIGN)
							g->vertex()->Vertex_Sign = PLUS;
						if (g->opposite()->vertex()->Vertex_Sign == NOSIGN)
							g->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					
					//Vertex Flags + Fill Halfedges queue.
					for (int i = 0; i < valence - 1; i++)
					{
						Border_edges[i]->vertex()->Vertex_Flag = CONQUERED;
						Border_edges[i]->opposite()->vertex()->Vertex_Flag = CONQUERED;

						if (i != (unsigned)Number_jump)
							Halfedges.push(Border_edges[i]);
					}
					
					if ((this->IsColored) && (!this->IsOneColor))
					{
						Color_Unit Average_color;

						Color_Unit Color_0 = Get_Vertex_Color(Border_edges[0]->opposite());
						Color_Unit Color_1 = Get_Vertex_Color(Border_edges[0]);
						Color_Unit Color_2 = Get_Vertex_Color(Border_edges[1]);
						Color_Unit Color_3 = Get_Vertex_Color(Border_edges[2]);
						Average_color.c0 = (int)(Color_0.c0 + Color_1.c0 + Color_2.c0 + Color_3.c0)/4.0;
						Average_color.c1 = (int)(Color_0.c1 + Color_1.c1 + Color_2.c1 + Color_3.c1)/4.0;
						Average_color.c2 = (int)(Color_0.c2 + Color_1.c2 + Color_2.c2 + Color_3.c2)/4.0;
						
						Color_Unit Color_diff = Removed_vertex_color - Average_color;						
						
						//#ifdef PREDICTION_METHOD
						this->InterVertexColor.push_front(Color_diff);
						//#endif
					}					
					this->InterGeometry.push_front(Vertex_position);					
				}					
			}

			else // Border vertex can not be removed
			{
				Number_symbol++;
				this->InterConnectivity.push_front(4);
				
				h->facet()->Facet_Flag = CONQUERED;				
				h->next()->vertex()->Vertex_Flag = CONQUERED;

				if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign == NOSIGN)
						h->next()->vertex()->Vertex_Sign = PLUS;
				}				
				
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());

				if (!h->prev()->is_border_edge())
					Halfedges.push(h->prev()->opposite());			
			}	
		}		

		// if its front face sholud be labelled as NULL PATCH		
		else
		{
			// Enter symbol 'NULL PATCH' into the list of symbols
			this->InterConnectivity.push_front(4);
			Number_symbol++;

			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;///////////////////////////

			if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}

			// two other output gates are pused into the fifo queue.
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)					
				Halfedges.push(h->prev()->opposite());
		}
	}	
	if ((this->IsColored) && (!this->IsOneColor))
	{
		//#ifdef PREDICTION_METHOD
		while(!this->InterVertexColor.empty())
		{
			Color_Unit Col = this->InterVertexColor.front();
			this->InterVertexColor.pop_front();
			this->VertexColor[Component_ID].push_front(Col);
		}
		//#endif
	}
	// InterConnectivity is the intermediate connectivity symbol container.
	// We put all symbols in the main container which is this->Connectivity
	while(!InterConnectivity.empty())
	{
		int Symbol = InterConnectivity.front();
		InterConnectivity.pop_front();
		this->Connectivity[Component_ID].push_front(Symbol);
	}
	
	// same operation than connectivity.
	while(!InterGeometry.empty())
	{
		Point_Int Geo = InterGeometry.front();
		InterGeometry.pop_front();
		this->Geometry[Component_ID].push_front(Geo);
	}
	
	this->NumberVertices[Component_ID].push_front(Number_vertices);
	this->NumberSymbol[Component_ID].push_front(Number_symbol);
	
	// if this decimation didn't remove any vertex, we should remove all connectivity symbols. for that we store number of symbols.
	this->DumpSymbolDecimation = Number_symbol;
	return Number_vertices;
}


// Description : Regulation conquest
int Compression_Valence_Component::JCW_Regulation_For_Segmentation(Polyhedron  & _pMesh,
																   const bool    Normal_flipping,
								                                   const bool    Use_metric,
											                       const float & Metric_thread,
											                       const bool    Use_forget_metric,
											                       const int   & Forget_value,
											                       const int   & Component_ID)
{    
	double Max_color, Mean_color;
	int Temp_NV = 0;
	int Number_facets;
	double Mean_area;
	if(this->IsColored)
	{
		this->Calculate_Edge_Color_Difference(_pMesh, Component_ID, Max_color, Mean_color, Temp_NV);
		this->Recalculate_Component_Area(_pMesh, Component_ID,Number_facets);
		Mean_area = (double)this->ComponentArea[Component_ID] / (double)Number_facets;
	}

	// Initialization
	Init(_pMesh);	
	
	// Number of removed vertices and number of connectivity symbols
	int Number_vertices = 0;
	int Number_symbol = 0;

	Halfedge_iterator hi = _pMesh.halfedges_begin();
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))
		hi++;	

	Halfedge_handle First_halfedge = &(*(hi));	

	First_halfedge->vertex()->Vertex_Flag = CONQUERED;
	First_halfedge->opposite()->vertex()->Vertex_Flag = CONQUERED;
		
	std::queue<Halfedge_handle> Halfedges;
	Halfedges.push(First_halfedge);

	Halfedge_handle h;

	while(!Halfedges.empty())
	{		
		h = Halfedges.front();
		Halfedges.pop();
		
		size_t valence = h->next()->vertex_degree();		

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;

		else if ((h->next()->vertex()->Vertex_Flag == FREE) && (valence == 3) && (Is_Border_Vertex(h->next()) == false)) // if valence is 3, remove the front vertex.
		{
			Halfedge_handle g = h;
			int type = 1; // ant type of valence 3
           
			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Geometric_metric_condition = false;
			if (Use_metric == true)
			{
				if (Use_forget_metric == true)
				{
					if (_pMesh.size_of_vertices() > (unsigned)Forget_value)
						Geometric_metric_condition = false;
					else
						Geometric_metric_condition = Is_Geometric_Metric_Violated(g, type, valence, Metric_thread);
				}

				else
					Geometric_metric_condition = Is_Geometric_Metric_Violated(g, type, valence, Metric_thread);
			}	
			bool Is_Color_Too_Important = false;	

			#ifdef USE_COLOR_METRIC
			if(this->IsColored)
				Is_Color_Too_Important = this->Error_Projected_Surface(_pMesh, g, Component_ID, Mean_color, Mean_area);
			#endif

			// remove the front vertex if its removal does not viloate the manifold property and some metrics			
			bool Check_all_condition = false;

			if ((!Geometric_metric_condition) && (!Is_Color_Too_Important))
				Check_all_condition = true;			
			
			// All conditions are goods. Remove the front vertex
			if (Check_all_condition)
			{
				Number_vertices++;
				Number_symbol++;
				
				// Insertion of a symbol 3 to the list
				this->InterConnectivity.push_front(0);
				
				Point3d Geo_info = g->next()->vertex()->point();
				Point_Int Geo = Change_Real_Int(Geo_info, Component_ID);

				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Removed_vertex_color;
					
					Removed_vertex_color = Get_Vertex_Color(g->next());

					Color_Unit Average_color;
					Average_color = Get_Average_Vertex_Color_Lee(g, valence);
					
					
					Color_Unit Color_diff = Removed_vertex_color - Average_color;					
					
					//#ifdef PREDICTION_METHOD
					this->InterVertexColor.push_front(Color_diff);
					//#endif
				}
				
				
				this->InterGeometry.push_front(Geo);
								
				g = h->next();
				g->vertex()->Vertex_Flag = TO_BE_REMOVED;
				g->facet()->Facet_Flag = TO_BE_REMOVED;

				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				
				if (!g->prev()->is_border_edge())
				{
					Halfedge_handle h1 = g->prev()->opposite();
					h1->facet()->Facet_Flag = CONQUERED;
					h1->next()->vertex()->Vertex_Flag = CONQUERED;
					if (!h1->next()->is_border_edge())
						Halfedges.push(h1->next()->opposite());
					if (!h1->prev()->is_border_edge())
						Halfedges.push(h1->prev()->opposite());					
				}
				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				if (!g->prev()->is_border_edge())
				{
					Halfedge_handle h2 = g->prev()->opposite();
					h2->facet()->Facet_Flag = CONQUERED;
					h2->next()->vertex()->Vertex_Flag = CONQUERED;
					if (!h2->next()->is_border_edge())
						Halfedges.push(h2->next()->opposite());
					if (!h2->prev()->is_border_edge())
						Halfedges.push(h2->prev()->opposite());
				}
			}

			else
			{
				this->InterConnectivity.push_front(1);
				Number_symbol++;
				h->facet()->Facet_Flag = CONQUERED;
				h->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());
				if (!h->prev()->is_border_edge())					
					Halfedges.push(h->prev()->opposite());				
			}
		}
		else  // NULL triangle
		{
			this->InterConnectivity.push_front(1);
			Number_symbol++;
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())					
				Halfedges.push(h->prev()->opposite());			
		}
	}
	
	
	// Removal of all vertices with TO_BE_REMOVED flag
	Vertex_iterator pVertex = NULL;
	for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end();)
	{
		Vertex_handle vh = pVertex;
		pVertex++;
		if (vh->Vertex_Flag == TO_BE_REMOVED)
		{
			Halfedge_handle temp = vh->halfedge();
			temp = _pMesh.erase_center_vertex(temp);
			temp->facet()->Component_Number = Component_ID;
		}
	}	

	if ((this->IsColored) && (!this->IsOneColor))
	{
		//#ifdef PREDICTION_METHOD
		while(!this->InterVertexColor.empty())
		{
			Color_Unit Col = this->InterVertexColor.front();
			this->InterVertexColor.pop_front();
			this->VertexColor[Component_ID].push_front(Col);
		}
		//#endif
	}
	while(!InterConnectivity.empty())
	{
		int Symbol = this->InterConnectivity.front();
		this->InterConnectivity.pop_front();
		this->Connectivity[Component_ID].push_front(Symbol);
	}
	
	while(!InterGeometry.empty())
	{
		Point_Int Geo = InterGeometry.front();
		InterGeometry.pop_front();
		this->Geometry[Component_ID].push_front(Geo);
	}
	
	this->NumberSymbol[Component_ID].push_front(Number_symbol);
	this->NumberVertices[Component_ID].push_front(Number_vertices);

	this->DumpSymbolRegulation = Number_symbol;
	
	

	return Number_vertices;
}

Point3d Compression_Valence_Component::JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence, const int & Direction)
{

	Halfedge_handle g = h;
	double x = 0., y = 0., z = 0.;

	for (int i = 0;i < valence; i++)
	{
		Point3d Pt = g->vertex()->point();
		double Spheric[3];
		this->Convert_To_Spherical(Pt, Spheric);
		if(Direction == 0)
			Spheric[0] -= this->m_Dist * this->m_EmbeddingStrength;
		else
			Spheric[0] += this->m_Dist * this->m_EmbeddingStrength;			
		
		g = g->next();

		double Cart[3];
	    this->Convert_To_Cartesian(Spheric, Cart);
		
		int Qx = (int)(ceil((Cart[0] - (double)this->xmin[0]) / (double)this->Quantization_Step[0])) - 1;
		if (Qx == -1)
			Qx = 0;
		int Qy = (int)(ceil((Cart[1] - (double)this->ymin[0]) / (double)this->Quantization_Step[0])) - 1;
		if (Qy == -1)
			Qy = 0;
		int Qz = (int)(ceil((Cart[2] - (double)this->zmin[0]) / (double)this->Quantization_Step[0])) - 1;
		if (Qz == -1)
			Qz = 0;		
		
		x += this->xmin[0] + (Qx + 0.5) * this->Quantization_Step[0];
		y += this->ymin[0] + (Qy + 0.5) * this->Quantization_Step[0];
		z += this->zmin[0] + (Qz + 0.5) * this->Quantization_Step[0];		
	}
	
	x = x / (double)valence;
	y = y / (double)valence;
	z = z / (double)valence;
	
	return Point3d(x,y,z);
}

Point3d Compression_Valence_Component::JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & valence, const int & Direction)
{

	Halfedge_handle g = h;
	double x = 0., y = 0., z = 0.;

	Halfedge_around_vertex_circulator Hvc = g->next()->vertex()->vertex_begin();
	Halfedge_around_vertex_circulator Hvc_end = Hvc;

	CGAL_For_all(Hvc, Hvc_end)
	{
		Point3d Pt = Hvc->opposite()->vertex()->point();
		double Spheric[3];
		this->Convert_To_Spherical(Pt, Spheric);
		if(Direction == 0)
			Spheric[0] -= this->m_Dist * this->m_EmbeddingStrength;
		else
			Spheric[0] += this->m_Dist * this->m_EmbeddingStrength;		

		double Cart[3];
	    this->Convert_To_Cartesian(Spheric, Cart);
		
		int Qx = (int)(ceil((Cart[0] - (double)this->xmin[0]) / (double)this->Quantization_Step[0])) - 1;
		if (Qx == -1)
			Qx = 0;
		int Qy = (int)(ceil((Cart[1] - (double)this->ymin[0]) / (double)this->Quantization_Step[0])) - 1;
		if (Qy == -1)
			Qy = 0;
		int Qz = (int)(ceil((Cart[2] - (double)this->zmin[0]) / (double)this->Quantization_Step[0])) - 1;
		if (Qz == -1)
			Qz = 0;		
		
		x += this->xmin[0] + (Qx + 0.5) * this->Quantization_Step[0];
		y += this->ymin[0] + (Qy + 0.5) * this->Quantization_Step[0];
		z += this->zmin[0] + (Qz + 0.5) * this->Quantization_Step[0];		
	}
	
	x = x / (double)valence;
	y = y / (double)valence;
	z = z / (double)valence;
	
	return Point3d(x,y,z);
}

// Description : Decoding function of decimation conquest
void Compression_Valence_Component::JCW_Un_Decimation_Conquest(Polyhedron       & _pMesh, 
															   Arithmetic_Codec & Decoder,
														       const int        & Component_ID)
{
	Init(_pMesh);	

	this->m_N_treated_vertices.clear();
	this->m_Rad_decision.clear();

	for(int i = 0; i < this->m_NumberRegion; i++)
	{		
		this->m_N_treated_vertices.push_back(0);
		this->m_Rad_decision.push_back(0.);		
	}
	
	int Number_connectivity_symbols;
	if (this->IsClosed[Component_ID])
		Number_connectivity_symbols = 5;
	else
		Number_connectivity_symbols = 7;

	Adaptive_Data_Model Connectivity(Number_connectivity_symbols);	
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;	

	unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];

	int Alpha_range = Decoder.get_bits(Qbit+1);
	int Alpha_offset = Decoder.get_bits(Qbit+1);
	int Gamma_range = Decoder.get_bits(Qbit+1);
	int Gamma_offset = Decoder.get_bits(Qbit+1);

	if (this->Smallest_Alpha < 0)
		Alpha_offset = Alpha_offset + this->Smallest_Alpha;
	
	if (this->Smallest_Gamma < 0)
		Gamma_offset = Gamma_offset + this->Smallest_Gamma;

	bool check_alpha = false;
	bool check_gamma = false;

	if ((Alpha_range == 0) || (Alpha_range == 1))
	{
		check_alpha = true;
		Alpha_range = 2;
	}

	if ((Gamma_range == 0) || (Gamma_range == 1))
	{
		check_gamma = true;
		Gamma_range = 2;
	}
	

	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);


	Adaptive_Data_Model alpha(Alpha_range);
	Adaptive_Data_Model gamma(Gamma_range);
	
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID+1))
		hi++;	

	// Two vertices of seed edges are flaged CONQUERED
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;	
		
	// These vertices are also flages with sign flags for retriangulation
	hi->vertex()->Vertex_Sign = PLUS;
	hi->opposite()->vertex()->Vertex_Sign = MINUS;
	
	// Two vertices of seed edges are flaged CONQUERED
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		unsigned int valence = 0,type = 0; // define type of retriangulation

		if (h->is_border() == true)
			continue;	
		
		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED
		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;
			
		valence = Decoder.decode(Connectivity) + 3;		
				
		// if its front vertex is free
		if ((valence >= 3) && (valence <= 6))
		{
			type = Find_Type(h, valence);
			
			// remove the front vertex if its removal does not viloate the manifold property
			Halfedge_handle g = h;
			
			Vector normal = Normal_Patch(g, valence);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(g, normal, T2);			
			
			if ((T1 == CGAL::NULL_VECTOR) || (T2 == CGAL::NULL_VECTOR) || (normal == CGAL::NULL_VECTOR))
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}		

			// remove edges to re_find polygon and (check. and attribute sign flag.)
			bool Check_Validity = Remove_Edges(_pMesh, g, type);
			Check_Validity = false;

			//Check_Validity = false;// * * * * * * * * * * *////
			if (Check_Validity == false)
			{
				g = h;
				//Halfedge_handle pass = h;			
				
				Point_Int Frenet;
				if (check_alpha == false)
				{
					Frenet.x = Decoder.decode(alpha);
					Frenet.y = Decoder.decode(alpha);
				}
				else if (check_alpha == true)
				{
					Frenet.x = 0;
					Frenet.y = 0;
				}
				if (check_gamma == false)
					Frenet.z = Decoder.decode(gamma);
				else if (check_gamma == true)
					Frenet.z = 0;

				Frenet.x -= Alpha_offset;
				Frenet.y -= Alpha_offset;
				Frenet.z -= Gamma_offset;

				Point_Int Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);

				// Calculate the position of barycenter.		
				Point3d Barycenter_real;
				Point3d Center_vertex;

				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				for(int i = 0; i < (int)valence; i++)
				{
					int N1 = reg->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}
					reg = reg->next();
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
				{
					Selected_region = T_possible_bin[0];
				}
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				int Region = Selected_region;

				if(this->m_N_remained_vertices[Region] < LIMIT_NUMBER)
				{
					Barycenter_real = Barycenter_Patch_After_Removal(h, valence);
					Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);
					Point_Int Center = BC + Diff;						
					Center_vertex = this->Change_Int_Real(Center, Component_ID);					
				}

				else if(this->m_N_treated_vertices[Region] < MINIMUM_PREDICTION_NUMBER)
				{
					Barycenter_real = Barycenter_Patch_After_Removal(h, valence);

					Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);
					Point_Int Center = BC + Diff;						
					Center_vertex = this->Change_Int_Real(Center, Component_ID);

					double BC_S[3];
					double CRV_S[3];

					this->Convert_To_Spherical(Center_vertex, CRV_S);
					this->Convert_To_Spherical(Barycenter_real, BC_S);
					
					double Diff_R = CRV_S[0] - BC_S[0];
					this->m_Rad_decision[Region] += Diff_R;
				}
				else
				{
					int Direction;

					if(this->m_Rad_decision[Region] > 0.0)
						Direction = 1;
					else
						Direction = 0;					

					Barycenter_real = this->JCW_Barycenter_Patch_After_Removal(h, valence, Direction);

					Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);

					Point_Int Center = BC + Diff;						
					Center_vertex = this->Change_Int_Real(Center, Component_ID);
				}
				this->m_N_treated_vertices[Region]++;				
				
				g = _pMesh.create_center_vertex(g);
				g->vertex()->point() = Center_vertex;

				g->vertex()->Region_Number = Selected_region;
				
				//if(Selected_region != -1)
				//	this->m_Number_Vertices_Per_Regions[Selected_region]++;
				
				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;

				g->vertex()->Vertex_Flag = CONQUERED;

				g = h;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				g->vertex()->Vertex_Flag = CONQUERED;
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;

				for (unsigned int k = 0; k < (valence - 1); k++)
				{
					g = g->next()->opposite()->next();
					g->facet()->Facet_Flag = TO_BE_REMOVED;
					g->vertex()->Vertex_Flag = CONQUERED;
					g->opposite()->vertex()->Vertex_Flag = CONQUERED;
					if (g->is_border_edge() == false)
						Halfedges.push(g->opposite());
				}
				
				g->next()->vertex()->Seed_Edge = -1;

				if ((this->IsColored) && (!this->IsOneColor))
				{				
					//#ifdef PREDICTION_METHOD
					g = h;				

					Color_Unit Predicted_color;

					Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
					
					Color_Unit Color_difference;
					Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
					Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
					Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;								
					
					Color_Unit CV = Predicted_color + Color_difference;					

					g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * Color_step;
					LAB[1] = this->C1_Min + CV.c1 * Color_step;
					LAB[2] = this->C2_Min + CV.c2 * Color_step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);
				
					//#endif
				}
				if ((this->IsColored) && (this->IsOneColor))
				{
					g = h->next();
					g->vertex()->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
				}

			}
		}
		// In case of border edge.
		else if ((valence == 8) || (valence == 9))
		{
			Halfedge_handle g = h;	

			type = Find_Type(g, valence - 5);				
			
			Vector normal = Normal_Patch(g, valence - 5);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(g, normal, T2);			
			
			//if ((T1 == CGAL::NULL_VECTOR) || (T2 == CGAL::NULL_VECTOR) || (normal == CGAL::NULL_VECTOR))
			{							
				T1 =     Vector(1.,0.,0.);
				T2 =     Vector(0.,1.,0.);	
				normal = Vector(0.,0.,1.);				
			}						
			
			Point_Int Frenet;
			if (check_alpha == false)
			{
				Frenet.x = Decoder.decode(alpha);
				Frenet.y = Decoder.decode(alpha);
			}
			else if (check_alpha == true)
			{
				Frenet.x = 0;
				Frenet.y = 0;
			}

			if (check_gamma == false)
				Frenet.z = Decoder.decode(gamma);
			else if (check_gamma == true)
				Frenet.z = 0;

			Frenet.x -= Alpha_offset;
			Frenet.y -= Alpha_offset;
			Frenet.z -= Gamma_offset;

			Point_Int Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);			
			
			//#ifdef PREDICTION_METHOD
			Color_Unit Predicted_color;
			if ((this->IsColored) && (!this->IsOneColor))
			{
				Predicted_color.c0 = Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
				Predicted_color.c1 = Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
				Predicted_color.c2 = Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
			}
			//#endif			
			// border edge with valence == 3
			if (valence == 8)
			{
				//#ifdef PREDICTION_METHOD
				Color_Unit Average_color;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Average_color = Get_Average_Vertex_Color_After_Removal(g, 3);
				}
				//#endif

				vector<Halfedge_handle> Border_edges;
				int Number_jump = 0;
				
				if (g->next()->is_border_edge())
					Number_jump = 0;
				if (g->prev()->is_border_edge())
				{
					Number_jump = 1;
					g = g->next();
				}
				
				Border_edges.push_back(g->opposite());
				
				g = g->prev();
				Border_edges.push_back(g->opposite());

				_pMesh.erase_facet(g);

				Halfedge_handle Prev_edge = Border_edges[1]->opposite()->prev();
					
				// g points the new vertex
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[1]->opposite());				
				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());

				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				Halfedge_around_vertex_circulator Hvc = h->next()->vertex()->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;
				
				CGAL_For_all(Hvc, Hvc_end)
				{
					int N1 = Hvc->opposite()->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}					
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
					Selected_region = T_possible_bin[0];
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				//int Region = Selected_region;				
				g->vertex()->Region_Number = Selected_region;

				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;
				g->vertex()->Vertex_Flag = CONQUERED;

				Point3d Barycenter_real;
				Point3d Center_vertex;
							
				Barycenter_real = Barycenter_Patch_Before_Removal(h);
				Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);
				Point_Int Center = BC + Diff;						
				Center_vertex = this->Change_Int_Real(Center, Component_ID);			

				//this->m_N_treated_vertices[Region]++;

				g->vertex()->point() = Center_vertex;
				

				//#ifdef PREDICTION_METHOD
				Color_Unit CV;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					CV = Average_color + Predicted_color;				

					g->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->vertex()->color(RGB[0], RGB[1], RGB[2]);
				}
				//#endif
				
				if ((this->IsColored) && (this->IsOneColor))
				{
					g->vertex()->color(this->OnlyColor[0],this->OnlyColor[1],this->OnlyColor[2]);
				}
				
				Halfedge_handle Tag_handle;

				if (Number_jump == 0)
					Tag_handle = Border_edges[1];
				if (Number_jump == 1)
					Tag_handle = Border_edges[0]->opposite();				
				Tag_handle->vertex()->Vertex_Flag = CONQUERED;

				if ((type == 1) || (type == 2) || (type == 4))
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = PLUS;
				}
				else if (type == 3)
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = MINUS;
				}
				for (int i=0; i <2; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}				
			}
			// border edge with valence == 4
			else if (valence == 9)
			{	
				g = h;
				int Number_jump = -1;				
				vector<Halfedge_handle> Border_edges;

				if ((type == 5) || (type == 8))
				{	
					// jump == 0;
					if (g->next()->is_border_edge())
					{
					
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());
						Border_edges.push_back(g->prev()->opposite()->next()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					// jump == 1;
					else if (g->prev()->opposite()->next()->is_border_edge())
					{						
						Number_jump = 1;
						
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[0]->opposite());
					}
					
					// jump == 2;
					else
					{						
						Number_jump = 2;

						Border_edges.push_back(g->prev()->opposite()->next()->opposite());
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
				}
				else
				{
					if (g->prev()->is_border_edge())
					{
						Number_jump = 2;

						Border_edges.push_back(g->next()->opposite()->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					else if (g->next()->opposite()->prev()->is_border_edge())
					{
						Number_jump = 1;

						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
					else
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

				}
				
				g = h;			
				
				// to create the new facets
				Halfedge_handle Prev_edge = Border_edges[2]->opposite()->prev();
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[2]->opposite());

				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());
					
				//#ifdef PREDICTION_METHOD
				Color_Unit CV;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Average_color;

					Color_Unit Color_0 = Get_Vertex_Color(Border_edges[0]->opposite());
					Color_Unit Color_1 = Get_Vertex_Color(Border_edges[0]);
					Color_Unit Color_2 = Get_Vertex_Color(Border_edges[1]);
					Color_Unit Color_3 = Get_Vertex_Color(Border_edges[2]);

					Average_color.c0 = (int)(Color_0.c0 + Color_1.c0 + Color_2.c0 + Color_3.c0)/4.0;
					Average_color.c1 = (int)(Color_0.c1 + Color_1.c1 + Color_2.c1 + Color_3.c1)/4.0;
					Average_color.c2 = (int)(Color_0.c2 + Color_1.c2 + Color_2.c2 + Color_3.c2)/4.0;					
					
					CV = Average_color + Predicted_color;					
			
					g->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->vertex()->color(RGB[0], RGB[1], RGB[2]);
				}
				//#endif
				if ((this->IsColored) && (this->IsOneColor))
				{
					g->vertex()->color(this->OnlyColor[0],this->OnlyColor[1],this->OnlyColor[2]);
				}
				

				////////////////////

				// Assign the region number to inserted vertex
				Halfedge_handle reg = g;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				Halfedge_around_vertex_circulator Hvc = reg->vertex()->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;
				
				CGAL_For_all(Hvc, Hvc_end)
				{
					int N1 = Hvc->opposite()->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}					
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
					Selected_region = T_possible_bin[0];
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				//int Region = Selected_region;

				Point3d Barycenter_real =  Barycenter_Patch_Before_Removal(h);
				Point3d Center_vertex;				
				
				Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);
				Point_Int Center = BC + Diff;						
				Center_vertex = this->Change_Int_Real(Center, Component_ID);				

				//this->m_N_treated_vertices[Region]++;
				
				g->vertex()->Region_Number = Selected_region;			
				g->vertex()->point() = Center_vertex;
				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;
				g->vertex()->Vertex_Flag = CONQUERED;

				////////////////////
				//vertex_tag
				if (Number_jump == 0)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}

				else if (Number_jump == 1)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}

				}
				else // jump == 2
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}
				
				for (int i=0; i < 3; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					
					Border_edges[i]->vertex()->Vertex_Flag = CONQUERED;
					Border_edges[i]->opposite()->vertex()->Vertex_Flag = CONQUERED;
					
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}					
			}
		}

			
		
		
		//  the symbol == N
		else if (valence == 7)
		{
			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;

			if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)
				Halfedges.push(h->prev()->opposite());	
		}
	}	

}

// Description : Decoding of the regulation conquest
void Compression_Valence_Component::JCW_Un_Regulation(Polyhedron &_pMesh, Arithmetic_Codec & Decoder, const int & Component_ID)
{	
	Init(_pMesh);

	this->m_N_remained_vertices.clear();
	this->m_N_treated_vertices.clear();
	this->m_Rad_decision.clear();

	for(int i = 0; i < this->m_NumberRegion; i++)
	{		
		this->m_N_treated_vertices.push_back(0);
		this->m_Rad_decision.push_back(0.);
		this->m_N_remained_vertices.push_back(0);
	}


	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		int Reg = pVert->Region_Number;
		this->m_N_remained_vertices[Reg]++;
	}

	Adaptive_Data_Model Connectivity(2);
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;

	unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];
	
	int Alpha_range = Decoder.get_bits(Qbit+1);
	int Alpha_offset = Decoder.get_bits(Qbit+1);
	int Gamma_range = Decoder.get_bits(Qbit+1);
	int Gamma_offset = Decoder.get_bits(Qbit+1);

	if (this->Smallest_Alpha < 0)
		Alpha_offset = Alpha_offset + this->Smallest_Alpha;

	if (this->Smallest_Gamma < 0)
		Gamma_offset = Gamma_offset + this->Smallest_Gamma;

	bool check_alpha = false;
	bool check_gamma = false;

	if ((Alpha_range == 0) || (Alpha_range == 1))
	{
		check_alpha = true;
		Alpha_range = 2;
	}

	if ((Gamma_range == 0) || (Gamma_range == 1))
	{
		check_gamma = true;
		Gamma_range = 2;
	}

	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);

	Adaptive_Data_Model alpha(Alpha_range);
	Adaptive_Data_Model gamma(Gamma_range);
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))	
		hi++;	
	
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;
	
	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))// already visited.
			continue;
		
		// read connectivity information
		int valence = Decoder.decode(Connectivity) + 3;
		
		// if valence is 3, remove the front vertex.
		if (valence == 3)
		{
			Halfedge_handle g = h;			

			// Insertion of a vertex
			Halfedge_handle pass = h;
			
			//vector<Point3d> Vertices; //contains 1-ring and 2-rings vertices		
			
			Vector normal = Triangle_Normal(pass);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(pass,normal,T2);

			if ((T1 == CGAL::NULL_VECTOR) || (T2 == CGAL::NULL_VECTOR) || (normal == CGAL::NULL_VECTOR))
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}

			Point_Int Frenet;
			
			if (check_alpha == false)
			{
				Frenet.x = Decoder.decode(alpha);
				Frenet.y = Decoder.decode(alpha);
			}
			else 
			{
				Frenet.x = 0;
				Frenet.y = 0;
			}

			if (check_gamma == false)
				Frenet.z = Decoder.decode(gamma);
			
			else 
				Frenet.z = 0;			

			Frenet.x -= Alpha_offset;
			Frenet.y -= Alpha_offset;
			Frenet.z -= Gamma_offset;

			Point_Int Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);			

			// Calculate the position of barycenter.		
			Point3d Barycenter_real;
			Point3d Center_vertex;

			// Assign the region number to inserted vertex
			Halfedge_handle reg = h;

			int Selected_region = 500000;
			//int Number_vertices = 500000;
			vector<int> T_Bin;
			vector<int> T_Number;

			for(int i = 0; i < valence; i++)
			{
				int N1 = reg->vertex()->Region_Number;
				bool Is_existed = false;
				for(unsigned int j = 0; j < T_Bin.size(); j++)
				{
					if(N1 == T_Bin[j])
					{
						T_Number[j]++;
						Is_existed = true;
					}
				}
				if(!Is_existed)
				{
					T_Bin.push_back(N1);
					T_Number.push_back(1);
				}
				reg = reg->next();
			}
			int Max = -5000;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] > Max)
					Max = T_Number[i];
			}
			vector<int> T_possible_bin;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] == Max)
					T_possible_bin.push_back(T_Bin[i]);
			}

			if(T_possible_bin.size() == 1)
			{
				Selected_region = T_possible_bin[0];
			}
			else
			{
				Selected_region = 5000;
				for(unsigned int i = 0; i < T_possible_bin.size(); i++)
				{
					if(T_possible_bin[i] < Selected_region)
						Selected_region = T_possible_bin[i];
				}
			}


			int Region = Selected_region;

			if(this->m_N_remained_vertices[Region] < LIMIT_NUMBER)
			{
				Barycenter_real = Barycenter_Patch_After_Removal(h, valence);
				
				Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);
				Point_Int Center = BC + Diff;						
				Center_vertex = this->Change_Int_Real(Center, Component_ID);
				
			}

			else if(this->m_N_treated_vertices[Region] < MINIMUM_PREDICTION_NUMBER)
			{
				Barycenter_real = Barycenter_Patch_After_Removal(h, valence);

				Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);
				Point_Int Center = BC + Diff;						
				Center_vertex = this->Change_Int_Real(Center, Component_ID);

				double BC_S[3];
				double CRV_S[3];

				this->Convert_To_Spherical(Center_vertex, CRV_S);
				this->Convert_To_Spherical(Barycenter_real, BC_S);

				this->m_Rad_decision[Region] += (CRV_S[0] - BC_S[0]);
			}
			else
			{
				int Direction;

				if(this->m_Rad_decision[Region] > 0.0)
					Direction = 1;
				else
					Direction = 0;		

				Barycenter_real = this->JCW_Barycenter_Patch_After_Removal(h, valence, Direction);
				Point_Int BC = Change_Real_Int(Barycenter_real, Component_ID);

				Point_Int Center = BC + Diff;						
				Center_vertex = this->Change_Int_Real(Center, Component_ID);
			}				
			this->m_N_treated_vertices[Region]++;

			// Vertex insertion
			g = _pMesh.create_center_vertex(g);
			
			g->vertex()->point() = Center_vertex;
			g->vertex()->Seed_Edge = -1;
			
			int RO = this->GlobalCountOperation - Decompress_count;
			g->vertex()->Removal_Order = RO;	
			g->vertex()->Region_Number = Selected_region;
			

			// Vertex flags
			g->vertex()->Vertex_Flag = CONQUERED;
			g = h;
			g->facet()->Facet_Flag = CONQUERED;
			g = g->next();
			g = g->prev_on_vertex();
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			g->facet()->Facet_Flag = CONQUERED;			
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h1 = g->prev()->opposite();
				h1->facet()->Facet_Flag = CONQUERED;
				h1->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h1->next()->is_border_edge())
					Halfedges.push(h1->next()->opposite());
				if (!h1->prev()->is_border_edge())
					Halfedges.push(h1->prev()->opposite());
			}
			g = g->prev_on_vertex();
			g->facet()->Facet_Flag = CONQUERED;
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h2 = g->prev()->opposite();
				h2->facet()->Facet_Flag = CONQUERED;
				h2->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h2->next()->is_border_edge())
					Halfedges.push(h2->next()->opposite());
				if (!h2->prev()->is_border_edge())
					Halfedges.push(h2->prev()->opposite());
			}
			

			if ((this->IsColored) && (!this->IsOneColor))
			{
				//#ifdef PREDICTION_METHOD
				g = h;				

				Color_Unit Predicted_color;

				Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
								
				Color_Unit Color_difference;
				Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
				Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
				Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;								
				
				Color_Unit CV = Predicted_color + Color_difference;				

				g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);

				float LAB[3];
				LAB[0] = this->C0_Min + CV.c0 * Color_step;
				LAB[1] = this->C1_Min + CV.c1 * Color_step;
				LAB[2] = this->C2_Min + CV.c2 * Color_step;

				float RGB[3];
				LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

				g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);			
				//#endif

			}
			
			if ((this->IsColored) && (this->IsOneColor))
			{
				g = h->next();
				g->vertex()->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
			}
		}

		else
		{
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())
				Halfedges.push(h->prev()->opposite());
		}
	}	
}

void Compression_Valence_Component::JCW_Code_Difference_Histogram_Shifting(Polyhedron &_pMesh, const int & Component_ID)
{
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin();pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
	}
	
	//int Count_treated_vertices = 0;
	std::queue<Vertex*> vertices;	
	
	// push a input gate to begin loop 
	// in order to treat all vertices.	
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));

	int Vertex_index = 0;
	int Count_treated_vertices = 0;
	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;			
			v->Vertex_Number = Vertex_index;

			Halfedge_around_vertex_circulator h;

			vector<int> Temp_error;
			
			int	RO = this->GlobalCountOperation - this->Decompress_count;
			
			if(v->Removal_Order == RO)
			{
				Temp_error.push_back(v->JCW_Move_Error[0]);
				Temp_error.push_back(v->JCW_Move_Error[1]);
				Temp_error.push_back(v->JCW_Move_Error[2]);

				this->m_JCW_Move_Error.push_back(Temp_error);
				Count_treated_vertices++;

				Point3d Pos_real = v->point();
				Point_Int Pos_int = this->Change_Real_Int(Pos_real, Component_ID);

				Pos_int.x = Pos_int.x - Temp_error[0];
				Pos_int.y = Pos_int.y - Temp_error[1];
				Pos_int.z = Pos_int.z - Temp_error[2];

				v->point() = this->Change_Int_Real(Pos_int, Component_ID);
			}

			// First_vertex -> To find the seed edge
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}

			else
			{
				int Comp_number = -2;
				h = v->vertex_begin();
				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}	
			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h, h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}
			Vertex_index++;
		}
	}


	this->m_N_Errors.push_back(Count_treated_vertices);	
}

void Compression_Valence_Component::JCW_Decode_Difference_Histogram_Shifting(Polyhedron &_pMesh, const int & Component_ID)
{
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin();pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
	}	

	std::queue<Vertex*> vertices;	
	
	// push a input gate to begin loop 
	// in order to treat all vertices.	
	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));

	int Vertex_index = 0;

	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;			
			v->Vertex_Number = Vertex_index;

			Halfedge_around_vertex_circulator h;

			vector<int> Temp_error;				
			
			int RO = this->GlobalCountOperation - this->Decompress_count;
			
			if (v->Removal_Order == RO)
			{			
				int T_x = Decoder.decode(this->DM_JCW_MOVE_ERROR);
				int T_y = Decoder.decode(this->DM_JCW_MOVE_ERROR);
				int T_z = Decoder.decode(this->DM_JCW_MOVE_ERROR);

				//int T_x = v->JCW_Move_Error[0];
				if(T_x == 2)
					T_x = -1;				
				//int T_y = v->JCW_Move_Error[1];
				if(T_y == 2)
					T_y = -1;
				//int T_z = v->JCW_Move_Error[2];
				if(T_z == 2)
					T_z = -1;


				Point3d V_real = v->point();
				Point_Int V_int = this->Change_Real_Int(V_real, Component_ID);

				V_int.x = V_int.x - T_x;
				V_int.y = V_int.y - T_y;
				V_int.z = V_int.z - T_z;

				Point3d Res = this->Change_Int_Real(V_int,Component_ID);

				v->point() = Res;				
			}

			// First_vertex -> To find the seed edge
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}

			else
			{
				int Comp_number = -2;
				h = v->vertex_begin();
				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}	
			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h, h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}
			Vertex_index++;
		}
	}
}

// Description : Decoding of the regulation conquest
void Compression_Valence_Component::JCW_Un_Regulation_For_Insertion(Polyhedron &_pMesh, 																	
																	const int & Component_ID,
																	list<int> & FP_Connectivity,
																	list<Point3d> & SP_Moved_Position,
																	list<Point3d> & SP_Original_Position,
																	list<Point_Int> & SP_Watermarked_Position,
                                                                    list<vector<int> > & JCW_ERROR)
{	
	Init(_pMesh);	

	this->m_N_remained_vertices.clear();
	this->m_N_treated_vertices.clear();
	this->m_Rad_decision.clear();

	for(int i = 0; i < this->m_NumberRegion; i++)
	{		
		this->m_N_treated_vertices.push_back(0);
		this->m_Rad_decision.push_back(0.);

		this->m_N_remained_vertices.push_back(0);
	}
	
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		int RN = pVert->Region_Number;
		this->m_N_remained_vertices[RN]++;
	}

	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;

	//unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];	
	
	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))	
		hi++;	
	
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))// already visited.
			continue;
		
		// read connectivity information																	
		int valence = FP_Connectivity.front();
		FP_Connectivity.pop_front();
		
		// Store connectivity for compression
		this->JCW_Connectivity.push_back(valence);		
		
		// if valence is 3, remove the front vertex.
		if (valence == 3)
		{
			Halfedge_handle g = h;	
			
			// Calculate 3 vectors of the encountered patch.
			Vector normal = Triangle_Normal(g);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(g, normal, T2);

			if ((T1 == CGAL::NULL_VECTOR) || (T2 == CGAL::NULL_VECTOR) || (normal == CGAL::NULL_VECTOR))
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}

			// Assign the region number to inserted vertex
			Halfedge_handle reg = h;

			int Selected_region = 500000;
			//int Number_vertices = 500000;
			vector<int> T_Bin;
			vector<int> T_Number;

			for(int i = 0; i < valence; i++)
			{
				int N1 = reg->vertex()->Region_Number;
				bool Is_existed = false;
				for(unsigned int j = 0; j < T_Bin.size(); j++)
				{
					if(N1 == T_Bin[j])
					{
						T_Number[j]++;
						Is_existed = true;
					}
				}
				if(!Is_existed)
				{
					T_Bin.push_back(N1);
					T_Number.push_back(1);
				}
				reg = reg->next();
			}
			int Max = -5000;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] > Max)
					Max = T_Number[i];
			}
			vector<int> T_possible_bin;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] == Max)
					T_possible_bin.push_back(T_Bin[i]);
			}

			if(T_possible_bin.size() == 1)
			{
				Selected_region = T_possible_bin[0];
			}
			else
			{
				Selected_region = 5000;
				for(unsigned int i = 0; i < T_possible_bin.size(); i++)
				{
					if(T_possible_bin[i] < Selected_region)
						Selected_region = T_possible_bin[i];
				}
			}

			int Region = Selected_region;
			

			// Calculate the position of barycenter.		
			Point3d Barycenter_real;// = Barycenter_Patch_After_Removal(h, valence);

			Point_Int CRV_int = SP_Watermarked_Position.front();
			SP_Watermarked_Position.pop_front();

			Point3d CRV_real = this->Change_Int_Real(CRV_int, Component_ID);		

			// Calculate position of new barycenter 
			if(this->m_N_remained_vertices[Region] < LIMIT_NUMBER)
			{
				Barycenter_real = Barycenter_Patch_After_Removal(h, valence);				
			}

			else if(this->m_N_treated_vertices[Region] < MINIMUM_PREDICTION_NUMBER)
			{
				double BC_S[3];
				double CRV_S[3];
				
				Barycenter_real = Barycenter_Patch_After_Removal(h, valence);

				this->Convert_To_Spherical(CRV_real, CRV_S);
				this->Convert_To_Spherical(Barycenter_real, BC_S);

				this->m_Rad_decision[Region] += (CRV_S[0] - BC_S[0]);				
			}
			else
			{
				int Direction;

				if(this->m_Rad_decision[Region] > 0.0)
					Direction = 1;
				else
					Direction = 0;
				Barycenter_real = this->JCW_Barycenter_Patch_After_Removal(h, valence, Direction);
			}
			Point_Int Barycenter_int = Change_Real_Int(Barycenter_real, Component_ID);

			this->m_N_treated_vertices[Region]++;		

			Point_Int Vec_diff = CRV_int - Barycenter_int;			

			Point_Int Frenet_coords = Frenet_Rotation(Vec_diff, T1, T2, normal);

			this->JCW_Geometry.push_back(Frenet_coords);
			
			Point3d Origin_pt = SP_Original_Position.front();
			SP_Original_Position.pop_front();

			Point3d Watermarked_pt = SP_Moved_Position.front();
			SP_Moved_Position.pop_front();

			// Vertex insertion
			g = _pMesh.create_center_vertex(g);			
			
			g->vertex()->point() = CRV_real;
			
			g->vertex()->Watermarked_Position[0] = Watermarked_pt.x();
			g->vertex()->Watermarked_Position[1] = Watermarked_pt.y();
			g->vertex()->Watermarked_Position[2] = Watermarked_pt.z();

			g->vertex()->Seed_Edge = -1;
			g->vertex()->Removal_Order = this->GlobalCountOperation - this->Decompress_count;
			g->vertex()->Region_Number = Selected_region;

			vector<int> Err = JCW_ERROR.front();
			JCW_ERROR.pop_front();
			g->vertex()->JCW_Move_Error[0] = Err[0];
			g->vertex()->JCW_Move_Error[1] = Err[1];
			g->vertex()->JCW_Move_Error[2] = Err[2];

			if((Err[0] != 0) || (Err[1] != 0) || (Err[2] != 0))
			{
				g->vertex()->color(1., 0., 0.);
				Number_non_reversible_vertices++;
			}

			// Vertex flags
			g->vertex()->Vertex_Flag = CONQUERED;
			g = h;
			g->facet()->Facet_Flag = CONQUERED;
			g = g->next();
			g = g->prev_on_vertex();
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			g->facet()->Facet_Flag = CONQUERED;			
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h1 = g->prev()->opposite();
				h1->facet()->Facet_Flag = CONQUERED;
				h1->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h1->next()->is_border_edge())
					Halfedges.push(h1->next()->opposite());
				if (!h1->prev()->is_border_edge())
					Halfedges.push(h1->prev()->opposite());
			}
			g = g->prev_on_vertex();
			g->facet()->Facet_Flag = CONQUERED;
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h2 = g->prev()->opposite();
				h2->facet()->Facet_Flag = CONQUERED;
				h2->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h2->next()->is_border_edge())
					Halfedges.push(h2->next()->opposite());
				if (!h2->prev()->is_border_edge())
					Halfedges.push(h2->prev()->opposite());
			}			
		}

		else
		{
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())
				Halfedges.push(h->prev()->opposite());
		}
	}

}



// Description : Decoding function of decimation conquest
void Compression_Valence_Component::JCW_Un_Decimation_For_Insertion(Polyhedron & _pMesh, 
																	const int  & Component_ID,
																	list<int> & FP_Connectivity,
																	list<Point3d> & SP_Moved_Position, 
																	list<Point3d> & SP_Original_Position,
																	list<Point_Int> & SP_Watermarked_Position,
                                                                                                                                        list<vector<int> > & JCW_ERROR)
{
	Init(_pMesh);
	
	this->m_N_treated_vertices.clear();
	this->m_Rad_decision.clear();

	for(int i = 0; i < this->m_NumberRegion; i++)
	{		
		this->m_N_treated_vertices.push_back(0);
		this->m_Rad_decision.push_back(0.);	
	}	

	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;	

	//unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];

	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);
	
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID+1))
		hi++;	

	// Two vertices of seed edges are flaged CONQUERED
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;	
		
	// These vertices are also flages with sign flags for retriangulation
	hi->vertex()->Vertex_Sign = PLUS;
	hi->opposite()->vertex()->Vertex_Sign = MINUS;
	
	// Two vertices of seed edges are flaged CONQUERED
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;	


	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		unsigned int valence = 0, type = 0; // define type of retriangulation

		if (h->is_border() == true)
			continue;	
		
		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED
		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;
			
		valence = FP_Connectivity.front();
		FP_Connectivity.pop_front();

		this->JCW_Connectivity.push_back(valence);
				
		// if its front vertex is free
		if ((valence >= 3) && (valence <= 6))
		{
			Halfedge_handle g = h;

			type = Find_Type(g, valence);	

			// remove the front vertex if its removal does not viloate the manifold property			
			Vector normal = Normal_Patch(g, valence);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(g, normal, T2);			
				
			if ((T1 == CGAL::NULL_VECTOR) || (T2 == CGAL::NULL_VECTOR) || (normal == CGAL::NULL_VECTOR))
			{				
				T1 = Vector(1., 0., 0.);
				T2 = Vector(0., 1., 0.);					
				normal = Vector(0., 0., 1.);	
			}			

			// remove edges to re_find polygon and (check. and attribute sign flag.)
			bool Check_Validity = Remove_Edges(_pMesh, g, type);
			Check_Validity = false;

			//Check_Validity = false;// * * * * * * * * * * *////
			if (Check_Validity == false)
			{
				g = h;				

				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				for(unsigned int i = 0; i < valence; i++)
				{
					int N1 = reg->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}
					reg = reg->next();
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
				{
					Selected_region = T_possible_bin[0];
				}
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				int Region = Selected_region;				
				Point3d Barycenter_real;				

				Point_Int CRV_int = SP_Watermarked_Position.front();
				SP_Watermarked_Position.pop_front();

				Point3d CRV_real = this->Change_Int_Real(CRV_int, Component_ID);		

				// Calculate position of new barycenter 
				if(this->m_N_remained_vertices[Region] < LIMIT_NUMBER)
				{
					Barycenter_real = Barycenter_Patch_After_Removal(h, valence);				
				}

				else if(this->m_N_treated_vertices[Region] < MINIMUM_PREDICTION_NUMBER)
				{
					Barycenter_real = Barycenter_Patch_After_Removal(h, valence);	

					double BC_S[3];
					double CRV_S[3];

					this->Convert_To_Spherical(CRV_real, CRV_S);
					this->Convert_To_Spherical(Barycenter_real, BC_S);

					this->m_Rad_decision[Region] += (CRV_S[0] - BC_S[0]);					
				}
				else
				{
					int Direction;

					if(this->m_Rad_decision[Region] > 0.0)
						Direction = 1;
					else
						Direction = 0;
					Barycenter_real = this->JCW_Barycenter_Patch_After_Removal(h, valence, Direction);
				}

				Point_Int Barycenter_int = Change_Real_Int(Barycenter_real, Component_ID);

				this->m_N_treated_vertices[Region]++;			
				
				Point_Int Vec_diff = CRV_int - Barycenter_int;							
				
				Point_Int Frenet_coords = Frenet_Rotation(Vec_diff, T1, T2, normal);
				this->JCW_Geometry.push_back(Frenet_coords);				
				
				Point3d WM_pt = SP_Moved_Position.front();
				SP_Moved_Position.pop_front();
				Point3d Original_pt = SP_Original_Position.front();
				SP_Original_Position.pop_front();

				g = _pMesh.create_center_vertex(g);
				g->vertex()->point() = Original_pt;
				
				g->vertex()->Watermarked_Position[0] = WM_pt.x();
				g->vertex()->Watermarked_Position[1] = WM_pt.y();
				g->vertex()->Watermarked_Position[2] = WM_pt.z();

				g->vertex()->Region_Number = Selected_region;			
				
				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;

				g->vertex()->Vertex_Flag = CONQUERED;

				vector<int> Err = JCW_ERROR.front();
				JCW_ERROR.pop_front();
				g->vertex()->JCW_Move_Error[0] = Err[0];
				g->vertex()->JCW_Move_Error[1] = Err[1];
				g->vertex()->JCW_Move_Error[2] = Err[2];

				if((Err[0] != 0) || (Err[1] != 0) || (Err[2] != 0))
				{
					g->vertex()->color(1., 0., 0.);
					Number_non_reversible_vertices++;
				}

				g = h;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				g->vertex()->Vertex_Flag = CONQUERED;
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;

				for (unsigned int k = 0; k < (valence - 1); k++)
				{
					g = g->next()->opposite()->next();
					g->facet()->Facet_Flag = TO_BE_REMOVED;
					g->vertex()->Vertex_Flag = CONQUERED;
					g->opposite()->vertex()->Vertex_Flag = CONQUERED;
					if (g->is_border_edge() == false)
						Halfedges.push(g->opposite());
				}
				
				g->next()->vertex()->Seed_Edge = -1;				
			}
		}
		// In case of border edge.
		else if ((valence == 8) || (valence == 9))
		{
			
			Halfedge_handle g = h;
			
			type = Find_Type(g, valence - 5);			
			
			Vector normal = Normal_Patch(g, valence - 5);
			Vector T2 = CGAL::NULL_VECTOR;
			Vector T1 = Calculate_T1_T2(g, normal, T2);			
			
			//if ((T1 == CGAL::NULL_VECTOR) || (T2 == CGAL::NULL_VECTOR) || (normal == CGAL::NULL_VECTOR))
			{							
				T1 =     Vector(1.,0.,0.);
				T2 =     Vector(0.,1.,0.);	
				normal = Vector(0.,0.,1.);				
			}				
					

			// border edge with valence == 3
			if (valence == 8)
			{
				g = h;
				vector<Halfedge_handle> Border_edges;
				int Number_jump = 0;
				
				if (g->next()->is_border_edge())
					Number_jump = 0;
				if (g->prev()->is_border_edge())
				{
					Number_jump = 1;
					g = g->next();
				}
				
				Border_edges.push_back(g->opposite());
				
				g = g->prev();
				Border_edges.push_back(g->opposite());

				_pMesh.erase_facet(g);

				Halfedge_handle Prev_edge = Border_edges[1]->opposite()->prev();
					
				// g points the new vertex
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[1]->opposite());				
				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());	

				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				Halfedge_around_vertex_circulator Hvc = g->vertex()->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;
				
				CGAL_For_all(Hvc, Hvc_end)
				{
					int N1 = Hvc->opposite()->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}					
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
					Selected_region = T_possible_bin[0];
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				//int Region = Selected_region;					

				g->vertex()->Region_Number = Selected_region;			
				
				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;
				g->vertex()->Vertex_Flag = CONQUERED;

				Point3d Barycenter_real;
				
				Point_Int CRV_int = SP_Watermarked_Position.front();
				SP_Watermarked_Position.pop_front();

				Point3d CRV_real = this->Change_Int_Real(CRV_int, Component_ID);	

				// Calculate position of new barycenter 									
				Barycenter_real = Barycenter_Patch_Before_Removal(h);
				
				//this->m_N_treated_vertices[Region]++;

				Point_Int Barycenter_int = Change_Real_Int(Barycenter_real, Component_ID);

				Point_Int Vec_diff = CRV_int - Barycenter_int;				
				
				Point_Int Frenet_coords = Frenet_Rotation(Vec_diff, T1, T2, normal);
				this->JCW_Geometry.push_back(Frenet_coords);				
				
				Point3d WM_pt = SP_Moved_Position.front();
				SP_Moved_Position.pop_front();
				Point3d Original_pt = SP_Original_Position.front();
				SP_Original_Position.pop_front();
				
				g->vertex()->point() = Original_pt;
				
				g->vertex()->Watermarked_Position[0] = WM_pt.x();
				g->vertex()->Watermarked_Position[1] = WM_pt.y();
				g->vertex()->Watermarked_Position[2] = WM_pt.z();			

				vector<int> Err = JCW_ERROR.front();
				JCW_ERROR.pop_front();
				g->vertex()->JCW_Move_Error[0] = Err[0];
				g->vertex()->JCW_Move_Error[1] = Err[1];
				g->vertex()->JCW_Move_Error[2] = Err[2];		

				if((Err[0] != 0) || (Err[1] != 0) || (Err[2] != 0))
				{
					g->vertex()->color(1., 0., 0.);
				}
		
				
				Halfedge_handle Tag_handle;

				if (Number_jump == 0)
					Tag_handle = Border_edges[1];
				if (Number_jump == 1)
					Tag_handle = Border_edges[0]->opposite();				
				Tag_handle->vertex()->Vertex_Flag = CONQUERED;

				if ((type == 1) || (type == 2) || (type == 4))
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = PLUS;
				}
				else if (type == 3)
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = MINUS;
				}
				for (int i=0; i <2; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}				
			}
			// border edge with valence == 4
			else if (valence == 9)
			{		
				g = h;
				int Number_jump = -1;				
				vector<Halfedge_handle> Border_edges;

				if ((type == 5) || (type == 8))
				{
					// jump == 0;
					if (g->next()->is_border_edge())
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());
						Border_edges.push_back(g->prev()->opposite()->next()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					// jump == 1;
					else if (g->prev()->opposite()->next()->is_border_edge())
					{
						Number_jump = 1;
						
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[0]->opposite());
					}
					
					// jump == 2;
					else
					{
						Number_jump = 2;

						Border_edges.push_back(g->prev()->opposite()->next()->opposite());
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
				}
				else
				{
					if (g->prev()->is_border_edge())
					{
						Number_jump = 2;

						Border_edges.push_back(g->next()->opposite()->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					else if (g->next()->opposite()->prev()->is_border_edge())
					{
						Number_jump = 1;

						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
					else
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

				}				
				g = h;											
				
				// to create the new facets
				Halfedge_handle Prev_edge = Border_edges[2]->opposite()->prev();
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[2]->opposite());								
				
				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());

				// Assign the region number to inserted vertex
				Halfedge_handle reg = g;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				Halfedge_around_vertex_circulator Hvc = reg->vertex()->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;
				
				CGAL_For_all(Hvc, Hvc_end)
				{
					int N1 = Hvc->opposite()->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}					
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
					Selected_region = T_possible_bin[0];
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				//int Region = Selected_region;

				Point_Int CRV_int = SP_Watermarked_Position.front();
				SP_Watermarked_Position.pop_front();

				Point3d CRV_real = this->Change_Int_Real(CRV_int, Component_ID);		
				
				// Calculate position of new barycenter 
				Point3d Barycenter_real = Barycenter_Patch_Before_Removal(h);
				Point_Int Barycenter_int = Change_Real_Int(Barycenter_real, Component_ID);

				//this->m_N_treated_vertices[Region]++;			

				Point_Int Vec_diff = CRV_int - Barycenter_int;				
				
				Point_Int Frenet_coords = Frenet_Rotation(Vec_diff, T1, T2, normal);
				this->JCW_Geometry.push_back(Frenet_coords);				
				
				Point3d WM_pt = SP_Moved_Position.front();
				SP_Moved_Position.pop_front();
				Point3d Original_pt = SP_Original_Position.front();
				SP_Original_Position.pop_front();
				
				g->vertex()->point() = Original_pt;
				
				g->vertex()->Watermarked_Position[0] = WM_pt.x();
				g->vertex()->Watermarked_Position[1] = WM_pt.y();
				g->vertex()->Watermarked_Position[2] = WM_pt.z();

				g->vertex()->Region_Number = Selected_region;			
				
				int RO = this->GlobalCountOperation - this->Decompress_count;
				g->vertex()->Removal_Order = RO;

				g->vertex()->Vertex_Flag = CONQUERED;

				vector<int> Err = JCW_ERROR.front();
				JCW_ERROR.pop_front();
				g->vertex()->JCW_Move_Error[0] = Err[0];
				g->vertex()->JCW_Move_Error[1] = Err[1];
				g->vertex()->JCW_Move_Error[2] = Err[2];

				if((Err[0] != 0) || (Err[1] != 0) || (Err[2] != 0))
				{
					g->vertex()->color(1., 0., 0.);
				}

				//vertex_tag
				if (Number_jump == 0)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}

				else if (Number_jump == 1)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}

				}
				else // jump == 2
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}
				
				for (int i=0; i < 3; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					
					Border_edges[i]->vertex()->Vertex_Flag = CONQUERED;
					Border_edges[i]->opposite()->vertex()->Vertex_Flag = CONQUERED;
					
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}					
			}
		}		
		
		
		//  the symbol == N
		else if (valence == 7)
		{
			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;

			if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)
				Halfedges.push(h->prev()->opposite());	
		}
	}
}


// Description : To obtain region number of each inserted vertices
void Compression_Valence_Component::JCW_Un_Regulation_For_Region_Detection(Polyhedron &_pMesh, 																	
																	       const int & Component_ID,
																		   list<int> & FP_Connectivity, 
																		   list<Point3d> & FP_Geometry, 
																		   list<int> & FP_Region_Number)
{	
	Init(_pMesh);
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;

	//unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];	
	
	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))	
		hi++;	
	
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))// already visited.
			continue;
		
		// read connectivity information
		int valence = this->Connectivity[Component_ID].front() + 3;
		this->Connectivity[Component_ID].pop_front();
		
		FP_Connectivity.push_back(valence);
		
		// if valence is 3, remove the front vertex.
		if (valence == 3)
		{
			Halfedge_handle g = h;			

			// Assign the region number to inserted vertex
			Halfedge_handle reg = h;

			int Selected_region = 500000;
			//int Number_vertices = 500000;
			vector<int> T_Bin;
			vector<int> T_Number;

			for(int i = 0; i < valence; i++)
			{
				int N1 = reg->vertex()->Region_Number;
				bool Is_existed = false;
				for(unsigned int j = 0; j < T_Bin.size(); j++)
				{
					if(N1 == T_Bin[j])
					{
						T_Number[j]++;
						Is_existed = true;
					}
				}
				if(!Is_existed)
				{
					T_Bin.push_back(N1);
					T_Number.push_back(1);
				}
				reg = reg->next();
			}
			int Max = -5000;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] > Max)
					Max = T_Number[i];
			}
			vector<int> T_possible_bin;
			for(unsigned int i = 0; i < T_Number.size(); i++)
			{
				if(T_Number[i] == Max)
					T_possible_bin.push_back(T_Bin[i]);
			}

			if(T_possible_bin.size() == 1)
			{
				Selected_region = T_possible_bin[0];
			}
			else
			{
				Selected_region = 5000;
				for(unsigned int i = 0; i < T_possible_bin.size(); i++)
				{
					if(T_possible_bin[i] < Selected_region)
						Selected_region = T_possible_bin[i];
				}
			}

			int Region = Selected_region;

			FP_Region_Number.push_back(Region);

			Point_Int CRV_int = this->Geometry[Component_ID].front();
			this->Geometry[Component_ID].pop_front();
			Point3d CRV_real = this->Change_Int_Real(CRV_int, Component_ID);

			FP_Geometry.push_back(CRV_real);

			// Vertex insertion
			g = _pMesh.create_center_vertex(g);			
			g->vertex()->point() = CRV_real;

			g->vertex()->Seed_Edge = -1;
			//g->vertex()->Removal_Order = this->Decompress_count;
			g->vertex()->Region_Number = Selected_region;			

			// Vertex flags
			g->vertex()->Vertex_Flag = CONQUERED;
			g = h;
			g->facet()->Facet_Flag = CONQUERED;
			g = g->next();
			g = g->prev_on_vertex();
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			g->facet()->Facet_Flag = CONQUERED;			
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h1 = g->prev()->opposite();
				h1->facet()->Facet_Flag = CONQUERED;
				h1->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h1->next()->is_border_edge())
					Halfedges.push(h1->next()->opposite());
				if (!h1->prev()->is_border_edge())
					Halfedges.push(h1->prev()->opposite());
			}
			g = g->prev_on_vertex();
			g->facet()->Facet_Flag = CONQUERED;
			g->opposite()->vertex()->Vertex_Flag = CONQUERED;
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h2 = g->prev()->opposite();
				h2->facet()->Facet_Flag = CONQUERED;
				h2->next()->vertex()->Vertex_Flag = CONQUERED;
				if (!h2->next()->is_border_edge())
					Halfedges.push(h2->next()->opposite());
				if (!h2->prev()->is_border_edge())
					Halfedges.push(h2->prev()->opposite());
			}			
		}

		else
		{
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())
				Halfedges.push(h->prev()->opposite());
		}
	}
}



// Description : Decoding function of decimation conquest
void Compression_Valence_Component::JCW_Un_Decimation_For_Region_Detection(Polyhedron & _pMesh, 
																	       const int  & Component_ID,
																		   list<int> & FP_Connectivity, 
																		   list<Point3d> & FP_Geometry, 
																		   list<int> & FP_Region_Number)
{
	Init(_pMesh);
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;	

	//unsigned Qbit = this->Qbit[Component_ID] + this->NumberChangeQuantization[Component_ID];	

	float Color_step = 0.0;
	if(this->NumberColorQuantization[Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, this->NumberColorQuantization[Component_ID]);
	
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID+1))
		hi++;	

	// Two vertices of seed edges are flaged CONQUERED
	hi->vertex()->Vertex_Flag = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag = CONQUERED;	
		
	// These vertices are also flages with sign flags for retriangulation
	hi->vertex()->Vertex_Sign = PLUS;
	hi->opposite()->vertex()->Vertex_Sign = MINUS;
	
	// Two vertices of seed edges are flaged CONQUERED
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;	

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		unsigned int valence = 0,type = 0; // define type of retriangulation

		if (h->is_border() == true)
			continue;	
		
		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED
		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;
			
		valence = this->Connectivity[Component_ID].front() + 3;
		this->Connectivity[Component_ID].pop_front();

		FP_Connectivity.push_back(valence);
				
		// if its front vertex is free
		if ((valence >= 3) && (valence <= 6))
		{
			Halfedge_handle g = h;

			type = Find_Type(g, valence);
						
			// remove edges to re_find polygon and (check. and attribute sign flag.)
			bool Check_Validity = Remove_Edges(_pMesh, g, type);
			Check_Validity = false;

			//Check_Validity = false;// * * * * * * * * * * *////
			if (Check_Validity == false)
			{
				g = h;				

				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				for(int i = 0; i < (int)valence; i++)
				{
					int N1 = reg->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}
					reg = reg->next();
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
				{
					Selected_region = T_possible_bin[0];
				}
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}

				Point_Int Center = this->Geometry[Component_ID].front();
				this->Geometry[Component_ID].pop_front();			

				Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);

				FP_Geometry.push_back(Center_vertex);
				
				g = _pMesh.create_center_vertex(g);
				g->vertex()->point() = Center_vertex;
				g->vertex()->Region_Number = Selected_region;				

				FP_Region_Number.push_back(Selected_region);

				g->vertex()->Vertex_Flag = CONQUERED;

				g = h;
				g->facet()->Facet_Flag = TO_BE_REMOVED;
				g->vertex()->Vertex_Flag = CONQUERED;
				g->opposite()->vertex()->Vertex_Flag = CONQUERED;

				for (unsigned int k = 0; k < (valence - 1); k++)
				{
					g = g->next()->opposite()->next();
					g->facet()->Facet_Flag = TO_BE_REMOVED;
					g->vertex()->Vertex_Flag = CONQUERED;
					g->opposite()->vertex()->Vertex_Flag = CONQUERED;
					if (g->is_border_edge() == false)
						Halfedges.push(g->opposite());
				}
				
				g->next()->vertex()->Seed_Edge = -1;
			}
		}

		// In case of border edge.
		else if ((valence == 8) || (valence == 9))
		{	
			Halfedge_handle g = h;	
			type = Find_Type(g, valence - 5);					

			// border edge with valence == 3
			if (valence == 8)
			{
				g = h;
				vector<Halfedge_handle> Border_edges;
				int Number_jump = 0;
				
				if (g->next()->is_border_edge())
					Number_jump = 0;
				if (g->prev()->is_border_edge())
				{
					Number_jump = 1;
					g = g->next();
				}
				
				Border_edges.push_back(g->opposite());
				
				g = g->prev();
				Border_edges.push_back(g->opposite());

				_pMesh.erase_facet(g);

				Halfedge_handle Prev_edge = Border_edges[1]->opposite()->prev();
					
				// g points the new vertex
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());	

				Point_Int Center = this->Geometry[Component_ID].front();
				this->Geometry[Component_ID].pop_front();			

				Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);

				FP_Geometry.push_back(Center_vertex);

				g->vertex()->point() = Center_vertex;

				// Assign the region number to inserted vertex
				Halfedge_handle reg = g;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				Halfedge_around_vertex_circulator Hvc = g->vertex()->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;
								
				
				
				CGAL_For_all(Hvc, Hvc_end)
				{
					int N1 = Hvc->opposite()->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}					
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
				{
					Selected_region = T_possible_bin[0];
				}
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}
				
				FP_Region_Number.push_back(Selected_region);
				g->vertex()->Region_Number = Selected_region;

				g->vertex()->Vertex_Flag = CONQUERED;
						
				
				Halfedge_handle Tag_handle;

				if (Number_jump == 0)
					Tag_handle = Border_edges[1];
				if (Number_jump == 1)
					Tag_handle = Border_edges[0]->opposite();				
				Tag_handle->vertex()->Vertex_Flag = CONQUERED;

				if ((type == 1) || (type == 2) || (type == 4))
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = PLUS;
				}
				else if (type == 3)
				{
					if (Tag_handle->vertex()->Vertex_Sign == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign = MINUS;
				}
				for (int i=0; i <2; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}				
			}
			// border edge with valence == 4
			else if (valence == 9)
			{			
				g = h;
				int Number_jump = -1;				
				vector<Halfedge_handle> Border_edges;

				if ((type == 5) || (type == 8))
				{
					// jump == 0;
					if (g->next()->is_border_edge())
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());
						Border_edges.push_back(g->prev()->opposite()->next()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					// jump == 1;
					else if (g->prev()->opposite()->next()->is_border_edge())
					{
						Number_jump = 1;
						
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[0]->opposite());
					}
					
					// jump == 2;
					else
					{
						Number_jump = 2;

						Border_edges.push_back(g->prev()->opposite()->next()->opposite());
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
				}
				else
				{
					if (g->prev()->is_border_edge())
					{
						Number_jump = 2;

						Border_edges.push_back(g->next()->opposite()->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

					else if (g->next()->opposite()->prev()->is_border_edge())
					{
						Number_jump = 1;

						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());

						_pMesh.erase_facet(Border_edges[0]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}
					else
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->prev()->opposite());

						_pMesh.erase_facet(Border_edges[2]->opposite());
						_pMesh.erase_facet(Border_edges[1]->opposite());
					}

				}
				
				g = h;
				
				// to create the new facets
				Halfedge_handle Prev_edge = Border_edges[2]->opposite()->prev();
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[2]->opposite());

				Prev_edge = Prev_edge->next();
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				_pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());

				
				Point_Int Center = this->Geometry[Component_ID].front();
				this->Geometry[Component_ID].pop_front();			

				Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);			
				FP_Geometry.push_back(Center_vertex);				
				g->vertex()->point() = Center_vertex;

				// Assign the region number to inserted vertex
				Halfedge_handle reg = g;

				int Selected_region = 500000;
				//int Number_vertices = 500000;
				vector<int> T_Bin;
				vector<int> T_Number;

				Halfedge_around_vertex_circulator Hvc = reg->vertex()->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;
				
				CGAL_For_all(Hvc, Hvc_end)
				{
					int N1 = Hvc->opposite()->vertex()->Region_Number;
					bool Is_existed = false;
					for(unsigned int j = 0; j < T_Bin.size(); j++)
					{
						if(N1 == T_Bin[j])
						{
							T_Number[j]++;
							Is_existed = true;
						}
					}
					if(!Is_existed)
					{
						T_Bin.push_back(N1);
						T_Number.push_back(1);
					}					
				}
				int Max = -5000;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] > Max)
						Max = T_Number[i];
				}
				vector<int> T_possible_bin;
				for(unsigned int i = 0; i < T_Number.size(); i++)
				{
					if(T_Number[i] == Max)
						T_possible_bin.push_back(T_Bin[i]);
				}

				if(T_possible_bin.size() == 1)
				{
					Selected_region = T_possible_bin[0];
				}
				else
				{
					Selected_region = 5000;
					for(unsigned int i = 0; i < T_possible_bin.size(); i++)
					{
						if(T_possible_bin[i] < Selected_region)
							Selected_region = T_possible_bin[i];
					}
				}
				
				FP_Region_Number.push_back(Selected_region);				
				g->vertex()->Region_Number = Selected_region;
				
				


				//vertex_tag
				if (Number_jump == 0)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign = PLUS;
					}

				}

				else if (Number_jump == 1)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}

				}
				else // jump == 2
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = MINUS;
					}
					else
					{
						if (Border_edges[0]->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign = PLUS;
					}
				}
				
				for (int i=0; i < 3; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag = CONQUERED;
					
					Border_edges[i]->vertex()->Vertex_Flag = CONQUERED;
					Border_edges[i]->opposite()->vertex()->Vertex_Flag = CONQUERED;
					
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}					
			}
		}		
		
		//  the symbol == N
		else if (valence == 7)
		{
			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag = CONQUERED;
			h->next()->vertex()->Vertex_Flag = CONQUERED;

			if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign == NOSIGN)
					h->next()->vertex()->Vertex_Sign = PLUS;
			}
			
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)
				Halfedges.push(h->prev()->opposite());	
		}
	}	
}

// Error metric which measures importance of color and geometry for each vertex.
// Used to prevent removal of the visually important vertex.
bool Compression_Valence_Component::Error_Projected_Surface(Polyhedron            &  _pMesh, 
														    const Halfedge_handle & _h, 
														    const int             & _Component_ID, 
														    const double          & Mean_color, 
														    const double          & Mean_area)
{
	Halfedge_handle g = _h;
	int Valence = (int)g->next()->vertex_degree();
	int Type = Find_Type(g, Valence);	

	float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;

	float Color_step = 0.0;

	if(this->NumberColorQuantization[_Component_ID] == 0)
		Color_step = this->Color_Quantization_Step;		
	else
		Color_step = this->Color_Quantization_Step * pow(2.0, (double)this->NumberColorQuantization[_Component_ID]);

	vector<float> Center_color;
	Center_color.push_back(C0_min + g->next()->vertex()->color_int(0)*Color_step);
	Center_color.push_back(C1_min + g->next()->vertex()->color_int(1)*Color_step);
	Center_color.push_back(C2_min + g->next()->vertex()->color_int(2)*Color_step);

	vector<Color_Unit> Neighbors_color;
	
	vector<float> Projected_color;
	
	double Patch_area = 0.0;
	
	// Calculate area of patch
	g = g->next();
	for(int i =0; i < Valence; i++)
	{
		g = g->opposite();
		Color_Unit c;

		c.c0 = g->vertex()->color_int(0);
		c.c1 = g->vertex()->color_int(1);
		c.c2 = g->vertex()->color_int(2);
		Neighbors_color.push_back(c);
		
		Patch_area += Area_Facet_Triangle(g);

		g = g->prev();
	}
	
	Color_Unit Average_color;

	if (Valence == 3)
	{		
		for(int i = 0; i < 3; i++)
			Average_color = Average_color + Neighbors_color[i];		
	}

	if (Valence == 4)
	{
		if ((Type == 5) || (Type == 8))
		{
			Average_color = Average_color + Neighbors_color[1];
			Average_color = Average_color + Neighbors_color[3];
		}

		if ((Type == 6) || (Type == 7))
		{
			Average_color = Average_color + Neighbors_color[0];
			Average_color = Average_color + Neighbors_color[2];
		}
	}

	if(Valence == 5)
	{
		if ((Type == 9) || (Type == 12))
		{
			Average_color = Average_color + Neighbors_color[1];
			Average_color = Average_color + Neighbors_color[3];
			Average_color = Average_color + Neighbors_color[4];
		}

		if (Type == 10)
		{
			Average_color = Average_color + Neighbors_color[0];
			Average_color = Average_color + Neighbors_color[1];
			Average_color = Average_color + Neighbors_color[3];
		}
		
		if (Type == 11)
		{
			Average_color = Average_color + Neighbors_color[0];
			Average_color = Average_color + Neighbors_color[2];
			Average_color = Average_color + Neighbors_color[4];
		}
	}
	

	if (Valence == 6)
	{
		if((Type == 13) || (Type == 16))
		{
			Average_color = Average_color + Neighbors_color[1];
			Average_color = Average_color + Neighbors_color[3];
			Average_color = Average_color + Neighbors_color[5];
		}
		else
		{		
			Average_color = Average_color + Neighbors_color[0];
			Average_color = Average_color + Neighbors_color[2];
			Average_color = Average_color + Neighbors_color[4];		
		}
	}
	
	float Divised_c0,Divised_c1,Divised_c2;
	if(Valence != 4)
	{
		Divised_c0 = (float)Average_color.c0 / 3.0;
		Divised_c1 = (float)Average_color.c1 / 3.0;
		Divised_c2 = (float)Average_color.c2 / 3.0;
	}
	else
	{
		Divised_c0 = (float)Average_color.c0 / 2.0;
		Divised_c1 = (float)Average_color.c1 / 2.0;
		Divised_c2 = (float)Average_color.c2 / 2.0;
	}

	Projected_color.push_back(C0_min + Divised_c0 * Color_step);
	Projected_color.push_back(C1_min + Divised_c1 * Color_step);
	Projected_color.push_back(C2_min + Divised_c2 * Color_step);
	

	// Color distance between the original color and estimated color
	double Color_distance = (Projected_color[0] - Center_color[0]) * (Projected_color[0] - Center_color[0]) + 
							(Projected_color[1] - Center_color[1]) * (Projected_color[1] - Center_color[1]) + 
							(Projected_color[2] - Center_color[2]) * (Projected_color[2] - Center_color[2]);

	Color_distance = sqrt(Color_distance);	
	
	// Mean color is obtained using number of vertices.
	// X 3.0 cause we should use the number of edges.
	double Relative_color_distance = Color_distance / Mean_color * 3.0;
	
	// Averaged area of triangles of patch
	double Area_per_triangle = Patch_area / double(Valence);
	Area_per_triangle *= pow(((double)10.0/(double)this->HighestLengthBB[_Component_ID]), 2.0);

	double Relative_geo_distance = Area_per_triangle / Mean_area; 

	double Global_distance = (Relative_color_distance) * (Relative_geo_distance);

	if(Global_distance > 0.5)
		return true;
	else
		return false;
}

int Compression_Valence_Component::JCW_Decompress_One_Level(Polyhedron &_pMesh, const char* File_Name, const int &Noise_number)
{

	int res = 1;

	if (this->Decompress_count == 0)
	{
		FILE * fp = fopen("Keys.txt", "rb");
		if(fp!=NULL)
		{		

			int division, reversibility;

            res = fread(&this->m_VC[0], sizeof(double), 1, fp);
            res = fread(&this->m_VC[1], sizeof(double), 1, fp);
            res = fread(&this->m_VC[2], sizeof(double), 1, fp);

            res = fread(&this->m_Dist, sizeof(double), 1, fp);

			res = fread(&this->m_NumberBin, sizeof(int), 1, fp);
			res = fread(&this->m_NumberRegion, sizeof(int), 1, fp);
			res = fread(&this->m_EmbeddingStrength, sizeof(int), 1, fp);
			res = fread(&division, sizeof(int), 1, fp);
			res = fread(&reversibility, sizeof(int), 1, fp);
			res = fread(&this->Division_Threshold, sizeof(int), 1, fp);

			if(division == 1)
				this->Is_Division_Big_Regions_Enabled = true;
			else
				this->Is_Division_Big_Regions_Enabled = false;

			if(reversibility == 1)
				this->Is_Complete_Reversibility_Enabled = true;
			else
				this->Is_Complete_Reversibility_Enabled = false;		
			
			fclose(fp);
		}	

		this->JCW_Generate_Regions_On_Base_Mesh(_pMesh);		
	}
	
	
	if(this->Decompress_count < this->GlobalCountOperation)
	{
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (this->Decompress_count < this->ComponentOperations[Component_ID])
			{	
				int Operation = Decoder.get_bits(2);
				if (Operation == 0)
				{					
					Processing_Component Process;					
					
					this->JCW_Un_Regulation(_pMesh, Decoder, Component_ID);					
					this->JCW_Un_Decimation_Conquest(_pMesh, Decoder, Component_ID);					

					this->Initialize_Spherical_Coordinates(_pMesh);										
					this->JCW_Region_Mass_Center_Extract_Watermark(_pMesh);				

					if (this->Is_Division_Big_Regions_Enabled)
						this->JCW_Divide_Big_Regions(_pMesh, this->Division_Threshold);				
					
					if(this->Is_Complete_Reversibility_Enabled)
						this->JCW_Decode_Difference_Histogram_Shifting(_pMesh, Component_ID);

				}
				else if (Operation == 1)
					this->Augment_Geometry_Quantization_Precision(_pMesh, Decoder, Component_ID);		
				else if (Operation == 2)
					this->Augment_Color_Quantization_Precision(_pMesh, Decoder, Component_ID);				
			}		
		}			
	}	

	JCW_Evaluate_Robustness();
		
	this->Decompress_count++;			
	//
	_pMesh.compute_normals();	
	
	return this->Decompress_count;
}

int Compression_Valence_Component::JCW_Decompress_One_Level_Without_Extraction(Polyhedron &_pMesh, const char* File_Name)
{
    int res;

	if (this->Decompress_count == 0)
	{
		FILE * fp = fopen("Keys.txt", "rb");
		if(fp != NULL)
		{		

			int division, reversibility;

			res = fread(&this->m_VC[0], sizeof(double), 1, fp);
			res = fread(&this->m_VC[1], sizeof(double), 1, fp);
			res = fread(&this->m_VC[2], sizeof(double), 1, fp);

			res = fread(&this->m_Dist, sizeof(double), 1, fp);

			res = fread(&this->m_NumberBin, sizeof(int), 1, fp);
			res = fread(&this->m_NumberRegion, sizeof(int), 1, fp);
			res = fread(&this->m_EmbeddingStrength, sizeof(int), 1, fp);
			res = fread(&division, sizeof(int), 1, fp);
			res = fread(&reversibility, sizeof(int), 1, fp);
			res = fread(&this->Division_Threshold, sizeof(int), 1, fp);

			if(division == 1)
				this->Is_Division_Big_Regions_Enabled = true;
			else
				this->Is_Division_Big_Regions_Enabled = false;

			if(reversibility == 1)
				this->Is_Complete_Reversibility_Enabled = true;
			else
				this->Is_Complete_Reversibility_Enabled = false;		
			fclose(fp);
		}		

		this->JCW_Generate_Regions_On_Base_Mesh(_pMesh);		
	}
	
	
	if (this->Decompress_count < this->GlobalCountOperation)
	{
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (this->Decompress_count < this->ComponentOperations[Component_ID])
			{	
				int Operation = Decoder.get_bits(2);
				if (Operation == 0)
				{
					this->JCW_Un_Regulation(_pMesh, Decoder, Component_ID);					
					this->JCW_Un_Decimation_Conquest(_pMesh, Decoder, Component_ID);	

					this->Initialize_Spherical_Coordinates(_pMesh);
					
					if (this->Is_Division_Big_Regions_Enabled)
						this->JCW_Divide_Big_Regions(_pMesh, this->Division_Threshold);				
					
					if(this->Is_Complete_Reversibility_Enabled)
						this->JCW_Decode_Difference_Histogram_Shifting(_pMesh, Component_ID);
				}
				else if (Operation == 1)
					this->Augment_Geometry_Quantization_Precision(_pMesh, Decoder, Component_ID);		
				else if (Operation == 2)
					this->Augment_Color_Quantization_Precision(_pMesh, Decoder, Component_ID);				
			}
		}			
	}	
			
	this->Decompress_count++;	
	_pMesh.compute_normals();
	
	
	return this->Decompress_count;
}

void Compression_Valence_Component::Read_Information_To_Hide(const char * Message)
{
   //char *Message;

	/*FILE * fp = fopen("Inserted_message.txt", "r");
	char buffer[200];
	if(fp != NULL)
	{
		res = fgets(buffer, 200, fp);
	}*/	

	int i = 0;	
	while (Message[i] != '\0')
	{
		int N_char = CHAR_BIT;
		char c = Message[i];
		while (N_char > 0)
		{
			--N_char;
			int bit = (c & (1 << N_char) ? 1 : 0);
			this->m_Watermarks.push_back(bit);
		}

		i++;
	}

	FILE * f_insert = fopen("Inserted_watermarks.txt", "w");	

	fclose(f_insert);	
}

QString Compression_Valence_Component::Write_Information_To_Hide()
{
	list<int>::iterator it = this->m_Watermarks.begin();

	unsigned int Size_watermarks = this->m_Watermarks.size();

	int N_char = (int)(Size_watermarks / 8);

	vector<char> ret;

	for(int i = 0; i < N_char; i++)
	{
		char c = 0;
		for(int j = 0; j < 8; j++)
		{
			if(*it == 1)
			{
				c += (char)pow(2.0, 7-j);
			}
			it++;
		}

		ret.push_back(c);
	}
	
	QString q;

	for(int i = 0; i <N_char; i++)
		q.append(ret[i]);	

	return q;
	
}

int Compression_Valence_Component::JCW_Divide_Big_Regions(Polyhedron &_pMesh, const int & _Thres_divide_regions)
{
	// To count number of vertices of each region.
	vector<int> Number_vertices_regions(this->m_NumberRegion, 0);
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		int Region_number = pVert->Region_Number;
		Number_vertices_regions[Region_number]++;		
	}	
	
	// Determine what region to divide.
	int Number_region_to_divide = 0;
	vector<int> Region_to_divide(this->m_NumberRegion, -1);
	vector<int> Matching_index_new_region(this->m_NumberRegion, -1);
	for(int i = 0; i < this->m_NumberRegion; i++)
	{
		if (Number_vertices_regions[i] >= _Thres_divide_regions)
		{			
			Region_to_divide[i] = 1;
			Matching_index_new_region[i] = this->m_NumberRegion + Number_region_to_divide;
			Number_region_to_divide++;
		}
	}

	if(Number_region_to_divide == 0)
		return 0;
	
	// Initialize vertex flags
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;		
	}
	
	// Upgrade region number of mesh
	this->m_NumberRegion += Number_region_to_divide;

	std::queue<Vertex*> vertices;

	std::queue<Vertex*> Vertices_region;

	vector<int> Verify_processed_region(this->m_NumberRegion, -1);

	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge != 0) || (hi->opposite()->vertex()->Seed_Edge != 1))
		hi++;	

	vertices.push(&(*(hi->vertex())));
	vertices.push(&(*(hi->opposite()->vertex())));

	int Vertex_index = 0;	

	while(!vertices.empty())
	{
		Vertex * v = vertices.front();
		vertices.pop();		
		
		if (v->Vertex_Flag == CONQUERED)
			continue;	
		
		else
		{
			v->Vertex_Flag = CONQUERED;		
			v->Vertex_Number = Vertex_index;

			int RN = v->Region_Number;

			if((Region_to_divide[RN] == 1) && (Verify_processed_region[RN] == -1))
			{
				Verify_processed_region[RN] = 1;
				Vertices_region.push(v);
			}

			Halfedge_around_vertex_circulator h;

			// First_vertex -> To find the seed edge
			if (Vertex_index == 0)
			{
				h = v->vertex_begin();
				do
				{
					h++;
				}while(&(*(h)) != &(*(hi)));
			}

			else
			{
				int Comp_number = -2;
				h = v->vertex_begin();
				Halfedge_around_vertex_circulator h2 = h;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number == Comp_number)
						break;
				}
			}	
			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
					vertices.push(&(*(h->opposite()->vertex())));
			}
			Vertex_index++;
		}
	}

	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
		pVert->Vertex_Flag = FREE;		

	vector<int> Number_New_Region(this->m_NumberRegion, 0);

	while(!Vertices_region.empty())
	{
		Vertex * v = Vertices_region.front();
		Vertices_region.pop();		
		
		int RN = v->Region_Number;
		
		if (v->Vertex_Flag == CONQUERED)
			continue;

		if (Region_to_divide[RN] == -1)
			continue;
		
		if (Number_New_Region[RN] > (Number_vertices_regions[RN] / 2))
			continue;
		
		else
		{
			v->Vertex_Flag = CONQUERED;
			Number_New_Region[RN]++;

			v->Region_Number = Matching_index_new_region[RN];

			Halfedge_around_vertex_circulator h;
			
			int Comp_number = -2;
			h = v->vertex_begin();
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Number > Comp_number)
					Comp_number = h->opposite()->vertex()->Vertex_Number;
			}

			h = h2;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Number == Comp_number)
					break;
			}				
			
			h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag == FREE)
				{
					if(h->opposite()->vertex()->Region_Number == RN)
					{
						Vertices_region.push(&(*(h->opposite()->vertex())));
					}
				}
			}			
		}
	}
	return Number_region_to_divide;
}


void Compression_Valence_Component::JCW_Colorify_Regions(Polyhedron &_pMesh)
{
	srand(time(NULL));
	//vector<vector<float>> Color;
	
	if ((unsigned)Region_Color.size() < (unsigned)this->m_NumberRegion)
	{
		int Number_Color_To_Insert = (unsigned)this->m_NumberRegion - (unsigned)Region_Color.size();

		for(int i = 0; i < Number_Color_To_Insert; i++)
		{				
			int R = rand() % 255;
			int G = rand() % 255;
			int B = rand() % 255;
			
			vector<float> fff;
			fff.push_back((float)R / 255.);
			fff.push_back((float)G / 255.);
			fff.push_back((float)B / 255.);

			this->Region_Color.push_back(fff);		
		}

	}

	/*for(int i = 0; i < this->m_NumberRegion; i++)
	{				
		int R = rand() % 255;
		int G = rand() % 255;
		int B = rand() % 255;
		
		vector<float> fff;
		fff.push_back((float)R / 255.);
		fff.push_back((float)G / 255.);
		fff.push_back((float)G / 255.);


		Color.push_back(fff);			
	}*/
	

	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		int NB = pVert->Region_Number;		
//		pVert->color(Color[NB][0], Color[NB][1], Color[NB][2]);
		pVert->color(this->Region_Color[NB][0], this->Region_Color[NB][1], this->Region_Color[NB][2]);
	}		
}
vector<double> Compression_Valence_Component::JCW_Evaluate_Robustness(void)
{
	FILE * insert = fopen("Inserted_watermarks.txt", "r");
	FILE * extract = fopen("Extracted_watermark.txt","r");

	vector<int> Inserted_bits;
	vector<int> Extracted_bits;
	
	char pLine[512];

	while(fgets(pLine,512,insert))
	{
		int w = 0, level = 0;
		if (sscanf(pLine,"%d %d",&w,&level) == 2)
		{
			Inserted_bits.push_back(w);			
		}
	}
	fclose(insert);

	while(fgets(pLine,512,extract))
	{
		int w = 0, level = 0;
		if (sscanf(pLine,"%d %d",&w,&level) == 2)
		{
			Extracted_bits.push_back(w);			
		}
	}
	fclose(extract);	
	
	FILE * f_res = fopen("Robustness_results.txt","a");

	int Count_correct_watermark = 0;
	int Count_incorrect_watermark = 0;
	
	double Average_insert = 0.;
	double Average_extract = 0.;
	unsigned int Number_watermarks = 0;
	
	
	if(Inserted_bits.size() < Extracted_bits.size())
		Number_watermarks = Inserted_bits.size();
	else
		Number_watermarks = Extracted_bits.size();



	for(unsigned int i = 0; i < Number_watermarks; i++)
	{
		if(Inserted_bits[i] == Extracted_bits[i])
			Count_correct_watermark++;
		else
			Count_incorrect_watermark++;
	}

	
	for(unsigned int i = 0; i < Inserted_bits.size(); i++)
	{
		Average_insert += (double)Inserted_bits[i];
	}

	for(unsigned int i = 0; i < Extracted_bits.size(); i++)
	{
		Average_extract += (double)Extracted_bits[i];
	}

	Average_insert /= (double)Inserted_bits.size();
	Average_extract /= (double)Extracted_bits.size();
	
		
	double Numerator = 0.0;
	for(unsigned int i = 0; i < Number_watermarks; i++)
	{
		Numerator += (Inserted_bits[i] - Average_insert) * (Extracted_bits[i] - Average_extract);
	}

	double denom1 = 0.0;
	double denom2 = 0.0;
	for(unsigned int i = 0; i < Inserted_bits.size(); i++)
	{
		denom1 += (Inserted_bits[i] - Average_insert) * (Inserted_bits[i] - Average_insert);
	}

	for(unsigned int i = 0; i < Extracted_bits.size(); i++)
	{
		denom2 += (Extracted_bits[i] - Average_extract) * (Extracted_bits[i] - Average_extract);
	}

	double denom = sqrt(denom1*denom2);

	double corr = Numerator/denom;

	

	double ber = (double)Count_correct_watermark/((double)Count_correct_watermark + (double)Count_incorrect_watermark);

	fprintf(f_res, "%lf \n", ber);
	fprintf(f_res, "%lf \n", corr);
	
	fclose(f_res);

	vector<double> res;
	res.push_back(ber);
	res.push_back(corr);

	return res;
}

void Compression_Valence_Component::Clear_After_Compression(void)
{
		
	for(unsigned i = 0; i < ListOperation.size(); i++)
		ListOperation.clear();
	for(unsigned i = 0; i < Connectivity.size(); i++)
		Connectivity.clear();
	for(unsigned i = 0; i < Geometry.size(); i++)
		Geometry.clear();
	for(unsigned i = 0; i < NumberSymbol.size(); i++)
		NumberSymbol.clear();
	for(unsigned i = 0; i < NumberVertices.size(); i++)
		NumberVertices.clear();
	for(unsigned i = 0; i < AlphaRange.size(); i++)
		AlphaRange.clear();
	for(unsigned i = 0; i < AlphaOffset.size(); i++)
		AlphaOffset.clear();
	for(unsigned i = 0; i < GammaRange.size(); i++)
		GammaRange.clear();
	for(unsigned i = 0; i < GammaOffset.size(); i++)
		GammaOffset.clear();
	for(unsigned i = 0; i < VertexColor.size(); i++)
		VertexColor.clear();
	for(unsigned i = 0; i < m_JCW_Move_Error.size(); i++)
		m_JCW_Move_Error.clear();    	
	 
	

// Colors
	m_JCW_Move_Error.clear();

	// Used for adatative quantization.				
       QuantizationCorrectVector.clear();
        NumberQuantizationLayer.clear();
		
		//for color
        NumberProcessedVertices.clear();
        ColorChildcellIndex.clear();
      ColorEncoderIndex.clear();

		ListOperation.clear();
       Connectivity.clear();
        Geometry.clear();
      NumberSymbol.clear();
     NumberVertices.clear();

AlphaRange.clear();
       AlphaOffset.clear();
         GammaRange.clear();
         GammaOffset.clear();

		// Colors
		VertexColor.clear(); // contain color error of all removed vertices
		InterVertexColor.clear();	
		IsClosed.clear();    // to know if the mesh is open or closed.				
       
			
	InterGeometry.clear();
	 InterConnectivity.clear();		

      
		
		// Quantization
		Qbit.clear(); 
			xmin.clear();
		ymin.clear();
		zmin.clear();		
		xmax.clear();
			ymax.clear();
		zmax.clear();		
			Quantization_Step.clear();
				
		
	 HighestLengthBB.clear();
	 ComponentVolume.clear();
		 ComponentArea.clear();
	   ComponentNumberVertices.clear();		
	 NumberColorQuantization.clear();			
		
	//	// mapping table
	// ColorArray.clear(); // Color table
 //PredictedColorArray.clear(); // Predicted values of colors in color table using prediction based on neighbors.		
	// DifferenceColor.clear(); // Difference between original and quantized vertex color. need to turn into lossless coding.
	//		   ColorIndex.clear();
	//   ReorderingColorIndex.clear();
	//	   InterColorIndex.clear();		
	//	   Number_color_index.clear(); // contain occurence of each initial color
	//   IsKnownIndex.clear();
		
		
	 NumberDecimation.clear(); // To stock number of Decimation.
	 NumberChangeQuantization.clear(); // to stock number of under_quantization.						
				
		Decompress_count = 0;

		 ComponentOperations.clear();

		
		 m_Number_Vertices_Per_Regions.clear();
		m_Watermarks.clear();

	

		 m_N_remained_vertices.clear();
		 m_N_treated_vertices.clear();
		 m_Rad_decision.clear();

        m_JCW_Move_Error.clear();
		m_N_Errors.clear();		
		JCW_Connectivity.clear();
		JCW_Geometry.clear();
}


//int numVertex = _pMesh.size_of_vertices();;
					//Vector centroid = Point3d(0,0,0) - CGAL::ORIGIN;
					//double distancetoCentroid = 0.0;
					//Vertex_iterator	pVertex;
					//for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
					//{
					//	Vector vectemp = pVertex->point() - CGAL::ORIGIN;
					//	centroid = centroid + vectemp;
					//}
					//centroid = centroid/numVertex;

					//// calculate the average distance from vertices to mesh centre
					//for (pVertex = _pMesh.vertices_begin(); pVertex!= _pMesh.vertices_end(); pVertex++)
					//{
					//	Vector vectemp = pVertex->point() - CGAL::ORIGIN;
					//	distancetoCentroid = distancetoCentroid + (double)std::sqrt((vectemp - centroid) * (vectemp - centroid));
					//}
					//distancetoCentroid = distancetoCentroid/numVertex;

					//// add random uniform-distributed (between [-noiseLevel, +noiseLevel])
					//srand((unsigned)time(NULL));
					//double noisex, noisey, noisez;
					//double noiseLevel = distancetoCentroid * 0.005;//noiseIntensity;
					//for (pVertex = _pMesh.vertices_begin(); pVertex!= _pMesh.vertices_end(); pVertex++)
					//{
					//	// keep boundaries untouched if demanded by user
					//	bool is_border_vertex = false;
					//	bool stopFlag = false;
					//	Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
					//	do
					//	{
					//		if (hav->is_border()==true)
					//		{
					//			is_border_vertex = true;
					//			stopFlag = true;
					//		}
					//		hav++;
					//	} while ((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

					//	/*if ((preserveBoundaries==true)&&(is_border_vertex==true))
					//		continue;*/

					//	noisex = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
					//	noisey = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
					//	noisez = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
					//	
					//	if(pVertex->Removal_Order == this->GlobalCountOperation - Decompress_count)
					//	{
					//		Vector temp = Point3d(noisex, noisey, noisez) - CGAL::ORIGIN;						
					//		pVertex->point() = pVertex->point() + temp;
					//	}
					//}

					//// for correct rendering, we need to update the mesh normals
					//_pMesh.compute_normals();


void Compression_Valence_Component::Decompression_From_File(Polyhedron &_pMesh)
{
	if (this->Current_level >= this->Total_layer)
		return;

	if (this->Process_level == 0)
		this->Write_Info(_pMesh);

	this->Current_level = this->Decompress_Each_Step(_pMesh, this->File_name.c_str());
	if (this->Current_level > this->Process_level)
	{
		this->Process_level++;
		this->Write_Info(_pMesh);
	}
}

void Compression_Valence_Component::Decompression_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh)
{
	if (this->Process_level == 0)
		this->Write_Info(pMesh);

	this->Copy_Polyhedron.copy(&(pMesh), &(New_mesh));
	this->Attibute_Seed_Gate_Flag(pMesh, New_mesh);
	New_mesh.compute_normals();			
	this->Process_level++;
	this->Visu_level++;

	this->Decompress_Each_Step(pMesh, this->File_name.c_str());

	Write_Info(pMesh);

	float prog = (float)this->Calculate_Current_File_Size() / this->Compressed_file_size * 100;
	float ratio = 1/((float)this->Calculate_Current_File_Size() / this->Initial_file_size);

	this->Prog.push_back(prog);
	this->Ratio.push_back(ratio);
}

void Compression_Valence_Component::Decompression_Specific_Level_From_File(Polyhedron &pMesh, const int & Wl)
{
	if (this->Process_level == 0)
		this->Write_Info(pMesh);

	int CLevel = this->Current_level;

	int Wanted_level = Wl;

	if (Wanted_level < 0)
		Wanted_level = 0;

	if (Wanted_level > this->Total_layer)
		Wanted_level = this->Total_layer;

	if (Wanted_level == CLevel)
		return;	

	if (Wanted_level > CLevel)
	{
		while(this->Current_level != Wanted_level)
		{
			this->Current_level = this->Decompress_Each_Step(pMesh, this->File_name.c_str());
		
			if (this->Current_level > this->Process_level)
			{
				this->Process_level++;
				this->Write_Info(pMesh);
			}
		}
	}
	else if (Wanted_level < CLevel)
	{
		this->Stop_Decoder();

		this->Decompress_Init(pMesh);
		
		while(this->Current_level != Wanted_level)
		{	
			this->Current_level = this->Decompress_Each_Step(pMesh, this->File_name.c_str());		
			this->Current_level++;
		}		
	}
}


void Compression_Valence_Component::Decompression_Coarser_From_File(Polyhedron &pMesh)
{
	if (this->Current_level <= 0)
		return;

	this->Stop_Decoder();
	
	int Temp_level = this->Current_level - 1;
	this->Decompress_Init(pMesh);
	

	this->Current_level = 0;
	while(this->Current_level != Temp_level)
		this->Current_level = this->Decompress_Each_Step(pMesh, this->File_name.c_str());
}

void Compression_Valence_Component::Decompression_All_From_File(Polyhedron &pMesh)
{
	if (this->Process_level == 0)
		Write_Info(pMesh);
	while(this->Current_level != this->Total_layer)
	{
		this->Current_level = this->Decompress_Each_Step(pMesh, this->File_name.c_str());
		this->Process_level++;
		Write_Info(pMesh);
	}
}

void Compression_Valence_Component::JCW_Decompression_From_File(Polyhedron &_pMesh)
{
	if (this->Current_level >= this->Total_layer)
		return;

	if (this->Process_level == 0)
		this->Write_Info(_pMesh);

	this->Current_level = this->JCW_Decompress_One_Level(_pMesh, this->File_name.c_str(), -1);

	this->Message = this->Write_Information_To_Hide();
	if (this->Current_level > this->Process_level)
	{
		this->Process_level++;
		this->Write_Info(_pMesh);
	}

}

void Compression_Valence_Component::JCW_Decompression_Without_Extraction_From_File(Polyhedron &_pMesh)
{
	if (this->Current_level >= this->Total_layer)
		return;

	if (this->Process_level == 0)
		this->Write_Info(_pMesh);

	this->Current_level = this->JCW_Decompress_One_Level_Without_Extraction(_pMesh, this->File_name.c_str());

	this->Message = this->Write_Information_To_Hide();
	if (this->Current_level > this->Process_level)
	{
		this->Process_level++;
		this->Write_Info(_pMesh);
	}
}
void Compression_Valence_Component::JCW_Decompression_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh)
{
	if (this->Process_level == 0)
		this->Write_Info(pMesh);

	this->Copy_Polyhedron.copy(&(pMesh), &(New_mesh));
	this->Attibute_Seed_Gate_Flag(pMesh, New_mesh);
	New_mesh.compute_normals();			
	this->Process_level++;
	this->Visu_level++;

	this->JCW_Decompress_One_Level(pMesh, this->File_name.c_str(),-1);
	this->Message = this->Write_Information_To_Hide();

	this->Write_Info(pMesh);

	float prog = (float)this->Calculate_Current_File_Size() / this->Compressed_file_size * 100;
	float ratio = 1/((float)this->Calculate_Current_File_Size() / this->Initial_file_size);

	this->Prog.push_back(prog);
	this->Ratio.push_back(ratio);
}
void Compression_Valence_Component::JCW_Decompression_Without_Extraction_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh)
{
	if (this->Process_level == 0)
		this->Write_Info(pMesh);

	this->Copy_Polyhedron.copy(&(pMesh), &(New_mesh));
	this->Attibute_Seed_Gate_Flag(pMesh, New_mesh);
	New_mesh.compute_normals();			
	this->Process_level++;
	this->Visu_level++;

	this->JCW_Decompress_One_Level_Without_Extraction(pMesh, this->File_name.c_str());
	this->Message = this->Write_Information_To_Hide();

	this->Write_Info(pMesh);

	float prog = (float)this->Calculate_Current_File_Size() / this->Compressed_file_size * 100;
	float ratio = 1/((float)this->Calculate_Current_File_Size() / this->Initial_file_size);

	this->Prog.push_back(prog);
	this->Ratio.push_back(ratio);
}

void Compression_Valence_Component::Write_Info(Polyhedron &_pMesh)
{
	
	if (this->Process_level == 0)
	{
		this->Dec_File_Info = this->File_name;
		size_t point = this->Dec_File_Info.find('.');
		this->Dec_File_Info.replace(point+1,3,"txt");

		this->Dec_Info = fopen(this->Dec_File_Info.c_str(),"w");
	}
	else
		this->Dec_Info = fopen(this->Dec_File_Info.c_str(),"a");

	int CLevel = 0, Number_vertices = 0;

	if (this->Sequence)
	{
		CLevel = this->Visu_level;
		Number_vertices = (int)_pMesh.size_of_vertices();
	}
	else
	{
		CLevel = this->Current_level;
		Number_vertices = (int)_pMesh.size_of_vertices();
	}

	unsigned Current_file_size = this->Calculate_Current_File_Size();

	float prog = (float)Current_file_size / this->Compressed_file_size * 100;			
			
	if (this->Process_level == this->Total_layer)
		prog = 100.0;

	float ratio = 1 / ((float)Current_file_size / this->Initial_file_size);

	fprintf(this->Dec_Info,"Level %2d   #v : %8d      %6u bytes     Prog : %7.3f %%    Ratio : %9.3f\n", CLevel, Number_vertices, Current_file_size, prog, ratio);
	fclose(this->Dec_Info);
}


QString Compression_Valence_Component::Show_Text(void)
{
	int CLevel = -1;
	if (this->Sequence)
		CLevel = this->Visu_level;
	else
		CLevel = this->Current_level;

	QString string = QString("Current level : %1/%2").arg(CLevel, 2).arg(this->Total_layer, 2);
			
	string += "   |   ";

	float prog =0, ratio = 0;

	if (this->Sequence)
	{
		prog = this->Prog[CLevel];
		ratio = this->Ratio[CLevel];
	}
	else
	{
		prog = (float)this->Calculate_Current_File_Size() / this->Compressed_file_size * 100;
		ratio = 1/((float)this->Calculate_Current_File_Size() / this->Initial_file_size);
	}

	if (CLevel != this->Total_layer)
		string += QString("Prog : %1 % ").arg(prog, 3, 'f', 3); 
	else
		string += "Prog : 100.000 %";

	string += "   |   ";
	string += QString("Ratio : %1\n").arg(ratio, 3, 'f', 3);
			
	string += this->Message;
	return string;
}

#endif
