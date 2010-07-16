#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

//////////////////////////////////////////////////////////////////////////////
//														 					//
// Author : Ho LEE													   		//
//														   					//
//////////////////////////////////////////////////////////////////////////////


#include "Compression_Valence_Component.h"
#include "Compression_Valence_Polyhedron.h"
#include "Matrix3X3.h"
#include "Compression_Valence_Common.h"
#include <CGAL/Timer.h>

#include "Color_Distance/VectorT.h"
#include "Color_Distance/Deviation.h"
#include <map>

#define PREDICTION_METHOD

// Choix de methodes de prediction.
#define USE_LEE_METHOD
//#define USE_SIMPLE_PREDICTION
//#define USE_YOON_METHOD


///!!! Do not use mapping table method for color !!!///
//#define MAPPING_TABLE_METHOD
#ifdef MAPPING_TABLE_METHOD
	#define COLOR_QUANTIZATION
#endif

#define NUMBER_SEEDS 256

#define USE_ESTIMATION_ADAPTIVE_QUANTIZATION

//#define MEASURE_COLOR_DEVIATION

//#define SAVE_INTERMEDIATE_MESHES

//#define GET_COLOR_TABLE_HISTOGRAM

//#define DEBUG_MODE


#define AC_BUFFER 1024 * 10000

//#define SPLIT_FILE_EACH_RESOLUTION

double Compression_Valence_Component::Main_Function(Polyhedron     & pMesh,
										const char*      File_Name,
										const int      & _Qbit,
										const int      & NVertices,
										const bool       Normal_flipping,
										const bool       Use_metric,
										const float    & Metric_thread,
										const bool       Use_forget_metric,
										const int      & Forget_value, 
										const bool       Compression_selected,
										const bool       Adaptive_quantization, 
										unsigned       & Number_layers, 
										unsigned       & Init_number_vertices,
										unsigned       & Final_number_vertices,
										unsigned       & Connectivity_size, 
										unsigned       & Color_size, 
										unsigned       & Total_size, 
										const unsigned & Initial_file_size)
{
	
	Timer timer;
	timer.start();	
	
	Init_number_vertices = (unsigned)pMesh.size_of_vertices();

	this->Global_Initialization(pMesh, _Qbit, File_Name);	

	if (Adaptive_quantization)
		this->Adaptive_Quantization(pMesh, NVertices, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value, _Qbit);
	else
		this->Simplification(pMesh, NVertices, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value);		
	
	// Compression
	if (Compression_selected)
		this->Compression(pMesh, File_Name, _Qbit, Connectivity_size, Color_size, Total_size, Initial_file_size);		

	Number_layers = this->GlobalCountOperation;
	Final_number_vertices = (unsigned)pMesh.size_of_vertices();

	timer.stop();
	return timer.time();	
}




// Description : To select the input gate.
void Compression_Valence_Component::Global_Initialization(Polyhedron &pMesh, 
											  const int & _Qbit,
											  const char * File_Name)
{

	#ifdef MEASURE_COLOR_DEVIATION
	
	// To compare with the original colored mesh
	this->Set_Original_Color_Mesh(pMesh, File_Name);
	#endif			
	
	#ifndef USE_ESTIMATION_ADAPTIVE_QUANTIZATION	
	// To compare with the original mesh
	if (Adaptive_quantization)
		pMesh.write_off("original_temp.off",false,false);
	#endif

	
	// (1) Determination if the mesh is colored.
	// (2) Conversion of color space.
	// (3) Quantization of converted colors into "vertex->color_int()"
	// (4) Re-paint mesh with re_calculated rgb colors from quantized converted color.
	// (5) Establish color palette.
	// (6) Color quantization - Descrease of possible color numbers 
	this->Color_Initialization(pMesh);
	
	
	for (Vertex_iterator pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
	{
		pVertex->Seed_Edge = OTHER_COORDINATE;	
		pVertex->Component_Number = -1;
	}

	
	// Several components
	// (1) To get the number of components;
	// (2) To get the volume, area and the number of vertices of each component;
	// (3) To know if each component is closed or not.
	pMesh.tag_facets(-1);
	int Component_index = 0;

	for (Facet_iterator pFacet = pMesh.facets_begin(); pFacet != pMesh.facets_end(); pFacet++)
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
				
				area += Area_Facet_Triangle(F->halfedge());

				F->tag(Component_index);

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
						if ((!pHalfedge->is_border_edge()) && (!Is_Border_Vertex(pHalfedge)) && (!Is_Border_Vertex(pHalfedge->opposite())) || (pHalfedge->next()->vertex_degree() != 6))
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
					}
				}
			}

			this->xmin.push_back(xmin);
			this->ymin.push_back(ymin);
			this->zmin.push_back(zmin);

			this->xmax.push_back(xmax);
			this->ymax.push_back(ymax);
			this->zmax.push_back(zmax);

			
			// volume, area, number of vertices of each component
			Vector e1(xmax-xmin, 0., 0.);
			Vector e2(0., ymax-ymin, 0.);

			Vector normal = CGAL::cross_product(e1,e2);
			double base_area = sqrt(normal * normal);
			double volume = base_area * (zmax - zmin);	
			
			this->ComponentVolume.push_back(volume);
			this->ComponentArea.push_back(area);
			this->ComponentNumberVertices.push_back(Number_vertices);

			// Pour savoir si ce composant est ouvert ou pas.
			this->IsClosed.push_back(Is_closed);

			Component_index++;
		}
	}

	this->NumberComponents = Component_index;	

	list<int> li;
	list<Point_Int> pi;
	list<Color_Unit> cu;
	
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

		// Qbit for each component
		this->Qbit.push_back(_Qbit);
		
		// Displacement vector for under quantization for each component
		this->QuantizationCorrectVector.push_back(li);

		// The number of vertices for each under_quantization
		this->NumberQuantizationLayer.push_back(li);
	}
	

	// Quantization of each component
	this->Quantization(pMesh);
}	


/*
	Description : Quantize all vertices so that the new positions 
	are reguliraly spaced in the 3D space. */
void Compression_Valence_Component::Quantization(Polyhedron & pMesh) 
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

		float Q_Pas = (float)max / (float)NbInteraval;
		this->QuantizationPas.push_back(Q_Pas);
				
	}	

	// Vertex quantization
	for (Vertex_iterator pVert = pMesh.vertices_begin();pVert != pMesh.vertices_end(); pVert++)
	{
		double x = pVert->point().x();
		double y = pVert->point().y();
		double z = pVert->point().z();

		int Component_ID = pVert->Component_Number;

		int Qx = (int)(ceil((x - this->xmin[Component_ID]) / this->QuantizationPas[Component_ID])) - 1;
		if (Qx == -1)
			Qx = 0;
		int Qy = (int)(ceil((y - this->ymin[Component_ID]) / this->QuantizationPas[Component_ID])) - 1;
		if (Qy == -1)
			Qy = 0;
		int Qz = (int)(ceil((z - this->zmin[Component_ID]) / this->QuantizationPas[Component_ID])) - 1;
		if (Qz == -1)
			Qz = 0;
		
		pVert->point() = Point3d(this->xmin[Component_ID] + (Qx + 0.5) * this->QuantizationPas[Component_ID],
							     this->ymin[Component_ID] + (Qy + 0.5) * this->QuantizationPas[Component_ID],
							     this->zmin[Component_ID] + (Qz + 0.5) * this->QuantizationPas[Component_ID]);		
	}
}



// this->ColorArray -> contains all initial colors present in the input mesh.
void Compression_Valence_Component::Color_Initialization(Polyhedron &pMesh)
{	
	
	// (1) To determine if the mesh is colored.
	//     We consider the mesh is colored if the number of colors are > 2.	
	Vertex_iterator pVertex = pMesh.vertices_begin();
	
	
	
	this->OnlyColor[0] = pVertex->color(0);
	this->OnlyColor[1] = pVertex->color(1);
	this->OnlyColor[2] = pVertex->color(2);
	
	for (; pVertex != pMesh.vertices_end(); pVertex++)
	{	
		if ((pVertex->color(0) != this->OnlyColor[0]) || (pVertex->color(1) != this->OnlyColor[1]) || (pVertex->color(2) != this->OnlyColor[2]))
		{
			this->IsColored = true;
			break;
		}		
	}
	
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
		for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
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
		
		const int Nb_interaval_c0 = (int)pow(2.0, C0_QUANTIZATION) - 1;
		//const int Nb_interaval_c1 = (int)pow(2.0, C1_QUANTIZATION) - 1;
		//const int Nb_interaval_c2 = (int)pow(2.0, C2_QUANTIZATION) - 1;		
		
		float Color_max = C0_max - C0_min;
		if ( C1_max - C1_min > Color_max)
			Color_max = C1_max - C1_min;
		if ( C2_max - C2_min > Color_max)
			Color_max = C2_max - C2_min;

		//Information needed to quantize converted color.
		this->C0_Min = C0_min;
		this->C1_Min = C1_min;
		this->C2_Min = C2_min;		

		this->Color_Quantization_Step = (float)((Color_max) / Nb_interaval_c0);
		//this->Color_Quantization_Step = (float)((C1_max - C1_min) / Nb_interaval_c1);
		//this->Color_Quantization_Step = (float)((C2_max - C2_min) / Nb_interaval_c2);				
		
		//Enter quantized color valued into "vertex->color_int" to use lated.
		//Also, recalculate vertex color using these quantized color values.
		Color_Unit Resulting_color;		
		float New_vertex_color[3];
		float Reconstructed_color[3];

		int Color_index = 0;
		for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
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
			

			#ifdef MAPPING_TABLE_METHOD
			/*
			 * To obtain color table. (chech if the color is already in the color table,
			 * If new color -> insert in the color table. */

			
			bool Is_existing_color = false;
			
			for (int i = 0; i < Color_index; i++)
			{
				if (Resulting_color == this->ColorArray[i])
				{
					pVertex->Vertex_Color_Index = i; // Assign a color index to vertex
					this->Number_color_index[i] += 1; // Increment number of this color.
					Is_existing_color = true;
					break;
				}
			}
		
			if (!Is_existing_color)	// new color
			{
				pVertex->Vertex_Color_Index = Color_index;

				Color_index++;			
				this->ColorArray.push_back(Resulting_color);
				this->Number_color_index.push_back(1);
			}
			#endif
		}

		#ifdef COLOR_QUANTIZATION
		this->Color_Quantization(pMesh);
		#endif

		#ifdef GET_COLOR_TABLE_HISTOGRAM
		// To observe histogram of color table.
		FILE *Hist_color_table;
		Hist_color_table = fopen("Histogram_Color_Table.txt", "w");
		
		if (Hist_color_table == NULL)
		{
			fprintf(stderr, "Can't open input file in.list!\n");
			exit(1);
		}
		
		for (int i=0; i< this->Number_color_index.size(); i++)
		{
			fprintf(Hist_color_table,"%d\n",Number_color_index[i]);		
		}
		
		//  Close of all FILE pointers.		
		fclose(Hist_color_table);
		#endif
	}

	
	
	#ifdef MAPPING_TABLE_METHOD
	Color_Unit TC;
			
	for (int i = 0; i < NUMBER_SEEDS; i++)
		this->PredictedColorArray.push_back(TC);
	this->Number_color_index.clear();
	this->ColorArray.clear();

	#endif

	
}


/* To quantize color */
#ifdef COLOR_QUANTIZATION
void Compression_Valence_Component::Color_Quantization(Polyhedron &pMesh)
{
	// Contains the most frequent colors(Seed colors).(index of this->ColorArray)
	vector<int> Seed_containers;		
	
	// select initial seeds =  K more frequent colors
	for (int i = 0; i < NUMBER_SEEDS; i++)
	{
		int Next_biggest_index = -1;
		int Actual_occurence = -1;
  
		for (unsigned j = 0; j < this->Number_color_index.size(); j++)
		{
			bool Check_existed = false;

			int Temp_occurence = this->Number_color_index[j];
			for (unsigned k = 0; k < Seed_containers.size(); k++)
			{
				if (j == Seed_containers[k])
				{
					Check_existed = true;
					break;
				}
			}
			if (Check_existed == false)
			{
				if (Temp_occurence > Actual_occurence)
				{
					Actual_occurence = Temp_occurence;
					Next_biggest_index = j;
				}
			}			
		}
		Seed_containers.push_back(Next_biggest_index);
	}

	vector<Color_Unit> Seeds;
	vector<Color_Unit> Centroids;

	vector<int> Number_colors_cluster; // number of vertices (colors) in a cluster.
	vector<int> Color_repartition; // mapping between ColorArray <-> QuantizedColorArray
	

	for (int i = 0; i < NUMBER_SEEDS; i++)
	{
		int Corresponding_color_index = Seed_containers[i];
		Centroids.push_back(this->ColorArray[Corresponding_color_index]);

		Number_colors_cluster.push_back(0);
	}

	for (unsigned i =0; i < this->ColorArray.size(); i++)	
		Color_repartition.push_back(-1);	
	
	int Count_loop = 0;
	bool Check_ending;
	
	do
	{
		Count_loop++;
		
		Seeds = Centroids;
		for (unsigned i = 0; i < Centroids.size(); i++)
		{
			Centroids[i].c0 = 0;
			Centroids[i].c1 = 0;
			Centroids[i].c2 = 0;
		}
		
		Number_colors_cluster.clear();

		for (int i = 0; i < NUMBER_SEEDS; i++)
			Number_colors_cluster.push_back(0);

		for (unsigned i = 0; i < this->ColorArray.size(); i++)
		{
			int Temp_seed_number = -1;
			double Temp_distance = 500000.0;

			for (unsigned j = 0; j < Seeds.size(); j++)
			{
				double distance = Sqrt_Color(this->ColorArray[i], Seeds[j]);
				if (distance < Temp_distance)
				{
					Temp_distance = distance;
					Temp_seed_number = j;
				}				
			}

			Color_repartition[i] = Temp_seed_number;
			Number_colors_cluster[Temp_seed_number] += 1;

			Centroids[Temp_seed_number].c0 += this->ColorArray[i].c0;
			Centroids[Temp_seed_number].c1 += this->ColorArray[i].c1;
			Centroids[Temp_seed_number].c2 += this->ColorArray[i].c2;
		}

		for (unsigned i = 0 ; i < Centroids.size(); i++)
		{
			if (Number_colors_cluster[i] != 0)
			{
				Centroids[i].c0 = (int)floor((double)Centroids[i].c0 / Number_colors_cluster[i] + 0.5);
				Centroids[i].c1 = (int)floor((double)Centroids[i].c1 / Number_colors_cluster[i] + 0.5);
				Centroids[i].c2 = (int)floor((double)Centroids[i].c2 / Number_colors_cluster[i] + 0.5);
			}
		}
		
		// Substitute found color by one of original color.
		for (unsigned i = 0; i < Centroids.size(); i++)
		{
			Color_Unit col;
			double Temp_dist = 500000.0;
			int Temp_index = -1;
			for (unsigned j = 0; j < this->ColorArray.size(); j++)
			{
				double distance = Sqrt_Color(Centroids[i], this->ColorArray[j]);
				if (distance < Temp_dist)
				{
					Temp_dist = distance;
					Temp_index = j;
				}
			}
			
			col = this->ColorArray[Temp_index];
			Centroids[i] = col;		
		}		
		Check_ending = Check_End_Clustering(Seeds, Centroids);
	}while((!Check_ending) && ( Count_loop < 200 ));
	
	float New_vertex_color[3];
	float Reconstructed_color[3];
	

	#ifdef MAPPING_TABLE_METHOD

	int Color_diff_max_c0 = -5000, Color_diff_max_c1 = -5000, Color_diff_max_c2 = -5000;
	int Color_diff_min_c0 =  5000, Color_diff_min_c1 =  5000, Color_diff_min_c2 =  5000;	
	
	for (Vertex_iterator pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
	{
		int Vertex_color_index = pVertex->Vertex_Color_Index; // original color index
		int New_index = Color_repartition[Vertex_color_index]; // corresponding quantized color index

		pVertex->Vertex_Color_Index = New_index;	
		
		Color_Unit col = Seeds[New_index];		
		Color_Unit Color_diff;
		 
		Color_diff.c0 = pVertex->color_int(0) - col.c0;
		Color_diff.c1 = pVertex->color_int(1) - col.c1;
		Color_diff.c2 = pVertex->color_int(2) - col.c2;
		
		if (Color_diff.c0 < Color_diff_min_c0)
			Color_diff_min_c0 = Color_diff.c0;		
		if (Color_diff.c1 < Color_diff_min_c1)
			Color_diff_min_c1 = Color_diff.c1;
		if (Color_diff.c2 < Color_diff_min_c2)
			Color_diff_min_c2 = Color_diff.c2;

		if (Color_diff.c0 > Color_diff_max_c0)
			Color_diff_max_c0 = Color_diff.c0;
		if (Color_diff.c1 > Color_diff_max_c1)
			Color_diff_max_c1 = Color_diff.c1;
		if (Color_diff.c2 > Color_diff_max_c2)
			Color_diff_max_c2 = Color_diff.c2;

		this->DifferenceColor.push_back(Color_diff); // difference between original and quantized color.

		pVertex->color_int(col.c0, col.c1, col.c2);		
		
		New_vertex_color[0] = this->C0_Min + (col.c0 * this->Color_Quantization_Step);
		New_vertex_color[1] = this->C1_Min + (col.c1 * this->Color_Quantization_Step);
		New_vertex_color[2] = this->C2_Min + (col.c2 * this->Color_Quantization_Step);
		
		pVertex->color_float(New_vertex_color[0], New_vertex_color[1], New_vertex_color[2]);		
		
		LAB_To_RGB(New_vertex_color[0], New_vertex_color[1], New_vertex_color[2], Reconstructed_color);			
		pVertex->color(Reconstructed_color[0], Reconstructed_color[1], Reconstructed_color[2]);		
	}
	
	// For reconstruction of colors. (At the end of the decompression stage),
	this->ColorDiffMinC0 = Color_diff_min_c0;
	this->ColorDiffMinC1 = Color_diff_min_c1;
	this->ColorDiffMinC2 = Color_diff_min_c2;

	this->ColorDiffRangeC0 = Color_diff_max_c0 - Color_diff_min_c0 + 1;
	this->ColorDiffRangeC1 = Color_diff_max_c1 - Color_diff_min_c1 + 1;
	this->ColorDiffRangeC2 = Color_diff_max_c2 - Color_diff_min_c2 + 1;	
	#endif
}
#endif





// Description : This function select a set of independent vertices to be removed
int Compression_Valence_Component::Decimation_Conquest(Polyhedron  & pMesh,
										   const bool    Normal_flipping,
										   const bool    Use_metric,
										   const float & Metric_thread,
										   const bool    Use_forget_metric,
										   const int   & Forget_value,
										   const int   & Component_ID)
{	
	
	#ifdef MAPPING_TABLE_METHOD
	vector<int> Verify_used_symbol;
	for (int i = 0; i < NUMBER_SEEDS; i++)
		Verify_used_symbol.push_back(-1);
	#endif


	// Initialize vertex and face flags.
	Init(pMesh);
		
	// to count number of independent vertices and number of symbol of connectivity.
	int Number_vertices = 0; 
	int Number_symbol = 0;	
	
	// To find first edge.
	Halfedge_iterator hi = pMesh.halfedges_begin();			
	
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
				  (Is_Border_Vertex(h->next()) == false))
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
					if (pMesh.size_of_vertices() > (unsigned)Forget_value)
						Geometric_metric_condition = false;
					else
						Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
				}
				else
					Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
			}					
			
			// remove the front vertex if its removal does not viloate the manifold property and some metrics			
			bool Check_all_condition = false;

			if ((!Manifold_condition) && (!Geometric_metric_condition) && (!Normal_flipping_condition))
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

					#ifdef USE_LEE_METHOD
						Average_color = Get_Average_Vertex_Color_Lee(g, valence);
					#endif

					#ifdef USE_SIMPLE_PREDICTION
						Average_color = Get_Average_Vertex_Color_Before_Removal(g, valence);
					#endif

					#ifdef USE_YOON_METHOD
						Average_color = Get_Average_Vertex_Color_Youn(g, valence);
					#endif					

					// Color difference from average color of neighbors					
					Color_Unit Color_diff = Removed_vertex_color - Average_color;					
					
					#ifdef PREDICTION_METHOD
					this->InterVertexColor.push_front(Color_diff);
					#endif
					

					#ifdef MAPPING_TABLE_METHOD

					int Vertex_color_index = g->next()->vertex()->Vertex_Color_Index;					
					// Color index is stored
					this->InterColorIndex.push_front(Vertex_color_index);

					if (Verify_used_symbol[Vertex_color_index] == -1)
					{
						Verify_used_symbol[Vertex_color_index] = 1;
						this->PredictedColorArray[Vertex_color_index] = Color_diff;						
					}
					#endif	
				}			
								

				// Enter symbol 'VALENCE CODE' into the list of symbols
				this->InterConnectivity.push_front(valence - 3);
				
				// Calculate the position of barycenter.		
				Point3d Barycenter = Barycenter_Patch_Before_Removal(g);								
				
				Point_Int BC = Change_Real_Int(Barycenter, Component_ID);

				// remove the front vertex
				pMesh.erase_center_vertex(g->next()); 

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
												
				Retriangulation(pMesh, g, valence, 0);
				
				Vector normal = Normal_Patch(g, valence);
				Vector T2 = CGAL::NULL_VECTOR;
				Vector T1 = Calculate_T1_T2(h, normal, T2);

				if (normal == CGAL::NULL_VECTOR)
				{					
					T1 = Vector(1,0,0);
					T2 = Vector(0,1,0);
					normal = Vector(0,0,1);
				}					
				
				Point_Int Dist = CRV - BC;

				//Bijection
				Point_Int Frenet_Coordinates = Frenet_Rotation(Dist, T1, T2,normal);				
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
			
			if (Is_Manifold_Property_Violated(g, type, valence))
				Check_border_structure = false;			
			
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
			
				Point3d Real_vertex_position = h->next()->vertex()->point();
				Point_Int Vertex_position = Change_Real_Int(Real_vertex_position, Component_ID);				
				
				Color_Unit Removed_vertex_color;
				int Vertex_color_index = -1;
				if ((this->IsColored) && (!this->IsOneColor))
				{	
					Removed_vertex_color = Get_Vertex_Color(h->next());
					
					#ifdef MAPPING_TABLE_METHOD
					Vertex_color_index = h->next()->vertex()->Vertex_Color_Index;
					this->InterColorIndex.push_front(Vertex_color_index);					
					#endif
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
					pMesh.erase_facet(Temp);
				}
								
				if (valence == 3)
				{
					Halfedge_handle Retriangulation_edge = Border_edges[valence - 2]->opposite()->prev();
					
					// One triangle has to be created to smooth the mesh boundary					
					pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());					
										
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

						#ifdef PREDICTION_METHOD
						this->InterVertexColor.push_front(Color_diff);
						#endif

						#ifdef MAPPING_TABLE_METHOD						
						if (Verify_used_symbol[Vertex_color_index] == -1)
						{
							this->PredictedColorArray[Vertex_color_index] = Color_diff;						
							Verify_used_symbol[Vertex_color_index] = 1;
						}						
						#endif						
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
					Point_Int Frenet_Coordinates = Frenet_Rotation(Dist,T1,T2,normal);
					
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
						pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[1]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag = CONQUERED;						

						pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());
						Border_edges[0]->opposite()->facet()->Facet_Flag = CONQUERED;						
					}

					else if (   ( (Number_jump == 0) && ((type == 6) || (type == 7)) ) ||
								( (Number_jump == 1) && ((type == 5) || (type == 8)) ) ||
								( (Number_jump == 2) && ((type == 6) || (type == 7)) )  )
					{												
						pMesh.add_facet_to_border(Border_edges[2]->opposite(), Border_edges[0]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag = CONQUERED;						
						Halfedge_handle Temp_border = Border_edges[2]->opposite()->next();
						
						pMesh.add_facet_to_border(Retriangulation_edge, Temp_border);
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
						
						#ifdef PREDICTION_METHOD
						this->InterVertexColor.push_front(Color_diff);
						#endif

						#ifdef MAPPING_TABLE_METHOD
						if (Verify_used_symbol[Vertex_color_index] == -1)
						{
							Verify_used_symbol[Vertex_color_index] = 1;
							this->PredictedColorArray[Vertex_color_index] = Color_diff;
						}
						#endif											
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
					Point_Int Frenet_Coordinates = Frenet_Rotation(Dist,T1,T2,normal);
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
		#ifdef PREDICTION_METHOD
		while(!this->InterVertexColor.empty())
		{
			Color_Unit Col = this->InterVertexColor.front();
			this->InterVertexColor.pop_front();
			this->VertexColor[Component_ID].push_front(Col);
		}
		#endif
		#ifdef MAPPING_TABLE_METHOD
		while(!this->InterColorIndex.empty())
		{
			int Color_index = this->InterColorIndex.front();
			this->InterColorIndex.pop_front();
			this->ColorIndex.push_front(Color_index);
		}
		#endif
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
	

	#ifdef SAVE_INTERMEDIATE_MESHES
	// Save intermediate meshes.
	wxString Outputfile = "output";
	Outputfile += wxString::Format("%d",this->GlobalCountOperation*2);
	Outputfile += ".off";
	pMesh.write_off(Outputfile.ToAscii(), true, false);		
	#endif

	return Number_vertices;
}


// Description : Regulation conquest
int Compression_Valence_Component::Regulation(Polyhedron  & pMesh,
								  const bool    Normal_flipping,
								  const bool    Use_metric,
								  const float & Metric_thread,
								  const bool    Use_forget_metric,
								  const int   & Forget_value,
								  const int   & Component_ID)
{    
	#ifdef MAPPING_TABLE_METHOD
	vector<int> Verify_used_symbol;
	for (int i = 0; i < NUMBER_SEEDS; i++)
		Verify_used_symbol.push_back(-1);
	#endif


	// Initialization
	Init(pMesh);
	
	
	
	// Number of removed vertices and number of connectivity symbols
	int Number_vertices = 0;
	int Number_symbol = 0;

	

	Halfedge_iterator hi = pMesh.halfedges_begin();
	
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
		int Component_ID = h->next()->vertex()->Component_Number;

		if ((h->facet()->Facet_Flag == CONQUERED) || (h->facet()->Facet_Flag == TO_BE_REMOVED))
			continue;

		else if ((h->next()->vertex()->Vertex_Flag == FREE) && (valence == 3) && (Is_Border_Vertex(h->next()) == false)) // if valence is 3, remove the front vertex.
		{
			Halfedge_handle g = h;
			int type = 1; // ant type of valence 3

			// Check if the manifold property is violated.
			bool Manifold_condition = Is_Manifold_Property_Violated(h, type, valence);			
           
			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Geometric_metric_condition = false;
			if (Use_metric == true)
			{
				if (Use_forget_metric == true)
				{
					if (pMesh.size_of_vertices() > (unsigned)Forget_value)
						Geometric_metric_condition = false;
					else
						Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
				}

				else
					Geometric_metric_condition = Is_Geometric_Metric_Violated(h, type, valence, Metric_thread);
			}			
			// remove the front vertex if its removal does not viloate the manifold property and some metrics			
			bool Check_all_condition = false;

			if ((Manifold_condition == false) && (Geometric_metric_condition == false))
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
					#ifdef USE_LEE_METHOD
						Average_color = Get_Average_Vertex_Color_Lee(g, valence);
					#endif
					#ifdef USE_YOON_METHOD
						Average_color = Get_Average_Vertex_Color_Youn(g, valence);
					#endif
					#ifdef USE_SIMPLE_PREDICTION
						Average_color = Get_Average_Vertex_Color_Before_Removal(g, valence);
					#endif
					
					Color_Unit Color_diff = Removed_vertex_color - Average_color;					

					#ifdef PREDICTION_METHOD
					this->InterVertexColor.push_front(Color_diff);
					#endif

					#ifdef MAPPING_TABLE_METHOD

					int Vertex_color_index = h->next()->vertex()->Vertex_Color_Index;					
					this->InterColorIndex.push_front(Vertex_color_index);

					if (Verify_used_symbol[Vertex_color_index] == -1)
					{
						Verify_used_symbol[Vertex_color_index] = 1;
						this->PredictedColorArray[Vertex_color_index] = Color_diff;
					}
					#endif					
				}


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
				
				Point_Int Dist = Geo - BC;					
				
				Point_Int Frenet_Coordinates = Frenet_Rotation(Dist,T1,T2,normal);
				
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
	for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); )
	{
		Vertex_handle vh = pVertex;
		pVertex++;
		Halfedge_handle del = vh->halfedge();

		if (vh->Vertex_Flag == TO_BE_REMOVED)
			pMesh.erase_center_vertex(del);
	}	


	if ((this->IsColored) && (!this->IsOneColor))
	{
		#ifdef PREDICTION_METHOD
		while(!this->InterVertexColor.empty())
		{
			Color_Unit Col = this->InterVertexColor.front();
			this->InterVertexColor.pop_front();
			this->VertexColor[Component_ID].push_front(Col);
		}
		#endif
		#ifdef MAPPING_TABLE_METHOD
		while(!this->InterColorIndex.empty())
		{
			int Color_index = this->InterColorIndex.front();
			this->InterColorIndex.pop_front();
			this->ColorIndex.push_front(Color_index);
		}
		#endif		
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
	
	
	#ifdef SAVE_INTERMEDIATE_MESHES
	// Save intermediate meshes.
	Outputfile = "output";
	Outputfile += wxString::Format("%d",this->GlobalCountOperation*2+1);
	Outputfile += ".off";
	pMesh.write_off(Outputfile.ToAscii(), true, false);		
	#endif

	return Number_vertices;
}



// Description : Decoding of the regulation conquest
void Compression_Valence_Component::Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID)
{	
	Init(pMesh);

	Adaptive_Data_Model Connectivity(2);
	
	Halfedge_iterator hi = pMesh.halfedges_begin();
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

			Point_Int Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);

			Point_Int Center = BC + Diff;			
			
			
			Point3d Center_vertex = this->Change_Int_Real(Center,Component_ID);			
			
			// Vertex insertion
			g = pMesh.create_center_vertex(g);
			
			g->vertex()->point() = Center_vertex;
			
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
				#ifdef PREDICTION_METHOD
				g = h;				

				Color_Unit Predicted_color;

				#ifdef USE_LEE_METHOD
					Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
				#endif
				#ifdef USE_YOON_METHOD
					Predicted_color = Get_Average_Vertex_Color_Youn(g, valence);
				#endif
				#ifdef USE_SIMPLE_PREDICTION
					Predicted_color = Get_Average_Vertex_Color_Before_Removal(g, valence);
				#endif												
				
				Color_Unit Color_difference;
				Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
				Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
				Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;								
				
				Color_Unit CV = Predicted_color + Color_difference;				

				g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);

				float LAB[3];
				LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
				LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
				LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

				float RGB[3];
				LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

				g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);			
				#endif


				#ifdef MAPPING_TABLE_METHOD
			
				g = h;
				Color_Unit CV;						

				int Color_index = this->Decoder.decode(this->Index_Model);
				g->next()->vertex()->Vertex_Color_Index = Color_index;

				if (this->IsKnownIndex[Color_index] == 1)
				{
					CV = this->ColorArray[Color_index];
				}

				// if this is a new color
				if (this->IsKnownIndex[Color_index] == -1)
				{
					this->IsKnownIndex[Color_index] = 1;

					Color_Unit Predicted_color;

					#ifdef USE_LEE_METHOD
						Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
					#endif
					#ifdef USE_YOON_METHOD
						Predicted_color = Get_Average_Vertex_Color_Youn(g, valence);
					#endif
					#ifdef USE_SIMPLE_PREDICTION
						Predicted_color = Get_Average_Vertex_Color_Before_Removal(g, valence);
					#endif

					Color_Unit Color_difference;
					Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
					Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
					Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
					
					CV = Predicted_color + Color_difference;
					//CV = Color_difference;					
					this->ColorArray[Color_index] = CV;
				}				
				
				g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);				

				float LAB[3];
				LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
				LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
				LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

				float RGB[3];
				LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

				g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);		

				#endif
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
void Compression_Valence_Component::Un_Decimation_Conquest(Polyhedron       & pMesh, 
											   Arithmetic_Codec & Decoder,
											   const int        & Component_ID)
{
	Init(pMesh);
	
	int Number_connectivity_symbols;
	if (this->IsClosed[Component_ID])
		Number_connectivity_symbols = 5;
	else
		Number_connectivity_symbols = 7;

	Adaptive_Data_Model Connectivity(Number_connectivity_symbols);	
	
	Halfedge_iterator hi = pMesh.halfedges_begin();
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
			
			if (normal == CGAL::NULL_VECTOR)
			{				
				T1 = Vector(1,0,0);
				T2 = Vector(0,1,0);					
				normal = Vector(0,0,1);	
			}						

			// remove edges to re_find polygon and (check. and attribute sign flag.)
			bool Check_Validity = Remove_Edges(pMesh, g, type);

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

				Point_Int Diff = Inverse_Frenet_Rotation(Frenet, T1, T2, normal);

				Point_Int Center = BC + Diff;				
				
				Point3d Center_vertex = this->Change_Int_Real(Center, Component_ID);
				
				g = pMesh.create_center_vertex(g);
				g->vertex()->point() = Center_vertex;
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


				if ((this->IsColored) && (!this->IsOneColor))
				{				
					#ifdef PREDICTION_METHOD
					g = h;				

					Color_Unit Predicted_color;

					#ifdef USE_LEE_METHOD
						Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
					#endif
					#ifdef USE_YOON_METHOD
						Predicted_color = Get_Average_Vertex_Color_Youn(g, valence);
					#endif
					#ifdef USE_SIMPLE_PREDICTION
						Predicted_color = Get_Average_Vertex_Color_Before_Removal(g, valence);
					#endif												
					
					Color_Unit Color_difference;
					Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
					Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
					Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;								
					
					Color_Unit CV = Predicted_color + Color_difference;					

					g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);

					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);
				
					#endif

					#ifdef MAPPING_TABLE_METHOD
					g = h;
					Color_Unit CV;						

					int Color_index = this->Decoder.decode(this->Index_Model);
					g->next()->vertex()->Vertex_Color_Index = Color_index;
					
					if (this->IsKnownIndex[Color_index] == 1)
					{
						CV = this->ColorArray[Color_index];
					}

					// if this is a new color
					if (this->IsKnownIndex[Color_index] == -1)
					{
						this->IsKnownIndex[Color_index] = 1;

						Color_Unit Predicted_color;

						#ifdef USE_LEE_METHOD
							Predicted_color = Get_Average_Vertex_Color_Lee(g, valence);
						#endif
						#ifdef USE_YOON_METHOD
							Predicted_color = Get_Average_Vertex_Color_Youn(g, valence);
						#endif
						#ifdef USE_SIMPLE_PREDICTION
							Predicted_color = Get_Average_Vertex_Color_Before_Removal(g, valence);
						#endif

						Color_Unit Color_difference;
						Color_difference.c0 = this->Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
						Color_difference.c1 = this->Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
						Color_difference.c2 = this->Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
						
						CV = Predicted_color + Color_difference;						
						this->ColorArray[Color_index] = CV;
					}					

					g->next()->vertex()->color_int(CV.c0, CV.c1, CV.c2);				

					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);
					#endif
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

			Point_Int Diff = Inverse_Frenet_Rotation(Frenet,T1,T2,normal);			
			Halfedge_handle g = h;
			
			
			#ifdef PREDICTION_METHOD
			Color_Unit Predicted_color;
			if ((this->IsColored) && (!this->IsOneColor))
			{
				Predicted_color.c0 = Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
				Predicted_color.c1 = Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
				Predicted_color.c2 = Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
			}
			#endif			
			
			#ifdef MAPPING_TABLE_METHOD
			if ((this->IsColored) && (!this->IsOneColor))
			{
				bool Is_new_color = false;
				Color_Unit Predicted_color;
				int Color_index;

				Color_index = Decoder.decode(this->Index_Model);
				if (this->IsKnownIndex[Color_index] == -1)
				{
					this->IsKnownIndex[Color_index] = 1;
					Is_new_color = true;

					Predicted_color.c0 = Decoder.decode(this->Color_0_Model) + this->Smallest_C0;
					Predicted_color.c1 = Decoder.decode(this->Color_1_Model) + this->Smallest_C1;
					Predicted_color.c2 = Decoder.decode(this->Color_2_Model) + this->Smallest_C2;
				}
			}
			#endif

			// border edge with valence == 3
			if (valence == 8)
			{				
				Point3d Barycenter = Barycenter_Patch_After_Removal(pass, 3);				
				Point_Int BC = Change_Real_Int(Barycenter,Component_ID);

				Point_Int Center = BC + Diff;				
				Point3d Center_vertex = this->Change_Int_Real(Center,Component_ID);
				
				#ifdef PREDICTION_METHOD
				Color_Unit Average_color;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Average_color = Get_Average_Vertex_Color_After_Removal(g, 3);
				}
				#endif

				#ifdef MAPPING_TABLE_METHOD
				Color_Unit Average_color;

				if ((this->IsColored) && (!this->IsOneColor))
				{					
					Average_color = Get_Average_Vertex_Color_After_Removal(g, 3);
				}
				#endif

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

				pMesh.erase_facet(g);

				Halfedge_handle Prev_edge = Border_edges[1]->opposite()->prev();
					
				// g points the new vertex
				g = pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				g->vertex()->point() = Center_vertex;
				g->vertex()->Vertex_Flag = CONQUERED;

				#ifdef PREDICTION_METHOD
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
				#endif

				#ifdef MAPPING_TABLE_METHOD
				Color_Unit CV;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					if (Is_new_color)
					{
						CV = Average_color + Predicted_color;						
						this->ColorArray[Color_index] = CV;
					}
					else
					{
						CV = this->ColorArray[Color_index];
					}

					if (( CV.c0 > 255) || (CV.c1 > 255) || ( CV.c2 > 255))
						int temp = 0;

					g->vertex()->color_int(CV.c0, CV.c1, CV.c2);

					g->vertex()->Vertex_Color_Index = Color_index;
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->vertex()->color(RGB[0], RGB[1], RGB[2]);					
				}
				#endif
				
				if ((this->IsColored) && (this->IsOneColor))
				{
					g->vertex()->color(this->OnlyColor[0],this->OnlyColor[1],this->OnlyColor[2]);
				}


				Prev_edge = Prev_edge->next();
				pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());			
				
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

						pMesh.erase_facet(Border_edges[0]->opposite());
						pMesh.erase_facet(Border_edges[1]->opposite());
					}

					// jump == 1;
					else if (g->prev()->opposite()->next()->is_border_edge())
					{
						Number_jump = 1;
						
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite()->prev()->opposite());

						pMesh.erase_facet(Border_edges[2]->opposite());
						pMesh.erase_facet(Border_edges[0]->opposite());
					}
					
					// jump == 2;
					else
					{
						Number_jump = 2;

						Border_edges.push_back(g->prev()->opposite()->next()->opposite());
						Border_edges.push_back(g->next()->opposite());
						Border_edges.push_back(g->opposite());

						pMesh.erase_facet(Border_edges[0]->opposite());
						pMesh.erase_facet(Border_edges[1]->opposite());
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

						pMesh.erase_facet(Border_edges[2]->opposite());
						pMesh.erase_facet(Border_edges[1]->opposite());
					}

					else if (g->next()->opposite()->prev()->is_border_edge())
					{
						Number_jump = 1;

						Border_edges.push_back(g->next()->opposite()->next()->opposite());
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());

						pMesh.erase_facet(Border_edges[0]->opposite());
						pMesh.erase_facet(Border_edges[1]->opposite());
					}
					else
					{
						Number_jump = 0;
						
						Border_edges.push_back(g->opposite());
						Border_edges.push_back(g->prev()->opposite());
						Border_edges.push_back(g->next()->opposite()->prev()->opposite());

						pMesh.erase_facet(Border_edges[2]->opposite());
						pMesh.erase_facet(Border_edges[1]->opposite());
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
				g = pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[2]->opposite());
				g->vertex()->point() = Center_vertex;				
					
				#ifdef PREDICTION_METHOD
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
				#endif

				#ifdef MAPPING_TABLE_METHOD
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

					if (Is_new_color)
					{
						CV = Average_color + Predicted_color;
						this->ColorArray[Color_index] = CV;
					}
					else
					{
						CV = this->ColorArray[Color_index];
					}

					g->vertex()->color_int(CV.c0, CV.c1, CV.c2);
					g->vertex()->Vertex_Color_Index = Color_index;
					
					float LAB[3];
					LAB[0] = this->C0_Min + CV.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + CV.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + CV.c2 * this->Color_Quantization_Step;

					float RGB[3];
					LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);

					g->vertex()->color(RGB[0], RGB[1], RGB[2]);					
				}
				#endif
				if ((this->IsColored) && (this->IsOneColor))
				{
					g->vertex()->color(this->OnlyColor[0],this->OnlyColor[1],this->OnlyColor[2]);
				}
				Prev_edge = Prev_edge->next();
				pMesh.add_facet_to_border(Prev_edge, Border_edges[1]->opposite());
				pMesh.add_facet_to_border(Prev_edge, Border_edges[0]->opposite());


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



int Compression_Valence_Component::Test_Next_Operation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread,
											const bool Use_forget_metric,const int &Forget_value,const int &Qbit)
{
		
	return 0;


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double mesh_area = 0;
	//for (Facet_iterator pFacet = pMesh.facets_begin(); pFacet != pMesh.facets_end(); pFacet++)
	//{
	//	Halfedge_handle h = pFacet->halfedge();
	//	mesh_area += Area_Facet_Triangle(h);
	//}
	//int Number_vvv = pMesh.size_of_vertices();

	//double Real_volume = 0;
	//for (Facet_iterator pFacet = pMesh.facets_begin(); pFacet != pMesh.facets_end(); pFacet++)
	//{
	//	Halfedge_handle h = pFacet->halfedge();
	//	Vector V0 = h->vertex()->point() - CGAL::ORIGIN;
	//	Vector V1 = h->next()->vertex()->point() - CGAL::ORIGIN;
	//	Vector V2 = h->prev()->vertex()->point() - CGAL::ORIGIN;
	//	
	//	Vector Temp = CGAL::cross_product(V1, V2);
	//	Real_volume += V0 * Temp;
	//}
	//Real_volume /= 6.0;
	//
	//double C = Real_volume / mesh_area / (double)Number_vvv;

	////////////////////////////////////////////////////////////////////////////////////////////////////////


	//this->CountOperation++;
	//
	//bool Is_decimation = true;
	//bool Is_limit_quantization = true;

	//double Deci_mrms = 0, Deci_mrmswrtBB = 0;
	//double Deci_hausdorff=0, Deci_hausdorffwrtBB = 0;	
	//
	//double Under_mrms = 0, Under_mrmswrtBB = 0;
	//double Under_hausdorff=0, Under_hausdorffwrtBB = 0;
	//
	//unsigned Under_number_vertices = 0;
	//unsigned Under_cost_bits = 0;
	//unsigned Decimation_cost_bits = 0;
	//int total_number_vertices = 0;

	//int Operation_choice = -1;		
	//
	//Polyhedron * Under_mesh = new Polyhedron;
	//Polyhedron * Deci_mesh = new Polyhedron;
	//
	//Copy_Polyhedron.copy(&(pMesh),Deci_mesh);
	//Attibute_Seed_Gate_Flag(pMesh,*Deci_mesh);	
	//
	//int v1 = this->Decimation_Conquest(*Deci_mesh,Normal_flipping,Use_metric,Metric_thread,Use_forget_metric,Forget_value);
	//int v2 = this->Regulation(*Deci_mesh,Normal_flipping,Use_metric,Metric_thread,Use_forget_metric,Forget_value);		
	//total_number_vertices = v1 + v2;		
	//Deci_mesh->write_off("DecimatedMesh.off", false, false);

	//
	//if (this->Qbit > LIMIT_QBIT)
	//{	
	//	Is_limit_quantization = false;		

	//	Copy_Polyhedron.copy(&(pMesh),Under_mesh);
	//	Attibute_Seed_Gate_Flag(pMesh,*Under_mesh);
	//	this->Under_Quantization(*Under_mesh);	
	//	Under_mesh->compute_normals();
	//	Under_mesh->write_off("UnderquantizedMesh.off", false, false);

	//	//calculate distorsion.
	//	this->Calculate_Distances("original_temp.off", "DecimatedMesh.off", Deci_mrms, Deci_mrmswrtBB, Deci_hausdorff, Deci_hausdorffwrtBB);	
	//	this->Calculate_Distances("original_temp.off","UnderquantizedMesh.off", Under_mrms, Under_mrmswrtBB, Under_hausdorff, Under_hausdorffwrtBB);

	//	
	//	
	//	////////    Calculates bits needs for decimation /////////
	//	Arithmetic_Codec TCodec(1024*300);	
	//	TCodec.start_encoder();	
	//
	//	list<int>::iterator it_number_vertices = this->NumberVertices.begin();	
	//	list<int>::iterator it_con = this->Connectivity.begin();		
	//	list<Point_Int>::iterator it_geo = this->Geometry.begin();


	//	int NV_regulation = *it_number_vertices;
	//	it_number_vertices++;
	//	int NV_decimation = *it_number_vertices;
	//	
	//	int amax, amin;
	//	int gmax, gmin;

	//	amax = -5000; amin = 5000;gmax = -5000;gmin = 5000;

	//	for (int i = 0; i< NV_regulation; i++)
	//	{
	//		Point_Int Temp = *it_geo;
	//		it_geo++;

	//		if (Temp.x > amax) amax = Temp.x;
	//		if (Temp.y > amax) amax = Temp.y;

	//		if (Temp.x < amin) amin = Temp.x;
	//		if (Temp.y < amin) amin = Temp.y;

	//		if (Temp.z > gmax) gmax = Temp.z;
	//		if (Temp.z < gmin) gmin = Temp.z;
	//	}
	//	
	//	int alpha_range1 = amax - amin + 1;
	//	int alpha_offset1 = -amin;
	//	
	//	int gamma_range1 = gmax - gmin + 1;
	//	int gamma_offset1 = - gmin;

	//	if (NV_regulation == 0)
	//	{
	//		alpha_range1 = 0;
	//		alpha_offset1 = 0;
	//		gamma_range1 = 0;
	//		gamma_offset1 = 0;
	//	}

	//	amax = -5000;amin = 5000;gmax = -5000;gmin = 5000;

	//	for (int i = 0; i< NV_decimation; i++)
	//	{
	//		Point_Int Temp = *it_geo;
	//		it_geo++;

	//		if (Temp.x > amax) amax = Temp.x;
	//		if (Temp.y > amax) amax = Temp.y;

	//		if (Temp.x < amin) amin = Temp.x;
	//		if (Temp.y < amin) amin = Temp.y;

	//		if (Temp.z > gmax) gmax = Temp.z;
	//		if (Temp.z < gmin) gmin = Temp.z;
	//	}
	//	
	//	int alpha_range2 = amax - amin + 1;
	//	int alpha_offset2 = -amin;
	//	
	//	int gamma_range2 = gmax - gmin + 1;
	//	int gamma_offset2 = -gmin;
	//	
	//	if (NV_decimation == 0)
	//	{
	//		alpha_range2 = 0;
	//		alpha_offset2 = 0;

	//		gamma_range2 = 0;
	//		gamma_offset2 = 0;
	//	}
	//	
	//	it_geo = this->Geometry.begin();

	//	if ((alpha_range1 == 0) || (alpha_range1 == 1)) alpha_range1 = 2;
	//	if ((alpha_range2 == 0) || (alpha_range2 == 1)) alpha_range2 = 2;
	//	if ((gamma_range1 == 0) || (gamma_range1 == 1)) gamma_range1 = 2;
	//	if ((gamma_range2 == 0) || (gamma_range2 == 1)) gamma_range2 = 2;
	//	
	//	Adaptive_Data_Model Connectivity1(2);
	//	Adaptive_Data_Model alpha1(alpha_range1); //regulation
	//	Adaptive_Data_Model gamma1(gamma_range1); // regulation
	//	
	//	int count_vertex = 0;
	//	int count_connec = 0;
	//	//regulation
	//	for (int i = 0; i<DumpSymbolRegulation;i++)
	//	{
	//		int symbol = *it_con;
	//		it_con++;
	//		TCodec.encode(symbol,Connectivity1);
	//		count_connec++;
	//		if (symbol != 1) // if the symbol is not a null code
	//		{
	//			count_vertex++;
	//			Point_Int Coeff = *it_geo;
	//			it_geo++;
	//			
	//			int a = Coeff.x + alpha_offset1;
	//			int b = Coeff.y + alpha_offset1;
	//			int g = Coeff.z + gamma_offset1;
	//			if (v2 != 0)
	//			{
	//				TCodec.encode(a,alpha1);
	//				TCodec.encode(b,alpha1);
	//				TCodec.encode(g,gamma1);	
	//			}					
	//		}		
	//	}
	//	Adaptive_Data_Model Connectivity2(this->NummberConnectivitySymbols);
	//	Adaptive_Data_Model alpha2(alpha_range2); //Decimation
	//	Adaptive_Data_Model gamma2(gamma_range2); //
	//	
	//	for (int i = 0; i<DumpSymbolDecimation;i++)
	//	{	
	//		count_connec++;
	//		int symbol = *it_con;
	//		it_con++;
	//		TCodec.encode(symbol,Connectivity2);
	//		
	//		if (symbol != 4) // if the symbol is not a null code
	//		{
	//			Point_Int Coeff = *it_geo;
	//			it_geo++;
	//			count_vertex++;
	//			int a = Coeff.x + alpha_offset2;
	//			int b = Coeff.y + alpha_offset2;
	//			int g = Coeff.z + gamma_offset2;
	//			if (v1 != 0)
	//			{
	//				TCodec.encode(a,alpha2);
	//				TCodec.encode(b,alpha2);
	//				TCodec.encode(g,gamma2);	
	//			}					
	//		}		
	//	}	
	//	unsigned Decimation_cost_bytes = TCodec.stop_encoder();	

	//	//change to bits and sum bits of the alpha, gamma range and offset information.
	//	unsigned Decimation_cost_bits = Decimation_cost_bytes * 8 + Qbit * 8;


	//	// bits need for under_quantization
	//	Under_number_vertices = pMesh.size_of_vertices();
	//	Under_cost_bits = 0;

	//	if (!Is_limit_quantization)
	//	{
	//		list<int>::iterator it = this->QuantizationCorrectVector.begin();
	//		
	//		Arithmetic_Codec TCodec2(1024*300);	
	//		TCodec2.start_encoder();	
	//		
	//		Adaptive_Data_Model Under_symbols(8);

	//		for (unsigned i = 0; i < Under_number_vertices; i++)
	//		{
	//			int symbol = *it;
	//			it++;
	//			TCodec2.encode(symbol,Under_symbols);
	//		}
	//		Under_cost_bits = TCodec2.stop_encoder() * 8;
	//	}			
	//	
	//	// calculate slopes to make a decision;
	//	double Under_slope = (Under_mrmswrtBB - this->OldDistortion) / (double)Under_cost_bits;
	//	double Deci_slope = (Deci_mrmswrtBB - this->OldDistortion) / (double)Decimation_cost_bits;

	//	if (Under_slope <= Deci_slope)
	//		Is_decimation = false;
	//}	

	//// Decimation is chosen to be applied.	
	//if (Is_decimation)
	//{		
	//	Operation_choice = 0;
	//	
	//	// Increase number of decimation.
	//	this->NumberDecimation++;

	//	// Change the actual mesh
	//	pMesh.clear();
	//	Copy_Polyhedron.copy(Deci_mesh, &(pMesh));
	//	Attibute_Seed_Gate_Flag(*Deci_mesh, pMesh);	
	//	pMesh.compute_normals();
	//	
	//	this->TotalBits += Decimation_cost_bits;
	//	//fprintf(this->LogFile,"%d \t %d \t %d \t %d \t %lf \t %lf \t %lf\n", Operation_choice, Number_vvv, pMesh.size_of_vertices(), this->Qbit, C, Real_volume, mesh_area);		
	//				
	//	//fprintf(this->RD_MRMS,"%d \t %f \t %f \n", this->CountOperation, Deci_mrms, Deci_mrmswrtBB);
	//	//fprintf(this->RD_HAUSDORFF,"%d \t %f \t %f \n", this->CountOperation, Deci_hausdorff, Deci_hausdorffwrtBB);
	//	
	//	this->OldDistortion = Deci_mrmswrtBB;
	//	if (!Is_limit_quantization)
	//	{
	//		this->NumberQuantizationLayer.pop_front();
	//		for (unsigned i = 0; i < Under_number_vertices; i++)
	//			this->QuantizationCorrectVector.pop_front();
	//	}
	//}
	//
	//// Under_quantization is chosen to be applied.
	//else
	//{
	//	Operation_choice = 1;

	//	this->Qbit--;
	//	this->NumberChangeQuantization++;

	//	pMesh.clear();
	//	Copy_Polyhedron.copy(Under_mesh,&(pMesh));
	//	Attibute_Seed_Gate_Flag(*Under_mesh,pMesh);
	//	pMesh.compute_normals();
	//	
	//	this->OldDistortion = Under_mrmswrtBB;
	//	
	//	// remove all symbols of decimation
	//	for (int i = 0; i < (this->DumpSymbolDecimation + this->DumpSymbolRegulation); i++)
	//		this->Connectivity.pop_front();
	//	
	//	for (int i = 0; i < total_number_vertices; i++)
	//		this->Geometry.pop_front();
	//	
	//	for (int i = 0; i < 2; i++)
	//	{				
	//		this->NumberSymbol.pop_front();
	//		this->NumberVertices.pop_front();			
	//	}


	//	this->TotalBits += Under_cost_bits;
	//	//fprintf(this->LogFile,"%d \t %d \t %d \t %d \t %lf \t %lf \t %lf\n", Operation_choice, Number_vvv, pMesh.size_of_vertices(), this->Qbit, C, Real_volume, mesh_area);		
	//	//fprintf(this->RD_MRMS,"%d\t%f\t%f\n", this->CountOperation, Under_mrms, Under_mrmswrtBB);
	//	//fprintf(this->RD_HAUSDORFF,"%d\t%f\t%f\n", this->CountOperation, Under_hausdorff, Under_hausdorffwrtBB);		
	//}
	//#ifdef SAVE_INTERMEDIATE_MESHES

	//// Save intermediate meshes.
	//wxString Outputfile = "output";
	//Outputfile += wxString::Format("%d",this->CountOperation);
	//Outputfile += ".off";
	//pMesh.write_off(Outputfile.ToAscii(), false, false);
	//#endif

	//delete Deci_mesh;
	//delete Under_mesh;

	//return Operation_choice;
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
void Compression_Valence_Component::Write_Base_Mesh(Polyhedron & pMesh, Arithmetic_Codec & enc, unsigned &Connectivity_size, unsigned & Color_size, const int & Num_color_base_mesh)
{

	unsigned int Max_Qbit = 0;
	for (int i =0; i< this->NumberComponents; i++)
	{
		enc.put_bits(this->ComponentOperations[i], 8);	 					  // number of operations < 256	
		enc.put_bits(this->Qbit[i]-4, 4); // quantization bit < 16
		enc.put_bits(this->NumberChangeQuantization[i], 4);		// number of decrease of quantization resolution

		

		if (this->Qbit[i] > Max_Qbit)
			Max_Qbit = this->Qbit[i];
	}
		
	enc.put_bits(pMesh.size_of_vertices(), 15);						  // number of vertices of base mesh < 4096
	enc.put_bits(pMesh.size_of_facets(), 16);						// number of facets of base mesh < 8192	

	int Base_color_index_bit = 0;

	#ifdef MAPPING_TABLE_METHOD
	Base_color_index_bit = (int)ceil(log((double)(Num_color_base_mesh+1))/log((double)2));	
	enc.put_bits(Base_color_index_bit, 12);	
	
	this->IsKnownIndex.clear();
	for (int i = 0; i < this->PredictedColorArray.size(); i++)
		this->IsKnownIndex.push_back(-1);
	
	#endif

	int Basemesh_vertex_number = 0;
	vector<int> Seed_edges(2 * this->NumberComponents, -1);

		
	// Encoding of vertex information of base mesh //
	for (Vertex_iterator pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); Basemesh_vertex_number++, pVertex++)
	{			
		pVertex->Vertex_Number = Basemesh_vertex_number;
		
		if (pVertex->Seed_Edge != OTHER_COORDINATE)
			Seed_edges[pVertex->Seed_Edge] = Basemesh_vertex_number;		
		
		int cid = pVertex->Component_Number;

		Point_Int Vertex = Change_Real_Int(pVertex->point(), pVertex->Component_Number);
		
		
		enc.put_bits(Vertex.x, Max_Qbit + 1);
		enc.put_bits(Vertex.y, Max_Qbit + 1);
		enc.put_bits(Vertex.z, Max_Qbit + 1);		
		
		if ((this->IsColored) && (!this->IsOneColor))
		{
			#ifdef PREDICTION_METHOD
			int C0 = pVertex->color_int(0);
			int C1 = pVertex->color_int(1);
			int C2 = pVertex->color_int(2);
			
			enc.put_bits(C0, C0_QUANTIZATION);
			enc.put_bits(C1, C1_QUANTIZATION);
			enc.put_bits(C2, C2_QUANTIZATION);
					
			Color_size += 3 * C0_QUANTIZATION;
			#endif

			#ifdef MAPPING_TABLE_METHOD

			// encode color index of vertex
			int Old_color_index = pVertex->Vertex_Color_Index;
			int New_color_index = this->ReorderingColorIndex[Old_color_index];
			enc.put_bits(New_color_index, Base_color_index_bit);				
			
			Color_size += Base_color_index_bit;
			
			// To know if it is a new color 
			if (this->IsKnownIndex[New_color_index] == -1)
			{
				this->IsKnownIndex[New_color_index] = 1;
				
				// If it's a new color, encode 3 color coordinates;
				int C0 = pVertex->color_int(0);
				int C1 = pVertex->color_int(1);
				int C2 = pVertex->color_int(2);
				
				enc.put_bits(C0, C0_QUANTIZATION);
				enc.put_bits(C1, C1_QUANTIZATION);
				enc.put_bits(C2, C2_QUANTIZATION);					
				
				Color_size += C0_QUANTIZATION + C1_QUANTIZATION + C2_QUANTIZATION;
			}			
			#endif
		}		
	}
	
	// Bits needed for each edge.
	int Facet_index_bit = (int)ceil(log((double)(pMesh.size_of_vertices()+1))/log((double)2));
	
	int Count_facet_index = 0;	 
	for (Facet_iterator pFacet = pMesh.facets_begin() ; pFacet != pMesh.facets_end() ; pFacet++)
	{		
		Halfedge_handle pHalfedge = pFacet->halfedge();				
		do
		{
			enc.put_bits(pHalfedge->vertex()->Vertex_Number, Facet_index_bit);			
			pHalfedge = pHalfedge->next();
			Count_facet_index++;

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
			#ifdef PREDICTION_METHOD
			list<Color_Unit>::iterator Vertex_color_iterator;
			for (Vertex_color_iterator = this->VertexColor[Component_ID].begin(); Vertex_color_iterator != this->VertexColor[Component_ID].end(); Vertex_color_iterator++)
			#endif
			#ifdef MAPPING_TABLE_METHOD
			vector<Color_Unit>::iterator Vertex_color_iterator;
			for (Vertex_color_iterator = this->PredictedColorArray[Component_ID].begin(); Vertex_color_iterator != this->PredictedColorArray[Component_ID].end(); Vertex_color_iterator++)
			#endif			
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


void Compression_Valence_Component::Color_Metric_Roy(Polyhedron &pMesh, double &min, double &max,double &mean, double &rms)
{
	Mesh_roy * Simplified = new Mesh_roy;
	Mesh_roy * Temp = new Mesh_roy;	
	
	Vector3d Pos, Color;
	Vector3i Face;

	for (int i = 0; i < this->Original->VertexNumber(); i++)
	{
		Pos = this->Original->Vertex(i);
		Temp->AddVertex(Pos);
		Color = this->Original->Color(i);
		Temp->AddColor(Color);
	}
	for (int i = 0; i < this->Original->FaceNumber(); i++)
	{
		Face = this->Original->Face(i);
		Temp->AddFace(Face);
	}


	for (Vertex_iterator pVert = pMesh.vertices_begin(); pVert != pMesh.vertices_end(); pVert++)
	{
		Point3d pt = pVert->point();

		Pos[0] = pt.x();
		Pos[1] = pt.y();
		Pos[2] = pt.z();

		Color[0] = pVert->color(0);
		Color[1] = pVert->color(1);
		Color[2] = pVert->color(2);

		Simplified->AddVertex(Pos);
		Simplified->AddColor(Color);
	}
	
	pMesh.set_index_vertices();

	for (Facet_iterator pFacet = pMesh.facets_begin(); pFacet != pMesh.facets_end(); pFacet++)
	{
		int count = 0;

		Halfedge_around_facet_circulator pH = pFacet->facet_begin();
		do
		{
			Face[count] = pH->vertex()->tag();
			count++;
		}while(++pH != pFacet->facet_begin());
		
		Simplified->AddFace(Face);
	}

	
	Deviation * dev = new Deviation;

	dev->Initialization(Temp, Simplified, 0, 0.5);
	dev->SetDeviationColorBound(0);

	bool IsOK = dev->Compute(COLOR_DEVIATION);

	min = dev->Min();
	max = dev->Max();
	mean = dev->Mean();
	rms = dev->Rms();	
	
	delete dev;
	delete Simplified;
	delete Temp;
}


void Compression_Valence_Component::Simplification(Polyhedron  & pMesh,
									   const int   & NVertices, 
									   const bool    Normal_flipping,
									   const bool    Use_metric,
									   const float & Metric_thread, 
									   const bool    Use_forget_metric,
									   const int   & Forget_value)
{			
	#ifdef MEASURE_COLOR_DEVIATION	

	FILE * temp = fopen("abc.txt","w");

	double Color_min, Color_max, Color_mean, Color_rms;
	if (this->IsColored)
	{
		this->Color_Metric_Roy(pMesh, Color_min, Color_max, Color_mean, Color_rms);
		//fprintf(this->LogColor," Q \t %d \t %lf \t %lf \t %lf \t %lf\n", pMesh.size_of_vertices(), Color_min, Color_max, Color_mean, Color_rms);
		fprintf(temp," Q \t %d \t %lf \t %lf \t %lf \t %lf\n", pMesh.size_of_vertices(), Color_min, Color_max, Color_mean, Color_rms);
		fclose(temp);

		//Write_SMF(pMesh, "Quantized.smf", true);
		//this->Calculate_Sqrt_Colors("Original.smf", "Quantized.smf", Color_min, Color_max, Color_mean);
		//fprintf(LogColor," Q \t %d \t %lf \t %lf \t %lf \n", pMesh.size_of_vertices(), Color_min, Color_max, Color_mean);
	}

	#endif

	bool Is_any_vertex_removed = true;
	
	unsigned Last_Number = 0;
	unsigned Current_Number = pMesh.size_of_vertices();

	int Operation_choice = -1;
	
	do
	{		
		Last_Number = Current_Number;
		
		// Simplify component by component
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{
			// if the ith component did not remove any vertex in last loop, it is not necessary to simplifiy it.
			if (this->ComponentOperations[Component_ID] == this->GlobalCountOperation)
			{
				unsigned Initial_number_vertices = pMesh.size_of_vertices();

				this->Decimation_Conquest(pMesh, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value, Component_ID);

				

				this->Regulation(pMesh, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value, Component_ID);


				int Diff_number_vertices = pMesh.size_of_vertices() - Initial_number_vertices;
				
				this->ComponentOperations[Component_ID] += 1;
				this->NumberDecimation[Component_ID] += 1;
				
				this->ListOperation[Component_ID].push_front(0);	
				
				if (Diff_number_vertices == 0)	
					this->Remove_Last_Phase_Elements(Component_ID);
			}
		}			
		
		#ifdef MEASURE_COLOR_DEVIATION
		if (this->IsColored)
		{	
			if (pMesh.size_of_vertices() > 30)
			{
				this->Color_Metric_Roy(pMesh,Color_min, Color_max, Color_mean, Color_rms);
				
				temp = fopen("abc.txt","a");
				//fprintf(this->LogColor,"%2d \t %d \t %lf \t %lf \t %lf \t %lf\n", this->CountOperation+1, pMesh.size_of_vertices(), Color_min, Color_max, Color_mean, Color_rms);		
				fprintf(temp,"%2d \t %d \t %lf \t %lf \t %lf \t %lf\n", this->GlobalCountOperation+1, pMesh.size_of_vertices(), Color_min, Color_max, Color_mean, Color_rms);		
				fclose(temp);

				//Write_SMF(pMesh, "Inter.smf", true);			
				//this->Calculate_Sqrt_Colors("Original.smf", "Inter.smf", Color_min, Color_max, Color_mean);
				//fprintf(LogColor,"%2d \t %d \t %lf \t %lf \t %lf \n", this->CountOperation+1, pMesh.size_of_vertices(), Color_min, Color_max, Color_mean);
			}
		}
		#endif

		Current_Number = pMesh.size_of_vertices();
		if (Current_Number != Last_Number)
			this->GlobalCountOperation++;

		if (Current_Number < (unsigned)NVertices) // MT
			break;
		
	}while((Current_Number != Last_Number));
	
	pMesh.compute_normals();	

#ifdef MEASURE_COLOR_DEVIATION
	fclose(this->LogColor);
#endif
}

void Compression_Valence_Component::Compression(Polyhedron     & pMesh, 
									const char     * File_Name, 
									const int      & _Qbit, 
									unsigned       & Connectivity_size, 
									unsigned       & Color_size, 
									unsigned       & Total_size, 
									const unsigned & Initial_file_size)
{
	// Calculate offset and range for the compression.
	this->Calculate_Geometry_Color_Offset_Range();	
	
	FILE * fp  = fopen(File_Name, "wb");												//Main FILE to save compression information.		
		
	fwrite(&this->Smallest_Alpha, sizeof(int), 1, fp);				  	 // smallest value of alpha (to save the negative value)
	fwrite(&this->Smallest_Gamma, sizeof(int), 1, fp);				 	 // smallest value of gamma (to save the negative value)
	fwrite(&Initial_file_size, sizeof(unsigned), 1, fp);    // Intial size of the input file (To visualize during decompression)	
	fwrite(&this->NumberComponents, sizeof(int), 1, fp);
	
	for (int i =0; i<this->NumberComponents; i++)
	{
		fwrite(&this->QuantizationPas[i], sizeof(float), 1, fp);							    // QuantizationPas(step of quantization)
		fwrite(&this->xmin[i], sizeof(float), 1, fp);																	   // xmin value
		fwrite(&this->ymin[i], sizeof(float), 1, fp);																	   // ymin value
		fwrite(&this->zmin[i], sizeof(float), 1, fp);																	   // zmin value
	}


#ifdef SPLIT_FILE_EACH_RESOLUTION

	std::string S_file_name(File_Name);
	size_t Pos_last_point = S_file_name.find_last_of('.');
	
	string Temp_split_file_name = S_file_name.substr(0, Pos_last_point);//,"_split.");
	Temp_split_file_name += "_split.";

	string Base_file_name = Temp_split_file_name + "p3d";


	FILE * fp_split  = fopen(Base_file_name.c_str(), "wb");												//Main FILE to save compression information.		
	
	fwrite(&this->Smallest_Alpha, sizeof(int), 1, fp_split);				  	 // smallest value of alpha (to save the negative value)
	fwrite(&this->Smallest_Gamma, sizeof(int), 1, fp_split);				 	 // smallest value of gamma (to save the negative value)
	fwrite(&Initial_file_size, sizeof(unsigned), 1, fp_split);    // Intial size of the input file (To visualize during decompression)	
	fwrite(&this->NumberComponents, sizeof(int), 1, fp_split);
	
	for (int i =0; i<this->NumberComponents; i++)
	{
		fwrite(&this->QuantizationPas[i], sizeof(float), 1, fp_split);							    // QuantizationPas(step of quantization)
		fwrite(&this->xmin[i], sizeof(float), 1, fp_split);																	   // xmin value
		fwrite(&this->ymin[i], sizeof(float), 1, fp_split);																	   // ymin value
		fwrite(&this->zmin[i], sizeof(float), 1, fp_split);
	}
#endif
	
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
	
	fwrite(&Colored, sizeof(char), 1, fp);
	if (this->IsColored)
		fwrite(&OneColor, sizeof(char), 1, fp);


#ifdef SPLIT_FILE_EACH_RESOLUTION
																		   
	fwrite(&Colored, sizeof(char), 1, fp_split);
	if (this->IsColored)
		fwrite(&OneColor, sizeof(char), 1, fp_split);

#endif


	int Num_color_base_mesh = 0;
	if ((this->IsColored) && (!this->IsOneColor))
	{
		#ifdef MAPPING_TABLE_METHOD
		// Reorganize the color index to follow the order of appearance
		Num_color_base_mesh = this->Mapping_Table_Index_Reordering(pMesh);		
		#endif				

		fwrite(&this->Color_Quantization_Step, sizeof(float), 1, fp);

		// En-tete pour la couleur
		fwrite(&this->C0_Min, sizeof(float), 1, fp); // smallest value of c0 
		fwrite(&this->C1_Min, sizeof(float), 1, fp); // smallest value of c1 
		fwrite(&this->C2_Min, sizeof(float), 1, fp); // smallest value of c2		
			
		fwrite(&this->Smallest_C0, sizeof(int), 1, fp);  
		fwrite(&this->Smallest_C1, sizeof(int), 1, fp); 
		fwrite(&this->Smallest_C2, sizeof(int), 1, fp);
		
		Color_size += sizeof(int) * 8 * 9;
		
		#ifdef SPLIT_FILE_EACH_RESOLUTION
		// En-tete pour la couleur
		fwrite(&this->Color_Quantization_Step, sizeof(float), 1, fp_split);

		fwrite(&this->C0_Min, sizeof(float), 1, fp_split); // smallest value of c0 
		fwrite(&this->C1_Min, sizeof(float), 1, fp_split); // smallest value of c1 
		fwrite(&this->C2_Min, sizeof(float), 1, fp_split); // smallest value of c2				
			
		fwrite(&this->Smallest_C0, sizeof(int), 1, fp_split);  
		fwrite(&this->Smallest_C1, sizeof(int), 1, fp_split); 
		fwrite(&this->Smallest_C2, sizeof(int), 1, fp_split);
		#endif
		
		#ifdef MAPPING_TABLE_METHOD

		// To encode the color difference of color quantization.
		fwrite(&this->ColorDiffMinC0, sizeof(int), 1, fp);
		fwrite(&this->ColorDiffMinC1, sizeof(int), 1, fp);
		fwrite(&this->ColorDiffMinC2, sizeof(int), 1, fp);
		Color_size += sizeof(int) * 8 * 3;

		#ifdef SPLIT_FILE_EACH_RESOLUTION

		fwrite(&this->ColorDiffMinC0, sizeof(int), 1, fp_split);
		fwrite(&this->ColorDiffMinC1, sizeof(int), 1, fp_split);
		fwrite(&this->ColorDiffMinC2, sizeof(int), 1, fp_split);

		#endif

		#endif			
	}	
	
	if ((this->IsColored) && (this->IsOneColor))
	{

		fwrite(&this->OnlyColor[0], sizeof(float), 1, fp); // smallest value of c0 
		fwrite(&this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1 
		fwrite(&this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2

		#ifdef SPLIT_FILE_EACH_RESOLUTION

		fwrite(&this->OnlyColor[0], sizeof(float), 1, fp_split); // smallest value of c0 
		fwrite(&this->OnlyColor[1], sizeof(float), 1, fp_split); // smallest value of c1 
		fwrite(&this->OnlyColor[2], sizeof(float), 1, fp_split); // smallest value of c2
		#endif
	}

	// Declaration du codeur.
	Arithmetic_Codec enc(AC_BUFFER); 
	enc.start_encoder(); 	

	// To calculate connectivity rate
	Arithmetic_Codec Connectivity_encoder(AC_BUFFER); 
	Connectivity_encoder.start_encoder(); 

	#ifdef SPLIT_FILE_EACH_RESOLUTION
	// To split compressed file.
	Arithmetic_Codec Split_encoder(AC_BUFFER);
	Split_encoder.start_encoder();
	#endif
		
	for (int i = 0; i < this->NumberComponents; i++)
	{
		if (this->IsClosed[i])
			enc.put_bits(0,1);			
					
		else
			enc.put_bits(1,1);			

		#ifdef SPLIT_FILE_EACH_RESOLUTION
		if (this->IsClosed[i])
			Split_encoder.put_bits(0,1);			
					
		else
			Split_encoder.put_bits(1,1);			
		#endif		
	}	

	/*	Write information of base mesh.
		geometry + connectivity + color information(if the mesh is colored). */
	this->Write_Base_Mesh(pMesh, enc, Connectivity_size, Color_size, Num_color_base_mesh);

	#ifdef SPLIT_FILE_EACH_RESOLUTION	
	unsigned int tx = 0, ty = 0, tz = 0;
	this->Write_Base_Mesh(pMesh, Split_encoder, tx, ty, tz);
	#endif

	// To calculate color rate
	Arithmetic_Codec Color_enc(AC_BUFFER); 
	Color_enc.start_encoder();


	#ifdef PREDICTION_METHOD

	Adaptive_Data_Model C0_model;
	Adaptive_Data_Model C1_model;
	Adaptive_Data_Model C2_model;
	
	//To calculate color rate	
	Adaptive_Data_Model PC_C0_model;
	Adaptive_Data_Model PC_C1_model;
	Adaptive_Data_Model PC_C2_model;
	
	#ifdef SPLIT_FILE_EACH_RESOLUTION

	Adaptive_Data_Model Split_C0_model;
	Adaptive_Data_Model Split_C1_model;
	Adaptive_Data_Model Split_C2_model;
	#endif

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

		#ifdef SPLIT_FILE_EACH_RESOLUTION
		Split_encoder.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
		Split_encoder.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
		Split_encoder.put_bits(this->C2_Range, C2_QUANTIZATION + 1);		

		Split_C0_model.set_alphabet(this->C0_Range);
		Split_C1_model.set_alphabet(this->C1_Range);
		Split_C2_model.set_alphabet(this->C2_Range);

		#endif

	}
	#endif		

	#ifdef MAPPING_TABLE_METHOD

	Adaptive_Data_Model Index_model;
	Adaptive_Data_Model C0_model;
	Adaptive_Data_Model C1_model;
	Adaptive_Data_Model C2_model;


	//To calculate color rate
	Adaptive_Data_Model Temp_index_model;

	Adaptive_Data_Model Temp_C0_model;
	Adaptive_Data_Model Temp_C1_model;
	Adaptive_Data_Model Temp_C2_model;

	#ifdef SPLIT_FILE_EACH_RESOLUTION

	Adaptive_Data_Model Split_index_model;	
	Adaptive_Data_Model Split_C0_model;
	Adaptive_Data_Model Split_C1_model;
	Adaptive_Data_Model Split_C2_model;
	#endif

	if ((this->IsColored) && (!this->IsOneColor))
	{
		enc.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
		enc.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
		enc.put_bits(this->C2_Range, C2_QUANTIZATION + 1);

		Color_size += 3 * (C0_QUANTIZATION + 1);		
		
		Index_model.set_alphabet(NUMBER_SEEDS);
		C0_model.set_alphabet(this->C0_Range);
		C1_model.set_alphabet(this->C1_Range);
		C2_model.set_alphabet(this->C2_Range);
		
		Temp_index_model.set_alphabet(NUMBER_SEEDS);
		Temp_C0_model.set_alphabet(this->C0_Range);
		Temp_C1_model.set_alphabet(this->C1_Range);
		Temp_C2_model.set_alphabet(this->C2_Range);

		#ifdef SPLIT_FILE_EACH_RESOLUTION
		Split_encoder.put_bits(this->C0_Range, C0_QUANTIZATION + 1);
		Split_encoder.put_bits(this->C1_Range, C1_QUANTIZATION + 1);
		Split_encoder.put_bits(this->C2_Range, C2_QUANTIZATION + 1);

		Split_index_model.set_alphabet(NUMBER_SEEDS);
		Split_C0_model.set_alphabet(this->C0_Range);
		Split_C1_model.set_alphabet(this->C1_Range);
		Split_C2_model.set_alphabet(this->C2_Range);

		#endif
	}
	#endif
	
	#ifdef SPLIT_FILE_EACH_RESOLUTION
	Split_encoder.write_to_file(fp_split);
	fclose(fp_split);
	#endif
	
	//char buffer[30];	// MT

	// Main loop of compression //
	for (int i = 0; i < this->GlobalCountOperation; i++)
	{
		#ifdef SPLIT_FILE_EACH_RESOLUTION
		Split_encoder.start_encoder();
		#endif

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
					enc.put_bits(0, 1);
					
					#ifdef SPLIT_FILE_EACH_RESOLUTION
					Split_encoder.put_bits(0,1);
					#endif
					
					for (int j = 0; j < 2 ; j++) //Decimation and regulation
					{
						Adaptive_Data_Model Connectivity;				
						
						if (j == 0)	Connectivity.set_alphabet(2);
						else		Connectivity.set_alphabet(Number_connectivity_symbols);
							
						Adaptive_Data_Model Temp_connectivity;
						if (j == 0)	Temp_connectivity.set_alphabet(2);
						else		Temp_connectivity.set_alphabet(Number_connectivity_symbols);

						#ifdef SPLIT_FILE_EACH_RESOLUTION
						Adaptive_Data_Model Split_connectivity;
						if (j == 0)	Split_connectivity.set_alphabet(2);
						else		Split_connectivity.set_alphabet(Number_connectivity_symbols);
						#endif
						
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
						
						#ifdef SPLIT_FILE_EACH_RESOLUTION
						Split_encoder.put_bits(Alpha_range, _Qbit+1);
						if (this->Smallest_Alpha < 0) 	Split_encoder.put_bits(Alpha_offset - this->Smallest_Alpha, _Qbit+1);
						else							Split_encoder.put_bits(Alpha_offset, _Qbit+1);

						Split_encoder.put_bits(Gamma_range, _Qbit+1);				
						if (this->Smallest_Gamma < 0)	Split_encoder.put_bits(Gamma_offset - this->Smallest_Gamma, _Qbit+1);
						else							Split_encoder.put_bits(Gamma_offset, _Qbit+1);				

						#endif

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

						#ifdef SPLIT_FILE_EACH_RESOLUTION

						Adaptive_Data_Model Split_alpha(Alpha_range);
						Adaptive_Data_Model Split_gamma(Gamma_range);

						#endif

						int Number_symbols = this->NumberSymbol[Component_ID].front();
						this->NumberSymbol[Component_ID].pop_front();

						for (unsigned k = 0; k < (unsigned)Number_symbols; k++)	// MT
						{
							unsigned symbol = this->Connectivity[Component_ID].front();
							this->Connectivity[Component_ID].pop_front();			

							enc.encode(symbol, Connectivity);

							// To calculare connectivity rate
							Connectivity_encoder.encode(symbol, Temp_connectivity);

							#ifdef SPLIT_FILE_EACH_RESOLUTION
							Split_encoder.encode(symbol, Split_connectivity);
							#endif
							
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
								
								#ifdef SPLIT_FILE_EACH_RESOLUTION
								if (check_alpha == false)
								{						
									Split_encoder.encode(x, Split_alpha);
									Split_encoder.encode(y, Split_alpha);						
								}
								if (check_gamma == false)
								{	
									Split_encoder.encode(z, Split_gamma);						
								}

								#endif
								
								#ifdef PREDICTION_METHOD
								if ((this->IsColored) && (!this->IsOneColor))
								{
									Color_Unit VC = this->VertexColor[Component_ID].front();
									this->VertexColor[Component_ID].pop_front();

									enc.encode(VC.c0 - this->Smallest_C0, C0_model);
									enc.encode(VC.c1 - this->Smallest_C1, C1_model);
									enc.encode(VC.c2 - this->Smallest_C2, C2_model);


									#ifdef SPLIT_FILE_EACH_RESOLUTION
									
									Split_encoder.encode(VC.c0 - this->Smallest_C0, Split_C0_model);
									Split_encoder.encode(VC.c1 - this->Smallest_C1, Split_C1_model);
									Split_encoder.encode(VC.c2 - this->Smallest_C2, Split_C2_model);
									
									#endif

									// To calculate color rate
									Color_enc.encode(VC.c0 - this->Smallest_C0, PC_C0_model);
									Color_enc.encode(VC.c1 - this->Smallest_C1, PC_C1_model);
									Color_enc.encode(VC.c2 - this->Smallest_C2, PC_C2_model);
								}
								#endif

								#ifdef MAPPING_TABLE_METHOD
								if ((this->IsColored) && (!this->IsOneColor))
								{
									int Old_index = this->ColorIndex.front();
									this->ColorIndex.pop_front();
									
									int New_index = this->ReorderingColorIndex[Old_index];

									enc.encode(New_index, Index_model);
									#ifdef SPLIT_FILE_EACH_RESOLUTION
									Split_encoder.encode(New_index, Split_index_model);
									#endif
									
									// To calculate color rate
									Color_enc.encode(New_index, Temp_index_model);

									if (this->IsKnownIndex[New_index] == -1)
									{
										this->IsKnownIndex[New_index] = 1;

										Color_Unit VC = this->PredictedColorArray[Old_index];
										enc.encode(VC.c0 - this->Smallest_C0, C0_model);
										enc.encode(VC.c1 - this->Smallest_C1, C1_model);
										enc.encode(VC.c2 - this->Smallest_C2, C2_model);

										#ifdef SPLIT_FILE_EACH_RESOLUTION
										Split_encoder.encode(VC.c0 - this->Smallest_C0, Split_C0_model);
										Split_encoder.encode(VC.c1 - this->Smallest_C1, Split_C1_model);
										Split_encoder.encode(VC.c2 - this->Smallest_C2, Split_C2_model);

										#endif
										
										// To calculate color rate
										Color_enc.encode(VC.c0 - this->Smallest_C0, Temp_C0_model);
										Color_enc.encode(VC.c1 - this->Smallest_C1, Temp_C1_model);
										Color_enc.encode(VC.c2 - this->Smallest_C2, Temp_C2_model);
									}
								}							
								#endif
							}
						}
						alpha.reset();
						gamma.reset();

						#ifdef SPLIT_FILE_EACH_RESOLUTION
						Split_alpha.reset();
						Split_gamma.reset();
						#endif
					}

				}
		
				// Decrease of quantization resolution.		
				else 
				{
					enc.put_bits(1,1);

					#ifdef SPLIT_FILE_EACH_RESOLUTION			
					
					Split_encoder.put_bits(1,1);
					Adaptive_Data_Model Split_under_quantization_model(8);
					#endif

					int Number_vertices = this->NumberQuantizationLayer[Component_ID].front();
					this->NumberQuantizationLayer[Component_ID].pop_front();
					
					Adaptive_Data_Model Under_quantization_model(8);
					
					for (int i = 0; i < Number_vertices; i++)
					{
						int Under_quantization_coeff = this->QuantizationCorrectVector[Component_ID].front();
						this->QuantizationCorrectVector[Component_ID].pop_front();

						enc.encode(Under_quantization_coeff, Under_quantization_model);

						#ifdef SPLIT_FILE_EACH_RESOLUTION
						Split_encoder.encode(Under_quantization_coeff, Split_under_quantization_model);
						#endif
					}		
				}
				
				
			}			
		}

		#ifdef SPLIT_FILE_EACH_RESOLUTION

		const char * Number = itoa(i, buffer, 10);
		string Number_iteration(Number);
		string New_file_name = Temp_split_file_name + "ps";
		New_file_name += Number_iteration;

		FILE * fp_split_iter = fopen(New_file_name.c_str(), "wb");

		Split_encoder.write_to_file(fp_split_iter);
		fclose(fp_split_iter);

		#endif


	}			

	#ifdef MAPPING_TABLE_METHOD		
	
	Adaptive_Data_Model Restorer_C0_model;
	Adaptive_Data_Model Restorer_C1_model;
	Adaptive_Data_Model Restorer_C2_model;

	Restorer_C0_model.set_alphabet(this->ColorDiffRangeC0);
	Restorer_C1_model.set_alphabet(this->ColorDiffRangeC1);
	Restorer_C2_model.set_alphabet(this->ColorDiffRangeC2);

	enc.put_bits(this->ColorDiffRangeC0, C0_QUANTIZATION + 1);
	enc.put_bits(this->ColorDiffRangeC1, C1_QUANTIZATION + 1);
	enc.put_bits(this->ColorDiffRangeC2, C2_QUANTIZATION + 1);


	

	// To calculate color rate
	Adaptive_Data_Model Temp_restorer_C0_model;
	Adaptive_Data_Model Temp_restorer_C1_model;
	Adaptive_Data_Model Temp_restorer_C2_model;	

	Temp_restorer_C0_model.set_alphabet(this->ColorDiffRangeC0);
	Temp_restorer_C1_model.set_alphabet(this->ColorDiffRangeC1);
	Temp_restorer_C2_model.set_alphabet(this->ColorDiffRangeC2);

	Color_enc.put_bits(this->ColorDiffRangeC0, C0_QUANTIZATION + 1);
	Color_enc.put_bits(this->ColorDiffRangeC1, C1_QUANTIZATION + 1);
	Color_enc.put_bits(this->ColorDiffRangeC2, C2_QUANTIZATION + 1);
	
	Color_size += (C0_QUANTIZATION + 1)*3;
	
	for (int i = 0; i < this->DifferenceColor.size(); i++)
	{
		Color_Unit Temp = this->DifferenceColor[i];

		enc.encode(Temp.c0 - this->ColorDiffMinC0, Restorer_C0_model);
		enc.encode(Temp.c1 - this->ColorDiffMinC1, Restorer_C1_model);
		enc.encode(Temp.c2 - this->ColorDiffMinC2, Restorer_C2_model);
		

		// To calculate color rate
		Color_enc.encode(Temp.c0 - this->ColorDiffMinC0, Temp_restorer_C0_model);
		Color_enc.encode(Temp.c1 - this->ColorDiffMinC1, Temp_restorer_C1_model);
		Color_enc.encode(Temp.c2 - this->ColorDiffMinC2, Temp_restorer_C2_model);
	}	
	#endif
	
	

	Connectivity_size += Connectivity_encoder.stop_encoder() * 8;
	Color_size += Color_enc.stop_encoder() * 8;

	enc.write_to_file(fp);	
	fclose(fp);

	FILE* f_size = fopen(File_Name,"rb");
	fseek(f_size, 0, SEEK_END);
	Total_size = ftell(f_size);
}

void Compression_Valence_Component::Adaptive_Quantization_Preparation(Polyhedron &pMesh, double & mrms, double &mrmswrtBB,double &hausdorff, double &hausdorffwrtBB)
{
	pMesh.write_off("quantized_temp.off",false,false);
		
	this->Calculate_Distances((char *)"original_temp.off",(char *)"quantized_temp.off",mrms,mrmswrtBB,hausdorff,hausdorffwrtBB);
	this->OldDistortion = mrmswrtBB; // reference distorsion (induced by quantization)
	
	/*fprintf(this->LogFile,"Input mesh : %d vertices \n",pMesh.size_of_vertices());
	fprintf(this->LogFile,"Quantization used : %d bits\n\n",Qbit);
	fprintf(this->LogFile,"******* Initial Quantization ******* \n");
	fprintf(this->LogFile,"MRMS : %f (%f wrt bounding box diagonal)\n",mrms,mrmswrtBB);
	fprintf(this->LogFile,"Hausdorff : %f (%f wrt bounding box diagonal)\n\n\n",hausdorff,hausdorffwrtBB);
		
	fprintf(this->RD_MRMS,"%d\t%f\t%f\n", this->CountOperation, mrms, mrmswrtBB);
	fprintf(this->RD_HAUSDORFF,"%d\t%f\t%f\n", this->CountOperation, hausdorff, hausdorffwrtBB);*/
}



void Compression_Valence_Component::Adaptive_Quantization(Polyhedron &pMesh, const int & NVertices, const bool Normal_flipping,const bool Use_metric,const float &Metric_thread,
											const bool Use_forget_metric,const int &Forget_value,const int &Qbit)
{


#ifdef USE_ESTIMATION_ADAPTIVE_QUANTIZATION

	//bool Is_any_vertex_removed = true;
	unsigned Last_Number = 0;
	unsigned Current_Number = pMesh.size_of_vertices();			
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


				int Q = Calulate_Proper_Quantization_Precision(pMesh, this->ComponentVolume[Component_ID],
																	  this->ComponentArea[Component_ID],
																	  this->ComponentNumberVertices[Component_ID]);

				if (Q < 4)
					Q = 4;
		
				if (Q < (int)this->Qbit[Component_ID])	// MT
				{
					Continue = true;
					Operation_choice = 1;
				
					this->Under_Quantization(pMesh, Component_ID);

					this->Qbit[Component_ID]--;
					this->NumberChangeQuantization[Component_ID]++;			

					this->ComponentOperations[Component_ID]++;
					this->ListOperation[Component_ID].push_front(Operation_choice);
				}		

				else
				{
					
					Operation_choice = 0;

					unsigned Initial_number_vertices = pMesh.size_of_vertices();
					
					this->Decimation_Conquest(pMesh, Normal_flipping, Use_metric,Metric_thread, Use_forget_metric, Forget_value,Component_ID);
					this->Regulation(pMesh, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value,Component_ID);

					this->NumberDecimation[Component_ID] += 1;

					unsigned Diff_number_vertices = pMesh.size_of_vertices() - Initial_number_vertices;
					
					this->ComponentOperations[Component_ID] += 1;
					this->ListOperation[Component_ID].push_front(Operation_choice);	
					
					this->ComponentNumberVertices[Component_ID] += Diff_number_vertices;

					if (Diff_number_vertices == 0)
						this->Remove_Last_Phase_Elements(Component_ID);

					else
						Continue = true;
						
					
				}				
			}
		}

		Current_Number = pMesh.size_of_vertices();	
		
		if (Continue)
			this->GlobalCountOperation++;
		
		#ifdef SAVE_INTERMEDIATE_MESHES
		// Save intermediate meshes.
		wxString Outputfile = "output";
		Outputfile += wxString::Format("%d",this->GlobalCountOperation);
		Outputfile += ".off";
		pMesh.write_off(Outputfile.ToAscii(), true, false);		
		#endif
		
		if (Current_Number < (unsigned)NVertices)	// MT
			break;
		
	}while((Current_Number != Last_Number) || (Continue));
	
	pMesh.compute_normals();


					
		

			

	
#else
	/////////////////////////////////////Original version///////////////////////////////////////////////////////////

	// Variable for adaptative quantization
	double mrms = 0, mrmswrtBB = 0;
	double hausdorff = 0, hausdorffwrtBB = 0;	

	// Calcul de distance avec quantization.
	this->Adaptive_Quantization_Preparation(pMesh, mrms, mrmswrtBB, hausdorff, hausdorffwrtBB);
	bool Is_any_vertex_removed = true;
	
	unsigned Last_Number = 0;
	unsigned Current_Number = 0;			
	int Operation_choice = -1;
	
	do
	{		
		Last_Number = Current_Number;		
		
		Operation_choice = this->Test_Next_Operation(pMesh, Normal_flipping, Use_metric, Metric_thread, Use_forget_metric, Forget_value, Qbit);
		this->ListOperation.push_front(Operation_choice);		
		
		// If no vertex has been decimated, stop the loop.
		Current_Number = pMesh.size_of_vertices();
		
		if (Current_Number <= (unsigned)NVertices)
		{
			Is_any_vertex_removed = false;
			break;
		}
		
	}while((Operation_choice == 1)||(Current_Number != Last_Number));
	
	pMesh.compute_normals();	

		// if last loop didn't remove any vertex, pop all element of the last loop.
	if (Is_any_vertex_removed)	
		this->Remove_Last_Phase_Elements();	

	remove("original_temp.off");
	remove("quantized_temp.off");
	remove("UnderquantizedMesh.off");
	remove("DecimatedMesh.off");
#endif
}



void Compression_Valence_Component::Up_Quantization(Polyhedron &pMesh,Arithmetic_Codec & Decoder, const int & Component_ID)
{
	Adaptive_Data_Model Under_quantization_model(8);	

	// Premier_Pas == distance d'un QuantizationPas de grill de quantification (Q)
	double Small_step = 0.0;
			
	if (this->NumberChangeQuantization[Component_ID] == 0)	
		Small_step = this->QuantizationPas[Component_ID];		
	
	else
		Small_step = this->QuantizationPas[Component_ID] * pow(2.0, this->NumberChangeQuantization[Component_ID] - 1);		
	
	//double Large_step = Small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = pMesh.halfedges_begin();			
	while((hi->vertex()->Seed_Edge != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_ID +1))
		hi++;	

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = pMesh.vertices_begin();pVert != pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;		
		pVert->Vertex_Number = -1;
		pVert->Component_Number = -1;
	}	
	//pMesh.compute_normals();	
	
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
			for (it_reorder = Priority_map.begin(); it_reorder != Priority_map.end(); it_reorder++)
			{
				Priority_reorder[it_reorder->second] = 7 - Temp_priority;
				Temp_priority++;
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
		}
	}
	Vertex_iterator pVertex = NULL;
	for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
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
void Compression_Valence_Component::Under_Quantization(Polyhedron &pMesh, const int & Component_ID)
{		

	// stock three mins for the reconstruction.
	float _xmin = this->xmin[Component_ID];
	float _ymin = this->ymin[Component_ID];
	float _zmin = this->zmin[Component_ID];

	// Premier_Pas == distance d'un QuantizationPas de grill de quantification (Q)
	double Small_step = 0.0;
			
	if (this->NumberChangeQuantization[Component_ID] == 0)		
		Small_step = this->QuantizationPas[Component_ID];		
	
	else
		Small_step = this->QuantizationPas[Component_ID] * pow(2.0, this->NumberChangeQuantization[Component_ID]);		
		
	// Large_step == distance d'un QuantizationPas de grille de quantification(Q - 1)
	double Large_step = Small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge != 2*Component_ID+1))
		hi++;

	// Vertex_Flag est donnee free a tous les sommets
	for (Vertex_iterator pVert = pMesh.vertices_begin();pVert != pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag = FREE;
		pVert->Vertex_Number = -1;
	}

	
	//premiere passe pour sous quantifie et donne l'indice de symbol a chaque sommet
	for (Vertex_iterator pVert = pMesh.vertices_begin(); pVert != pMesh.vertices_end(); pVert++)
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
	
	/////// appliquer calcul des normales?? ou QuantizationPas???	
	//pMesh.compute_normals();	

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
			for (it_reorder = Priority_map.begin(); it_reorder != Priority_map.end(); it_reorder++)
			{
				Priority_reorder[it_reorder->second] = 7 - Temp_priority;
				Temp_priority++;
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

void Compression_Valence_Component::Set_Original_Color_Mesh(Polyhedron &pMesh, const char * Filename)
{
	this->Original = new Mesh_roy;
	
	string File = Filename;
	size_t dot = File.find_last_of('.');
	File.erase(dot,4);
	File += "_color_distortion.txt";
	this->LogColor = fopen(File.c_str(), "w");

	Vector3d Pos, Color;
	Vector3i Face;
	for (Vertex_iterator pVert = pMesh.vertices_begin(); pVert != pMesh.vertices_end(); pVert++)
	{
		Point3d pt = pVert->point();

		Pos[0] = pt.x();
		Pos[1] = pt.y();
		Pos[2] = pt.z();

		Color[0] = pVert->color(0);
		Color[1] = pVert->color(1);
		Color[2] = pVert->color(2);

		this->Original->AddVertex(Pos);
		this->Original->AddColor(Color);
	}
	
	pMesh.set_index_vertices();

	for (Facet_iterator pFacet = pMesh.facets_begin(); pFacet != pMesh.facets_end(); pFacet++)
	{
		int count = 0;

		Halfedge_around_facet_circulator pH = pFacet->facet_begin();
		do
		{
			Face[count] = pH->vertex()->tag();
			count++;
		}while(++pH != pFacet->facet_begin());
		
		this->Original->AddFace(Face);
	}
	

}




int Compression_Valence_Component::Decompress_Init(Polyhedron &pMesh, unsigned & Initial_file_size, const char* File_Name)
{	
	pMesh.clear();
	
	this->Decompress_count = 0;	
	
	FILE* fp = fopen(File_Name,"rb");

	int res;		
			
	res=fread(&this->Smallest_Alpha, sizeof(int), 1, fp);
	res=fread(&this->Smallest_Gamma, sizeof(int), 1, fp);
	res=fread(&Initial_file_size, sizeof(unsigned), 1, fp);	 
	res=fread(&this->NumberComponents, sizeof(int), 1, fp);
	
	float Qpas;
	float t_xmin, t_ymin, t_zmin;

	for (int i =0; i < this->NumberComponents; i++)
	{
		res=fread(&Qpas, sizeof(float), 1, fp);
		res=fread(&t_xmin, sizeof(float), 1, fp);
		res=fread(&t_ymin, sizeof(float), 1, fp);
		res=fread(&t_zmin, sizeof(float), 1, fp);

		this->QuantizationPas.push_back(Qpas);
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

	if ((this->IsColored) && (!this->IsOneColor))
	{	
		res=fread(&this->Color_Quantization_Step, sizeof(float), 1, fp);	
		
		res=fread(&this->C0_Min, sizeof(float), 1, fp); // smallest value of c0 
		res=fread(&this->C1_Min, sizeof(float), 1, fp); // smallest value of c1 
		res=fread(&this->C2_Min, sizeof(float), 1, fp); // smallest value of c2			
			
		res=fread(&this->Smallest_C0, sizeof(int), 1, fp);  
		res=fread(&this->Smallest_C1, sizeof(int), 1, fp); 
		res=fread(&this->Smallest_C2, sizeof(int), 1, fp);	

		#ifdef MAPPING_TABLE_METHOD		
		res=fread(&this->ColorDiffMinC0, sizeof(int), 1, fp);
		res=fread(&this->ColorDiffMinC1, sizeof(int), 1, fp);
		res=fread(&this->ColorDiffMinC2, sizeof(int), 1, fp);
		#endif
	}
	
	if ((this->IsColored) && (this->IsOneColor))
	{
		res=fread(&this->OnlyColor[0], sizeof(float), 1, fp); // smallest value of c0 
		res=fread(&this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1 
		res=fread(&this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2
	}
	this->Decoder.set_buffer(AC_BUFFER);
	this->Decoder.read_from_file(fp);
	
	
	for (int i = 0; i < this->NumberComponents; i++)
	{
		if (Decoder.get_bits(1) == 0)
			this->IsClosed.push_back(true);					
		else
			this->IsClosed.push_back(false);		
	}
	
	this->GlobalCountOperation = -1;
	
	unsigned Max_Qbit = 0;
	for (int i = 0; i < this->NumberComponents; i++)
	{
		int Number_operation = Decoder.get_bits(8);
		this->ComponentOperations.push_back(Number_operation);
		
		if (Number_operation > this->GlobalCountOperation)
			this->GlobalCountOperation = Number_operation;
			
		int qbit = Decoder.get_bits(4);
		qbit += 4;

		this->Qbit.push_back(qbit);
		int NCQ = Decoder.get_bits(4);

		this->NumberChangeQuantization.push_back(NCQ);

		if (this->Qbit[i] > Max_Qbit)
			Max_Qbit = this->Qbit[i];
	}
	
	int Number_basemesh_vertex = Decoder.get_bits(15);
	int Number_basemesh_facet  = Decoder.get_bits(16);
	
	
	// Creation of base mesh
	vector<Point3d> vlist;
	vector<int>     flist;
	vector<float>   clist;
	vector<int>     Color_index_list;

	#ifdef MAPPING_TABLE_METHOD
	int Base_color_index_bit = Decoder.get_bits(12);
	
	this->IsKnownIndex.clear();	
	this->ColorArray.clear();
	
	Color_Unit Temp_color;

	for (int i = 0; i < NUMBER_SEEDS; i++)
	{
		this->IsKnownIndex.push_back(-1);
		this->ColorArray.push_back(Temp_color);
	}	
	#endif


	for (int i = 0; i < Number_basemesh_vertex; i++)
	{
		Point_Int Pt_int;
		Pt_int.x = Decoder.get_bits(Max_Qbit + 1);
		Pt_int.y = Decoder.get_bits(Max_Qbit + 1);
		Pt_int.z = Decoder.get_bits(Max_Qbit + 1);
		
		Point3d Pt_real = Change_Int_Real(Pt_int, 0);
		vlist.push_back(Pt_real);
		

		if ((this->IsColored) && (!this->IsOneColor))
		{
			#ifdef PREDICTION_METHOD
		
			Color_Unit TC;
			TC.c0 = Decoder.get_bits(C0_QUANTIZATION);
			TC.c1 = Decoder.get_bits(C1_QUANTIZATION);
			TC.c2 = Decoder.get_bits(C2_QUANTIZATION);

			float L = this->C0_Min + TC.c0 * this->Color_Quantization_Step;
			float a = this->C1_Min + TC.c1 * this->Color_Quantization_Step;
			float b = this->C2_Min + TC.c2 * this->Color_Quantization_Step;
			
			clist.push_back(L);
			clist.push_back(a);
			clist.push_back(b);							
			#endif

			
			#ifdef MAPPING_TABLE_METHOD
		
			int Color_index = Decoder.get_bits(Base_color_index_bit);
			Color_index_list.push_back(Color_index);

			Color_Unit TC;
			
			if (this->IsKnownIndex[Color_index] == 1)
				TC = this->ColorArray[Color_index];

			if (this->IsKnownIndex[Color_index] == -1)
			{
				this->IsKnownIndex[Color_index] = 1;				
				
				TC.c0 = Decoder.get_bits(C0_QUANTIZATION);
				TC.c1 = Decoder.get_bits(C1_QUANTIZATION);
				TC.c2 = Decoder.get_bits(C2_QUANTIZATION);
				
				this->ColorArray[Color_index] = TC;
			}			

			float L = this->C0_Min + TC.c0 * this->Color_Quantization_Step;
			float a = this->C1_Min + TC.c1 * this->Color_Quantization_Step;
			float b = this->C2_Min + TC.c2 * this->Color_Quantization_Step;
			
			clist.push_back(L);
			clist.push_back(a);
			clist.push_back(b);
			#endif
		}
		
		
	}	

	int Facet_index_bit = (int)ceil(log((double)(Number_basemesh_vertex + 1)) / log((double)2));
	
	for (int i = 0; i < (Number_basemesh_facet * 3); i++)
	{
		int v = Decoder.get_bits(Facet_index_bit);
		flist.push_back(v);
	}

	CModifyBasemeshBuilder<HalfedgeDS, Polyhedron, Enriched_kernel> builder(&vlist, &flist, &clist, &Color_index_list);
	pMesh.delegate(builder);
	
	pMesh.compute_normals();
	


	// Seed Edges;
	map<int, int> Seed_Edges;

	for (int i =0; i<2*this->NumberComponents; i++)
	{
		int Vertex_number = Decoder.get_bits(Facet_index_bit);
		Seed_Edges.insert(pair<int,int>(Vertex_number, i));
	}	

	int Basemesh_vertex_number = 0;
	
	Vertex_iterator pVertex = NULL;
	
	map<int,int>::iterator Seed_edge_iterator = Seed_Edges.begin();
	
	int Count_detected_vertices = 0;
	for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); Basemesh_vertex_number++, pVertex++)
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
	


	if ((this->IsColored) && (!this->IsOneColor))
	{
		for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
		{
			float LAB[3];
			LAB[0] = pVertex->color(0);
			LAB[1] = pVertex->color(1);
			LAB[2] = pVertex->color(2);
			
			float RGB[3];
			LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
			pVertex->color(RGB[0], RGB[1], RGB[2]);		
			
			int L = (int)(floor((LAB[0] - this->C0_Min) / this->Color_Quantization_Step + 0.5));
			int a = (int)(floor((LAB[1] - this->C1_Min) / this->Color_Quantization_Step + 0.5));
			int b = (int)(floor((LAB[2] - this->C2_Min) / this->Color_Quantization_Step + 0.5));
			
			pVertex->color_int(L,a,b);			
		}
	
		this->C0_Range = Decoder.get_bits(C0_QUANTIZATION + 1);
		this->Color_0_Model.set_alphabet(this->C0_Range);		
		this->C1_Range = Decoder.get_bits(C1_QUANTIZATION + 1);
		this->Color_1_Model.set_alphabet(this->C1_Range);		
		this->C2_Range = Decoder.get_bits(C2_QUANTIZATION + 1);
		this->Color_2_Model.set_alphabet(this->C2_Range);

		#ifdef MAPPING_TABLE_METHOD
		this->Index_Model.set_alphabet(NUMBER_SEEDS);
		#endif
	}
	if ((this->IsColored) && (this->IsOneColor))
	{
		for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
		{
			pVertex->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
		}		
	}

	// To get number of vertices of each component and restore the real position of vertices
	if (this->NumberComponents != 1)
	{
		pMesh.tag_facets(-1);
		int Component_index = 0;
		
		for (int Component_number = 1; Component_number < this->NumberComponents; Component_number++)
		{
			Halfedge_iterator hi = pMesh.halfedges_begin();			
	
			while((hi->vertex()->Seed_Edge != 2 * Component_number) || (hi->opposite()->vertex()->Seed_Edge != 2 * Component_number + 1))
				hi++;			

			//int Number_vertices = 0;

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

						Point_Int Temp_pos = Change_Real_Int(Wrong_position, 0);
						
						Point3d Real_position = Change_Int_Real(Temp_pos, Component_number);

						pHalfedge->vertex()->point() = Real_position;

						//Number_vertices++;

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

	return this->GlobalCountOperation;
}


// Description : To decode step by step - show intermediate meshes
int Compression_Valence_Component::Decompress_Each_Step(Polyhedron &pMesh, const char* File_Name)
{	
	#ifdef SPLIT_FILE_EACH_RESOLUTION
	this->Decoder.stop_decoder();
	
	char buffer[30];

	std::string S_file_name(File_Name);
	size_t Pos_last_point = S_file_name.find_last_of('.');
	
	string Temp_split_file_name = S_file_name.substr(0, Pos_last_point);//,"_split.");	

	const char * Number = itoa(this->Decompress_count, buffer, 10);
	string Number_iteration(Number);
	string New_file_name = Temp_split_file_name + ".ps";
	New_file_name += Number_iteration;

	FILE * fp_split_iter = fopen(New_file_name.c_str(), "rb");
	this->Decoder.read_from_file(fp_split_iter);
	#endif
	
	if (this->Decompress_count < this->GlobalCountOperation)
	{
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (this->Decompress_count < this->ComponentOperations[Component_ID])
			{				
				if (Decoder.get_bits(1) == 0)
				{
					this->Un_Regulation(pMesh, Decoder, Component_ID);					
					this->Un_Decimation_Conquest(pMesh, Decoder, Component_ID);				
				}
				else
					this->Up_Quantization(pMesh, Decoder,Component_ID);		
				
			}
		}			
	}	
		
	this->Decompress_count++;	
	pMesh.compute_normals();
	
	#ifdef SPLIT_FILE_EACH_RESOLUTION
	fclose(fp_split_iter);
	#endif
	return this->Decompress_count;
	
}




int Compression_Valence_Component::Decompress_To_Level(Polyhedron &pMesh, const int Wanted_level)
{
	#ifdef SPLIT_FILE_EACH_RESOLUTION
	//this->Decoder.stop_decoder();
	//
	//char buffer[30];

	//std::string S_file_name(File_Name);
	//size_t Pos_last_point = S_file_name.find_last_of('.');
	//
	//string Temp_split_file_name = S_file_name.substr(0, Pos_last_point);//,"_split.");	

	//const char * Number = itoa(this->Decompress_count, buffer, 10);
	//string Number_iteration(Number);
	//string New_file_name = Temp_split_file_name + ".ps";
	//New_file_name += Number_iteration;

	//FILE * fp_split_iter = fopen(New_file_name.c_str(), "rb");
	//this->Decoder.read_from_file(fp_split_iter);
	#endif
	
	
	/*while (this->Decompress_count < Wanted_level)
	{
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (this->Decompress_count < this->ComponentOperations[Component_ID])
			{				
				
				if (Decoder.get_bits(1) == 0)
				{
					this->Un_Regulation(pMesh, Decoder, Component_ID);
					this->Un_Decimation_Conquest(pMesh, Decoder, Component_ID);		
				}
				else
					this->Up_Quantization(pMesh, Decoder);		
				
			}
		}			

		this->Decompress_count++;	
		pMesh.compute_normals();
	}
	
	return this->Decompress_count;*/
	return 0;

}


// Description :: Decompressing function. Just give the final mesh.
double Compression_Valence_Component::Decompress_All(Polyhedron &pMesh)
{	
	/*
	Timer timer;
	timer.start();	
	
	#ifdef MESURE_DECOMPRESSION_TIME
	FILE * time = fopen("decompress_time.txt","w");
	#endif
	
	for (; this->Decompress_count < this->CountOperation; this->Decompress_count++)
	{
		if (Decoder.get_bits(1) == 0)
		{			
			this->Un_Regulation(pMesh, Decoder);
			this->Un_Decimation_Conquest(pMesh, Decoder);							
		}
		else
			this->Up_Quantization(pMesh, Decoder);	
		//temp
		#ifdef MESURE_DECOMPRESSION_TIME
		fprintf(time,"%d\t%f\n",pMesh.size_of_vertices(), timer.time());
		#endif
	}

	#ifdef MESURE_DECOMPRESSION_TIME
	fclose(time);
	#endif

	pMesh.compute_normals();

	timer.stop();
	return timer.time();
	*/
	return 0;
}




// Description : Change a point coordinates in real to integer coordinates
Point_Int Compression_Valence_Component::Change_Real_Int(const Point3d &pt, const int & Component_ID)
{
	Point_Int Point;

	double Quantization_step = 0.0;

	// If the quantization resolution is decreased, 
	// we increase the step of quantization by a power of two.
	if (this->NumberChangeQuantization[Component_ID] == 0)
		Quantization_step = this->QuantizationPas[Component_ID];
	else
		Quantization_step = this->QuantizationPas[Component_ID] * pow(2.0, (int)this->NumberChangeQuantization[Component_ID]);	
	
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
		Quantization_step = this->QuantizationPas[Component_ID];
	else
		Quantization_step = this->QuantizationPas[Component_ID] * pow(2.0,(int)this->NumberChangeQuantization[Component_ID]);		
		
	float xmin = this->xmin[Component_ID];
	float ymin = this->ymin[Component_ID];
	float zmin = this->zmin[Component_ID];	

	double x = xmin + (Point.x + 0.5) * Quantization_step;
	double y = ymin + (Point.y + 0.5) * Quantization_step;
	double z = zmin + (Point.z + 0.5) * Quantization_step;

	Point3d pt(x, y, z);

	return pt;
}


	/*
	Description : calculate the geometric distance(error) between two meshes for adaptative quantization.
	Result : Maximum of RMS and Hausdorff distance. */
void Compression_Valence_Component::Calculate_Distances(char * originalMesh,char * attackedMesh,double & mrms,double & mrmswrtBB,double & hausdorff,double & hausdorffwrtBB)
{
#if (0)	// MT
	STARTUPINFO si;
	PROCESS_INFORMATION pi;

	SECURITY_ATTRIBUTES sa = {0};
	sa.nLength = sizeof(SECURITY_ATTRIBUTES);
	sa.bInheritHandle = TRUE;

	HANDLE fp = CreateFile(_T("tempfileformetro.txt"),GENERIC_WRITE,0,&sa,CREATE_ALWAYS,FILE_ATTRIBUTE_NORMAL,NULL);

	ZeroMemory(&si,sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi,sizeof(pi));

	si.dwFlags = STARTF_USESHOWWINDOW|STARTF_USESTDHANDLES;
	si.wShowWindow = SW_HIDE;
	si.hStdOutput = fp;
	si.hStdError = fp;

	wxString inFilename = originalMesh;
	wxString outFilename = attackedMesh;
	wxString commandLine = " " + inFilename + " " + outFilename;
	LPTSTR buff = new TCHAR[512];
	_tcscpy(buff,commandLine);
	LPTSTR process_metro = "metro.exe";
	if (!CreateProcess(process_metro,buff,NULL,NULL,TRUE,0,NULL,NULL,&si,&pi))
	{
		fprintf(LogFile,"CreateProcess failed for metro.exe (%d)\n",GetLastError());
	}

	WaitForSingleObject(pi.hProcess,INFINITE);

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	CloseHandle(fp);

	delete [] buff;
	buff = 0;

	mrms = -1.0;
	FILE * tempMetroFile = fopen("tempfileformetro.txt","r");
	char * pLine = new char[512];
	while(fgets(pLine,512,tempMetroFile))
	{
		wxString currentLine = pLine;
		currentLine.Trim();

		if (currentLine.Left(5)=="  RMS")
		{
			int stringLength = currentLine.Length();
			int commaPosition = currentLine.Find(':');
			wxString stringTemp = currentLine.Mid(commaPosition+2,stringLength-commaPosition-2);
			double rmsTemp;
			stringTemp.ToDouble(&rmsTemp);
			if (rmsTemp>mrms)
				mrms = rmsTemp;
		}

		if (currentLine.Left(9)=="Hausdorff")
		{
			int spacePosition1 = currentLine.Find(':') + 1;
			int spacePosition2 = currentLine.Find(' ',spacePosition1+1);
			int leftParenPosition = currentLine.Find('(');
			int spacePosition3 = currentLine.Find(' ',leftParenPosition+1);
			wxString stringTemp1 = currentLine.Mid(spacePosition1+1,spacePosition2-spacePosition1-1);
			stringTemp1.ToDouble(&hausdorff);
			wxString stringTemp2 = currentLine.Mid(leftParenPosition+1,spacePosition3-leftParenPosition-1);
			stringTemp2.ToDouble(&hausdorffwrtBB);
		}
	}

	fclose(tempMetroFile);

	delete [] pLine;
	pLine = 0;

	remove("tempfileformetro.txt");

	mrmswrtBB = hausdorffwrtBB / hausdorff * mrms;
	if (mrms<=0.00000000001)
		mrmswrtBB = 0.0;
#endif
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


// To reordering the color index for mapping table method.
int Compression_Valence_Component::Mapping_Table_Index_Reordering(Polyhedron &pMesh)
{
	int			     New_index = 0;
	int				 Num_color_base_mesh;
	vector<int>	     Reordered_color_index;
	Vertex_iterator  pVertex;
	
	this->ReorderingColorIndex.clear();

	for (int i = 0; i < NUMBER_SEEDS; i++)
		this->ReorderingColorIndex.push_back(-1);
	
	for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); pVertex++)
	{
		int Current_vertex_index = pVertex->Vertex_Color_Index;
		
		if (this->ReorderingColorIndex[Current_vertex_index] == -1)
		{
			this->ReorderingColorIndex[Current_vertex_index] = New_index;
			New_index++;
		}
	}	
	Num_color_base_mesh = New_index;


	list<int>::iterator Color_iter = this->ColorIndex.begin();

	for (;Color_iter != this->ColorIndex.end(); Color_iter++)
	{
		int Current_vertex_index = *Color_iter;
		
		if (this->ReorderingColorIndex[Current_vertex_index] == -1)
		{
			this->ReorderingColorIndex[Current_vertex_index] = New_index;
			New_index++;
		}
	}	

	return Num_color_base_mesh;
}



void Compression_Valence_Component::Separate_Components(Polyhedron &pMesh)
{
	CCopyPoly<Polyhedron, Enriched_kernel> Copy_Polyhedron;

	int Number_components = pMesh.nb_components();
	
	if (Number_components == 1)
		return;

	for (int i = 0; i < Number_components; i++)
	{
		Polyhedron * New_mesh = new Polyhedron;		
		Copy_Polyhedron.copy(&(pMesh), New_mesh);
		New_mesh->tag_facets(-1);
		int Component_index = 0;

		for (Facet_iterator pFacet = New_mesh->facets_begin(); pFacet != New_mesh->facets_end(); pFacet++)
		{
			if (pFacet->tag() == -1)
			{				
				New_mesh->tag_component(pFacet, -1, Component_index);
				Component_index++;
			}
		}
		
		for (int j = 0; j < Number_components - 1; j++)
		{

			for (Halfedge_iterator pHedge = New_mesh->halfedges_begin(); pHedge != New_mesh->halfedges_end(); )
			{
				Halfedge_handle h = pHedge;
				pHedge++;							
				
				if (!h->is_border())
				{
					int Facet_tag = h->facet()->tag();

					if (Facet_tag != i)
					{
						New_mesh->erase_connected_component(h);							
						break;
					}
				}
			}
		}

		QString Outputfile = QString("Separate_%1").arg(i);	//Outputfile += wxString::Format("%d",i);
		Outputfile += ".off";
		New_mesh->write_off(Outputfile.toStdString(), true, false);
	}
}

int Compression_Valence_Component::GetResolutionChange(Polyhedron *pMesh, float Prec)
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

	for (Vertex_iterator pVertex = pMesh->vertices_begin(); pVertex!= pMesh->vertices_end(); pVertex++)
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

	////~la moitie des triangle sont affichés on choisi n pixel/triangle
	
	return 2*floor(aire/Prec);
}
#endif
