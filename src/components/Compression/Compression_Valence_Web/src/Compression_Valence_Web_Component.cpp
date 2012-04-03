///////////////////////////////////////////////////////////////////////////
// Author: Ho LEE
// Year: 2011
// Month: MAY
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence_Web

#include "Compression_Valence_Web_Component.h"
#include "Compression_Valence_Web_Polyhedron.h"

#include "Compression_Valence_Common.h"


#include <CGAL/Timer.h>
#include <iostream>


#include <map>
#include <set>
#include <bitset>

#define COLOR_NUMBER 10000
#define USE_COLOR_METRIC
#define AC_BUFFER 1024 * 10000

const int MINIMUM_PREDICTION_NUMBER = 3;
const int LIMIT_NUMBER = 50;

QString Compression_Valence_Web_Component::Main_Function(Polyhedron     & _pMesh,
													const char *     _Input_File_Name,
													const char*      _File_Name,
													const int      & _Qbit,
													const int      & _NVertices,
													/*const bool       _Normal_flipping,
													const bool       _Use_metric,
													const float    & _Metric_thread,
													const bool       _Use_forget_metric,
													const int      & _Forget_value, */
													const bool       _Compression_selected)
													///const bool       _Adaptive_quantization,
													//const bool       _Is_bijection_selected)												
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


	unsigned Init_number_vertices = (unsigned)_pMesh.size_of_vertices();
	
	// Initialization - Quantization, Color, Multiple components
	this->Global_Initialization(_pMesh, _Qbit, _File_Name);
		
	// When use of adaptive quantization is selected

	this->Simplification(_pMesh, _NVertices, false, false, false, false, false);		
	
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
void Compression_Valence_Web_Component::Global_Initialization(Polyhedron & _pMesh, 
														  const int  & _Qbit,
														  const char * _File_Name)
{
	// (1) Determination if the mesh is colored.
	// (2) Conversion of color space.
	// (3) Quantization of converted colors into "vertex->color_int_web()"
	// (4) Re-paint mesh with re_calculated rgb colors from quantized converted color.
	// (5) Establish color palette.
	// (6) Color quantization - Descrease of possible color numbers 
	this->Color_Initialization(_pMesh);	
	
	// Initialization of multiple components (quantization is performed separately for each component)
	this->Multiple_Components_Initialization(_pMesh, _Qbit);		

	// Quantization of each component
	this->Quantization(_pMesh);
}	

void Compression_Valence_Web_Component::Multiple_Components_Initialization(Polyhedron & _pMesh, 
																	   int const  & _Qbit)
{
	// Initialize vertex flags
	for (Vertex_iterator pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
	{
		pVertex->Seed_Edge_web = OTHER_COORDINATE;	
		pVertex->Component_Number_web = -1;
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
			bool Is_Seed_Edge_web_found = false;

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
				F->Component_Number_web = Component_index;
				
				area += Area_Facet_Triangle(F->halfedge());			

				Halfedge_around_facet_circulator pHalfedge = F->facet_begin();
				Halfedge_around_facet_circulator end = pHalfedge;
				
				CGAL_For_all(pHalfedge, end)
				{
					// tag the vertex to its corresponding component number
					if (pHalfedge->vertex()->Component_Number_web == -1)
					{
						pHalfedge->vertex()->Component_Number_web = Component_index;
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
					
					if (!Is_Seed_Edge_web_found)
					{
						if ( ( (!pHalfedge->is_border_edge()) && (!Is_Border_Vertex(pHalfedge)) && (!Is_Border_Vertex(pHalfedge->opposite())) ) || (pHalfedge->next()->vertex_degree() != 6) )
						{
							Is_Seed_Edge_web_found = true;
							
							//Seed edge of each component
							pHalfedge->vertex()->Seed_Edge_web = 2 * Component_index;
							pHalfedge->opposite()->vertex()->Seed_Edge_web = 2 * Component_index + 1;							
						}
					}

					Facet_handle pNFacet = pHalfedge->opposite()->facet();
					if (pNFacet != NULL && pNFacet->tag() == -1)
					{
						facets.push_front(pNFacet);
						pNFacet->tag(Component_index);
						pNFacet->Component_Number_web = Component_index;
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

		//connectivity
		this->Connectivity.push_back(li);

		//geometry
		this->Geometry.push_back(pi);

		//vertex color
		this->VertexColor.push_back(cu);
		
		// Number of connectivity and geometry symbols
		this->NumberSymbol.push_back(li);
		this->NumberVertices.push_back(li);
			
		// Number of operations for each component
		this->ComponentOperations.push_back(0);

		// Number of decimations for each component
		this->NumberDecimation.push_back(0);
		

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
void Compression_Valence_Web_Component::Quantization(Polyhedron & _pMesh) 
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

		int Component_ID = pVert->Component_Number_web;
		int NbInteraval = pow(2., (int)this->Qbit[Component_ID]);

		int Qx = (int)(ceil((x - (double)this->xmin[Component_ID]) / (double)this->Quantization_Step[Component_ID])) - 1;
		if (Qx <= -1)
			Qx = 0;
		else if (Qx>=NbInteraval)
			Qx=NbInteraval-1;
		int Qy = (int)(ceil((y - (double)this->ymin[Component_ID]) / (double)this->Quantization_Step[Component_ID])) - 1;
		if (Qy <= -1)
			Qy = 0;
		else if (Qy>=NbInteraval)
			Qy=NbInteraval-1;
		int Qz = (int)(ceil((z - (double)this->zmin[Component_ID]) / (double)this->Quantization_Step[Component_ID])) - 1;
		if (Qz <= -1)
			Qz = 0;
		else if (Qz>=NbInteraval)
			Qz=NbInteraval-1;
		
		pVert->point() = Point3d(this->xmin[Component_ID] + (Qx + 0.5) * this->Quantization_Step[Component_ID],
							     this->ymin[Component_ID] + (Qy + 0.5) * this->Quantization_Step[Component_ID],
							     this->zmin[Component_ID] + (Qz + 0.5) * this->Quantization_Step[Component_ID]);		
	}
}



// this->ColorArray -> contains all initial colors present in the input mesh.
void Compression_Valence_Web_Component::Color_Initialization(Polyhedron &_pMesh)
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
		
		//Enter quantized color valued into "vertex->color_int_web" to use lated.
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
			pVertex->color_int_web(Resulting_color.c0, Resulting_color.c1, Resulting_color.c2);			

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




// Description : This function select a set of independent vertices to be removed
int Compression_Valence_Web_Component::Decimation_Conquest(Polyhedron  & _pMesh,
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
	
	while((hi->vertex()->Seed_Edge_web != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2 * Component_ID+1))
		hi++;

	Halfedge_handle First_halfedge = &(*(hi));
	
	// Two vertices of seed edge are flaged CONQUERED
	First_halfedge->vertex()->Vertex_Flag_web = CONQUERED;
	First_halfedge->opposite()->vertex()->Vertex_Flag_web = CONQUERED;

	// These vertices are also flaged with sign flags for retriangulation
	First_halfedge->vertex()->Vertex_Sign_web = PLUS;
	First_halfedge->opposite()->vertex()->Vertex_Sign_web = MINUS;	
	
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
		if ((h->facet()->Facet_Flag_web == CONQUERED) || (h->facet()->Facet_Flag_web == TO_BE_REMOVED))
			continue;

		// if its front vertex is free and has a valence <= 6 and it is not a border vertex.
		else if ( (h->next()->vertex()->Vertex_Flag_web == FREE) && 
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
				this->InterGeometry.push_front(CRV);	
				// encoding of vertex color
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Removed_vertex_color;
					Removed_vertex_color = Get_Vertex_Color(g->next());					
					this->InterVertexColor.push_front(Removed_vertex_color);
				}			
								

				// Enter symbol 'VALENCE CODE' into the list of symbols
				this->InterConnectivity.push_front(valence - 3);

				// remove the front vertex
				_pMesh.erase_center_vertex(g->next()); 

				g = h;
				g->facet()->Facet_Flag_web = TEMP_FLAG;				
				
				// Flags and fill queue
				for (int j = 0; j < (valence-1); j++)
				{
					g = g->next();
					g->vertex()->Vertex_Flag_web = CONQUERED;
					
					g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
					
					if (!g->is_border_edge())
						Halfedges.push(g->opposite());
				}
				g = h;							
												
				Retriangulation(_pMesh, g, valence, 0, Component_ID);

			}
			
			// Conditions are not satisfied -> NULL CODE
			else 
			{
				// Enter symbol 'NULL PATCH' into the list of symbols
				this->InterConnectivity.push_front(4);

				Number_symbol++;
				h->facet()->Facet_Flag_web = CONQUERED;
				h->next()->vertex()->Vertex_Flag_web = CONQUERED;
				
				if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = PLUS;
				}
				
				if (!h->next()->is_border_edge())
					Halfedges.push(h->next()->opposite());				
				
				if (!h->prev()->is_border_edge())					
					Halfedges.push(h->prev()->opposite());				
			}
		}
		// border edge.
		else if ( (h->next()->vertex()->Vertex_Flag_web == FREE) && 
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
				this->InterGeometry.push_front(Vertex_position);

				Color_Unit Removed_vertex_color;
				//int Vertex_color_index = -1;
				if ((this->IsColored) && (!this->IsOneColor))
				{	
					Removed_vertex_color = Get_Vertex_Color(h->next());
 					this->InterVertexColor.push_front(Removed_vertex_color);	
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
					Input_gate->facet()->Facet_Flag_web = CONQUERED;					
					Input_gate->next()->vertex()->Vertex_Flag_web = CONQUERED;

					if ((type == 1) || (type == 2) || (type == 4))
					{
						if (Input_gate->next()->vertex()->Vertex_Sign_web == NOSIGN)
							Input_gate->next()->vertex()->Vertex_Sign_web = PLUS;
					}
					else if (type == 3)
					{
						if (Input_gate->next()->vertex()->Vertex_Sign_web == NOSIGN)
							Input_gate->next()->vertex()->Vertex_Sign_web = MINUS;
					}
					

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
						Border_edges[1]->opposite()->facet()->Facet_Flag_web = CONQUERED;						

						_pMesh.add_facet_to_border(Retriangulation_edge, Border_edges[0]->opposite());
						Border_edges[0]->opposite()->facet()->Facet_Flag_web = CONQUERED;						
					}

					else if (   ( (Number_jump == 0) && ((type == 6) || (type == 7)) ) ||
								( (Number_jump == 1) && ((type == 5) || (type == 8)) ) ||
								( (Number_jump == 2) && ((type == 6) || (type == 7)) )  )
					{												
						_pMesh.add_facet_to_border(Border_edges[2]->opposite(), Border_edges[0]->opposite());
						Border_edges[1]->opposite()->facet()->Facet_Flag_web = CONQUERED;						
						Halfedge_handle Temp_border = Border_edges[2]->opposite()->next();
						
						_pMesh.add_facet_to_border(Retriangulation_edge, Temp_border);
						Border_edges[2]->opposite()->facet()->Facet_Flag_web = CONQUERED;						
					}					

					Halfedge_handle Input_gate = Border_edges[Number_jump]->opposite();

					
					// Vertex Signs
					if ((type == 5) || (type == 8))
					{
						Halfedge_handle g = Input_gate;
						g = g->prev()->opposite()->next();
						if (g->vertex()->Vertex_Sign_web == NOSIGN)
							g->vertex()->Vertex_Sign_web = MINUS;
						if (g->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							g->opposite()->vertex()->Vertex_Sign_web = PLUS;
					}
					else if ((type == 6) || (type == 7))
					{
						Halfedge_handle g = Input_gate;
						g = g->next()->opposite()->prev();
						if (g->vertex()->Vertex_Sign_web == NOSIGN)
							g->vertex()->Vertex_Sign_web = PLUS;
						if (g->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							g->opposite()->vertex()->Vertex_Sign_web = MINUS;
					}
					
					//Vertex Flags + Fill Halfedges queue.
					for (int i = 0; i < valence - 1; i++)
					{
						Border_edges[i]->vertex()->Vertex_Flag_web = CONQUERED;
						Border_edges[i]->opposite()->vertex()->Vertex_Flag_web = CONQUERED;

						if (i != (unsigned)Number_jump)
							Halfedges.push(Border_edges[i]);
					}

				}					
			}

			else // Border vertex can not be removed
			{
				Number_symbol++;
				this->InterConnectivity.push_front(4);
				
				h->facet()->Facet_Flag_web = CONQUERED;				
				h->next()->vertex()->Vertex_Flag_web = CONQUERED;

				if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
						h->next()->vertex()->Vertex_Sign_web = PLUS;
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
			h->facet()->Facet_Flag_web = CONQUERED;
			h->next()->vertex()->Vertex_Flag_web = CONQUERED;///////////////////////////

			if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = PLUS;
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
int Compression_Valence_Web_Component::Regulation(Polyhedron  & _pMesh,
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
	
	while((hi->vertex()->Seed_Edge_web != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2*Component_ID+1))
		hi++;	

	Halfedge_handle First_halfedge = &(*(hi));	

	First_halfedge->vertex()->Vertex_Flag_web = CONQUERED;
	First_halfedge->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
		
	std::queue<Halfedge_handle> Halfedges;
	Halfedges.push(First_halfedge);

	Halfedge_handle h;

	while(!Halfedges.empty())
	{		
		h = Halfedges.front();
		Halfedges.pop();
		
		size_t valence = h->next()->vertex_degree();
		//int Component_ID = h->next()->vertex()->Component_Number_web;

		if ((h->facet()->Facet_Flag_web == CONQUERED) || (h->facet()->Facet_Flag_web == TO_BE_REMOVED))
			continue;

		else if ((h->next()->vertex()->Vertex_Flag_web == FREE) && (valence == 3) && (Is_Border_Vertex(h->next()) == false)) // if valence is 3, remove the front vertex.
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
//cata new
				this->InterGeometry.push_front(Geo);
//cata new end
				
//cata
//				Point3d Barycenter = Barycenter_Patch_Before_Removal(g);								
//				Point_Int BC = Change_Real_Int(Barycenter, Component_ID);		
// cata end				
				
				if ((this->IsColored) && (!this->IsOneColor))
				{
					Color_Unit Removed_vertex_color;
					
					Removed_vertex_color = Get_Vertex_Color(h->next());
//cata new
					this->InterVertexColor.push_front(Removed_vertex_color);
//cata new end

				}

								
				g = h->next();
				g->vertex()->Vertex_Flag_web = TO_BE_REMOVED;
				g->facet()->Facet_Flag_web = TO_BE_REMOVED;

				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
				g->facet()->Facet_Flag_web = TO_BE_REMOVED;
				
				if (!g->prev()->is_border_edge())
				{
					Halfedge_handle h1 = g->prev()->opposite();
					h1->facet()->Facet_Flag_web = CONQUERED;
					h1->next()->vertex()->Vertex_Flag_web = CONQUERED;
					if (!h1->next()->is_border_edge())
						Halfedges.push(h1->next()->opposite());
					if (!h1->prev()->is_border_edge())
						Halfedges.push(h1->prev()->opposite());					
				}
				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
				g->facet()->Facet_Flag_web = TO_BE_REMOVED;
				if (!g->prev()->is_border_edge())
				{
					Halfedge_handle h2 = g->prev()->opposite();
					h2->facet()->Facet_Flag_web = CONQUERED;
					h2->next()->vertex()->Vertex_Flag_web = CONQUERED;
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
				h->facet()->Facet_Flag_web = CONQUERED;
				h->next()->vertex()->Vertex_Flag_web = CONQUERED;
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
			h->facet()->Facet_Flag_web = CONQUERED;
			h->next()->vertex()->Vertex_Flag_web = CONQUERED;
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

		if (vh->Vertex_Flag_web == TO_BE_REMOVED)
		{
			Halfedge_handle temp = vh->halfedge();
			temp = _pMesh.erase_center_vertex(temp);
			temp->facet()->Component_Number_web = Component_ID;
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
void Compression_Valence_Web_Component::Un_Regulation(Polyhedron &_pMesh, const int & Component_ID)
{	
	Init(_pMesh);
	int Connnectivity_Range=1;

	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;

	unsigned Qbit = this->Qbit[Component_ID];
	float Color_step = 0.0;
	Color_step = this->Color_Quantization_Step;

	
	while((hi->vertex()->Seed_Edge_web != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2*Component_ID+1))	
		hi++;	
	
	hi->vertex()->Vertex_Flag_web = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		if ((h->facet()->Facet_Flag_web == CONQUERED) || (h->facet()->Facet_Flag_web == TO_BE_REMOVED))// already visited.
			continue;
		
		// read connectivity information
		int valence = Decoder.decode_bits(Connnectivity_Range) + 3;
		
		// if valence is 3, remove the front vertex.
		if (valence == 3)
		{
			Halfedge_handle g = h;			

			// Insertion of a vertex
			Halfedge_handle pass = h;
			
			vector<Point3d> Vertices; //contains 1-ring and 2-rings vertices



//cata new		
			Point_Int geo;
			geo.x = Decoder.decode_bits(Qbit);
			geo.y = Decoder.decode_bits(Qbit);
			geo.z = Decoder.decode_bits(Qbit);
			Point3d Center_vertex = this->Change_Int_Real(geo,Component_ID);	
//cata new end

			// Vertex insertion

			g = _pMesh.create_center_vertex(g);
			
			g->vertex()->point() = Center_vertex;
			g->vertex()->Seed_Edge_web = -1;		
			

			// Vertex flags
			g->vertex()->Vertex_Flag_web = CONQUERED;
			g = h;
			g->facet()->Facet_Flag_web = CONQUERED;
			g = g->next();
			g = g->prev_on_vertex();
			g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
			g->facet()->Facet_Flag_web = CONQUERED;			
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h1 = g->prev()->opposite();
				h1->facet()->Facet_Flag_web = CONQUERED;
				h1->next()->vertex()->Vertex_Flag_web = CONQUERED;
				if (!h1->next()->is_border_edge())
					Halfedges.push(h1->next()->opposite());
				if (!h1->prev()->is_border_edge())
					Halfedges.push(h1->prev()->opposite());
			}
			g = g->prev_on_vertex();
			g->facet()->Facet_Flag_web = CONQUERED;
			g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
			if (!g->prev()->is_border_edge())
			{
				Halfedge_handle h2 = g->prev()->opposite();
				h2->facet()->Facet_Flag_web = CONQUERED;
				h2->next()->vertex()->Vertex_Flag_web = CONQUERED;
				if (!h2->next()->is_border_edge())
					Halfedges.push(h2->next()->opposite());
				if (!h2->prev()->is_border_edge())
					Halfedges.push(h2->prev()->opposite());
			}
			
			
			if ((this->IsColored) && (!this->IsOneColor))
			{				
				g = h;				

//cata new		
			Color_Unit color;
			color.c0 = Decoder.decode_bits(HEADER_L_COLOR_QUANT); // Read color info.
			color.c1 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);
			color.c2 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);

			g->next()->vertex()->color_int_web(color.c0, color.c1, color.c2);
			float LAB[3];
			LAB[0] = this->C0_Min + color.c0 * Color_step;
			LAB[1] = this->C1_Min + color.c1 * Color_step;
			LAB[2] = this->C2_Min + color.c2 * Color_step;

			float RGB[3];
			LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
			g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);

//cata new end

			}
			
			if ((this->IsColored) && (this->IsOneColor))
			{
				g = h->next();
				g->vertex()->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
			}
		}

		else
		{
			h->facet()->Facet_Flag_web = CONQUERED;
			h->next()->vertex()->Vertex_Flag_web = CONQUERED;
			if (!h->next()->is_border_edge())
				Halfedges.push(h->next()->opposite());
			if (!h->prev()->is_border_edge())
				Halfedges.push(h->prev()->opposite());
		}
	}
}



// Description : Decoding function of decimation conquest
void Compression_Valence_Web_Component::Un_Decimation_Conquest(Polyhedron       & _pMesh, 
														 const int        & Component_ID)
{
	Init(_pMesh);
	
	int Number_connectivity_symbols;
	if (this->IsClosed[Component_ID])
		Number_connectivity_symbols = 5;
	else
		Number_connectivity_symbols = 7;
	int Connectivity_Range = 3;

//	Adaptive_Data_Model Connectivity(Number_connectivity_symbols);	
	
	Halfedge_iterator hi = _pMesh.halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;	

	unsigned Qbit = this->Qbit[Component_ID];


	float Color_step = 0.0;
	Color_step = this->Color_Quantization_Step;


	while((hi->vertex()->Seed_Edge_web != 2 * Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2 * Component_ID+1))
		hi++;	

	// Two vertices of seed edges are flaged CONQUERED
	hi->vertex()->Vertex_Flag_web = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag_web = CONQUERED;	
		
	// These vertices are also flages with sign flags for retriangulation
	hi->vertex()->Vertex_Sign_web = PLUS;
	hi->opposite()->vertex()->Vertex_Sign_web = MINUS;
	
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
		if ((h->facet()->Facet_Flag_web == CONQUERED) || (h->facet()->Facet_Flag_web == TO_BE_REMOVED))
			continue;
			
		valence = Decoder.decode_bits(Connectivity_Range) + 3;		
				
		// if its front vertex is free
		if ((valence >= 3) && (valence <= 6))
		{
			type = Find_Type(h, valence);
			
			// remove the front vertex if its removal does not viloate the manifold property
			Halfedge_handle g = h;

			// remove edges to re_find polygon and (check. and attribute sign flag.)
			bool Check_Validity = Remove_Edges(_pMesh, g, type);
			Check_Validity = false;

			//Check_Validity = false;// * * * * * * * * * * *////
			if (Check_Validity == false)
			{
				g = h;
				Halfedge_handle pass = h;
				
				vector<Point3d> Vertices; //contains 1-ring and 2-ring vertices;


//cata new		
			Point_Int geo;
			geo.x = Decoder.decode_bits(Qbit);
			geo.y = Decoder.decode_bits(Qbit);
			geo.z = Decoder.decode_bits(Qbit);
			Point3d Center_vertex = this->Change_Int_Real(geo,Component_ID);	
//cata new end
				// Assign the region number to inserted vertex
				Halfedge_handle reg = h;

				g = _pMesh.create_center_vertex(g);
				g->vertex()->point() = Center_vertex;

				g->vertex()->Vertex_Flag_web = CONQUERED;

				g = h;
				g->facet()->Facet_Flag_web = TO_BE_REMOVED;
				g->vertex()->Vertex_Flag_web = CONQUERED;
				g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;

				for (unsigned int k = 0; k < (valence - 1); k++)
				{
					g = g->next()->opposite()->next();
					g->facet()->Facet_Flag_web = TO_BE_REMOVED;
					g->vertex()->Vertex_Flag_web = CONQUERED;
					g->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
					if (g->is_border_edge() == false)
						Halfedges.push(g->opposite());
				}
				
				g->next()->vertex()->Seed_Edge_web = -1;

			if ((this->IsColored) && (!this->IsOneColor))
			{				
				g = h;				

//cata new		
			Color_Unit color;
			color.c0 = Decoder.decode_bits(HEADER_L_COLOR_QUANT); // Read color info.
			color.c1 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);
			color.c2 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);

			g->next()->vertex()->color_int_web(color.c0, color.c1, color.c2);
			float LAB[3];
			LAB[0] = this->C0_Min + color.c0 * Color_step;
			LAB[1] = this->C1_Min + color.c1 * Color_step;
			LAB[2] = this->C2_Min + color.c2 * Color_step;

			float RGB[3];
			LAB_To_RGB(LAB[0], LAB[1], LAB[2], RGB);
			g->next()->vertex()->color(RGB[0], RGB[1], RGB[2]);

//cata new end

					
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

//cata new		
			Point_Int geo;
			geo.x = Decoder.decode_bits(Qbit);
			geo.y = Decoder.decode_bits(Qbit);
			geo.z = Decoder.decode_bits(Qbit);
			Point3d Center_vertex = this->Change_Int_Real(geo,Component_ID);	
//cata new end
		
			Halfedge_handle g = h;			
			Color_Unit color;
						
			if ((this->IsColored) && (!this->IsOneColor))
			{
//cata new		
			color.c0 = Decoder.decode_bits(HEADER_L_COLOR_QUANT); // Read color info.
			color.c1 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);
			color.c2 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);

//			g->next()->vertex()->color_int_web(color.c0, color.c1, color.c2);
//cata new end
			}			

			// border edge with valence == 3
			if (valence == 8)
			{				
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
				g->vertex()->Vertex_Flag_web = CONQUERED;

				//#ifdef PREDICTION_METHOD
				if ((this->IsColored) && (!this->IsOneColor))
				{
					g->vertex()->color_int_web(color.c0, color.c1, color.c2);				
					float LAB[3];
					LAB[0] = this->C0_Min + color.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + color.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + color.c2 * this->Color_Quantization_Step;

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
				Tag_handle->vertex()->Vertex_Flag_web = CONQUERED;

				if ((type == 1) || (type == 2) || (type == 4))
				{
					if (Tag_handle->vertex()->Vertex_Sign_web == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign_web = PLUS;
				}
				else if (type == 3)
				{
					if (Tag_handle->vertex()->Vertex_Sign_web == NOSIGN)
						Tag_handle->vertex()->Vertex_Sign_web = MINUS;
				}
				for (int i=0; i <2; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag_web = CONQUERED;
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

				// to create the new facets
				Halfedge_handle Prev_edge = Border_edges[2]->opposite()->prev();
				g = _pMesh.add_vertex_and_facet_to_border(Prev_edge, Border_edges[2]->opposite());
				g->vertex()->point() = Center_vertex;				
					
				//#ifdef PREDICTION_METHOD
				Color_Unit CV;
				if ((this->IsColored) && (!this->IsOneColor))
				{
					g->vertex()->color_int_web(color.c0, color.c1, color.c2);				
					float LAB[3];
					LAB[0] = this->C0_Min + color.c0 * this->Color_Quantization_Step;
					LAB[1] = this->C1_Min + color.c1 * this->Color_Quantization_Step;
					LAB[2] = this->C2_Min + color.c2 * this->Color_Quantization_Step;

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
						if (Border_edges[2]->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign_web = PLUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign_web = MINUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign_web = MINUS;
						if (Border_edges[2]->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[2]->opposite()->vertex()->Vertex_Sign_web = PLUS;
					}

				}

				else if (Number_jump == 1)
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[2]->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign_web = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign_web = PLUS;
					}
					else
					{
						if (Border_edges[2]->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[2]->vertex()->Vertex_Sign_web = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign_web = MINUS;
					}

				}
				else // jump == 2
				{
					if ((type == 5) || (type == 8))
					{
						if (Border_edges[0]->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign_web = PLUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign_web = MINUS;
					}
					else
					{
						if (Border_edges[0]->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[0]->vertex()->Vertex_Sign_web = MINUS;
						if (Border_edges[0]->opposite()->vertex()->Vertex_Sign_web == NOSIGN)
							Border_edges[0]->opposite()->vertex()->Vertex_Sign_web = PLUS;
					}

				}
				
				for (int i=0; i < 3; i++)
				{
					Border_edges[i]->opposite()->facet()->Facet_Flag_web = CONQUERED;
					
					Border_edges[i]->vertex()->Vertex_Flag_web = CONQUERED;
					Border_edges[i]->opposite()->vertex()->Vertex_Flag_web = CONQUERED;
					
					if (i != Number_jump)
						Halfedges.push(Border_edges[i]);
				}					
			}
		}

			
		
		
		//  the symbol == N
		else if (valence == 7)
		{
			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag_web = CONQUERED;
			h->next()->vertex()->Vertex_Flag_web = CONQUERED;

			if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign_web == PLUS) && (h->opposite()->vertex()->Vertex_Sign_web == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign_web == MINUS) && (h->opposite()->vertex()->Vertex_Sign_web == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign_web == NOSIGN)
					h->next()->vertex()->Vertex_Sign_web = PLUS;
			}
			
			if (h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());

			if (h->prev()->is_border_edge() == false)
				Halfedges.push(h->prev()->opposite());	
		}
	}	
}


// When the last loop doen not remove any vertex, all information are deleted
void Compression_Valence_Web_Component::Remove_Last_Phase_Elements(const int & Component_ID)
{
	this->NumberDecimation[Component_ID] -= 1; //this->NumberDecimation[Component_ID] - 1;
	this->ComponentOperations[Component_ID] -= 1;	
//cata
//	this->ListOperation[Component_ID].pop_front();
//cata end
	for (int i = 0; i < (this->DumpSymbolDecimation + this->DumpSymbolRegulation);i++)	
		this->Connectivity[Component_ID].pop_front();		

	for (int i = 0; i < 2; i++)
	{			
		this->NumberSymbol[Component_ID].pop_front();
		this->NumberVertices[Component_ID].pop_front();
	}	
}



// To store information needed for the reconstruction of the base mesh.
void Compression_Valence_Web_Component::Write_Base_Mesh(Polyhedron & _pMesh, bit_file_c & enc, unsigned &Connectivity_size, unsigned & Color_size, const int & Num_color_base_mesh)
{

	unsigned int Max_Qbit = 0;
	for (int i =0; i< this->NumberComponents; i++)
	{

		if (this->Qbit[i] > Max_Qbit)
			Max_Qbit = this->Qbit[i];
	}
		
	enc.encode_bits(_pMesh.size_of_vertices(), HEADER_L_NUM_VERTICES);									    // number of vertices of base mesh < 4096
	enc.encode_bits(_pMesh.size_of_facets(), HEADER_L_NUM_FACETS);										  // number of facets of base mesh < 8192	

	//int Base_color_index_bit = 0;

	int Basemesh_Vertex_Number_web = 0;
	vector<int> Seed_Edge_webs(2 * this->NumberComponents, -1);

	// Encoding of vertex information of base mesh //
	int counter=0;
	for (Vertex_iterator pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); Basemesh_Vertex_Number_web++, pVertex++)
	{			
		pVertex->Vertex_Number_web = Basemesh_Vertex_Number_web;
		
		//int SEDD = pVertex->Seed_Edge_web;

		if (pVertex->Seed_Edge_web != OTHER_COORDINATE)
			Seed_Edge_webs[pVertex->Seed_Edge_web] = Basemesh_Vertex_Number_web;		
		
		//int cid = pVertex->Component_Number_web;

		Point_Int Vertex = Change_Real_Int(pVertex->point(), pVertex->Component_Number_web);
		
////cata
//		std::ofstream file("base.txt",fstream::app);
//		file<<Vertex;
//		file.close();
////cata end

		enc.encode_bits(Vertex.x, Max_Qbit);
		enc.encode_bits(Vertex.y, Max_Qbit);
		enc.encode_bits(Vertex.z, Max_Qbit);
		counter++;

		if ((this->IsColored) && (!this->IsOneColor))
		{
			//#ifdef PREDICTION_METHOD
			int C0 = pVertex->color_int_web(0);
			int C1 = pVertex->color_int_web(1);
			int C2 = pVertex->color_int_web(2);
			
			enc.encode_bits(C0, HEADER_L_COLOR_QUANT);
			enc.encode_bits(C1, HEADER_L_COLOR_QUANT);
			enc.encode_bits(C2, HEADER_L_COLOR_QUANT);
					
			Color_size += 3 * (HEADER_L_COLOR_QUANT);
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
			enc.encode_bits(pHalfedge->vertex()->Vertex_Number_web, Facet_index_bit);			
			pHalfedge = pHalfedge->next();
			

		} while (pHalfedge != pFacet->halfedge());		
	}
	
	// Store seed edge information.
	for (int i = 0 ; i < (int)Seed_Edge_webs.size(); i++)	// MT
		enc.encode_bits(Seed_Edge_webs[i], Facet_index_bit);	
	
	Connectivity_size += Facet_index_bit * (Count_facet_index + Seed_Edge_webs.size());	
}


void Compression_Valence_Web_Component::Simplification(Polyhedron  & _pMesh,
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
////cata new
//		char file_q[70]="down_";
//		char file_number[70]="";
//		strcpy(file_number,"");
//		sprintf(file_number,"%.3d",this->GlobalCountOperation);
//		_pMesh.write_off(strcat(strcat(file_q,file_number),".off"),false,false);
////cata new end

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

////cata
//				this->ListOperation[Component_ID].push_front(0);	
//cata end				
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
//cata new	
//	_pMesh.write_off("down_base.off",false,false);
//cata new end

	_pMesh.compute_normals();
}

void Compression_Valence_Web_Component::Compression(Polyhedron     & _pMesh, 
												const char     * File_Name, 
												const int      & _Qbit, 
												unsigned       & Connectivity_size, 
												unsigned       & Color_size, 
												unsigned       & Total_size)
												//const unsigned & Initial_file_size)
{

	FILE * fp  = fopen(File_Name, "wb");													//Main FILE to save compression information.		
	
	int res;
	res=fwrite(&this->Initial_file_size, sizeof(unsigned), 1, fp);			// Intial size of the input file (To visualize during decompression)	
	res=fwrite(&this->NumberComponents, sizeof(int), 1, fp);
	
		// Color Type of mesh
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

////cata
//		res=fwrite(&this->Smallest_C0, sizeof(int), 1, fp);  
//		res=fwrite(&this->Smallest_C1, sizeof(int), 1, fp); 
//		res=fwrite(&this->Smallest_C2, sizeof(int), 1, fp);
////cata end		
		Color_size += sizeof(float) * 8 * 4;		
	
	}	
	
	if ((this->IsColored) && (this->IsOneColor))
	{

		res=fwrite(&this->OnlyColor[0], sizeof(float), 1, fp); // smallest value of c0 
		res=fwrite(&this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1 
		res=fwrite(&this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2		
	}


	for (int i =0; i<this->NumberComponents; i++)
	{
		res=fwrite(&this->Quantization_Step[i], sizeof(float), 1, fp);							    // Quantization_Step(step of quantization)
		res=fwrite(&this->xmin[i], sizeof(float), 1, fp);																	   // xmin value
		res=fwrite(&this->ymin[i], sizeof(float), 1, fp);																	   // ymin value
		res=fwrite(&this->zmin[i], sizeof(float), 1, fp);																	   // zmin value
	}
	
	


	
	fclose(fp);
	//close the header file.. reopen it as bitfile for coding

	res+=0; // just to avoid gcc 4.6 warning


	//test of header size
	FILE* f_size = fopen(File_Name,"rb");
	fseek(f_size, 0, SEEK_END);
	Total_size = ftell(f_size);
	fclose(f_size);
	// Declaration du codeur.
   // open bit file for writing
	bit_file_c enc;
	enc.Open(File_Name, BF_APPEND);
   
	//CBitStream enc;
   
	for (int i = 0; i < this->NumberComponents; i++)
	{
		if (this->IsClosed[i])
			enc.encode_bits(0,HEADER_L_CLOSED);								
		else		
			 enc.encode_bits(1,HEADER_L_CLOSED);	
		enc.encode_bits(this->ComponentOperations[i], HEADER_L_COMPONENT_OPERATION);	 									   // number of operations < 256	
		enc.encode_bits(this->Qbit[i], HEADER_L_QBIT);															// quantization bit < 16

		
	}



	/*	Write information of base mesh.
		geometry + connectivity + color information(if the mesh is colored). */
	this->Write_Base_Mesh(_pMesh, enc, Connectivity_size, Color_size, Num_color_base_mesh);	


////cata
//	if ((this->IsColored) && (!this->IsOneColor))
//	{
//		enc.WriteBits(this->C0_Range, C0_QUANTIZATION);
//		enc.WriteBits(this->C1_Range, C1_QUANTIZATION);
//		enc.WriteBits(this->C2_Range, C2_QUANTIZATION);		
//	
//	}
////cata end
		
	int connectivity_max = *std::max_element(this->Connectivity[0].begin(), this->Connectivity[0].end());

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

////cata				int Type_operation = this->ListOperation[Component_ID].front();
//				this->ListOperation[Component_ID].pop_front();				
//
//				// decimation is chosen to be applied.
//				if (Type_operation == 0) 
//				{			
//					enc.WriteBits(0, 2);					
////cata end					
					for (int j = 0; j < 2 ; j++) //Decimation and regulation
					{
						int Number_symbols = this->NumberSymbol[Component_ID].front();
						this->NumberSymbol[Component_ID].pop_front();

						int Connectivity_Range;
						if (j == 0)	Connectivity_Range=1; 
						else Connectivity_Range=3;

						for (unsigned k = 0; k < (unsigned)Number_symbols; k++)	// MT
						{


							unsigned symbol = this->Connectivity[Component_ID].front();
							this->Connectivity[Component_ID].pop_front();			

							enc.encode_bits(symbol, Connectivity_Range);
							Connectivity_size += Connectivity_Range;

//							// To calculare connectivity rate
//							Connectivity_encoder.encode(symbol, Temp_connectivity);
														
							if (((j == 0) && (symbol != 1)) || ((j == 1) && (symbol != 4)))
							{
								Point_Int Coeff = this->Geometry[Component_ID].front();
								this->Geometry[Component_ID].pop_front();
												
								enc.encode_bits(Coeff.x, this->Qbit[Component_ID]);
								enc.encode_bits(Coeff.y, this->Qbit[Component_ID]);						
								enc.encode_bits(Coeff.z, this->Qbit[Component_ID]);						
																
							
								
								//#ifdef PREDICTION_METHOD
								if ((this->IsColored) && (!this->IsOneColor))
								{
									Color_Unit VC = this->VertexColor[Component_ID].front();
									this->VertexColor[Component_ID].pop_front();

									enc.encode_bits(VC.c0, HEADER_L_COLOR_QUANT);
									enc.encode_bits(VC.c1, HEADER_L_COLOR_QUANT);
									enc.encode_bits(VC.c2, HEADER_L_COLOR_QUANT);
									Color_size += (3 * (HEADER_L_COLOR_QUANT));

								}
								//#endif							
							}
						}
									
					}
					

//cata					
				//}
//cata end

			}			
		}

	}			

	//enc.SaveStream((LPSTR)File_Name);
	
	enc.Close();
	f_size = fopen(File_Name,"rb");
	fseek(f_size, 0, SEEK_END);
	Total_size = ftell(f_size);
	fclose(f_size);
}



void Compression_Valence_Web_Component::Calculate_Edge_Color_Difference(Polyhedron & _pMesh, 
																	const int & _Component_ID, 
																	double & _Max_color, 
																	double & _Mean_color,
																	int & Number_of_vertices)
{	

	float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;

	float Color_small_step = 0.0;

	Color_small_step = this->Color_Quantization_Step;

	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge_web != 2*_Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2*_Component_ID+1))
		hi++;

	// Vertex_Flag_web est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag_web = FREE;
		pVert->Vertex_Number_web = -1;
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
		
		if (v->Vertex_Flag_web == CONQUERED)
			continue;		
		else
		{			
			v->Vertex_Flag_web = CONQUERED;	
			v->Vertex_Number_web = Count_vertices;			
						
			Halfedge_around_vertex_circulator hvc = v->vertex_begin();
			Halfedge_around_vertex_circulator phvc = hvc;			
			
			double Mean = 0.0;
			int Count = 0;
			CGAL_For_all(hvc, phvc)
			{
				if(hvc->opposite()->vertex()->Vertex_Flag_web == FREE)
				{
					Color_Unit Color_0, Color_1;
					Color_0.c0 = hvc->vertex()->color_int_web(0);
					Color_0.c1 = hvc->vertex()->color_int_web(1);
					Color_0.c2 = hvc->vertex()->color_int_web(2);

					Color_1.c0 = hvc->opposite()->vertex()->color_int_web(0);
					Color_1.c1 = hvc->opposite()->vertex()->color_int_web(1);
					Color_1.c2 = hvc->opposite()->vertex()->color_int_web(2);

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
				if (h->opposite()->vertex()->Vertex_Flag_web == FREE)
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
void Compression_Valence_Web_Component::Diminush_Color_Quantization_Precision(Polyhedron &_pMesh, const int _Component_ID)
{

	//int N = 16; // Number of neighboring colors to use;
	//float Color_diff_seuil = 30.0;

	/*float C0_min = this->C0_Min;
	float C1_min = this->C1_Min;
	float C2_min = this->C2_Min;*/

	float Color_small_step = 0.0;

		Color_small_step = this->Color_Quantization_Step;		

	double Color_large_step = Color_small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();		
	
	while((hi->vertex()->Seed_Edge_web != 2*_Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2*_Component_ID+1))
		hi++;

	// Vertex_Flag_web est donnee free a tous les sommets
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag_web = FREE;
		pVert->Vertex_Number_web = -1;
	}	

	//premiere passe pour sous quantifier et donne l'indice de symbol a chaque sommet
	for(Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		if (pVert->Component_Number_web == _Component_ID)
		{
			int LAB_init_C0 = pVert->color_int_web(0);
			int LAB_init_C1 = pVert->color_int_web(1);
			int LAB_init_C2 = pVert->color_int_web(2);				
			
			// Les coordonnees apres under-quantization
			int LAB_final_C0 = LAB_init_C0 / 2;
			int LAB_final_C1 = LAB_init_C1 / 2;
			int LAB_final_C2 = LAB_init_C2 / 2;
			
			int i = LAB_init_C0 % 2;
			int j = LAB_init_C1 % 2;
			int k = LAB_init_C2 % 2;
			
			int Q_index = Get_Correct_Vector(i, j, k);
			
			pVert->color_int_web(LAB_final_C0, LAB_final_C1, LAB_final_C2);
			pVert->Q_Index_web = Q_index;

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
		
		if (v->Vertex_Flag_web == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag_web = CONQUERED;			
			v->Vertex_Number_web = Vertex_index;

			int Q_index = v->Q_Index_web;		
			vector<int> LAB;
		
			LAB.push_back(v->color_int_web(0));
			LAB.push_back(v->color_int_web(1));
			LAB.push_back(v->color_int_web(2));

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
					if(h->opposite()->vertex()->Vertex_Number_web > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number_web;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if(h->opposite()->vertex()->Vertex_Number_web == Comp_number)
						break;
				}
			}

			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag_web == FREE)
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

void Compression_Valence_Web_Component::Recalculate_Component_Area(Polyhedron & _pMesh, const int & _Component_ID, int & Number_facets)
{
	double Area = 0.0;	
	Number_facets = 0;
	for(Facet_iterator pFacet = _pMesh.facets_begin(); pFacet != _pMesh.facets_end(); pFacet++)
	{
		if(pFacet->Component_Number_web == _Component_ID)
		{
			Number_facets += 1;
			Area += Area_Facet_Triangle(pFacet->halfedge());
		}
	}
	Area *= pow(((double)10.0/(double)this->HighestLengthBB[_Component_ID]), 2.0);
	this->ComponentArea[_Component_ID] = Area;
}


double Compression_Valence_Web_Component::Calculate_Area(Polyhedron & _pMesh)
{	
	double mesh_area = 0;
	for(Facet_iterator pFacet = _pMesh.facets_begin(); pFacet != _pMesh.facets_end(); pFacet++)
	{
		Halfedge_handle h = pFacet->halfedge();
		mesh_area += Area_Facet_Triangle(h);
	}

	return mesh_area;
}


/*  Description : ADAPTIVE_QUANTIZATION
Decreasing of quantization resolution based on the prediction of PENG.
Opposite function is up_quantization. */
void Compression_Valence_Web_Component::Diminush_Geometry_Quantization_Precision(Polyhedron &_pMesh, const int & Component_ID)
{		

	// stock three mins for the reconstruction.
	float _xmin = this->xmin[Component_ID];
	float _ymin = this->ymin[Component_ID];
	float _zmin = this->zmin[Component_ID];

	// Premier_Pas == distance d'un Quantization_Step de grill de quantification (Q)
	double Small_step = 0.0;
			
		Small_step = this->Quantization_Step[Component_ID];		
		
	// Large_step == distance d'un Quantization_Step de grille de quantification(Q - 1)
	double Large_step = Small_step * 2;	
	
	// To find first points to start the conquest.	
	Halfedge_iterator hi = _pMesh.halfedges_begin();	
	
	while((hi->vertex()->Seed_Edge_web != 2*Component_ID) || (hi->opposite()->vertex()->Seed_Edge_web != 2*Component_ID+1))
		hi++;

	// Vertex_Flag_web est donnee free a tous les sommets
	for (Vertex_iterator pVert = _pMesh.vertices_begin();pVert != _pMesh.vertices_end(); pVert++)
	{
		pVert->Vertex_Flag_web = FREE;
		pVert->Vertex_Number_web = -1;
	}

	
	//premiere passe pour sous quantifie et donne l'indice de symbol a chaque sommet
	for (Vertex_iterator pVert = _pMesh.vertices_begin(); pVert != _pMesh.vertices_end(); pVert++)
	{
		if (pVert->Component_Number_web == Component_ID)
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
			pVert->Q_Index_web = Q_index;

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
		
		if (v->Vertex_Flag_web == CONQUERED)
			continue;
		
		else
		{
			v->Vertex_Flag_web = CONQUERED;			
			v->Vertex_Number_web = Vertex_index;			

			int Q_index = v->Q_Index_web;			
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
					if (h->opposite()->vertex()->Vertex_Number_web > Comp_number)
						Comp_number = h->opposite()->vertex()->Vertex_Number_web;
				}

				h = h2;
				CGAL_For_all(h,h2)
				{
					if (h->opposite()->vertex()->Vertex_Number_web == Comp_number)
						break;
				}
			}	
			
			Halfedge_around_vertex_circulator h2 = h;
			CGAL_For_all(h,h2)
			{
				if (h->opposite()->vertex()->Vertex_Flag_web == FREE)
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


QString Compression_Valence_Web_Component::Decompress_Init(Polyhedron &_pMesh)//, unsigned & _Initial_file_size, const char* File_Name)
{

	if (FILE *file2 = fopen(this->File_name.c_str(), "r"))
	{
		fseek(file2,0,SEEK_END);
		this->Compressed_file_size = ftell(file2)*8;
		fclose(file2);
	}	
	//calculate header size
	this->header_size = 0; 

	_pMesh.clear();
	this->IsClosed.clear();
	this->Qbit.clear();
	this->xmin.clear();
	this->ymin.clear();
	this->zmin.clear();
	this->Quantization_Step.clear();
	this->ComponentOperations.clear();
	
	this->Decompress_count = 0;	
	
	FILE* fp = fopen(this->File_name.c_str(),"rb");
	

	int res;		


	res=fread(&this->Initial_file_size, sizeof(unsigned), 1, fp);	 
	res=fread(&this->NumberComponents, sizeof(int), 1, fp);
	this->header_size+=sizeof(unsigned);
	this->header_size+= (sizeof(int));

	char Colored, One_color;
	res=fread(&Colored, sizeof(char), 1 , fp);		
	this->header_size+= sizeof(char);
	
	this->IsOneColor = false;
	if (Colored == 1)
	{
		this->IsColored = true;	
		res=fread(&One_color, sizeof(char), 1 , fp);
		this->header_size+= sizeof(char);

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
		this->header_size+= (sizeof(float)*4);
			
////cata
//		res=fread(&this->Smallest_C0, sizeof(int), 1, fp);  //smallest quantized positions
//		res=fread(&this->Smallest_C1, sizeof(int), 1, fp); 
//		res=fread(&this->Smallest_C2, sizeof(int), 1, fp);		
//		this->header_size+= (sizeof(int)*3);
////cata end
	}
	
	if ((this->IsColored) && (this->IsOneColor))
	{
		res=fread(&this->OnlyColor[0], sizeof(float), 1, fp); // smallest absolute position of c0 
		res=fread(&this->OnlyColor[1], sizeof(float), 1, fp); // smallest value of c1 
		res=fread(&this->OnlyColor[2], sizeof(float), 1, fp); // smallest value of c2
		this->header_size+= (sizeof(float)*3);
}	
	float Qpas;
	float t_xmin, t_ymin, t_zmin;
	
	// Read geometry information for each component;
	for (int i =0; i < this->NumberComponents; i++)
	{
		res=fread(&Qpas, sizeof(float), 1, fp); // quantization step
		res=fread(&t_xmin, sizeof(float), 1, fp); // x_min
		res=fread(&t_ymin, sizeof(float), 1, fp); // y_min
		res=fread(&t_zmin, sizeof(float), 1, fp); // z_min
		this->header_size+= (sizeof(float)*4);

		this->Quantization_Step.push_back(Qpas);
		this->xmin.push_back(t_xmin);
		this->ymin.push_back(t_ymin);
		this->zmin.push_back(t_zmin);
	}
	

	fclose(fp);
	
	Decoder.Open(this->File_name.c_str(),BF_READ);
//		LoadStream(this->File_name.c_str(),);
	int headers;
	for(unsigned int h=0;h<this->header_size;h++)
		headers= Decoder.decode_bits(8);
	
//cata
// //	this->Decoder.set_buffer(AC_BUFFER);
//	this->Decoder.read_from_file(fp);
// cata end
// 
	res+=0; // just to avoid gcc 4.6 warning
	
	this->GlobalCountOperation = -1;
	
	unsigned Max_Qbit = 0;

	// To know if each component is colored or not, and closed of not.
	for (int i = 0; i < this->NumberComponents; i++)
	{
		if (Decoder.decode_bits(HEADER_L_CLOSED) == 0)
			this->IsClosed.push_back(true);					
		else
			this->IsClosed.push_back(false);	
				
		int Number_operation = Decoder.decode_bits(HEADER_L_COMPONENT_OPERATION);  // Number of total operations.
		
		this->ComponentOperations.push_back(Number_operation);
		
		if (Number_operation > this->GlobalCountOperation)
			this->GlobalCountOperation = Number_operation;
			
		int Qbit = Decoder.decode_bits(HEADER_L_QBIT); // Initial quantization bit of geometry
		this->Qbit.push_back(Qbit);
		
		if (this->Qbit[i] > Max_Qbit)
			Max_Qbit = this->Qbit[i];
	}

	
	int Number_basemesh_vertex = Decoder.decode_bits(HEADER_L_NUM_VERTICES); // Number of vertices of base mesh
	int Number_basemesh_facet  = Decoder.decode_bits(HEADER_L_NUM_FACETS); // Number of facets of base mesh
		
	// Vectors for generation of base mesh
	vector<Point3d> vlist;
	vector<int>     flist;
	vector<float>   clist;
	vector<int>     Color_index_list;


	int minval=2000;
	int maxval=-200;



	for (int i = 0; i < Number_basemesh_vertex; i++)
	{

		Point_Int Pt_int;
		Pt_int.x = Decoder.decode_bits(Max_Qbit); // Read geometry info
		Pt_int.y = Decoder.decode_bits(Max_Qbit);
		Pt_int.z = Decoder.decode_bits(Max_Qbit);

	std::ofstream file("base_up.txt",fstream::app);
	file<<Pt_int;
	file.close();

	
	 
		// All vertices have quantization precision of component 0
		// That'll be corrected below.
		Point3d Pt_real = Change_Int_Real(Pt_int, 0); 
		vlist.push_back(Pt_real);		

		if ((this->IsColored) && (!this->IsOneColor))
		{
			//#ifdef PREDICTION_METHOD
		
			Color_Unit TC;
			TC.c0 = Decoder.decode_bits(HEADER_L_COLOR_QUANT); // Read color info.
			TC.c1 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);
			TC.c2 = Decoder.decode_bits(HEADER_L_COLOR_QUANT);

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
		int v = Decoder.decode_bits(Facet_index_bit);
		flist.push_back(v);
	}
	
	// Generation of base mesh using builder()
	CModifyBasemeshBuilder<HalfedgeDS, Polyhedron, Enriched_kernel> builder(&vlist, &flist, &clist, &Color_index_list);
	_pMesh.delegate(builder);
	
	_pMesh.compute_normals();

	// Seed Edges;
	map<int, int> Seed_Edge_webs;

	for (int i = 0; i < 2*this->NumberComponents; i++)
	{
		int Vertex_Number_web = Decoder.decode_bits(Facet_index_bit);
		Seed_Edge_webs.insert(pair<int,int>(Vertex_Number_web, i));
	}	

	int Basemesh_Vertex_Number_web = 0;
	
	Vertex_iterator pVertex = NULL;
	
	map<int,int>::iterator Seed_Edge_web_iterator = Seed_Edge_webs.begin();
	
	int Count_detected_vertices = 0;
	for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); Basemesh_Vertex_Number_web++, pVertex++)
	{	
		if (Count_detected_vertices < this->NumberComponents * 2)
		{
			if (Basemesh_Vertex_Number_web == Seed_Edge_web_iterator->first) 
			{
				pVertex->Seed_Edge_web = Seed_Edge_web_iterator->second;
				Seed_Edge_web_iterator++;
				Count_detected_vertices++;
			}	
		}
		else
			pVertex->Seed_Edge_web = OTHER_COORDINATE;

		pVertex->Component_Number_web = -1;
	}


	if ((this->IsColored) && (this->IsOneColor))
	{
		for (pVertex = _pMesh.vertices_begin(); pVertex != _pMesh.vertices_end(); pVertex++)
		{
			pVertex->color(this->OnlyColor[0], this->OnlyColor[1], this->OnlyColor[2]);
		}		
	}

	float Color_small_step = 0.0;	

	// To get number of vertices of each component and restore the real position of vertices
	//if (this->NumberComponents != 1)
	{
		_pMesh.tag_facets(-1);
		int Component_index = 0;
		
		for (int Component_Number_web = 0; Component_Number_web < this->NumberComponents; Component_Number_web++)
		{
			Halfedge_iterator hi = _pMesh.halfedges_begin();			
	
			while((hi->vertex()->Seed_Edge_web != 2 * Component_Number_web) || (hi->opposite()->vertex()->Seed_Edge_web != 2 * Component_Number_web + 1))
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
					if (pHalfedge->vertex()->Component_Number_web == -1)
					{
						pHalfedge->vertex()->Component_Number_web = Component_index;
						Point3d Wrong_position = pHalfedge->vertex()->point();
						
						// The correct position of vertex is restored.
						Point_Int Temp_pos = Change_Real_Int(Wrong_position, 0); 
						Point3d Real_position = Change_Int_Real(Temp_pos, Component_Number_web);

						pHalfedge->vertex()->point() = Real_position;
						
						float Wrong_LAB[3];
						Wrong_LAB[0] = pHalfedge->vertex()->color(0);
						Wrong_LAB[1] = pHalfedge->vertex()->color(1);
						Wrong_LAB[2] = pHalfedge->vertex()->color(2);												

						int Original_LAB[3];
						Original_LAB[0] = (int)floor((Wrong_LAB[0] - this->C0_Min)/this->Color_Quantization_Step + 0.5);
						Original_LAB[1] = (int)floor((Wrong_LAB[1] - this->C1_Min)/this->Color_Quantization_Step + 0.5);
						Original_LAB[2] = (int)floor((Wrong_LAB[2] - this->C2_Min)/this->Color_Quantization_Step + 0.5);
												
						pHalfedge->vertex()->color_int_web(Original_LAB[0], Original_LAB[1], Original_LAB[2]);						

							Color_small_step = this->Color_Quantization_Step;		

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


//cata new	
//		_pMesh.write_off("up_base.off",false,false);
//cata new end

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
int Compression_Valence_Web_Component::Decompress_Each_Step(Polyhedron &_pMesh, const char* File_Name)
{		
	if (this->Decompress_count < this->GlobalCountOperation)
	{
		for (int Component_ID = 0; Component_ID < this->NumberComponents; Component_ID++)
		{			
			if (this->Decompress_count < this->ComponentOperations[Component_ID])
			{	

	              this->Un_Regulation(_pMesh, Component_ID);					
				this->Un_Decimation_Conquest(_pMesh, Component_ID);
			}
		}			
	}	
		
			
	this->Decompress_count++;	
	_pMesh.compute_normals();

//	//cata new	
//		char file_q[70]="up_";
//		char file_number[70]="";
//		strcpy(file_number,"");
//		sprintf(file_number,"%.3d",this->Decompress_count);
//		_pMesh.write_off(strcat(strcat(file_q,file_number),".off"),false,false);
////cata new end


	return this->Decompress_count;
	
}





// Description : Change a point coordinates in real to integer coordinates
Point_Int Compression_Valence_Web_Component::Change_Real_Int(const Point3d &pt, const int & Component_ID)
{
	Point_Int Point;

	double Quantization_step = 0.0;

	Quantization_step = this->Quantization_Step[Component_ID];
	
	float xmin = this->xmin[Component_ID];
	float ymin = this->ymin[Component_ID];
	float zmin = this->zmin[Component_ID];
	
	double x = pt.x();
	double y = pt.y();
	double z = pt.z();

	Point.x = (int)(ceil((x-xmin)/Quantization_step)) - 1;
	if (Point.x == -1)
		Point.x = 0;
	else if (Point.x ==  (1<< this->Qbit[Component_ID]))
		Point.x--;

	Point.y = (int)(ceil((y-ymin)/Quantization_step)) - 1;
	if (Point.y == -1)
		Point.y = 0;
	else if (Point.y ==  (1<< this->Qbit[Component_ID]))
		Point.y--;

	Point.z = (int)(ceil((z-zmin)/Quantization_step)) - 1;
	if (Point.z == -1)
		Point.z = 0;
	else if (Point.z ==  (1<< this->Qbit[Component_ID]))
		Point.z--;

	return Point;
}



// Description : Change a point coordinates in integer to real coordinates
Point3d Compression_Valence_Web_Component::Change_Int_Real(const Point_Int &Point, const int & Component_ID)
{
	float Quantization_step = 0;
	
	Quantization_step = this->Quantization_Step[Component_ID];
		
	float xmin = this->xmin[Component_ID];
	float ymin = this->ymin[Component_ID];
	float zmin = this->zmin[Component_ID];	

	if (Point.y == 1023)
		int k=1;
	double x = xmin + (Point.x + 0.5) * Quantization_step;
	double y = ymin + (Point.y + 0.5) * Quantization_step;
	double z = zmin + (Point.z + 0.5) * Quantization_step;

	Point3d pt(x, y, z);

	return pt;
}




// To match first vertices flags between two meshes. Must be used when copying a mesh.
void Compression_Valence_Web_Component::Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh)
{
	Vertex_iterator pVertex = NULL;
	Vertex_iterator pVertex2 = New_mesh.vertices_begin();

	for (pVertex = Original_mesh.vertices_begin(); pVertex != Original_mesh.vertices_end(); pVertex++, pVertex2++)
	{
		//pVertex2->Seed_Edge_web = pVertex->Seed_Edge_web;
		
		int LAB[3];
		LAB[0] = pVertex->color_int_web(0);
		LAB[1] = pVertex->color_int_web(1);
		LAB[2] = pVertex->color_int_web(2);
		pVertex2->color_int_web(LAB[0], LAB[1], LAB[2]);

		float RGB[3];
		RGB[0] = pVertex->color(0);
		RGB[1] = pVertex->color(1);
		RGB[2] = pVertex->color(2);

		pVertex2->color(RGB[0], RGB[1], RGB[2]);
	}	
}


//// To reordering the color index for mapping table method.
//int Compression_Valence_Web_Component::Mapping_Table_Index_Reordering(Polyhedron &_pMesh)
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


//void Compression_Valence_Web_Component::Separate_Components(Polyhedron &_pMesh)
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

int Compression_Valence_Web_Component::GetResolutionChange(Polyhedron *_pMesh, float Prec)
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



// Error metric which measures importance of color and geometry for each vertex.
// Used to prevent removal of the visually important vertex.
bool Compression_Valence_Web_Component::Error_Projected_Surface(Polyhedron            &  _pMesh, 
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

		Color_step = this->Color_Quantization_Step;		

	vector<float> Center_color;
	Center_color.push_back(C0_min + g->next()->vertex()->color_int_web(0)*Color_step);
	Center_color.push_back(C1_min + g->next()->vertex()->color_int_web(1)*Color_step);
	Center_color.push_back(C2_min + g->next()->vertex()->color_int_web(2)*Color_step);

	vector<Color_Unit> Neighbors_color;
	
	vector<float> Projected_color;
	
	double Patch_area = 0.0;
	
	// Calculate area of patch
	g = g->next();
	for(int i =0; i < Valence; i++)
	{
		g = g->opposite();
		Color_Unit c;

		c.c0 = g->vertex()->color_int_web(0);
		c.c1 = g->vertex()->color_int_web(1);
		c.c2 = g->vertex()->color_int_web(2);
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

void Compression_Valence_Web_Component::Clear_After_Compression(void)
{
		

	for(unsigned i = 0; i < Connectivity.size(); i++)
		Connectivity.clear();
	for(unsigned i = 0; i < Geometry.size(); i++)
		Geometry.clear();
	for(unsigned i = 0; i < NumberSymbol.size(); i++)
		NumberSymbol.clear();
	for(unsigned i = 0; i < NumberVertices.size(); i++)
		NumberVertices.clear();
	for(unsigned i = 0; i < VertexColor.size(); i++)
		VertexColor.clear();
 
	

// Colors

	// Used for adatative quantization.				
	QuantizationCorrectVector.clear();
	NumberQuantizationLayer.clear();
		
		//for color
    NumberProcessedVertices.clear();
    ColorChildcellIndex.clear();
    ColorEncoderIndex.clear();

    Connectivity.clear();
    Geometry.clear();
    NumberSymbol.clear();
    NumberVertices.clear();

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
				
		Decompress_count = 0;

		 ComponentOperations.clear();

}




void Compression_Valence_Web_Component::Decompression_From_File(Polyhedron &_pMesh)
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

void Compression_Valence_Web_Component::Decompression_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh)
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

void Compression_Valence_Web_Component::Decompression_Specific_Level_From_File(Polyhedron &pMesh, const int & Wl)
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


void Compression_Valence_Web_Component::Decompression_Coarser_From_File(Polyhedron &pMesh)
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

void Compression_Valence_Web_Component::Decompression_All_From_File(Polyhedron &pMesh)
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


void Compression_Valence_Web_Component::Write_Info(Polyhedron &_pMesh)
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


QString Compression_Valence_Web_Component::Show_Text(void)
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
