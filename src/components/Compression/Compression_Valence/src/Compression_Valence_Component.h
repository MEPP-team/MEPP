#ifndef Compression_Valence_COMPONENT_H
#define Compression_Valence_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "../../../../mepp/mepp_component.h"

#include "Compression_Valence_Polyhedron.h"
#include "arithmetic_codec.h"
#include "Compression_Valence_Basemesh_builder.h"

#include <queue>
#include <list>

// struct of integer coordinates

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \struct	Point_Int
///
/// \brief	Point.
///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Point_Int
{
	
	int x;///< The x coordinate
	int y;///< The y coordinate
	int z;///< The z coordinate

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	Point_Int()
	///
	/// \brief	Default constructor.
	///
	////////////////////////////////////////////////////////////////////////////////////////////////////

	Point_Int()
	{
		x=y=z=0;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	const Point_Int operator+(const Point_Int &Pt) const
	///
	/// \brief	Addition operator.
	///
	/// \param	Pt	The point.
	///
	/// \return	The result of the operation.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	const Point_Int operator+(const Point_Int &Pt) const
	{
		Point_Int Temp;
		Temp.x = x + Pt.x;
		Temp.y = y + Pt.y;
		Temp.z = z + Pt.z;

		return Temp;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	const Point_Int operator-(const Point_Int &Pt) const
	///
	/// \brief	Negation operator.
	///
	/// \param	Pt	The point.
	///
	/// \return	The result of the operation.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	const Point_Int operator-(const Point_Int &Pt) const
	{
		Point_Int Temp;
		Temp.x = x - Pt.x;
		Temp.y = y - Pt.y;
		Temp.z = z - Pt.z;

		return Temp;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator==(const Point_Int &Pt) const
	///
	/// \brief	Equality operator.
	///
	/// \param	Pt	The point.
	///
	/// \return	true if the parameters are considered equivalent.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator ==(const Point_Int &Pt) const
	{
		return (x == Pt.x && y == Pt.y && z == Pt.z);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator!=(const Point_Int &Pt) const
	///
	/// \brief	Finaliser.
	///
	///
	/// \param	Pt	The point.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator !=(const Point_Int &Pt) const
	{
		return !(*this == Pt);
	}
};

 

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \struct	Color_Unit
///
/// \brief	Struct of integer color components. 
///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Color_Unit
{
	int c0;
	int c1;
	int c2;
	Color_Unit()
	{
		c0=c1=c2=0;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	const Color_Unit operator+(const Color_Unit &Col) const
	///
	/// \brief	Addition operator.
	///
	///
	/// \param	Col	The col.
	///
	/// \return	The result of the operation.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	const Color_Unit operator+(const Color_Unit &Col) const
	{
		Color_Unit Temp;
		Temp.c0 = c0 + Col.c0;
		Temp.c1 = c1 + Col.c1;
		Temp.c2 = c2 + Col.c2;
		
		return Temp;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	const Color_Unit operator-(const Color_Unit &Col) const
	///
	/// \brief	Negation operator.
	///
	/// \param	Col	The col.
	///
	/// \return	The result of the operation.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	const Color_Unit operator-(const Color_Unit &Col) const
	{
		Color_Unit Temp;
		Temp.c0 = c0 - Col.c0;
		Temp.c1 = c1 - Col.c1;
		Temp.c2 = c2 - Col.c2;
		
		return Temp;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator==(const Color_Unit &Col) const
	///
	/// \brief	Equality operator.
	///
	///
	/// \param	Col	The col.
	///
	/// \return	true if the parameters are considered equivalent.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator ==(const Color_Unit &Col) const
	{
		return( c0 == Col.c0 && c1 == Col.c1 && c2 == Col.c2);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator!=(const Color_Unit &Col) const
	///
	/// \brief	Finaliser.
	///
	/// \param	Col	The col.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator !=(const Color_Unit &Col) const
	{
		return!(*this == Col);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \class	Compression_Valence_Component
///
/// \brief	Compression valence component. 
///
////////////////////////////////////////////////////////////////////////////////////////////////////

class Compression_Valence_Component : 

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \property	public mepp_component
  ///
  /// \brief	Gets the mepp component.
  ///
  /// \value	.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  public mepp_component
{
	public:
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Compression_Valence_Component(Viewer* v, PolyhedronPtr p)
		///
		/// \brief	Default Constructor.
		///
		/// \param [in,out]	v	If non-null, the v.
		/// \param	p		 	The.
		////////////////////////////////////////////////////////////////////////////////////////////////////

		Compression_Valence_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
		{					
			
			
			IsOneColor = false;  
			Decompress_count = 0; 			

			Smallest_Alpha = 5000;
			Smallest_Gamma = 5000;
			
			Smallest_C0 = 5000;
			Smallest_C1 = 5000;
			Smallest_C2 = 5000;

			ColorDiffMinC0 = 5000;
			ColorDiffMinC1 = 5000;
			ColorDiffMinC2 = 5000;
			
			//IsClosed = false;
			IsColored = false;			
			
			GlobalCountOperation = 0;

			TotalBits = 0;		

			// from IHM
			IsCompressed = false;
			IsDecompress = false;
			Afficher = false;

			Possible_change_sequence = true;
			Sequence = true;
			Visu_level = 0;
			Process_level = 0;
			Writing_level = -1;

			m_Mode = 0;

			N_Inserted_Watermarks = 0;

			m_VC[0] = 0.;
			m_VC[1] = 0.;
			m_VC[2] = 0.;
			m_Dist = 0.0005;

			Number_non_reversible_vertices = 0;

			// MEPP 2
			componentName = "Compression_Valence_Component";
			init = 1;
		}

		/**
		 \fn	~Compression_Valence_Component()
		
		 \brief	Finaliser.
		
		 */

		~Compression_Valence_Component() 
		{}
		

		// Main Function

		/**
		 \fn	double Main_Function(Polyhedron &pMesh, const char* File_Name,const int &Qbit,
		 		const int &NVertices, const bool Normal_flipping, const bool Use_metric,
		 		const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,
		 		const bool Compression_selected, const bool Adaptive_quantization,unsigned &Number_layers,
		 		unsigned &Init_number_vertices,unsigned &Final_number_vertices, unsigned &Connectivity_size,
		 		unsigned & Color_size, unsigned & Total_size, const unsigned & Initial_file_size);
		
		 \brief	Main function.
		
		
		 \param [in,out]	pMesh				 	The mesh.
		 \param	File_Name						 	Filename of the file.
		 \param	Qbit							 	The qbit.
		 \param	NVertices						 	The vertices.
		 \param	Normal_flipping					 	The normal flipping.
		 \param	Use_metric						 	The use metric.
		 \param	Metric_thread					 	The metric thread.
		 \param	Use_forget_metric				 	The use forget metric.
		 \param	Forget_value					 	The forget value.
		 \param	Compression_selected			 	The compression selected.
		 \param	Adaptive_quantization			 	The adaptive quantization.
		 \param [in,out]	Number_layers		 	Number of layers.
		 \param [in,out]	Init_number_vertices 	The initialise number vertices.
		 \param [in,out]	Final_number_vertices	The final number vertices.
		 \param [in,out]	Connectivity_size	 	Size of the connectivity.
		 \param [in,out]	Color_size			 	Size of the color.
		 \param [in,out]	Total_size			 	Size of the total.
		 \param	Initial_file_size				 	Size of the initial file.
		
		 \return	.
		 */

		double Main_Function(Polyhedron &pMesh, const char* File_Name,const int &Qbit,const int &NVertices, const bool Normal_flipping,
			                 const bool Use_metric, const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const bool Compression_selected,
						     const bool Adaptive_quantization,unsigned &Number_layers, unsigned &Init_number_vertices,unsigned &Final_number_vertices, 
						     unsigned &Connectivity_size, unsigned & Color_size, unsigned & Total_size, const unsigned & Initial_file_size);


		//Initialization

		/**
		 \fn	void Global_Initialization(Polyhedron &pMesh, const int & qbit, const char * File_name);
		
		 \brief	Global initialization to select the input gate.
		
		\param [in,out]	pMesh	The mesh.
		 \param	qbit			 	The qbit.
		 \param	File_name		 	Filename of the file.
		 */

		void Global_Initialization(Polyhedron &pMesh, const int & qbit, const char * File_name);

			/**
			 \fn	void Quantization(Polyhedron &pMesh);
			
			 \brief	Quantize all vertices so that the new positions are reguliraly spaced in the 3D space..
			
			 \param [in,out]	pMesh	The mesh.
			 */

			void Quantization(Polyhedron &pMesh);

			/**
			 \fn	void Color_Initialization(Polyhedron &pMesh);
			
			 \brief	Color initialization.
			
			 \param [in,out]	pMesh	The mesh.
			 */

			void Color_Initialization(Polyhedron &pMesh);

			/**
			 \fn	void Color_Quantization(Polyhedron &pMesh);
			
			 \brief	Color quantization.
			
			 \param [in,out]	pMesh	The mesh.
			 */

			void Color_Quantization(Polyhedron &pMesh);		

			/**
			 \fn	void Multiple_Components_Initialization(Polyhedron & pMesh, const int & Qbit);
			
			 \brief	Multiple components initialization.
			
			 \param [in,out]	pMesh	The mesh.
			 \param	Qbit			 	The qbit.
			 */

			void Multiple_Components_Initialization(Polyhedron & pMesh, const int & Qbit);
		
		
		/**
		 \fn	void Simplification(Polyhedron &pMesh, const int & NVertices, const bool Normal_flipping,
		 		const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,
		 		const int &Forget_value);
		
		 \brief	Mesh Simplifications.
		
		 \param [in,out]	pMesh	The mesh.
		 \param	NVertices		 	The vertices.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 */

		void Simplification(Polyhedron &pMesh, const int & NVertices, const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value);		

			/**
			 \fn	int Decimation_Conquest(Polyhedron &pMesh,const bool Normal_flipping,
			 		const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,
			 		const int &Forget_value, const int & Component_ID);
			
			 \brief	selects a set of independent vertices to be removed using decimation conquest.
					
			 \param [in,out]	pMesh	The mesh.
			 \param	Normal_flipping  	The normal flipping.
			 \param	Use_metric		 	The use metric.
			 \param	Metric_thread	 	The metric thread.
			 \param	Use_forget_metric	The use forget metric.
			 \param	Forget_value	 	The forget value.
			 \param	Component_ID	 	Identifier for the component.
			
			 \return	.
			 */

			int Decimation_Conquest(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);													

			/**
			 \fn	int Regulation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,
			 		const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,
			 		const int & Component_ID);
			
			 \brief	regulation conquest.
						
			 \param [in,out]	pMesh	The mesh.
			 \param	Normal_flipping  	The normal flipping.
			 \param	Use_metric		 	The use metric.
			 \param	Metric_thread	 	The metric thread.
			 \param	Use_forget_metric	The use forget metric.
			 \param	Forget_value	 	The forget value.
			 \param	Component_ID	 	Identifier for the component.
			
			 \return	.
			 */
			int Regulation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);			
		

		// Adaptive Quantization		

		/**
		 \fn	void Adaptive_Quantization(Polyhedron &pMesh, const int & NVertices,
		 		const bool Normal_flipping,const bool Use_metric,const float &Metric_thread,
		 		const bool Use_forget_metric,const int &Forget_value,const int &Qbit);
		
		 \brief	Adaptive quantization. 
		
		 \param [in,out]	pMesh	The mesh.
		 \param	NVertices		 	The vertices.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 \param	Qbit			 	The qbit.
		 */

		void Adaptive_Quantization(Polyhedron &pMesh, const int & NVertices, const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,const int &Qbit);

		/**
		 \fn	int Test_Next_Operation(Polyhedron &pMesh,const bool Normal_flipping,
		 		const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,
		 		const int &Forget_value,const int &Qbit);
		
		 \brief	Tests next operation.
				
		 \param [in,out]	pMesh	The mesh.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 \param	Qbit			 	The qbit.
		
		 \return	.
		 */

		int  Test_Next_Operation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,const int &Qbit);	

		/**
		 \fn	void Up_Quantization(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & a);
		
		 \brief	Up quantization.
		
		 \param [in,out]	pMesh  	The mesh.
		 \param [in,out]	Decoder	The decoder.
		 \param	a				   	a.
		 */

		void Up_Quantization(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & a);

			/**
			 \fn	void Under_Quantization(Polyhedron &pMesh, const int & a);
			
			 \brief	Under quantization. Decreasing of quantization resolution based on the prediction of PENG.
Opposite function is up_quantization
					
			 \param [in,out]	pMesh	The mesh.
			 \param	a				 	a.
			 */


		void Under_Quantization(Polyhedron &pMesh, const int & a);

		/**
		 \fn	void Calculate_Edge_Color_Difference(Polyhedron & pMesh, const int & Component_ID,
		 		double & Max_color, double & Mean_color, int & Number_of_vertices);
		
		 \brief	Calculates the edge color difference.
		
		 \param [in,out]	pMesh			  	The mesh.
		 \param	Component_ID				  	Identifier for the component.
		 \param [in,out]	Max_color		  	The maximum color.
		 \param [in,out]	Mean_color		  	The mean color.
		 \param [in,out]	Number_of_vertices	Number of vertices.
		 */

		void Calculate_Edge_Color_Difference(Polyhedron & pMesh, const int & Component_ID, double & Max_color, double & Mean_color, int & Number_of_vertices);

		/**
		 \fn	void Color_Under_Quantization(Polyhedron &pMesh, const int Component_ID);
		
		 \brief	Color under quantization.
				
		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Identifier for the component.
		 */

		void Color_Under_Quantization(Polyhedron &pMesh, const int Component_ID);

		/**
		 \fn	void Color_Up_Quantization(Polyhedron &pMesh, Arithmetic_Codec & Decoder,
		 		const int & Component_ID);
		
		 \brief	Color up quantization.
		
		 \param [in,out]	pMesh  	The mesh.
		 \param [in,out]	Decoder	The decoder.
		 \param	Component_ID	   	Identifier for the component.
		 */

		void Color_Up_Quantization(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
			// Compression 

		/**
		 \fn	void Compression(Polyhedron &pMesh,const char* File_Name, const int &Qbit,
		 		unsigned &Connectivity_size, unsigned & Color_size, unsigned & Total_size,
		 		const unsigned & Initial_file_size);
		
		 \brief	Compressions.
			
		 \param [in,out]	pMesh			 	The mesh.
		 \param	File_Name					 	Filename of the file.
		 \param	Qbit						 	The qbit.
		 \param [in,out]	Connectivity_size	Size of the connectivity.
		 \param [in,out]	Color_size		 	Size of the color.
		 \param [in,out]	Total_size		 	Size of the total.
		 \param	Initial_file_size			 	Size of the initial file.
		 */

		void Compression(Polyhedron &pMesh,const char* File_Name, const int &Qbit, unsigned &Connectivity_size, unsigned & Color_size, unsigned & Total_size, const unsigned & Initial_file_size);				

		/**
		 \fn	int Calculate_Connectivity_Rate(void);
		
		 \brief	Calculates the connectivity rate.
		

		 \return	The calculated connectivity rate.
		 */

		int  Calculate_Connectivity_Rate(void);

		/**
		 \fn	void Calculate_Geometry_Color_Offset_Range();
		
		 \brief	Calculates the geometry color offset range.
		
		 */
		void Calculate_Geometry_Color_Offset_Range();				

		/**
		 \fn	void Remove_Last_Phase_Elements(const int & Component_ID);
		
		 \brief	Removes the last phase elements described by Component_ID.
			
		 \param	Component_ID	Identifier for the component.
		 */

		void Remove_Last_Phase_Elements(const int & Component_ID);

		/**
		 \fn	void Write_Base_Mesh(Polyhedron &pMesh, Arithmetic_Codec & Enc,
		 		unsigned &Connectivity_size, unsigned &Color_size, const int &Num_color_base_mesh);
		
		 \brief	Writes a base mesh.
			
		 \param [in,out]	pMesh			 	The mesh.
		 \param [in,out]	Enc				 	The encode.
		 \param [in,out]	Connectivity_size	Size of the connectivity.
		 \param [in,out]	Color_size		 	Size of the color.
		 \param	Num_color_base_mesh			 	Number of color base meshs.
		 */

		void Write_Base_Mesh(Polyhedron &pMesh, Arithmetic_Codec & Enc, unsigned &Connectivity_size, unsigned &Color_size, const int &Num_color_base_mesh);
		

		// Decompression		

		/**
		 \fn	int Decompress_Init(Polyhedron &pMesh,unsigned & Initial_file_size,
		 		const char* File_Name);
		
		 \brief	Initialize the Decompression by loading the base mesh into pMesh.
				
		 \param [in,out]	pMesh			 	The mesh.
		 \param [in,out]	Initial_file_size	Size of the initial file.
		 \param	File_Name					 	Filename of the file.
		
		 \return	.
		 */

		int    Decompress_Init(Polyhedron &pMesh,unsigned & Initial_file_size, const char* File_Name);

		/**
		 \fn	double Decompress_All(Polyhedron &pMesh);
		
		 \brief	Decompresses all the progressive levels described by pMesh. Just gives the final mesh.
		
		 \param [in,out]	pMesh	The mesh.
		
		 \return	.
		 */

		double Decompress_All(Polyhedron &pMesh);				

		/**
		 \fn	int Decompress_Each_Step(Polyhedron &pMesh, const char* File_Name);
		
		 \brief	Decompress the each step to visualize intermediate meshes
		
		 \param [in,out]	pMesh	The mesh.
		 \param	File_Name		 	Filename of the file.
		
		 \return	.
		 */

		int    Decompress_Each_Step(Polyhedron &pMesh, const char* File_Name);

		/**
		 \fn	int Decompress_To_Level(Polyhedron &pMesh,const int Wanted_level);
		
		 \brief	Decompress to level.
		
		
		 \param [in,out]	pMesh	The mesh.
		 \param	Wanted_level	 	The wanted level.
		
		 \return	.
		 */

		int    Decompress_To_Level(Polyhedron &pMesh,const int Wanted_level);

		/**
		 \fn	void Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & dec, const int & id);
		
		 \brief	Decoding of the decimation conquest.
			
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	dec  	The decrement.
		 \param	id				 	The identifier.
		 */

		void   Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & dec, const int & id);		

		/**
		 \fn	void Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);
		
		 \brief	Decoding of the regulation conquest.
				
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	dec  	The decrement.
		 \param	id				 	The identifier.
		 */

		void   Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);		

		// Other functions
		void   Separate_Components(Polyhedron &pMesh);		

		/**
		 \fn	double Calculate_Area(Polyhedron & pMesh);
		
		 \brief	Calculates the area of pMesh.
				
		 \param [in,out]	pMesh	The mesh.
		
		 \return	The calculated area.
		 */

		double Calculate_Area(Polyhedron & pMesh);
		bool   Error_Projected_Surface(Polyhedron & pMesh, const Halfedge_handle & _h, const int & _Component_ID, const double & Mean_color, const double & Mean_area);
		void   Recalculate_Component_Area(Polyhedron & pMesh, const int & Component_ID, int & Number_facets);

		// 

		/**
		 \fn	Point_Int Change_Real_Int(const Point3d &pt, const int & Component_ID);
		
		 \brief	Change from real to int point(related to quantization bit)
		
		 
		 \param	pt				The point.
		 \param	Component_ID	Identifier for the component.
		
		 \return	.
		 */

		Point_Int Change_Real_Int(const Point3d &pt, const int & Component_ID);		

		/**
		 \fn	Point3d Change_Int_Real(const Point_Int &vec, const int & Component_ID);
		
		 \brief	Change from int to real point(related to quantization bit)
		
		 \param	vec				The vector.
		 \param	Component_ID	Identifier for the component.
		
		 \return	.
		 */

		Point3d   Change_Int_Real(const Point_Int &vec, const int & Component_ID);
		
		
		/**
		 \fn	void Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh);
		
		 \brief	To add flags to copied mesh

		
		 \param [in,out]	Original_mesh	The original mesh.
		 \param [in,out]	New_mesh	 	The new mesh.
		 */

		void Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh);
		


		//REMOVE THESE Functions and variables
//REMOVE for cleaning
//		void Color_Metric_Roy(Polyhedron &pMesh, double &min, double &max, double &mean, double &rms);
//		Mesh_roy *Original;
//		void Set_Original_Color_Mesh(Polyhedron &pMesh, const char * Filename);
//		#ifdef MEASURE_COLOR_DEVIATION		
		// Distance calculation
//		void Calculate_Distances(char * originalMesh, char * attackedMesh, double & mrms, double & mrmswrtBB, double & hausdorff, double & hausdorffwrtBB);		
//		void Calculate_Color_Distances(char * originalMesh, char * attackedMesh, double & min, double & max, double & mean);		

		int Mapping_Table_Index_Reordering(Polyhedron &pMesh);



		/**
		 \fn	int Coding_Cost_Color_Simple_Prediction();
		
		 \brief	Color coding cost - prediction method.
		
		
		 \return	.
		 */

		int Coding_Cost_Color_Simple_Prediction();	
		
		
		/**
		 \fn	int Coding_Cost_Each_Layer_Color_Clustering(Polyhedron &pMesh);
		
		 \brief	Color coding cost - prediction method

		 \param [in,out]	pMesh	The mesh.
		
		 \return	.
		 */
		int Coding_Cost_Each_Layer_Color_Clustering(Polyhedron &pMesh);	
		
	 
		/**
		 \fn	unsigned Calculate_Current_File_Size(void)
		
		 \brief	Calculates the current file size. Measure bits used for decompression.
		
	 	
		 \return	The calculated current file size.
		 */

		unsigned Calculate_Current_File_Size(void) 
		{
			// To measure exact quantity of bits used for decompression.
			unsigned Adjust_value; 
			if (this->IsColored)
				Adjust_value = 25 + 9 * 4;
			else
				Adjust_value = 25;

			return this->Decoder.calculate_current_decoded_size() + Adjust_value;
		}

		/**
		 \fn	int GetResolutionChange(Polyhedron *pMesh, float Prec);
		
		 \brief	Gets a resolution change.
				
		 \param [in,out]	pMesh	If non-null, the mesh.
		 \param	Prec			 	The prec.
		
		 \return	The resolution change.
		 */

		int GetResolutionChange(Polyhedron *pMesh, float Prec);

		/**
		 \fn	void Stop_Decoder(void)
		
		 \brief	To stop the decoder (Used "Previous level" while decoding).
		
		 */

		void Stop_Decoder(void) {this->Decoder.stop_decoder(); }		

		//////////////////////////////////////////
		// Joint Compression Watermarking (JCW) //
		//////////////////////////////////////////

		/**
		 \fn	double Joint_Compression_Watermarking(Polyhedron &pMesh, const int &NVertices,
		 		const bool Normal_flipping, const bool Use_metric, const float & Metric_thread,
		 		const bool Use_forget_metric, const int &Forget_value, const int &Qbit, int & Number_bits,
		 		unsigned int & Connectivity_size,unsigned int & Color_size, unsigned int & Total_size,
		 		unsigned int & Initial_file_size);
		
		 \brief	Joint Compression Watermarking (JCW) 
		
		 \param [in,out]	pMesh			 	The mesh.
		 \param	NVertices					 	The vertices.
		 \param	Normal_flipping				 	The normal flipping.
		 \param	Use_metric					 	The use metric.
		 \param	Metric_thread				 	The metric thread.
		 \param	Use_forget_metric			 	The use forget metric.
		 \param	Forget_value				 	The forget value.
		 \param	Qbit						 	The qbit.
		 \param [in,out]	Number_bits		 	Number of bits.
		 \param [in,out]	Connectivity_size	Size of the connectivity.
		 \param [in,out]	Color_size		 	Size of the color.
		 \param [in,out]	Total_size		 	Size of the total.
		 \param [in,out]	Initial_file_size	Size of the initial file.
		
		 \return	.
		 */

		double Joint_Compression_Watermarking(Polyhedron &pMesh, const int &NVertices, const bool Normal_flipping, const bool Use_metric, const float & Metric_thread, const bool Use_forget_metric, const int &Forget_value, const int &Qbit, int & Number_bits,
											  unsigned int & Connectivity_size,unsigned int & Color_size, unsigned int & Total_size,unsigned int & Initial_file_size);

		/**
		 \fn	int JCW_Decimation_For_Segmentation(Polyhedron &pMesh,const bool Normal_flipping,
		 		const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,
		 		const int &Forget_value, const int & Component_ID);
		
		 \brief	JCW decimation for segmentation.
				
		 \param [in,out]	pMesh	The mesh.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 \param	Component_ID	 	Identifier for the component.
		
		 \return	.
		 */

		int JCW_Decimation_For_Segmentation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);													

		/**
		 \fn	int JCW_Regulation_For_Segmentation(Polyhedron &pMesh,const bool Normal_flipping,
		 		const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,
		 		const int &Forget_value, const int & Component_ID);
		
		 \brief	JCW regulation for segmentation.
		
		
		 \param [in,out]	pMesh	The mesh.
		 \param	Normal_flipping  	The normal flipping.
		 \param	Use_metric		 	The use metric.
		 \param	Metric_thread	 	The metric thread.
		 \param	Use_forget_metric	The use forget metric.
		 \param	Forget_value	 	The forget value.
		 \param	Component_ID	 	Identifier for the component.
		
		 \return	.
		 */

		int JCW_Regulation_For_Segmentation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);					

        /**
         \fn	void JCW_Un_Regulation_For_Insertion(Polyhedron &pMesh, const int & Component_ID,
         		list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position,
         		list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position,
         		list<vector<int> > & JCW_ERROR);
        
         \brief	JCW un regulation for insertion.
        
         \param [in,out]	pMesh				   	The mesh.
         \param	Component_ID					   	Identifier for the component.
         \param [in,out]	FP_Connectivity		   	The fp connectivity.
         \param [in,out]	SP_Moved_Position	   	The sp moved position.
         \param [in,out]	SP_Original_Position   	The sp original position.
         \param [in,out]	SP_Watermarked_Position	The sp watermarked position.
         \param [in,out]	JCW_ERROR			   	The jcw error.
         */

        void JCW_Un_Regulation_For_Insertion(Polyhedron &pMesh, const int & Component_ID, list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position, list<vector<int> > & JCW_ERROR);

		/**
		 \fn	void JCW_Un_Decimation_For_Insertion(Polyhedron &pMesh, const int & Component_ID,
		 		list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position,
		 		list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position,
		 		list<vector<int> > & JCW_ERROR);
		
		 \brief	JCW un decimation for insertion.
		
		 	
		 \param [in,out]	pMesh				   	The mesh.
		 \param	Component_ID					   	Identifier for the component.
		 \param [in,out]	FP_Connectivity		   	The fp connectivity.
		 \param [in,out]	SP_Moved_Position	   	The sp moved position.
		 \param [in,out]	SP_Original_Position   	The sp original position.
		 \param [in,out]	SP_Watermarked_Position	The sp watermarked position.
		 \param [in,out]	JCW_ERROR			   	The jcw error.
		 */

		void JCW_Un_Decimation_For_Insertion(Polyhedron &pMesh, const int & Component_ID, list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position, list<vector<int> > & JCW_ERROR);

		/**
		 \fn	Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h,
		 		const int & Direction);
		
		 \brief	JCW barycenter patch before removal.
			
		 \param	h		 	The.
		 \param	Direction	The direction.
		
		 \return	.
		 */

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & Direction);

		/**
		 \fn	Point3d JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence,
		 		const int & Direction);
		
		 \brief	JCW barycenter patch after removal.
		
		
		 \param	h		 	The.
		 \param	valence  	The valence.
		 \param	Direction	The direction.
		
		 \return	.
		 */

		Point3d JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);

		/**
		 \fn	Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h,
		 		const int & valence, const int & Direction);
		
		 \brief	JCW barycenter patch before removal.
		 
		
		 \param	h		 	The.
		 \param	valence  	The valence.
		 \param	Direction	The direction.
		
		 \return	.
		 */

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);

		/**
		 \fn	void JCW_Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & dec, const int & id);
		

		 \brief	JCW un decimation conquest.
		 	
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	dec  	The decrement.
		 \param	id				 	The identifier.
		 */

		void JCW_Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & dec, const int & id);		

		/**
		 \fn	void JCW_Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);
		
		 \brief	JCW un regulation.
		
	 	
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	dec  	The decrement.
		 \param	id				 	The identifier.
		 */

		void JCW_Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);		

		/**
		 \fn	void JCW_Un_Regulation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID,
		 		list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);
		
		 \brief	JCW un regulation for region detection, in order to obtain region number of each inserted vertices
		
 		 \param [in,out]	pMesh	  	The mesh.
		 \param	Component_ID		  	Identifier for the component.
		 \param [in,out]	FP_connect	The fp connect.
		 \param [in,out]	FP_Geo	  	The fp geo.
		 \param [in,out]	FP_RN	  	The fp rn.
		 */

		void JCW_Un_Regulation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID, list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);

		/**
		 \fn	void JCW_Un_Decimation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID,
		 		list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);
		
		 \brief	JCW un decimation for region detection.
		
	 	
		 \param [in,out]	pMesh	  	The mesh.
		 \param	Component_ID		  	Identifier for the component.
		 \param [in,out]	FP_connect	The fp connect.
		 \param [in,out]	FP_Geo	  	The fp geo.
		 \param [in,out]	FP_RN	  	The fp rn.
		 */

		void JCW_Un_Decimation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID, list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);

		/**
		 \fn	vector<double> JCW_Evaluate_Robustness(void);
		
		 \brief	Evaluates the JCW robustness.
		
		
		 \return	.
		 */

		vector<double> JCW_Evaluate_Robustness(void);

		/**
		 \fn	void Set_Number_Bin(const int & NB)
		
		 \brief	Sets the bin number to NB.
		
		 \param	NB	The nb.
		 */

		void Set_Number_Bin(const int & NB)
		{
			this->m_NumberBin = NB;
		}

		/**
		 \fn	int Get_Number_Bin(void)
		
		 \brief	Gets the bin number.
				
		 \return	The number bin.
		 */

		int Get_Number_Bin(void)
		{
			return this->m_NumberBin;
		}		

		/**
		 \fn	void Set_Embedding_Level(const int &EL)
		
		 \brief	Sets the embedding level.
		
	 	
		 \param	EL	The el.
		 */

		void Set_Embedding_Level(const int &EL)
		{
			this->m_EmbeddingLevel = EL;
		}

		/**
		 \fn	int Get_Embedding_Level(void)
		
		 \brief	Gets the embedding level.
		
	 	
		 \return	The embedding level.
		 */

		int Get_Embedding_Level(void)
		{
			return this->m_EmbeddingLevel;
		}

		/**
		 \fn	void Set_Number_Region(const int &NR)
		
		 \brief	Sets the region number.
		
			 \param	NR	The nr.
		 */

		void Set_Number_Region(const int &NR)
		{
			this->m_NumberRegion = NR;
		}

		/**
		 \fn	int Get_Number_Region(void)
		
		 \brief	Gets the region number.
		
		 \return	The number region.
		 */

		int Get_Number_Region(void)
		{
			return this->m_NumberRegion;
		}

		/**
		 \fn	void JCW_Calculate_Mesh_Center(Polyhedron &pMesh);
		
		 \brief	 calculates the mesh center for JCW.
				
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Calculate_Mesh_Center(Polyhedron &pMesh);		

		/**
		 \fn	void JCW_Calculate_Radius(Polyhedron &pMesh);
		
		 \brief	Calculateq the radius for JCW.
			
		 \param [in,out]	pMesh	The mesh.
		 */
		void JCW_Calculate_Radius(Polyhedron &pMesh);

		/**
		 \fn	void JCW_Expand_Mesh(Polyhedron &pMesh);
		
		 \brief	JCW expand mesh.
		
		 
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Expand_Mesh(Polyhedron &pMesh);

		/**
		 \fn	void JCW_Quantization(Polyhedron &pMesh);
		
		 \brief	quantization for JCW.
		
	 	
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Quantization(Polyhedron &pMesh);

		/**
		 \fn	void JCW_Generate_Regions_On_Base_Mesh(Polyhedron &pMesh);
		
		 \brief	Generate regions on base mesh for JCW.
		
	 	
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Generate_Regions_On_Base_Mesh(Polyhedron &pMesh);

        /**
         \fn	void JCW_Region_Mass_Center_Insert_Watermark(Polyhedron & pMesh,
         		list<Point3d> & FP_Geometry, list<int> & FP_Region_Number,
         		list<Point_Int> & SP_Watermarked_Position, list<Point3d> & SP_Moved_Position,
         		list<Point3d> & SP_Original_Position, list<vector<int> > & JCW_ERROR);
        
         \brief	JCW region mass center insert watermark.
             
         \param [in,out]	pMesh				   	The mesh.
         \param [in,out]	FP_Geometry			   	The fp geometry.
         \param [in,out]	FP_Region_Number	   	The fp region number.
         \param [in,out]	SP_Watermarked_Position	The sp watermarked position.
         \param [in,out]	SP_Moved_Position	   	The sp moved position.
         \param [in,out]	SP_Original_Position   	The sp original position.
         \param [in,out]	JCW_ERROR			   	The jcw error.
         */

        void JCW_Region_Mass_Center_Insert_Watermark(Polyhedron & pMesh, list<Point3d> & FP_Geometry, list<int> & FP_Region_Number, list<Point_Int> & SP_Watermarked_Position, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<vector<int> > & JCW_ERROR);

		/**
		 \fn	void JCW_Choose_Valid_Vertices(Polyhedron & pMesh);
		
		 \brief	JCW choose valid vertices.
	
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Choose_Valid_Vertices(Polyhedron & pMesh);

		/**
		 \fn	vector<int> JCW_Region_Mass_Center_Extract_Watermark(Polyhedron & pMesh);
		
		 \brief	JCW region mass center extract watermark.
		
		 \param [in,out]	pMesh	The mesh.
		
		 \return	.
		 */

		vector<int> JCW_Region_Mass_Center_Extract_Watermark(Polyhedron & pMesh);

		/**
		 \fn	void JCW_Code_Difference_Histogram_Shifting(Polyhedron &pMesh,const int & Component_ID);
		
		 \brief	JCW code difference histogram shifting.
				
		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Identifier for the component.
		 */

		void JCW_Code_Difference_Histogram_Shifting(Polyhedron &pMesh,const int & Component_ID);

		/**
		 \fn	void JCW_Decode_Difference_Histogram_Shifting(Polyhedron &pMesh,
		 		const int & Component_ID);
		
		 \brief	JCW decode difference histogram shifting.
		
		
		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Identifier for the component.
		 */

		void JCW_Decode_Difference_Histogram_Shifting(Polyhedron &pMesh, const int & Component_ID);

		/**
		 \fn	void Initialize_Spherical_Coordinates(Polyhedron &pMesh);
		
		 \brief	Initializes the spherical coordinates.
				
		 \param [in,out]	pMesh	The mesh.
		 */

		void Initialize_Spherical_Coordinates(Polyhedron &pMesh);

		/**
		 \fn	void Convert_To_Spherical(const Point3d & Pt, double * Spheric);
		
		 \brief	Converts this object to a spherical.
			
		 \param	Pt				   	The point.
		 \param [in,out]	Spheric	If non-null, the spheric.
		 */

		void Convert_To_Spherical(const Point3d & Pt, double * Spheric);

		/**
		 \fn	void Convert_To_Cartesian(const double * Spheric, double * Cartesian);
		
		 \brief	Converts this object to a cartesian.
		
		
		 \param	Spheric				 	The spheric.
		 \param [in,out]	Cartesian	If non-null, the cartesian.
		 */

		void Convert_To_Cartesian(const double * Spheric, double * Cartesian);	

		/**
		 \fn	int JCW_Decompress_One_Level(Polyhedron &pMesh, const char* File_Name,
		 		const int & Noise_mode);
		
		 \brief	Jcw decompress one level.
		
		 \param [in,out]	pMesh	The mesh.
		 \param	File_Name		 	Filename of the file.
		 \param	Noise_mode		 	The noise mode.
		
		 \return	.
		 */

		int  JCW_Decompress_One_Level(Polyhedron &pMesh, const char* File_Name, const int & Noise_mode);

		/**
		 \fn	int JCW_Decompress_One_Level_Without_Extraction(Polyhedron &pMesh, const char* File_Name);
		
		 \brief	Jcw decompress one level without extraction.
			
		 \param [in,out]	pMesh	The mesh.
		 \param	File_Name		 	Filename of the file.
		
		 \return	.
		 */

		int  JCW_Decompress_One_Level_Without_Extraction(Polyhedron &pMesh, const char* File_Name);

		/**
		 \fn	void JCW_Run(void);
		
		 \brief	Jcw run.
		
		 */

		void JCW_Run(void);

		/**
		 \fn	int JCW_Divide_Big_Regions(Polyhedron &pMesh);
		
		 \brief	Jcw divide big regions.
		
		
		 \param [in,out]	pMesh	The mesh.
		
		 \return	.
		 */

		int  JCW_Divide_Big_Regions(Polyhedron &pMesh);

		/**
		 \fn	void JCW_Colorify_Regions(Polyhedron & pMesh);
		
		 \brief	Jcw colorify regions.
		
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Colorify_Regions(Polyhedron & pMesh);

		/**
		 \fn	void Read_Information_To_Hide();
		
		 \brief	Reads the information to hide.
		*/

		void    Read_Information_To_Hide();

		/**
		 \fn	QString Write_Information_To_Hide();
		
		 \brief	Writes the information to hide.
				
		 \return	.
		 */

		QString Write_Information_To_Hide();

		/**
		 \fn	void Clear_After_Compression();
		
		 \brief	Clears after the compression.

		 */

		void Clear_After_Compression();
	private:
		
	
		Polyhedron * UnderquantizedMesh;	///< The underquantized mesh


		
		Polyhedron * DecimatedMesh;		///< The decimated mesh

		
		FILE * LogColor;///< The log color
		
		vector<bool> IsClosed;    ///< The is closed. To know if the mesh is open or closed.
		
		
		
		bool IsColored;///< true if is colored
		
		bool IsOneColor;///< true if is one color
		
		float OnlyColor[3];///< The only color
		// Number of connectivity symbol types. Without boundary = 5, with boundary = 7;
		//int NummberConnectivitySymbols;

		//CCopyPoly<Polyhedron, Enriched_kernel> Copy_Polyhedron;	//MT
		
		// To encode each type of operation between decimation and increase of quantization resolution
        vector<list<int> >		 ListOperation;///< The list operation
        vector<list<int> >		 Connectivity;///< The connectivity
        vector<list<Point_Int> > Geometry;///< The geometry
        vector<list<int> >		 NumberSymbol;///< Number of symbols
        vector<list<int> >		 NumberVertices;///< Number of vertices
        
			
		list<Point_Int>          InterGeometry;///< The inter geometry
		list<int>		         InterConnectivity;		///< The inter connectivity
		

        vector<list<int> > AlphaRange;///< The alpha range
        vector<list<int> > AlphaOffset;///< The alpha offset
        vector<list<int> > GammaRange;///< The gamma range
        vector<list<int> > GammaOffset;///< The gamma offset
        
		
		// Quantization
		
		vector<unsigned>	Qbit; ///< The qbit
		vector<float>		xmin;///< The xmin
		vector<float>		ymin;///< The ymin
		vector<float>		zmin;		///< The zmin
		vector<float>		xmax;///< The xmax
		vector<float>		ymax;///< The ymax
		vector<float>		zmax;		///< The zmax
		vector<float>		QuantizationPas;///< The quantization pas
		
		
		int			  TotalBits;		///< The total bits
		double		  OldDistortion;	///< The old distortion
		int			  Smallest_Alpha;		///< The smallest alpha
		int			  Smallest_Gamma;			///< The smallest gamma
			
		
		vector<double> HighestLengthBB;///< The highest length bb
		vector<double> ComponentVolume;///< The component volume
		vector<double> ComponentArea;///< The component area
		vector<int>    ComponentNumberVertices;///< The component number vertices
		
		
		// Used for adatative quantization.				
        vector<list<int> > QuantizationCorrectVector; ///< The quantization correct vector
        vector<list<int> > NumberQuantizationLayer; ///< Number of quantization layers
       
		
		//for color
        vector<list<int> > NumberProcessedVertices;///< Number of processed vertices
        vector<list<int> > ColorChildcellIndex;///< the color childcell index
        vector<list<int> > ColorEncoderIndex;///< the color encoder index
        

		// Colors
		vector<list<Color_Unit> > VertexColor; ///< Contain color error of all removed vertices
		list<Color_Unit> InterVertexColor;		///< The inter vertex color
		
		
		
		float C0_Min;///< The C0 minimum
		float C1_Min;///< The C1 minimum
		float C2_Min;///< The C2 minimum
		


		float Color_Quantization_Step;		///< The color quantization step
		vector<int> NumberColorQuantization;///< Number of color quantizations
		

		int Smallest_C0; ///< The smallest value of C0 used for preventing from negative sylbols
		int Smallest_C1; ///< The smallest value of C1 used for preventing from negative sylbols
		int Smallest_C2; ///< The smallest value of C2 used for preventing from negative sylbols
		
		
		int C0_Range;///< The C0 range
		int C1_Range;///< The C1 range
		int C2_Range;///< The C2 range
		
		int ColorDiffMinC0;///< The color difference minimum c0 for restoration of colors in Mapping table
		int ColorDiffMinC1;///< The first color difference minimum c1 for restoration of colors in Mapping table
		int ColorDiffMinC2;///< The second color difference minimum c2 for restoration of colors in Mapping table
		
		
		
		int ColorDiffRangeC0;///< The color difference range c 0 for restoration of colors in Mapping table
		int ColorDiffRangeC1;///< The first color difference range c1 for restoration of colors in Mapping table
		int ColorDiffRangeC2;///< The second color difference range c2 for restoration of colors in Mapping table
		
		
		// mapping table
		vector<Color_Unit> ColorArray; ///< Color table
		vector<Color_Unit> PredictedColorArray; ///< Predicted values of colors in color table using prediction based on neighbors.		
		vector<Color_Unit> DifferenceColor; ///< Difference between original and quantized vertex color. need to turn into lossless coding.
		
		list<int>		   ColorIndex;///< List of color index
		vector<int>		   ReorderingColorIndex;///< Vector of reordering color index
		list<int>		   InterColorIndex;		///< Vector of inter color index
		vector<int>		   Number_color_index; ///< contain occurence of each initial color
		vector<int>		   IsKnownIndex;///< is known index
		

		
		// Decoder
		
		Arithmetic_Codec Decoder;		///< The arithmetic decoder
		
		Adaptive_Data_Model Color_0_Model;///< The color 0 model
		Adaptive_Data_Model Color_1_Model;///< The color 1 model
		Adaptive_Data_Model Color_2_Model;	///< The color 2 model
		

		Adaptive_Data_Model Index_Model;///< The index model
		


		vector<int> NumberDecimation; ///< To stock number of Decimation.
		vector<int> NumberChangeQuantization; ///< to stock number of under_quantization.						
				
		int DumpSymbolDecimation;///< The dump symbol decimation in order to eliminate symbols of the last conquest if it is useless.
		int DumpSymbolRegulation;///< The dump symbol regulation in order to eliminate symbols of the last conquest if it is useless.
		

		int Decompress_count;///< Number of decompress in order to know how many steps to go for the decompression step by step.
		int NumberComponents;///< Number of components in order to know how many steps to go for the decompression step by step.
		

		int GlobalCountOperation;///< The global count operation
	
		vector<int> ComponentOperations;///< The component operations
		

		// JCW
		int   m_Mode;
		int   m_NumberBin;
		int   m_EmbeddingLevel;
		int   m_NumberRegion;
		double m_VC[3];
		double m_Rmin;
		double m_Rmax;
		double m_Dist;
		vector<int> m_Number_Vertices_Per_Regions;
		list<int> m_Watermarks;

		int N_Inserted_Watermarks;


		vector<int> m_N_remained_vertices;
		vector<int> m_N_treated_vertices;
		vector<double> m_Rad_decision;

        list<vector<int> > m_JCW_Move_Error;
		list<int> m_N_Errors;

		Adaptive_Data_Model DM_JCW_MOVE_ERROR;
		
		list<int> JCW_Connectivity;
		list<Point_Int> JCW_Geometry;

		//double LUT_CourbureClust[3*256];
		vector< vector<float> > Region_Color;
		int Number_non_reversible_vertices;
	// from IHM
	public:
		bool IsDecompress;///< true if is decompress
		bool IsCompressed;///< true if is compressed
		int Current_level;///< The current level
		int Total_layer;///< The total number of layer
		unsigned Initial_file_size;///< Size of the initial file
		unsigned Compressed_file_size;///< Size of the compressed file
		string File_name;///< Filename of the file
		
		vector<float> Prog;
		vector<float> Ratio;
		QString Message;

	public://private:
		CCopyPoly<Polyhedron, Enriched_kernel> Copy_Polyhedron;
		bool Afficher;
		bool Sequence;
		bool Possible_change_sequence;
		
		int Visu_level;
		int Process_level;		

		FILE *Dec_Info;
		int Writing_level;
		string Dec_File_Info;
};

#endif

#endif // Compression_Valence_COMPONENT_H
