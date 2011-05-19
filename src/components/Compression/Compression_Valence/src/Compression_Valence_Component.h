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
		Point_Int Res;
		Res.x = x + Pt.x;
		Res.y = y + Pt.y;
		Res.z = z + Pt.z;

		return Res;
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
		Point_Int Res;
		Res.x = x - Pt.x;
		Res.y = y - Pt.y;
		Res.z = z - Pt.z;

		return Res;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator==(const Point_Int &Pt) const
	///
	/// \brief	Comparison operatior ==.
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
	/// \brief	Comparison operatior !=.
	///
	///
	/// \param	Pt	The point.
	/// \return true if the parameters are not identical.
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
	/// \fn	const Color_Unit operator+(const Color_Unit & m_color) const
	///
	/// \brief	Addition operator.	
	///
	/// \param	m_color	The color.
	///
	/// \return	The result of the operation.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	const Color_Unit operator+(const Color_Unit & m_color) const
	{
		Color_Unit Res;
		Res.c0 = c0 + m_color.c0;
		Res.c1 = c1 + m_color.c1;
		Res.c2 = c2 + m_color.c2;
		
		return Res;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	const Color_Unit operator-(const Color_Unit &m_color) const
	///
	/// \brief	Negation operator.
	///
	/// \param	m_color	The color.
	///
	/// \return	The result of the operation.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	const Color_Unit operator-(const Color_Unit & m_color) const
	{
		Color_Unit Res;
		Res.c0 = c0 - m_color.c0;
		Res.c1 = c1 - m_color.c1;
		Res.c2 = c2 - m_color.c2;
		
		return Res;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator==(const Color_Unit & m_color) const
	///
	/// \brief	Comparison operatior ==.
	///	
	/// \param	Col	The col.
	///
	/// \return	true if the parameters are considered equivalent.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator ==(const Color_Unit &m_color) const
	{
		return( c0 == m_color.c0 && c1 == m_color.c1 && c2 == m_color.c2);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	bool operator!=(const Color_Unit & m_color) const
	///
	/// \brief	Comparison operatior !=.
	///
	/// \param	Col	The col.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator !=(const Color_Unit & m_color) const
	{
		return!(*this == m_color);
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

			/*ColorDiffMinC0 = 5000;
			ColorDiffMinC1 = 5000;
			ColorDiffMinC2 = 5000;*/
			
			IsColored = false;			
			GlobalCountOperation = 0;

			// from IHM
			IsCompressed = false;
			IsDecompress = false;
			//Afficher = false;

			Possible_change_sequence = true;
			Sequence = true;
			Visu_level = 0;
			Process_level = 0;
			//Writing_level = -1;
			
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
		 \brief	Destructor.		
		 */

		~Compression_Valence_Component() 
		{}
		

		// Main Function

		/**
		 \fn	double Main_Function(Polyhedron &pMesh, 
							 const char * Input_File_Name, 
							 const char* File_Name,
							 const int &Qbit,
							 const int &NVertices, 
							 const bool Is_normal_flipping_selected,
			                 const bool Is_use_metric_selected, 
							 const float &Metric_thread, 
							 const bool Is_use_forget_metric_selected,
							 const int &Forget_value, 
							 const bool Is_compression_selected,
							 const bool Is_adaptive_quantization_selected, 
							 const bool Is_bijection_selected);
		
		 \brief	Main function of compression.		
		
		 \param [in,out]	pMesh				 	The mesh.
		 \param	Input_File_Name						Filename of the input file.
		 \param	File_Name						 	Filename of the compressed(output) file.
		 \param	Qbit							 	The qbit.
		 \param	NVertices						 	The vertices.
		 \param	Is_normal_flipping_selected					 	The normal flipping.
		 \param	Is_use_metric_selected						 	The use metric.
		 \param	Metric_thread					 	The metric thread.
		 \param	Is_use_forget_metric_selected				 	The use forget metric.
		 \param	Forget_value					 	The forget value.
		 \param	Is_compression_selected			 	The compression selected.
		 \param	Is_adaptive_quantization_selected			 	The adaptive quantization.
		 \param	Is_bijection_selected			 	The bijection of the geometry coding.
		 \return	Information of compression results.
		 */

		QString Main_Function(Polyhedron &pMesh, 
							 const char * Input_File_Name, 
							 const char* File_Name,
							 const int &Qbit,
							 const int &NVertices, 
							 const bool Is_normal_flipping_selected,
			                 const bool Is_use_metric_selected, 
							 const float &Metric_thread, 
							 const bool Is_use_forget_metric_selected,
							 const int &Forget_value, 
							 const bool Is_compression_selected,
							 const bool Is_adaptive_quantization_selected, 
							 const bool Is_bijection_selected); 


		//Initialization

		/**
		 \fn	void Global_Initialization(Polyhedron &pMesh, const int & Quantization_bit, const char * File_name);
		
		 \brief	Global initialization to select the input gate.
		
		\param [in,out]	pMesh		The mesh.
		 \param	Quantization_bit	Number of bits used for geometry quantization.
		 \param	File_name		 	Filename of the file.
		 */

		void Global_Initialization(Polyhedron & pMesh, 
								   const int  & Quantization_bit, 
								   const char * File_name);

		/**
			\fn	void Quantization(Polyhedron &pMesh);
			
			\brief	Quantize all vertices so that the new positions are reguliraly spaced in the 3D space.
			
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
			
			\brief	Initialization to deal with a mesh composed of multiple separated components.
			
			\param [in,out]	pMesh	The mesh.
			\param	Quantization_bit		The number of quantization for geometry.
			*/

		void Multiple_Components_Initialization(Polyhedron & pMesh,
												const int & Quantization_bit);
		
		
		/**
		 \fn	void Simplification(Polyhedron &pMesh, 
									const int & NVertices, 
									const bool Is_normal_flipping_selected,
		 							const bool Is_use_metric_selected,
									const float &Metric_thread, 
									const bool Is_use_forget_metric_selected,
		 							const int &Forget_value);
		
		 \brief	Mesh Simplification which applies iteratively decimation and regulation in pair.
		
		 \param [in,out]	pMesh				The mesh.
		 \param	NVertices		 				The vertices.
		 \param	Is_normal_flipping_selected  	The normal flipping.
		 \param	Is_use_metric_selected		 	The use metric.
		 \param	Metric_thread	 				The metric thread.
		 \param	Is_use_forget_metric_selected	The use forget metric.
		 \param	Forget_value	 				The forget value.
		 */

		void Simplification(Polyhedron &pMesh,
							const int & NVertices,
							const bool Is_normal_flipping_selected,
							const bool Is_use_metric_selected,
							const float &Metric_thread, 
							const bool Is_use_forget_metric_selected,
							const int &Forget_value);		

		/**
			\fn	int Decimation_Conquest(Polyhedron &pMesh,const bool Normal_flipping,
			 	const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,
			 	const int &Forget_value, const int & Component_ID);
			
			\brief	Removal of a set of independent vertices.
					
			\param [in,out]	pMesh	The mesh.
			\param	Is_normal_flipping_selected  	The normal flipping.
			\param	Is_use_metric_selected		 	The use metric.
			\param	Metric_thread	 	The metric thread.
			\param	Is_use_forget_metric_selected	The use forget metric.
			\param	Forget_value	 	The forget value.
			\param	Component_ID	 	Identifier for the component.
			
			\return	Number of decimated vertices.
			*/

		int Decimation_Conquest(Polyhedron &pMesh,
								const bool Is_normal_flipping_selected,
								const bool Is_use_metric_selected,
								const float &Metric_thread, 
								const bool Is_use_forget_metric_selected,
								const int &Forget_value, 
								const int & Component_ID);													

		/**
			\fn	int Regulation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,
			 	const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,
			 	const int & Component_ID);
			
			\brief	Removal of a set of independent vertices.
						
			\param [in,out]	pMesh	The mesh.
			\param	Is_normal_flipping_selected  	The normal flipping.
			\param	Is_use_metric_selected		 	The use metric.
			\param	Metric_thread	 	The metric thread.
			\param	Is_use_forget_metric_selected	The use forget metric.
			\param	Forget_value	 	The forget value.
			\param	Component_ID	 	Identifier for the component.
			
			\return	Number of decimated vertices.
			*/
		int Regulation(Polyhedron &pMesh,
						const bool Is_normal_flipping_selected,
						const bool Is_use_metric_selected,
						const float &Metric_thread, 
						const bool Is_use_forget_metric_selected,
						const int &Forget_value, 
						const int & Component_ID);			
		

		// Adaptive Quantization		

		/**
		 \fn	void Adaptive_Quantization(Polyhedron &pMesh, const int & NVertices,
		 		const bool Normal_flipping,const bool Use_metric,const float &Metric_thread,
		 		const bool Use_forget_metric,const int &Forget_value,const int &Qbit);
		
		 \brief	Adaptive quantization which not only decimates a input mesh but also adapts quantization precision for all intermediate meshes.
		
		 \param [in,out]	pMesh				The mesh.
		 \param	NVertices		 				The vertices.
		 \param	Is_normal_flipping_selected  	The normal flipping.
		 \param	Is_use_metric_selected		 	The use metric.
		 \param	Metric_thread	 				The metric thread.
		 \param	Is_use_forget_metric_selected	The use forget metric.
		 \param	Forget_value	 				The forget value.
		 \param	Qbit			 				The qbit.
		 */

		void Adaptive_Quantization(Polyhedron &pMesh, 
								   const int & NVertices, 
								   const bool Is_normal_flipping_selected,
								   const bool Is_use_metric_selected,
								   const float &Metric_thread, 
								   const bool Is_use_forget_metric_selected,
								   const int &Forget_value,
								   const int &Qbit);


		/**
		 \fn	void Augment_Geometry_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
		 \brief	Refines quantization precision of mesh geometry.
		
		 \param [in,out]						pMesh  	The mesh.
		 \param [in,out]						Decoder	The decoder.
		 \param	Component_ID				   	Component ID.
		 */

		void Augment_Geometry_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);

			/**
			 \fn	void Diminush_Geometry_Quantization_Precision(Polyhedron &pMesh, const int & Component_ID);
			
			 \brief	Decreasing of quantization resolution based on the prediction of PENG.
					Opposite function is Augment_Geometry_Quantization_Precision.
					
			 \param [in,out]						pMesh	The mesh.
			 \param	Component_ID				 	Component ID.
			 */


		void Diminush_Geometry_Quantization_Precision(Polyhedron &pMesh, const int & Component_ID);

		/**
		 \fn	void Calculate_Edge_Color_Difference(Polyhedron & pMesh, const int & Component_ID,
		 		double & Max_color, double & Mean_color, int & Number_of_vertices);
		
		 \brief	Calculates the edge color difference.
		
		 \param [in,out]	pMesh			  	The mesh.
		 \param	Component_ID				  	Identifier for the component.
		 \param [in,out]	Max_color		  	The maximum of color difference.
		 \param [in,out]	Mean_color		  	The mean of color difference.
		 \param [in,out]	Number_of_vertices	Number of vertices.
		 */

		void Calculate_Edge_Color_Difference(Polyhedron & pMesh, const int & Component_ID, double & Max_color, double & Mean_color, int & Number_of_vertices);

		/**
		 \fn	void Diminush_Color_Quantization_Precision(Polyhedron &pMesh, const int Component_ID);
		
		 \brief Reduces quantization precision of color coordinates.
				
		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Identifier for the component.
		 */

		void Diminush_Color_Quantization_Precision(Polyhedron &pMesh, const int Component_ID);

		/**
		 \fn	void Augment_Color_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder,
		 		const int & Component_ID);
		
		 \brief	Refine quantization precision of mesh color coordinates.
		
		 \param [in,out]	pMesh  	The mesh.
		 \param [in,out]	Decoder	The decoder.
		 \param	Component_ID	   	Identifier for the component.
		 */

		void Augment_Color_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
			// Compression 

		/**
		 \fn	void Compression(Polyhedron &pMesh,
								 const char* File_Name, 
								 const int &Qbit,
		 						 unsigned &Connectivity_size, 
								 unsigned & Color_size, 
								 unsigned & Total_size,
		 						 const unsigned & Initial_file_size);
		
		 \brief	Compressions.
			
		 \param [in,out]	pMesh			 	The mesh.
		 \param	File_Name					 	Filename of the output file.
		 \param	Qbit						 	The quantization of geometry.
		 \param [in,out]	Connectivity_size	Size of the connectivity.
		 \param [in,out]	Color_size		 	Size of the color.
		 \param [in,out]	Total_size		 	Size of the total.
		 \param	Initial_file_size			 	Size of the initial file.
		 */

		void Compression(Polyhedron &pMesh,const char* File_Name, const int &Qbit, unsigned &Connectivity_size, unsigned & Color_size, unsigned & Total_size);

		/**
		 \fn	int Calculate_Connectivity_Rate(void);
		
		 \brief	Calculates the connectivity rate.
		

		 \return	The calculated connectivity rate.
		 */

		int  Calculate_Connectivity_Rate(void);

		/**
		 \fn	void Calculate_Geometry_Color_Offset_Range();
		
		 \brief	Calculates the geometry color offset range.
		        This function is needed since the coder/decoder deals only with the positive symbol numbers.
		
		 */
		void Calculate_Geometry_Color_Offset_Range();				

		/**
		 \fn	void Remove_Last_Phase_Elements(const int & Component_ID);
		
		 \brief	When the last simplification operation has to be cancelled,
		        we need to remove its elements in order not to compress these elements. 
			
		 \param	Component_ID	Identifier for the component.
		 */

		void Remove_Last_Phase_Elements(const int & Component_ID);

		/**
		 \fn	void Write_Base_Mesh(Polyhedron &pMesh, Arithmetic_Codec & Enc,
		 		unsigned &Connectivity_size, unsigned &Color_size, const int &Num_color_base_mesh);
		
		 \brief	Writes a base mesh.
			
		 \param [in,out]	pMesh			 	The mesh.
		 \param [in,out]	Enc				 	The encoder.
		 \param [in,out]	Connectivity_size	Size of the connectivity.
		 \param [in,out]	Color_size		 	Size of the color.
		 \param	Num_color_base_mesh			 	Number of color in the base mesh.
		 */

		void Write_Base_Mesh(Polyhedron &pMesh, Arithmetic_Codec & Enc, unsigned &Connectivity_size, unsigned &Color_size, const int &Num_color_base_mesh);
		

		// Decompression		

		/**
		 \fn	int Decompress_Init(Polyhedron &pMesh,unsigned & Initial_file_size,
		 		const char* File_Name);
		
		 \brief	Initialize the Decompression by loading the base mesh into pMesh.
				
		 \param [in,out]	pMesh			 	The mesh.
		 \param [in,out]	Initial_file_size	Size of the initial file.
		 \param	File_Name					 	Filename of the compressed file.
		
		 \return	Information related to the decompression (number of LoDs, size of the input file, etc).
		 */

		QString    Decompress_Init(Polyhedron &pMesh);//,unsigned & Initial_file_size, const char* File_Name);		
		
		/**
		 \fn	void Decompression_From_File(Polyhedron &pMesh);
		
		 \brief	Decompression from file (No creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */

		void Decompression_From_File(Polyhedron &pMesh);

		/**
		 \fn	void Write_Info(Polyhedron &pMesh);
		
		 \brief	Write information of each LoD at decompression in a file.
		 \param [in,out]	pMesh	The mesh.		 
		 		 
		 */
		void Write_Info(Polyhedron &pMesh);

		/**
		 \fn	void Decompression_From_Sequence(Polyhedron &pMesh);
		
		 \brief	Decompression from sequence (Creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void Decompression_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh);


		/**
		 \fn	void Decompression_Coarser_From_File(Polyhedron &pMesh);
		
		 \brief	Obtain the previous LoD in the "from file" mode (no mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void Decompression_Coarser_From_File(Polyhedron &pMesh);


		/**
		 \fn	void Decompression_All_From_File(Polyhedron &pMesh);
		
		 \brief	Decompression of all LoDs from file.
		 The finest LoD is visualized without creating mesh sequence.
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void Decompression_All_From_File(Polyhedron &pMesh);


		/**
		 \fn	void Decompression_Specific_Level_From_File(Polyhedron &pMesh);
		
		 \brief	To obtain a user-desired LoD from file (No creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void Decompression_Specific_Level_From_File(Polyhedron &pMesh, const int & WL);

		/**
		 \fn	void JCW_Decompression_From_File(Polyhedron &pMesh);
		
		 \brief	Decompression from file for JCW (No creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void JCW_Decompression_From_File(Polyhedron &pMesh);


		/**
		 \fn	void JCW_Decompression_Without_Extraction_From_File(Polyhedron &pMesh);
		
		 \brief	Decompression from file without watermark extraction and geometry correction (No creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void JCW_Decompression_Without_Extraction_From_File(Polyhedron &pMesh);


		/**
		 \fn	void JCW_Decompression_From_Sequence(Polyhedron &pMesh);
		
		 \brief	Decompression from sequence for JCW (Creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */
		void JCW_Decompression_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh);


		/**
		 \fn	void JCW_Decompression_Without_Extraction_From_Sequence(Polyhedron &pMesh);
		
		 \brief	Decompression from sequence without watermark extraction (Creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 */

		void JCW_Decompression_Without_Extraction_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh);

		/**
		 \fn	void Show_Text(void);
		
		 \brief	Obtain information to the main window.

		 \return Information of each LoD at decompression.
		 		 
		 */
		QString Show_Text(void);


		/**
		 \fn	int Decompress_Each_Step(Polyhedron &pMesh, const char* File_Name);
		
		 \brief	Decompress the each step to visualize intermediate meshes
		
		 \param [in,out]	pMesh	The mesh.
		 \param	File_Name	Filename of the file.
		
		 \return Index of the current LoD.
		 */

		int    Decompress_Each_Step(Polyhedron &pMesh, const char* File_Name);

		/**
		 \fn	int Decompress_To_Level(Polyhedron &pMesh,const int Wanted_level);
		
		 \brief	Decompress to level.
		
		
		 \param [in,out]	pMesh	The mesh.
		 \param	Wanted_level	 	The wanted level.
		
		 \return Index of the current LoD.
		 */

		//int    Decompress_To_Level(Polyhedron &pMesh,const int Wanted_level);

		/**
		 \fn	void Un_Decimation_Conquest(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
		 \brief	Decoding of the decimation conquest.
			
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decoder.
		 \param	Component_ID			Component ID.
		 */

		void   Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & Decoder, const int & Component_ID);		

		/**
		 \fn	void Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);
		
		 \brief	Decoding of the regulation conquest.
				
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decoder.
		 \param	Component_ID				 	Component ID.
		 */

		void   Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);		

		// Other functions
		//void   Separate_Components(Polyhedron &pMesh);		

		/**
		 \fn	double Calculate_Area(Polyhedron & pMesh);
		
		 \brief	Calculates the area of pMesh.
				
		 \param [in,out]	pMesh	The mesh.
		
		 \return	The calculated area.
		 */

		double Calculate_Area(Polyhedron & pMesh);

		/**
		 \fn	bool Error_Projected_Surface(Polyhedron & pMesh, const Halfedge_handle & _h, const int & _Component_ID, const double & Mean_color, const double & Mean_area);
		
		 \brief	Calculate an error cause by a vertex removal and decide in order not to remove an important vertex.
				
		 \param [in,out]	pMesh	The mesh.
		 \param h				Input gate.
		 \param Component_ID    Component ID.
		 \param Mean_color      Mean color value.
		 \param Mean_area       Average of facets areas.
		 \return	Decision of vertex removal.
		 */

		bool   Error_Projected_Surface(Polyhedron & pMesh, const Halfedge_handle & h, const int & Component_ID, const double & Mean_color, const double & Mean_area);

		/**
		 \fn	void Recalculate_Component_Area(Polyhedron & pMesh, const int & Component_ID, int & Number_facets);
		
		 \brief	Update the area of each component.
				
		 \param [in,out]	pMesh	The mesh.
		 \param Component_ID  Component ID.
		 \param Number_facets Number of facets
		
		 \return	Decision of vertex removal.
		 */
		void   Recalculate_Component_Area(Polyhedron & pMesh, const int & Component_ID, int & Number_facets);		

		/**
		 \fn	Point_Int Change_Real_Int(const Point3d &pt, const int & Component_ID);
		
		 \brief	Change from real to int point(related to quantization bit)		
		 
		 \param	pt				The point.
		 \param	Component_ID	Component ID.
		
		 \return Integer coordinates.
		 */

		inline Point_Int Change_Real_Int(const Point3d &pt, const int & Component_ID);		

		/**
		 \fn	Point3d Change_Int_Real(const Point_Int &vec, const int & Component_ID);
		
		 \brief	Change from int to real point(related to quantization bit)
		
		 \param	pt				The point(integers).
		 \param	Component_ID	Component ID.
		
		 \return	.
		 */

		inline Point3d   Change_Int_Real(const Point_Int & pt, const int & Component_ID);
		
		
		/**
		 \fn	void Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh);
		
		 \brief	To add flags to copied mesh

		
		 \param [in,out]	Original_mesh	The original mesh.
		 \param [in,out]	New_mesh	 	The new mesh.
		 */

		void Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh);			
		
		/**
		 \fn	int Coding_Cost_Each_Layer_Color_Clustering(Polyhedron &pMesh);
		
		 \brief	Color coding cost - prediction method

		 \param [in,out]	pMesh	The mesh.
		
		 \return	.
		 */
		//int Coding_Cost_Each_Layer_Color_Clustering(Polyhedron &pMesh);	
		
	 
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
		 \fn	double QString Joint_Compression_Watermarking(Polyhedron &pMesh, 
											  const char * Input_File_Name, 
											  const char * Output_File_Name, 
											  const int & Number_bins, 
											  const int & Number_regions,
											  const int & Embedding_strength, 
											  const char * Embedding_message, 
											  const bool Is_complete_reversibility_selected, 
											  const bool Is_divide_regions_selected,
											  const int & Thres_divide_regions, 
											  const int &Qbit, 
											  const int  & NVertices, 
											  const bool Normal_flipping, 
											  const bool Use_metric, 
											  const float & Metric_thread, 
											  const bool Use_forget_metric, 
											  const int &Forget_value);
		
		 \brief	Joint Compression Watermarking (JCW) 
		
		 \param [in,out]	pMesh			 		The mesh.
		 \param Input_File_Name						The name of the input file.
		 \param Output_File_Name					The name of the output file.
		 \param Number_bins							The number of bins.
		 \param Number_regions						The number of regions.
		 \param Embedding_strength				    The strengh of embedding (number of shifted bins).
		 \param Embedding_message                   The message of watermarking
		 \param Is_complete_reversibility_selected  The selection of complete reversibility. 
		 \param Is_divide_regions_selected          The selection of division of regions.
		 \param Thres_divide_regions                The threshold to divide big regions.
		 \param Qbit                                The geometry quantization.
		 \param NVertices                           The wanted number of base mesh.
		 \param Normal_flipping                     The selection of normal flipping.
		 \param Use_metric                          The selection of geometric metric use
		 \param Metric_thread                       The threshold of metric.
		 \param Use_forget_metric                   The selection of use of "forget metric"
		 \param Forget_value                        The threshold of "Use_forget_metric"
		
		 \return	Information of JCW.
		 */

		QString Joint_Compression_Watermarking(Polyhedron &pMesh, 
											  const char * Input_File_Name, 
											  const char * Output_File_Name, 
											  const int & Number_bins, 
											  const int & Number_regions,
											  const int & Embedding_strength, 
											  const char * Embedding_message, 
											  const bool Is_complete_reversibility_selected, 
											  const bool Is_divide_regions_selected,
											  const int & Thres_divide_regions, 
											  const int &Qbit, 
											  const int  & NVertices, 
											  const bool Normal_flipping, 
											  const bool Use_metric, 
											  const float & Metric_thread, 
											  const bool Use_forget_metric, 
											  const int &Forget_value);


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
         \param	Component_ID					   	Component ID.
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
		 \param	Component_ID					   	Component ID
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
			
		 \param	h		 	The input gate.
		 \param	Direction	The direction.
		
		 \return Predicted position.
		 */

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & Direction);

		/**
		 \fn	Point3d JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence,
		 		const int & Direction);
		
		 \brief	Caclulates the barycentric position for JCW after removal.
		
		
		 \param	h		 	The input gate.
		 \param	valence  	The valence.
		 \param	Direction	The direction.
		
		 \return Predicted position.
		 */

		Point3d JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);

		/**
		 \fn	Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h,
		 		const int & valence, const int & Direction);
		
		 \brief	JCW barycenter patch before removal.
		 
		
		 \param	h		 	The input gate.
		 \param	valence  	The valence.
		 \param	Direction	The direction.
		
		 \return	.
		 */

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);

		/**
		 \fn	void JCW_Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & Decoder, const int & Component_ID);
		

		 \brief	Undecimation conquest for JCW.
		 	
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decrement.
		 \param	Component_ID				 	Component ID.
		 */

		void JCW_Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & Decoder, const int & Component_ID);		

		/**
		 \fn	void JCW_Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);
		
		 \brief	JCW unregulation.
		
	 	
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decrement.
		 \param	Component_ID				 	Component ID.
		 */

		void JCW_Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);		

		/**
		 \fn	void JCW_Un_Regulation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID,
		 		list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);
		
		 \brief	JCW unregulation for region detection, in order to obtain region number of each inserted vertices
		
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
		
		 \brief	JCW undecimation for region detection.
		
	 	
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
		
		
		 \return Results of robustness evaluation.
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
			this->m_EmbeddingStrength = EL;
		}

		/**
		 \fn	int Get_Embedding_Level(void)
		
		 \brief	Gets the embedding level.
		
	 	
		 \return	The embedding level.
		 */

		int Get_Embedding_Level(void)
		{
			return this->m_EmbeddingStrength;
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

		//void JCW_Choose_Valid_Vertices(Polyhedron & pMesh);

		/**
		 \fn	vector<int> JCW_Region_Mass_Center_Extract_Watermark(Polyhedron & pMesh);
		
		 \brief	JCW region mass center extract watermark.
		
		 \param [in,out]	pMesh	The mesh.		
		 
		 */

		void JCW_Region_Mass_Center_Extract_Watermark(Polyhedron & pMesh);

		/**
		 \fn	void JCW_Code_Difference_Histogram_Shifting(Polyhedron &pMesh,const int & Component_ID);
		
		 \brief	JCW code difference histogram shifting.
				
		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Component ID.
		 */

		void JCW_Code_Difference_Histogram_Shifting(Polyhedron &pMesh,const int & Component_ID);

		/**
		 \fn	void JCW_Decode_Difference_Histogram_Shifting(Polyhedron &pMesh,
		 		const int & Component_ID);
		
		 \brief	JCW decode difference histogram shifting.
		
		
		 \param [in,out]	pMesh	The mesh.
		 \param	Component_ID	 	Component ID.
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
		
		 \brief	Converts a cartesian coordinates into spherical coordinates.
			
		 \param	Pt				   	The point.
		 \param [in,out]	Spheric	If non-null, the spheric.
		 */

		void Convert_To_Spherical(const Point3d & Pt, double * Spheric);

		/**
		 \fn	void Convert_To_Cartesian(const double * Spheric, double * Cartesian);
		
		 \brief	Converts spherical coordinates into cartesian coordinates.
		
		
		 \param	Spheric				 	The spheric.
		 \param [in,out]	Cartesian	If non-null, the cartesian.
		 */

		void Convert_To_Cartesian(const double * Spheric, double * Cartesian);	

		/**
		 \fn	int JCW_Decompress_One_Level(Polyhedron &pMesh, const char* File_Name,
		 		const int & Noise_mode);
		
		 \brief	JCW decompress one level.
		
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
		 \fn	int JCW_Divide_Big_Regions(Polyhedron &pMesh);
		
		 \brief	Jcw divide big regions.
		
		
		 \param [in,out]	pMesh	The mesh.
		
		 \return Number of divided regions.
		 */

		int  JCW_Divide_Big_Regions(Polyhedron &pMesh, const int & _Thres_divide_regions);

		/**
		 \fn	void JCW_Colorify_Regions(Polyhedron & pMesh);
		
		 \brief	Colorify regions for JCW.
		
		 \param [in,out]	pMesh	The mesh.
		 */

		void JCW_Colorify_Regions(Polyhedron & pMesh);

		/**
		 \fn	void Read_Information_To_Hide();
		
		 \brief	Reads the watermarked information.
		*/

		void    Read_Information_To_Hide(const char * Embedding_message);

		/**
		 \fn	QString Write_Information_To_Hide();
		
		 \brief	Writes the information to hide.
				
		 \return Inserted to watermarking message.
		 */

		QString Write_Information_To_Hide();

		/**
		 \fn	void Clear_After_Compression();
		
		 \brief	Clears after the compression.

		 */

		void Clear_After_Compression();
	private:		
		
		vector<bool> IsClosed;    ///< The is closed. To know if the mesh is open or closed.		
		
		bool IsColored;///< true if is colored		
		bool IsOneColor;///< true if is one color

		float OnlyColor[3];///< The coordinates of color when there is only one color
		// Number of connectivity symbol types. Without boundary = 5, with boundary = 7;
		//int NummberConnectivitySymbols;		
		
		// To encode each type of operation between decimation and increase of quantization resolution
        vector<list<int> >		 ListOperation;///< The list of operation
        vector<list<int> >		 Connectivity;///< The information of connectivity to compress
        vector<list<Point_Int> > Geometry;///< The geometry information to compress
        vector<list<int> >		 NumberSymbol;///< Number of symbols of each simplification step
        vector<list<int> >		 NumberVertices;///< Number of vertices of each simplification step
        			
		list<Point_Int>          InterGeometry;///< The intermediate information of geometry
		list<int>		         InterConnectivity;		///< The intermediate information of connectivity		

        vector<list<int> > AlphaRange;///< The range of alpha (Frenet coordinates) of each LoD
        vector<list<int> > AlphaOffset;///< The offset of alpha
        vector<list<int> > GammaRange;///< The range of gamma (Frenet coordinates) of each LoD
        vector<list<int> > GammaOffset;///< The offset of gamma
        
		
		// Quantization		
		vector<unsigned>	Qbit; ///< The Quantization bits
		vector<float>		xmin; ///< The xmin
		vector<float>		ymin; ///< The ymin
		vector<float>		zmin; ///< The zmin
		vector<float>		xmax; ///< The xmax
		vector<float>		ymax; ///< The ymax
		vector<float>		zmax; ///< The zmax
		vector<float>		Quantization_Step; ///< The quantization step				
		
		int			  Smallest_Alpha;		    ///< The smallest alpha
		int			  Smallest_Gamma;			///< The smallest gamma
					
		vector<double> HighestLengthBB; ///< The highest length bb
		vector<double> ComponentVolume; ///< The volume of each component
		vector<double> ComponentArea; ///< The area of each component
		vector<int>    ComponentNumberVertices;///< The number of vertices of each components
				
		// Used for adatative quantization.				
        vector<list<int> > QuantizationCorrectVector; ///< The quantization correct vector
        vector<list<int> > NumberQuantizationLayer; ///< Number of quantization layers
       		
		//for color
        vector<list<int> > NumberProcessedVertices;///< Number of processed vertices
        vector<list<int> > ColorChildcellIndex;///< the color childcell index
        vector<list<int> > ColorEncoderIndex;///< the color encoder index        

		// Colors
		vector<list<Color_Unit> > VertexColor; ///< Contain color error of all removed vertices
		list<Color_Unit> InterVertexColor;		///< The intermediate information vertex color
		
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
		
		/*int ColorDiffMinC0;///< The color difference minimum c0 for restoration of colors in Mapping table
		int ColorDiffMinC1;///< The first color difference minimum c1 for restoration of colors in Mapping table
		int ColorDiffMinC2;///< The second color difference minimum c2 for restoration of colors in Mapping table
				
		int ColorDiffRangeC0;///< The color difference range c 0 for restoration of colors in Mapping table
		int ColorDiffRangeC1;///< The first color difference range c1 for restoration of colors in Mapping table
		int ColorDiffRangeC2;///< The second color difference range c2 for restoration of colors in Mapping table*/
				
		//// mapping table
		//vector<Color_Unit> ColorArray; ///< Color table
		//vector<Color_Unit> PredictedColorArray; ///< Predicted values of colors in color table using prediction based on neighbors.		
		//vector<Color_Unit> DifferenceColor; ///< Difference between original and quantized vertex color. need to turn into lossless coding.
		//
		//list<int>		   ColorIndex;///< List of color index
		//vector<int>		   ReorderingColorIndex;///< Vector of reordering color index
		//list<int>		   InterColorIndex;		///< Vector of inter color index
		//vector<int>		   Number_color_index; ///< contain occurence of each initial color
		//vector<int>		   IsKnownIndex;///< is known index
		

		
		// Decoder
		
		Arithmetic_Codec Decoder;		///< The arithmetic decoder
		
		Adaptive_Data_Model Color_0_Model;///< The statistical model used for color_0
		Adaptive_Data_Model Color_1_Model;///< The statistical model used for color_1
		Adaptive_Data_Model Color_2_Model;///< The statistical model used for color_2

		Adaptive_Data_Model Index_Model;///< The index model

		vector<int> NumberDecimation; ///< To stock number of Decimation operation.
		vector<int> NumberChangeQuantization; ///< to stock number of diminution of quantization.
				
		int DumpSymbolDecimation;///< The dump symbol decimation in order to eliminate symbols of the last conquest if it is useless.
		int DumpSymbolRegulation;///< The dump symbol regulation in order to eliminate symbols of the last conquest if it is useless.
		
		int Decompress_count;///< Number of decompress in order to know how many steps to go for the decompression step by step.
		int NumberComponents;///< Number of components in order to know how many steps to go for the decompression step by step.
		
		int GlobalCountOperation;///< The global count operation
	
		vector<int> ComponentOperations;///< The operations of each component
		
		// JCW
		int   m_NumberBin;///< Number of bins for JCW
		int   m_EmbeddingStrength; ///< Embedding strength for JCW
		int   m_NumberRegion; ///< Number of regions for JCW
		double m_VC[3]; ///< Mesh center position for JCW
		double m_Rmin; ///< Distance of farthest vertex from mesh center for JCW
		double m_Rmax; ///< Distance of nearst vertex from mesh center for JCW
		double m_Dist; ///< Distance of each bin
		vector<int> m_Number_Vertices_Per_Regions; ///< Number of vertices in each region
		list<int> m_Watermarks; ///< watermark

		int N_Inserted_Watermarks; ///< Number of inserted watermarks

		vector<int> m_N_remained_vertices; ///< Number of remained vertices in each region
		vector<int> m_N_treated_vertices; ///< Number of treated vertices in each region
		vector<double> m_Rad_decision; 

        list<vector<int> > m_JCW_Move_Error; ///< Stock difference related to complete reversibility
		list<int> m_N_Errors; ///< Index of error related to complete reversibility.

		Adaptive_Data_Model DM_JCW_MOVE_ERROR;
		
		list<int> JCW_Connectivity; ///< Stock connectivity information for JCW
		list<Point_Int> JCW_Geometry; ///< Stock geometry information for JCW

		//double LUT_CourbureClust[3*256];
		vector< vector<float> > Region_Color; ///< Color of each region
		int Number_non_reversible_vertices; ///< Number of vertices which violate complte reversibility
		int Number_Save_Over_bins; ///< Number of empty bins to shift the current bins
		bool Is_Division_Big_Regions_Enabled; ///< true if "division_big_regions" option is seleted
		bool Is_Complete_Reversibility_Enabled; ///< true if "complete reversibility" option is seleted
		bool Is_Bijection_Enabled; ///< true if "bijection" option is seleted
		int Division_Threshold; ///< Threshold of region division
		
	// from IHM
	public:
		bool IsDecompress;///< true if is decompress
		bool IsCompressed;///< true if is compressed
		int Current_level;///< The current level
		int Total_layer;///< The total number of layer
		unsigned Initial_file_size;///< Size of the initial file
		unsigned Compressed_file_size;///< Size of the compressed file
		string File_name;///< Filename of the file
		
		vector<float> Prog; ///< Stock information of progression of decompression in pourcentage (100 = total decompression)
		vector<float> Ratio;  ///< Stock information of ration regarding size of the input file
		QString Message; ///< Message to be shown in the main window

	public://private:
		CCopyPoly<Polyhedron, Enriched_kernel> Copy_Polyhedron;
		
		//bool Afficher; ///<
		bool Sequence; ///< Decompression mode (sequence mode or file mode) for IHM
		bool Possible_change_sequence; ///< To disable the mode change during decompression for IHM
		
		int Visu_level; ///< Level of visualized LoD for IHM
		int Process_level; ///< Level of processed LoD for IHM

		FILE *Dec_Info; ///< File to write information of decompression for IHM
		//int Writing_level; ///< Level of 
		string Dec_File_Info; ///< File name to write decompression information for IHM
};

#endif

/*! \mainpage component_Compression_Valence Documentation
 *
 * \section auth Authors
 * H. Lee, C. Dikici, G. Lavoué and F. Dupont \n
 * M2Disco Team, LIRIS, CNRS, Université Lyon 1/INSA-Lyon, France. \n
 * \n
 * Please contact Ho Lee (hosmail@gmail.com) in case you have any question, comment or a bug to report.
 * \n
 * \n
 *
 * \section paper_reference Related publications
 * Here is the list of related papers if you want to cite our methods. \n
 * It is also strongly recommended that you read these papers
 * in order to know the parameters used in our methods. 
 * \n 
 * 
 * \subsection compression_paper Progressive compression of colored meshes
 * H. Lee, G. Lavoué, F. Dupont, \n
 * Rate-distortion optimization for progressive compression of 3D mesh with color attributes, \n
 * The Visual Computer, 2011, (accepted for publication).
 * \n
 * 
 * \subsection jcw_paper Joint watermarking and progressive compression of 3D meshes
 * H. Lee, C. Dikici, G. Lavoué, F. Dupont, \n
 * Joint reversible watermarking and progressive compression of 3D meshes, \n
 * The Visual Computer, 2011, (accepted for publication).
 * \n 
 * \n 
 * 
 * \section overview Overview
 * This Compression Valence component has two main functions:
 * (1) a progressive compresion and (2) a joint progressive compression and reversible watermarking for 3D meshes. \n
 * The progressive compression method simplifies iteratively an input mesh to generate differents levels of details (LoDs). \n
 * These LoDs are then transmitted progressively in a coarse-to-fine way at the decompression. \n
 * In particular, our method adapts the quantization precision (both geometry and color) to each LoD in order to optimize the rate-distortion (R-D) trade-off. \n \n
 * Our joint progressive compression and reversible watermarking method offers a possibility to embed an watermark information
 * in order to protect the ownership of the input mesh. \n
 * Hence, at each decompression step, the inserted message can be extracted also progressively. \n
 * The embedded watermarks are reversible, meaning that the deformation caused by watermarking embedding step can be removed when extracting the watermark.
 * \n
 * \n
 * 
 * \section howto_use How to use this component
 * \subsection howto_progressive_compression Progressive compression
 * First of all you have to load a mesh. \n
 * To compress the input mesh, choose "Compression" and a windows appears. \n
 * (1) File name              : to choose a name for the compressed file. Please be aware that the file extension should be ".p3d". \n
 * (2) Mode                   : the mode "Simplification" does not generate a compressed file. \n
 * (3) Compression mode       : to enable or disable the use of the adaptation of quantization precision. \n
 * (4) Quantization bits      : the number of bits for the geometry quantization. \n
 * (5) # Vertices wanted      : the number of vertices of the base (the coarsest) mesh. \n
 * (6) Use bijection          : to enable or disable the use of the bijection for the geometry encoding. This bijection reduces the coding rates but it needs a longer processing time. \n
 * (7) Forbid normal flipping : this option is used to forbid a normal flipping when removing a vertex. \n
 * (8) Use Metric             : this option is used to forbid a vertex removal if it induces a significant deformation. The threshold value is initially set to 0.25. \n
                                The use of this metric can be "forgetted" if the number of vertices of the current intermediate mesh is superior to an user-defined threshold. \n
 * \n
 * For the decompression, first you have to load a .p3d file. Then, the base mesh is rendered. \n
 * To obtain higher LoDs, you can use : \n
 * (1) "Decompression : all LoDs", to decompress all LoDs and the finest intermediate is visualized, \n
 * (2) "Decompression : next LoD", to decompress one level and the next LoD is rendered (ALT + left mouse button), \n
 * (3) "Decompression : go to specific LoD", to reach the desired LoD. \n 
 * You can also visualize the previous LoDs by selecting "Decompression : previous LoD" (ALT + right mouse button). \n \n 
 * After loading the .p3d file, you can enable or disable the option of mesh sequence generation with the menu "Activate/Desactivate mesh sequence generation". \n
 * When this option is enabled, all LoDs are stored in the memory, so that the navigation of differents LoDs can be performed quickly. \n
 * You can disable this option in order to save the memory.
 * In this case, when you want to visualize the previous LoD, the decompression is performed again until getting the previous LoD.
 * The navigation takes a longer time. \n
 * Note that this option can be modified only after the rendering of the base mesh and it cannot be modified when any decompression operation is achieved. \n 
 * \n
 * \subsection howto_jcw Joint compression and watermarking
 * First of all you have to load a mesh. \n
 * You can apply our joint method by selecting "JCW - Compression and Watermarking embedding". \n
 * (1) File name              : to choose a name for the compressed file (.p3d). \n
 * (2) Q bits                 : the number of bits for the geometry quantization. \n
 * (3) # vertices             : the number of the base mesh after an iterative simplification. \n
 * (4) # Bins				  : the number of the histogram bins. \n
 * (5) # Regions              : the number of regions. One watermark bit is embedded in each region. \n
 * (6) Embedding Strength     : the number of shifted bins when embedding/extracting the watermark. \n
 * (7) Embedding Message      : the message to insert. When this field is empty or the length of the message is shorter than necessary, the message is generated randomly. \n
 * (8) Divide Regions         : when this option is selected, the big region is divided in two in order to insert more watermark bits. \n
 * (9) Complete Reversibility : When this option is checked, the initial positions of all vertices are exactly restored. Some extra coding bits are necessary.
 * \n
 * For the decompression and the watermark extraction,
 * you have to load a .p3d file. \n
 * To obtain higher LoDs, you can select "JCW - Decompression and Watermark extraction : next LoD". \n
 * The extracted message is shown in the status bar of the main window. \n
 * To visualize the results without watermark extraction (non authorized users), you can use "JCW - Decompression without extraction : next LoD". \n
 * \n
 * 
 * \section last_updated Last updated
 * 18 May, 2011
 * 
 */

#endif // Compression_Valence_COMPONENT_H
