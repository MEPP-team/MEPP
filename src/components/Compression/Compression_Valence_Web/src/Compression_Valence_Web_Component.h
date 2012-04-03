///////////////////////////////////////////////////////////////////////////
// Author: Ho LEE
// Year: 2011
// Month: MAY
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

#ifndef Compression_Valence_Web_COMPONENT_H
#define Compression_Valence_Web_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence_Web

#include "../../../../mepp/mepp_component.h"

#include "Compression_Valence_Web_Polyhedron.h"
#include "Compression_Valence_Basemesh_builder.h"
#include "bitfile.h"
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

 
	inline std::ostream& operator<<(std::ostream& os, const Point_Int & bb)
{
	
	os<<bb.x<<","<<bb.y<<","<<bb.z<<"\n";
	
	return os;
}
	inline std::ostream& operator<<(std::ostream& os, const vector<Point_Int>& bb)
	{
		for (unsigned int i=0;i<bb.size();i++)
		{
			os<<bb[i];
		}
		return os;
	}
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
	/// \param	m_color	The col.
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
	/// \param	m_color	The col.
	////////////////////////////////////////////////////////////////////////////////////////////////////

	bool operator !=(const Color_Unit & m_color) const
	{
		return!(*this == m_color);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \class	Compression_Valence_Web_Component
///
/// \brief	Compression valence web component. 
///
////////////////////////////////////////////////////////////////////////////////////////////////////

class Compression_Valence_Web_Component : 

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
		/// \fn	Compression_Valence_Web_Component(Viewer* v, PolyhedronPtr p)
		///
		/// \brief	Default Constructor.
		///
		/// \param [in,out]	v	If non-null, the v.
		/// \param	p		 	The.
		////////////////////////////////////////////////////////////////////////////////////////////////////

		Compression_Valence_Web_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
		{			
			IsOneColor = false;  
			Decompress_count = 0; 			
					
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
			

			// MEPP 2
			componentName = "Compression_Valence_Web_Component";
			init = 1;
		}

		/**
		 \fn	~Compression_Valence_Web_Component()		
		 \brief	Destructor.		
		 */

		~Compression_Valence_Web_Component() 
		{}
		

		// Main Function

		/**
		 \fn	double Main_Function(Polyhedron &pMesh, 
							 const char * Input_File_Name, 
							 const char* File_Name,
							 const int &Qbit,
							 const int &NVertices, 
							 const bool Is_compression_selected);
		
		 \brief	Main function of compression.		
		
		 \param [in,out]	pMesh				 	The mesh.
		 \param	Input_File_Name						Filename of the input file.
		 \param	File_Name						 	Filename of the compressed(output) file.
		 \param	Qbit							 	The qbit.
		 \param	NVertices						 	The vertices.
		 \param	Is_compression_selected			 	The compression selected.
		 \return	Information of compression results.
		 */

		QString Main_Function(Polyhedron &pMesh, 
							 const char * Input_File_Name, 
							 const char* File_Name,
							 const int &Qbit,
							 const int &NVertices, 
							const bool Is_compression_selected);
							 


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
			\fn	void Multiple_Components_Initialization(Polyhedron & pMesh, const int & Quantization_bit);
			
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
			\fn	int Decimation_Conquest(Polyhedron &pMesh,const bool Is_normal_flipping_selected,
			 	const bool Is_use_metric_selected,const float &Metric_thread, const bool Is_use_forget_metric_selected,
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
			\fn	int Regulation(Polyhedron &pMesh,const bool Is_normal_flipping_selected,const bool Is_use_metric_selected,
			 	const float &Metric_thread, const bool Is_use_forget_metric_selected,const int &Forget_value,
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
		


		/**
		 \fn	void Augment_Geometry_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
		 \brief	Refines quantization precision of mesh geometry.
		
		 \param [in,out]						pMesh  	The mesh.
		 \param [in,out]						Decoder	The decoder.
		 \param	Component_ID				   	Component ID.
		 */

//		void Augment_Geometry_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);

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

//		void Augment_Color_Quantization_Precision(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
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

		void Write_Base_Mesh(Polyhedron &pMesh, bit_file_c & Enc, unsigned &Connectivity_size, unsigned &Color_size, const int &Num_color_base_mesh);
		

		// Decompression		

		/**
		 \fn	int Decompress_Init(Polyhedron &pMesh,unsigned & Initial_file_size,
		 		const char* File_Name);
		
		 \brief	Initialize the Decompression by loading the base mesh into pMesh.
				
		 \param [in,out]	pMesh			 	The mesh.		
		
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
		 \fn	void Decompression_From_Sequence(Polyhedron &pMesh, Polyhedron &New_mesh);
		
		 \brief	Decompression from sequence (Creation of mesh sequence).
		
		 \param [in,out]	pMesh	The mesh.		 
		 \param [in,out]	New_mesh	The new copied mesh.
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
		 \param             WL      The wanted level.
		 */
		void Decompression_Specific_Level_From_File(Polyhedron &pMesh, const int & WL);

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
		 \fn	void Un_Decimation_Conquest(Polyhedron &pMesh, bit_file_c & Decoder, const int & Component_ID);
		
		 \brief	Decoding of the decimation conquest.
			
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decoder.
		 \param	Component_ID			Component ID.
		 */

		void   Un_Decimation_Conquest(Polyhedron &pMesh, const int & Component_ID);		

		/**
		 \fn	void Un_Regulation(Polyhedron &pMesh, bit_file_c & Decoder, const int & Component_ID);
		
		 \brief	Decoding of the regulation conquest.
				
		 \param [in,out]	pMesh	The mesh.
		 \param [in,out]	Decoder  	The decoder.
		 \param	Component_ID				 	Component ID.
		 */

		void   Un_Regulation(Polyhedron &pMesh, const int & Component_ID);		

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
		 \fn	bool Error_Projected_Surface(Polyhedron & pMesh, const Halfedge_handle & h, const int & Component_ID, const double & Mean_color, const double & Mean_area);
		
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
		 \fn	Point3d Change_Int_Real(const Point_Int &pt, const int & Component_ID);
		
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
		 \fn	unsigned Calculate_Current_File_Size(void)
		
		 \brief	Calculates the current file size. Measure bits used for decompression.		
	 	
		 \return	The calculated current file size.
		 */

		unsigned Calculate_Current_File_Size(void) 
		{
			// To measure exact quantity of bits used for decompression.

			return this->Decoder.calculate_current_decoded_size();
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

		void Stop_Decoder(void) {

		}
		void Clear_After_Compression();
	private:		
		
		vector<bool> IsClosed;    ///< The is closed. To know if the mesh is open or closed.		
		
		bool IsColored;///< true if is colored		
		bool IsOneColor;///< true if is one color

		float OnlyColor[3];///< The coordinates of color when there is only one color
		// Number of connectivity symbol types. Without boundary = 5, with boundary = 7;
		//int NummberConnectivitySymbols;		
		
		// To encode each type of operation between decimation and increase of quantization resolution
              
        vector<list<int> >		 Connectivity;///< The information of connectivity to compress
        vector<list<Point_Int> > Geometry;///< The geometry information to compress
        vector<list<int> >		 NumberSymbol;///< Number of symbols of each simplification step
        vector<list<int> >		 NumberVertices;///< Number of vertices of each simplification step
        			
		list<Point_Int>          InterGeometry;///< The intermediate information of geometry
		list<int>		         InterConnectivity;		///< The intermediate information of connectivity		      
		
		// Quantization		
		vector<unsigned>	Qbit; ///< The Quantization bits
		vector<float>		xmin; ///< The xmin
		vector<float>		ymin; ///< The ymin
		vector<float>		zmin; ///< The zmin
		vector<float>		xmax; ///< The xmax
		vector<float>		ymax; ///< The ymax
		vector<float>		zmax; ///< The zmax
		vector<float>		Quantization_Step; ///< The quantization step				
							
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
		
		// Decoder
		
		//Arithmetic_Codec Decoder;		///< The arithmetic decoder
		bit_file_c Decoder;
		unsigned int header_size;///< The header size of p3dw file, used in decoding


		vector<int> NumberDecimation; ///< To stock number of Decimation operation.
				
		int DumpSymbolDecimation;///< The dump symbol decimation in order to eliminate symbols of the last conquest if it is useless.
		int DumpSymbolRegulation;///< The dump symbol regulation in order to eliminate symbols of the last conquest if it is useless.
		
		int Decompress_count;///< Number of decompress in order to know how many steps to go for the decompression step by step.
		int NumberComponents;///< Number of components in order to know how many steps to go for the decompression step by step.
		
		int GlobalCountOperation;///< The global count operation
	
		vector<int> ComponentOperations;///< The operations of each component
		
	
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

/*! \mainpage component_Compression_Valence_Web Documentation
 *
 * \section auth Authors
 * H. LEE, C. DIKICI, G. LAVOUE and F. DUPONT \n
 * M2DisCo Team, LIRIS, CNRS, Universite Lyon 1/INSA-Lyon, France. \n
 * \n
 * Please contact the authors \n H. LEE (hosmail@gmail.com), C. DIKICI (cagataydikici@yahoo.com),
 * G. LAVOUE (glavoue@liris.cnrs.fr) and F. DUPONT (fdupont@liris.cnrs.fr) \n
 * in case you have any question, comment or a bug to report.
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
 * H. LEE, G. LAVOUE, F. DUPONT, \n
 * Rate-distortion optimization for progressive compression of 3D mesh with color attributes, \n
 * The Visual Computer, 2011.
 * \n
 * 
 * \subsection jcw_paper Joint watermarking and progressive compression of 3D meshes
 * H. LEE, C. DIKICI, G. LAVOUE, F. DUPONT, \n
 * Joint reversible watermarking and progressive compression of 3D meshes, \n
 * The Visual Computer 27(6-8):781-792, 2011.
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
 * (3) Quantization bits      : the number of bits for the geometry quantization. \n
 * (4) # Vertices wanted      : the number of vertices of the base (the coarsest) mesh. \n
 * \n
 * For the decompression, first you have to load a .p3dw file. Then, the base mesh is rendered. \n
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
 * \section last_updated Last updated
 * 18 May, 2011
 * 
 */

#endif // Compression_Valence_Web_COMPONENT_H
