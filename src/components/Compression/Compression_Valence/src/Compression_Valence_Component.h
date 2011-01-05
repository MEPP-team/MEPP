#ifndef Compression_Valence_COMPONENT_H
#define Compression_Valence_COMPONENT_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "../../../../mepp/mepp_component.h"

#include "Compression_Valence_Polyhedron.h"
#include "arithmetic_codec.h"
#include "Compression_Valence_Basemesh_builder.h"
#include "Color_Distance/ColorMesh.h"

#include <queue>
#include <list>

// struct of integer coordinates
struct Point_Int
{
	int x;
	int y;
	int z;
	Point_Int()
	{
		x=y=z=0;
	}
	const Point_Int operator+(const Point_Int &Pt) const
	{
		Point_Int Temp;
		Temp.x = x + Pt.x;
		Temp.y = y + Pt.y;
		Temp.z = z + Pt.z;

		return Temp;
	}
	const Point_Int operator-(const Point_Int &Pt) const
	{
		Point_Int Temp;
		Temp.x = x - Pt.x;
		Temp.y = y - Pt.y;
		Temp.z = z - Pt.z;

		return Temp;
	}
	bool operator ==(const Point_Int &Pt) const
	{
		return (x == Pt.x && y == Pt.y && z == Pt.z);
	}
	bool operator !=(const Point_Int &Pt) const
	{
		return !(*this == Pt);
	}
};

// struct of integer color components.
struct Color_Unit
{
	int c0;
	int c1;
	int c2;
	Color_Unit()
	{
		c0=c1=c2=0;
	}

	const Color_Unit operator+(const Color_Unit &Col) const
	{
		Color_Unit Temp;
		Temp.c0 = c0 + Col.c0;
		Temp.c1 = c1 + Col.c1;
		Temp.c2 = c2 + Col.c2;
		
		return Temp;
	}
	const Color_Unit operator-(const Color_Unit &Col) const
	{
		Color_Unit Temp;
		Temp.c0 = c0 - Col.c0;
		Temp.c1 = c1 - Col.c1;
		Temp.c2 = c2 - Col.c2;
		
		return Temp;
	}
	bool operator ==(const Color_Unit &Col) const
	{
		return( c0 == Col.c0 && c1 == Col.c1 && c2 == Col.c2);
	}
	bool operator !=(const Color_Unit &Col) const
	{
		return!(*this == Col);
	}
};


class Compression_Valence_Component : 
  public mepp_component
{
	public:
		//default constructor
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

			m_VC[0] = 0.;
			m_VC[1] = 0.;
			m_VC[2] = 0.;
			m_Dist = 0.5;

			// MEPP 2
			componentName = "Compression_Valence_Component";
			init = 1;		
		}

		~Compression_Valence_Component() 
		{}
		

		// Main Function
		double Main_Function(Polyhedron &pMesh, const char* File_Name,const int &Qbit,const int &NVertices, const bool Normal_flipping,
			                 const bool Use_metric, const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const bool Compression_selected,
						     const bool Adaptive_quantization,unsigned &Number_layers, unsigned &Init_number_vertices,unsigned &Final_number_vertices, 
						     unsigned &Connectivity_size, unsigned & Color_size, unsigned & Total_size, const unsigned & Initial_file_size);


		//Initialization
		void Global_Initialization(Polyhedron &pMesh, const int & qbit, const char * File_name);
			void Quantization(Polyhedron &pMesh);
			void Color_Initialization(Polyhedron &pMesh);
			void Color_Quantization(Polyhedron &pMesh);		
			void Multiple_Components_Initialization(Polyhedron & pMesh, const int & Qbit);
		
		
		// Mesh Simplification
		void Simplification(Polyhedron &pMesh, const int & NVertices, const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value);		
			int Decimation_Conquest(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);													
			int Regulation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);			
		

		// Adaptive Quantization		
		void Adaptive_Quantization(Polyhedron &pMesh, const int & NVertices, const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,const int &Qbit);
			void Adaptive_Quantization_Preparation(Polyhedron &pMesh, double & mrms, double &mrmswrtBB,double &hausdorff, double &hausdorffwrtBB);
			int  Test_Next_Operation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value,const int &Qbit);	
			void Up_Quantization(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & a);
			void Under_Quantization(Polyhedron &pMesh, const int & a);
			void Calculate_Edge_Color_Difference(Polyhedron & pMesh, const int & Component_ID, double & Max_color, double & Mean_color, int & Number_of_vertices);
			void Color_Under_Quantization(Polyhedron &pMesh, const int Component_ID);
			void Color_Up_Quantization(Polyhedron &pMesh, Arithmetic_Codec & Decoder, const int & Component_ID);
		
			// Compression 
		void Compression(Polyhedron &pMesh,const char* File_Name, const int &Qbit, unsigned &Connectivity_size, unsigned & Color_size, unsigned & Total_size, const unsigned & Initial_file_size);				
		    int  Calculate_Connectivity_Rate(void);
			void Calculate_Geometry_Color_Offset_Range();				
			void Remove_Last_Phase_Elements(const int & Component_ID);
			void Write_Base_Mesh(Polyhedron &pMesh, Arithmetic_Codec & Enc, unsigned &Connectivity_size, unsigned &Color_size, const int &Num_color_base_mesh);
		

		// Decompression		
		int    Decompress_Init(Polyhedron &pMesh,unsigned & Initial_file_size, const char* File_Name);
		double Decompress_All(Polyhedron &pMesh);				
		int    Decompress_Each_Step(Polyhedron &pMesh, const char* File_Name);
		int    Decompress_To_Level(Polyhedron &pMesh,const int Wanted_level);
		void   Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & dec, const int & id);		
		void   Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);		

		// Other functions
		void   Separate_Components(Polyhedron &pMesh);		
		double Calculate_Area(Polyhedron & pMesh);
		bool   Error_Projected_Surface(Polyhedron & pMesh, const Halfedge_handle & _h, const int & _Component_ID, const double & Mean_color, const double & Mean_area);
		void   Recalculate_Component_Area(Polyhedron & pMesh, const int & Component_ID, int & Number_facets);

		// Change between real and int (related to quantization bit)
		Point_Int Change_Real_Int(const Point3d &pt, const int & Component_ID);		
		Point3d   Change_Int_Real(const Point_Int &vec, const int & Component_ID);
		
		// To add flags to copied mesh
		void Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh);
		
		// Distance calculation
		void Calculate_Distances(char * originalMesh, char * attackedMesh, double & mrms, double & mrmswrtBB, double & hausdorff, double & hausdorffwrtBB);		
		void Calculate_Color_Distances(char * originalMesh, char * attackedMesh, double & min, double & max, double & mean);		

		void Set_Original_Color_Mesh(Polyhedron &pMesh, const char * Filename);
		void Color_Metric_Roy(Polyhedron &pMesh, double &min, double &max, double &mean, double &rms);
		
		int Mapping_Table_Index_Reordering(Polyhedron &pMesh);

		// Color coding cost - prediction method
		int Coding_Cost_Color_Simple_Prediction();	
		
		// Color coding cost - prediction method
		int Coding_Cost_Each_Layer_Color_Clustering(Polyhedron &pMesh);	
		
		// Measure bits used for decompression
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

		
		int GetResolutionChange(Polyhedron *pMesh, float Prec);

		// To stop the decoder (Used "Previous level" while decoding).
		void Stop_Decoder(void) {this->Decoder.stop_decoder(); }		

		//////////////////////////////////////////
		// Joint Compression Watermarking (JCW) //
		//////////////////////////////////////////

		double Joint_Compression_Watermarking(Polyhedron &pMesh, const int &NVertices, const bool Normal_flipping, const bool Use_metric, const float & Metric_thread, const bool Use_forget_metric, const int &Forget_value, const int &Qbit, int & Number_bits,
											  unsigned int & Connectivity_size,unsigned int & Color_size, unsigned int & Total_size,unsigned int & Initial_file_size);
		int JCW_Decimation_For_Segmentation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);													
		int JCW_Regulation_For_Segmentation(Polyhedron &pMesh,const bool Normal_flipping,const bool Use_metric,const float &Metric_thread, const bool Use_forget_metric,const int &Forget_value, const int & Component_ID);			
		
		//list<Point3d> SP_Moved_Position;
		
		void JCW_Un_Regulation_For_Insertion(Polyhedron &pMesh, const int & Component_ID, list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position, list<vector<int>> & JCW_ERROR);
		void JCW_Un_Decimation_For_Insertion(Polyhedron &pMesh, const int & Component_ID, list<int> & FP_Connectivity, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<Point_Int> & SP_Watermarked_Position, list<vector<int>> & JCW_ERROR);

		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & Direction);
		Point3d JCW_Barycenter_Patch_After_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);
		Point3d JCW_Barycenter_Patch_Before_Removal(const Halfedge_handle & h, const int & valence, const int & Direction);
		void JCW_Un_Decimation_Conquest(Polyhedron &pMesh,Arithmetic_Codec & dec, const int & id);		
		void JCW_Un_Regulation(Polyhedron &pMesh, Arithmetic_Codec & dec, const int & id);		
		void JCW_Un_Regulation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID, list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);
		void JCW_Un_Decimation_For_Region_Detection(Polyhedron & pMesh, const int & Component_ID, list<int> & FP_connect, list<Point3d> & FP_Geo, list<int> & FP_RN);
		
		void Set_Number_Bin(const int & NB)
		{
			this->m_NumberBin = NB;
		}
		int Get_Number_Bin(void)
		{
			return this->m_NumberBin;
		}		
		void Set_Embedding_Level(const int &EL)
		{
			this->m_EmbeddingLevel = EL;
		}
		int Get_Embedding_Level(void)
		{
			return this->m_EmbeddingLevel;
		}
		void Set_Number_Region(const int &NR)
		{
			this->m_NumberRegion = NR;
		}
		int Get_Number_Region(void)
		{
			return this->m_NumberRegion;
		}


		void JCW_Calculate_Mesh_Center(Polyhedron &pMesh);		
		void JCW_Calculate_Radius(Polyhedron &pMesh);
		void JCW_Expand_Mesh(Polyhedron &pMesh);
		void JCW_Quantization(Polyhedron &pMesh);
		void JCW_Generate_Regions_On_Base_Mesh(Polyhedron &pMesh);
		void JCW_Region_Mass_Center_Insert_Watermark(Polyhedron & pMesh, list<Point3d> & FP_Geometry, list<int> & FP_Region_Number, list<Point_Int> & SP_Watermarked_Position, list<Point3d> & SP_Moved_Position, list<Point3d> & SP_Original_Position, list<vector<int>> & JCW_ERROR);
		void JCW_Choose_Valid_Vertices(Polyhedron & pMesh);
		vector<int> JCW_Region_Mass_Center_Extract_Watermark(Polyhedron & pMesh);

		void JCW_Code_Difference_Histogram_Shifting(Polyhedron &pMesh,const int & Component_ID);
		void JCW_Decode_Difference_Histogram_Shifting(Polyhedron &pMesh, const int & Component_ID);

		void Initialize_Spherical_Coordinates(Polyhedron &pMesh);

		void Convert_To_Spherical(const Point3d & Pt, double * Spheric);
		void Convert_To_Cartesian(const double * Spheric, double * Cartesian);	

		int  JCW_Decompress_One_Level(Polyhedron &pMesh, const char* File_Name);
		int  JCW_Decompress_One_Level_Without_Extraction(Polyhedron &pMesh, const char* File_Name);
		
		void Read_Information_To_Hide();
	private:
		
		Polyhedron * UnderquantizedMesh;
		Polyhedron * DecimatedMesh;		

		FILE * LogColor;
		
		vector<bool> IsClosed;    // to know if the mesh is open or closed.
		bool IsColored;
		bool IsOneColor;
		float OnlyColor[3];
		// Number of connectivity symbol types. Without boundary = 5, with boundary = 7;
		//int NummberConnectivitySymbols;

		//CCopyPoly<Polyhedron, Enriched_kernel> Copy_Polyhedron;	//MT
		
		// To encode each type of operation between decimation and increase of quantization resolution
		vector<list<int>>		 ListOperation;
		vector<list<int>>		 Connectivity;
		vector<list<Point_Int>>  Geometry;
		vector<list<int>>		 NumberSymbol;
		vector<list<int>>		 NumberVertices;
			
		list<Point_Int>          InterGeometry;
		list<int>		         InterConnectivity;		

		vector<list<int>> AlphaRange;
		vector<list<int>> AlphaOffset;
		vector<list<int>> GammaRange;
		vector<list<int>> GammaOffset;
		
		// Quantization
		vector<unsigned>	Qbit; 
		vector<float>		xmin;
		vector<float>		ymin;
		vector<float>		zmin;		
		vector<float>		xmax;
		vector<float>		ymax;
		vector<float>		zmax;		
		vector<float>		QuantizationPas;
		
		int			  TotalBits;
		double		  OldDistortion;		
		int			  Smallest_Alpha;
		int			  Smallest_Gamma;				
		
		vector<double> HighestLengthBB;
		vector<double> ComponentVolume;
		vector<double> ComponentArea;
		vector<int>    ComponentNumberVertices;
		
		// Used for adatative quantization.				
		vector<list<int>> QuantizationCorrectVector;
		vector<list<int>> NumberQuantizationLayer;
		
		//for color
		vector<list<int>> NumberProcessedVertices;
		vector<list<int>> ColorChildcellIndex;
		vector<list<int>> ColorEncoderIndex;

		// Colors
		vector<list<Color_Unit> > VertexColor; // contain color error of all removed vertices
		list<Color_Unit> InterVertexColor;		
		
		float C0_Min;
		float C1_Min;
		float C2_Min;

		float Color_Quantization_Step;
		vector<int> NumberColorQuantization;

		int Smallest_C0; // the smallest value of C0 used for preventing from negative sylbols
		int Smallest_C1; // the smallest value of C1 used for preventing from negative sylbols
		int Smallest_C2; // the smallest value of C2 used for preventing from negative sylbols
		
		int C0_Range;
		int C1_Range;
		int C2_Range;
		
		//Restoration des couleurs pour Mapping table
		int ColorDiffMinC0;
		int ColorDiffMinC1;
		int ColorDiffMinC2;

		int ColorDiffRangeC0;
		int ColorDiffRangeC1;
		int ColorDiffRangeC2;			
		
		// mapping table
		vector<Color_Unit> ColorArray; // Color table
		vector<Color_Unit> PredictedColorArray; // Predicted values of colors in color table using prediction based on neighbors.		
		vector<Color_Unit> DifferenceColor; // Difference between original and quantized vertex color. need to turn into lossless coding.
		list<int>		   ColorIndex;
		vector<int>		   ReorderingColorIndex;
		list<int>		   InterColorIndex;		
		vector<int>		   Number_color_index; // contain occurence of each initial color
		vector<int>		   IsKnownIndex;

		
		// Decoder
		Arithmetic_Codec Decoder;		
		
		Adaptive_Data_Model Color_0_Model;
		Adaptive_Data_Model Color_1_Model;
		Adaptive_Data_Model Color_2_Model;	

		Adaptive_Data_Model Index_Model;

		Mesh_roy *Original;

		vector<int> NumberDecimation; // To stock number of Decimation.
		vector<int> NumberChangeQuantization; // to stock number of under_quantization.						
				
		// to eliminate symbols of the last conquest if it is useless.
		int DumpSymbolDecimation;
		int DumpSymbolRegulation;			

		// to know how many steps to go for the decompression step by step.
		int Decompress_count;
		int NumberComponents;

		int GlobalCountOperation;

		vector<int> ComponentOperations;

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

		vector<int> m_N_remained_vertices;
		vector<int> m_N_treated_vertices;
		vector<double> m_Rad_decision;

		list<vector<int>> m_JCW_Move_Error;
		list<int> m_N_Errors;

		Adaptive_Data_Model DM_JCW_MOVE_ERROR;
		
		list<int> JCW_Connectivity;
		list<Point_Int> JCW_Geometry;

		

	// from IHM
	public:
		bool IsDecompress;
		bool IsCompressed;
		int Current_level;
		int Total_layer;
		unsigned Initial_file_size;
		unsigned Compressed_file_size;
		string File_name;
		
		vector<float> Prog;
		vector<float> Ratio;

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
