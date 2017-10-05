///////////////////////////////////////////////////////////////////////////
// Author: Ho LEE
// Year: 2011
// Month: MAY
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

#ifndef Compression_Valence_COMMON_H
#define Compression_Valence_COMMON_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

typedef CGAL::MatrixC33<Enriched_kernel> Matrix;

const int RGB_Range = pow((int)2, 4.0) - 1;
/* Quantization bits used.*/

const int RGB_QUANTIZATION = 8;
const int C0_QUANTIZATION = 8;
const int C1_QUANTIZATION = 8;
const int C2_QUANTIZATION = 8;

const int FREE = -1;
const int CONQUERED = 0;
const int TO_BE_REMOVED = 1;
const int TEMP_FLAG = 2;

const int PLUS = 1;///< The plus tag for retriangulation 
const int MINUS = -1;///< The minus tag for retriangulation 
const int NOSIGN = 0;///< The nosign tag for retriangulation 

const int FIRST_COORDINATE = 1;///< The first coordinate tag for seed gate 
const int SECOND_COORDINATE = 2;///< The second coordinate tag for seed gate
const int OTHER_COORDINATE = -1;///< The other coordinate tag for seed gate



const int DECIMATION_CONNECTIVITY_SYMBOL = 5;///< The decimation connectivity symbol
const int LIMIT_QBIT = 4;///< The limit qbit
const double PI = 3.14159265358979323846264338327950288419716939937510582097494459;///< The pi



 
/**
 \fn	void Init(Polyhedron &pMesh)

 \brief	Initialize all flags of verticeces and facets to FREE and give order to vertices (manifold property check).
This function is called every conquest.

 \param [in,out]	pMesh	The mesh.
 */

void Init(Polyhedron &pMesh)
{		
	int i = 0;	
	
	// vertices flags initialization
	Vertex_iterator pVertex = NULL;
	for (pVertex = pMesh.vertices_begin(); pVertex != pMesh.vertices_end(); i++,pVertex++)
	{
		pVertex->Vertex_Flag = FREE;
		pVertex->Vertex_Number = i;
		pVertex->Vertex_Sign = NOSIGN;
	}

	// facets flag initialization.
	Facet_iterator pFace = NULL;
	for (pFace = pMesh.facets_begin(); pFace != pMesh.facets_end(); pFace++)
		pFace->Facet_Flag = FREE;		
}

/**
 \fn	void Write_SMF(Polyhedron &pMesh, const char *_Name, bool Is_colored)

 \brief	Writes a smf.

 \param [in,out]	pMesh	The mesh.
 \param	_Name			 	The name.
 \param	Is_colored		 	true if is colored.
 */

void Write_SMF(Polyhedron &pMesh, const char *_Name, bool Is_colored)
{	
	string Filename = _Name;
	std::ofstream file(Filename.c_str());
	vector<float> Color_container;

	file << "begin" << endl;

	// output vertices
	for (Vertex_iterator pVert = pMesh.vertices_begin(); pVert != pMesh.vertices_end(); pVert++)
	{
		file << "v "<< pVert->point().x() << " " << pVert->point().y() << " " << pVert->point().z();
		file << endl;

		if (Is_colored)
		{
			Color_container.push_back(pVert->color(0));
			Color_container.push_back(pVert->color(1));
			Color_container.push_back(pVert->color(2));
		}		
	}

	// precompute vertex indices
	pMesh.set_index_vertices();

	// output facets
	for (Facet_iterator pFacet = pMesh.facets_begin(); pFacet != pMesh.facets_end(); pFacet++)
	{
		vector<int> indices;
		Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();

		do
		{
			indices.push_back(pHalfedge->vertex()->tag());
		}
		while(++pHalfedge != pFacet->facet_begin());

		file << "f";
		for (unsigned int i=0; i<indices.size(); i++)
		{
			file << " " << indices[i] + 1;
		}

		file << endl;
	}
	if (Is_colored)
	{
		file<<"bind c vertex"<<endl;
		for (unsigned i = 0; i < Color_container.size()/3; i++)
		{
			file<<"c"<<" ";
			file<<Color_container[3*i]<<" ";			
			file<<Color_container[3*i+1]<<" ";			
			file<<Color_container[3*i+2];
			file<<endl;
		}
	}
	file<<"end";

	file.flush();
	file.close();
}

inline int Signe(const double & x)
{
  return (x < 0.) ? -1 : 1 ;
}


/**
 \fn	double Area_Facet_Triangle(const Halfedge_handle &h)

 \brief	Gives an area of the triangle which contain the halfedge_handle h.

 \param	h	The halfedge_handle.

 \return	.
 */

double Area_Facet_Triangle(const Halfedge_handle &h)
{

	Point3d P = h->vertex()->point();
	Point3d Q = h->next()->vertex()->point();
	Point3d R = h->next()->next()->vertex()->point();

	Vector PQ = Q-P;
	Vector QR = R-Q;

	Vector normal =	CGAL::cross_product(PQ, QR);
	double area = 0.5 * sqrt(normal * normal);

	return area;
}


/**
 \fn	double Area_Facet_Triangle(const Point3d &P,const Point3d &Q, const Point3d &R)

 \brief	to calculate the area of a triangle.

 \param	P	The.
 \param	Q	The.
 \param	R	The.

 \return	.
 */

double Area_Facet_Triangle(const Point3d &P,const Point3d &Q, const Point3d &R)
{
	Vector PQ = Q - P;
        //Vector PR = R - P; // MT
	Vector QR = R - Q;

	Vector normal = CGAL::cross_product(PQ,QR);
	double area = 0.5 * sqrt(normal*normal);

	return area;
}

/**
 \fn	int Estimate_Geometry_Quantization(Polyhedron & pMesh, const double & volume,
 		const double & area, const int & number_vertices)

 \brief	Estimate geometry quantization.

 \param [in,out]	pMesh	The mesh.
 \param	volume			 	The volume.
 \param	area			 	The area.
 \param	number_vertices  	Number of vertices.

 \return	.
 */

int Estimate_Geometry_Quantization(Polyhedron & pMesh, const double & volume, const double & area, const int & number_vertices)
{
	double C = (double)volume / (double)area / number_vertices;

	double a = -1.248;
	double b = -0.954;

	int Q = floor(a * log(C) + b + 0.5);

	return Q;
}


/**
 \fn	Vector Triangle_Normal(const Halfedge_handle &h)

 \brief	Gives a normal vector of the triangle containing the halfedge_handle h.

 
 \param	h	The.

 \return	.
 */

Vector Triangle_Normal(const Halfedge_handle &h)
{
	Point3d P = h->vertex()->point();
	Point3d Q = h->next()->vertex()->point();
	Point3d R = h->next()->next()->vertex()->point();

	Vector PQ = Q-P;
        //Vector PR = R-P; // MT
	Vector QR = R-Q;

	Vector normal = CGAL::cross_product(PQ,QR);
	double length = std::sqrt(normal*normal);
	if (length != 0.0)
		normal = normal / length;

	return normal;
}

 
/**
 \fn	Vector Triangle_Normal(const Point3d & P,const Point3d & Q,const Point3d &R)

 \brief	Gives a normal vector of the triangle formed by three points P Q R in the counterclockwise way..


 \param	P	The.
 \param	Q	The.
 \param	R	The.

 \return	.
 */

Vector Triangle_Normal(const Point3d & P,const Point3d & Q,const Point3d &R)
{
	Vector PQ = Q - P;
        //Vector PR = R - P; // MT
	Vector QR = R - Q;

	Vector normal = CGAL::cross_product(PQ,QR);
	double length = std::sqrt(normal*normal);
	if (length != 0.0)
		normal = normal / length;

	return normal;
}

 
/**
 \fn	int Find_Type(const Halfedge_handle &h, const int & valence)

 \brief	To find a correspondent type to retriangulate.
 
 \param	h	   	The.
 \param	valence	The valence.

 \return	The found type.
 */

int Find_Type(const Halfedge_handle &h, const int & valence)
{
	int type = 0;	

	if (valence == 3)
	{
		if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 1;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 2;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 3;
		else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 4;
	}

	else if (valence == 4)
	{
		if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 5;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 6;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 7;
		else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 8;
	}

	else if (valence == 5)
	{
		if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 9;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 10;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 11;
		else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 12;
	}

	else if (valence == 6)
	{
		if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 13;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 14;
		else if ((h->vertex()->Vertex_Sign == PLUS) && (h->opposite()->vertex()->Vertex_Sign == PLUS))
			type = 15;
		else if ((h->vertex()->Vertex_Sign == MINUS) && (h->opposite()->vertex()->Vertex_Sign == MINUS))
			type = 16;
	}
	return type;
}

// ADPATIVE_QUANTIZATION 

/**
 \fn	void Get_Coefficient_Up_Quantization(const int & Correct_symbol, int coeff[3])

 \brief	Gets a coefficient up quantization.

 
 \param	Correct_symbol	The correct symbol.
 \param	coeff		  	The coeff.
 */

void Get_Coefficient_Up_Quantization(const int & Correct_symbol, int coeff[3])
{
	if (Correct_symbol == 0)
	{
		coeff[0] = -1;
		coeff[1] = -1;
		coeff[2] = -1;
	}
	else if (Correct_symbol == 1)
	{
		coeff[0] = 1;
		coeff[1] = -1;
		coeff[2] = -1;
	}
	else if (Correct_symbol == 2)
	{
		coeff[0] = -1;
		coeff[1] = 1;
		coeff[2] = -1;
	}
	else if (Correct_symbol == 3)
	{
		coeff[0] = 1;
		coeff[1] = 1;
		coeff[2] = -1;
	}
	else if (Correct_symbol == 4)
	{
		coeff[0] = -1;
		coeff[1] = -1;
		coeff[2] = 1;
	}
	else if (Correct_symbol == 5)
	{
		coeff[0] = 1;
		coeff[1] = -1;
		coeff[2] = 1;
	}
	else if (Correct_symbol == 6)
	{
		coeff[0] = -1;
		coeff[1] = 1;
		coeff[2] = 1;
	}
	else if (Correct_symbol == 7)
	{
		coeff[0] = 1;
		coeff[1] = 1;
		coeff[2] = 1;
	}
}
 
/**
 \fn	int Get_Correct_Vector(const int & i, const int & j, const int & k)

 \brief	ADAPTIVE_QUANTIZATION : gets symbols to correct Vector of under_quantization.

 
 \param	i	The.
 \param	j	The.
 \param	k	The.

 \return	The correct vector.
 */

int Get_Correct_Vector(const int & i, const int & j, const int & k)
{
	int Correct_symbol = -1;
	
	if ((i == 0) && (j == 0) && (k == 0))
		Correct_symbol = 0;

	else if ((i == 1) && (j == 0) && (k == 0))
		Correct_symbol = 1;

	else if ((i == 0) && (j == 1) && (k == 0))
		Correct_symbol = 2;

	else if ((i == 1) && (j == 1) && (k == 0))
		Correct_symbol = 3;

	else if ((i == 0) && (j == 0) && (k == 1))
		Correct_symbol = 4;

	else if ((i == 1) && (j == 0) && (k == 1))
		Correct_symbol = 5;

	else if ((i == 0) && (j == 1) && (k == 1))
		Correct_symbol = 6;
	else 
		Correct_symbol = 7;
		
	return Correct_symbol;
}



/*
bool Is_Same_Color(const Color_Unit & Color0, const Color_Unit & Color1)
{
	bool check = false;

	if ((Color0.c0 == Color1.c0) && (Color0.c1 == Color1.c1) && (Color0.c2 == Color1.c2))
		check = true;
	return check;
}
*/


 
/**
 \fn	bool Check_End_Clustering(vector<Color_Unit> v1, vector<Color_Unit> v2)

 \brief	checks if all centroids are unchanged after one iteration.

 
 \param	v1	The first vector<Color_Unit>
 \param	v2	The second vector<Color_Unit>

 \return	true if it succeeds, false if it fails.
 */

bool Check_End_Clustering(vector<Color_Unit> v1, vector<Color_Unit> v2)
{
	double Err = 0;
	for (unsigned i = 0; i < v1.size(); i++)
	{
		double Err_c0 = v1[i].c0 - v2[i].c0;
		double Err_c1 = v1[i].c1 - v2[i].c1;
		double Err_c2 = v1[i].c2 - v2[i].c2;

		Err += sqrt(Err_c0 * Err_c0 + Err_c1 * Err_c1 + Err_c2 * Err_c2);
	}

	Err = Err / (double)v1.size();

	if (Err < 0.0001)
		return true;
	else
		return false;
}

/* .*/

/**
 \fn	Color_Unit Get_Vertex_Color(const Halfedge_handle & h)

 \brief	gets a color of front vertex of h.

 
 \param	h	The.

 \return	The vertex color.
 */

Color_Unit Get_Vertex_Color(const Halfedge_handle & h)
{
	Color_Unit Resulting_color;
	Resulting_color.c0 = h->vertex()->color_int(0);
	Resulting_color.c1 = h->vertex()->color_int(1);
	Resulting_color.c2 = h->vertex()->color_int(2);
	
	return Resulting_color;
	
}


 
/**
 \fn	Color_Unit Get_Vertex_Color(Vertex v)

 \brief	Gets color of vertex v.

 
 \param	v	The vertex v.

 \return	The vertex color.
 */

Color_Unit Get_Vertex_Color(Vertex v)
{	
	Color_Unit Resulting_color;
	Resulting_color.c0 = v.color_int(0);
	Resulting_color.c1 = v.color_int(1);
	Resulting_color.c2 = v.color_int(2);
	
	return Resulting_color;	
}

 
/**
 \fn	Vector Calculate_T1_T2(const Halfedge_handle &h, const Vector & normal, Vector & T2)

 \brief	Calculates the base vectors of new coordinates system which is frenet system.
 
 \param	h			  	The.
 \param	normal		  	The normal.
 \param [in,out]	T2	The second Vector &.

 \return	The calculated t 1 t 2.
 */

Vector Calculate_T1_T2(const Halfedge_handle &h, const Vector & normal, Vector & T2)
{
	Point3d P = h->vertex()->point();
	Point3d Q = h->opposite()->vertex()->point();

	Vector T1 = CGAL::NULL_VECTOR;
	Vector u = P - Q;

	double length = std::sqrt(u*u);
	if ( length != 0)
        u = u/length;

	double product = std::sqrt(u*u) * std::sqrt(normal*normal);

	// cosine
	double dot = u * normal;
	double cosine = 0;
	if (product != 0)
		cosine = dot / product;
	
	if (cosine > 1.0)
		cosine = 1.0;
	if (cosine < -1.0)
		cosine = -1.0;
	
	double cosine_rad = std::acos(cosine);

	double beta_rad = 0;
	if (cosine_rad <= PI/2)
		beta_rad = PI/2 - cosine_rad;
	else
		beta_rad = cosine_rad - PI/2;
	double beta = std::cos(beta_rad);

	if (beta != 0)
		T1 = (u - cosine * normal) / beta;
	else
		T1 = CGAL::NULL_VECTOR;
		
	T2 = CGAL::cross_product(normal,T1);
	
	return T1;
}

 
/**
 \fn	Point_Int Frenet_Rotation(const Point_Int &Dist, const Vector &T1,const Vector &T2,
 		const Vector &normal)

 \brief	 finds a bijection through a rotation transformation in 3D with only integer coordinates.



 \param	Dist  	The distance.
 \param	T1	  	The first const Vector &.
 \param	T2	  	The second const Vector &.
 \param	normal	The normal.

 \return	.
 */

inline Point_Int Frenet_Rotation(const Point_Int &Dist, const Vector &T1,const Vector &T2,const Vector &normal)
{
	Matrix R(T1.x(),T2.x(),normal.x(),T1.y(),T2.y(),normal.y(),T1.z(),T2.z(),normal.z());
	Matrix M = R;
	/*Matrix Rt = */R.transpose(); // MT
	Vector Dist1(Dist.x, Dist.y, Dist.z);
	//Vector Dist2 = Rt * Dist1; // MT
	Matrix D1(1,0,0,  0,1,0,  0,0,1);
	Matrix D2(0,1,0,  1,0,0,  0,0,1);
	Matrix D3(-1,0,0,  0,-1,0,  0,0,1);
	Matrix D4(1,0,0,  0,0,1,  0,1,0);
	Matrix D5(1,0,0,  0,-1,0,  0,0,-1);
	Matrix *S = new Matrix[16];	
		
	// Verify in order to find the smallest rotation angle.
	if (std::abs(M.Get(0,2)) > std::abs(M.Get(1,2)))	
		S[0] = D2;
	else	
		S[0] = D1;	
	M = product_matrices(S[0], M);
	
	if (M.Get(1,2) < 0)
		S[1] = D3;	
	else
		S[1] = D1;
	
	M = product_matrices(S[1],M);
	
	/// first rotation angle : phi;
	double phi = -100;

	if (CGAL::square(M.Get(0,2)) + CGAL::square(M.Get(1,2)) == 0)
		phi = 0;
	else
		phi = Signe(-1 * M.Get(0,2)) * std::acos(M.Get(1,2) / std::sqrt(CGAL::square(M.Get(0,2)) + CGAL::square(M.Get(1,2))));
	
	Matrix R1(std::cos(phi),-std::sin(phi),0,  std::sin(phi),std::cos(phi),0,  0,0,1);

	S[2] = Matrix(1,-std::tan(phi/2),0,  0,1,0,  0,0,1);
	S[3] = Matrix(1,0,0,  std::sin(phi),1,0,  0,0,1);
	S[4] = S[2];

	Matrix R1inv(std::cos(phi),std::sin(phi),0,  -std::sin(phi),std::cos(phi),0,  0,0,1);

	M = product_matrices(R1inv,M);
	
	if (std::abs(M.Get(1,2)) > std::abs(M.Get(2,2)))
		S[5] = D4;
	else
		S[5] = D1;

	M = product_matrices(S[5], M);

	if (M.Get(2,2) < 0)
		S[6] = D5;
	else
		S[6] = D1;
	
	M = product_matrices(S[6], M);
	double psi = -100;
	
	/// Second rotation angle psi.
	if (CGAL::square(M.Get(1,2)) + CGAL::square(M.Get(2,2)) == 0)
		psi = 0;
	else
		psi = Signe(-1*M.Get(1,2)) * std::acos(M.Get(2,2) / std::sqrt(CGAL::square(M.Get(1,2)) + CGAL::square(M.Get(2,2))));

	Matrix R2(1,0,0,  0,std::cos(psi),-std::sin(psi),  0,std::sin(psi),std::cos(psi));
	S[7] = Matrix(1,0,0,  0,1,-std::tan(psi/2),  0,0,1);
	S[8] = Matrix(1,0,0,  0,1,0,  0,std::sin(psi),1);
	S[9] = S[7];

	Matrix R2inv(1,0,0,  0,std::cos(psi),std::sin(psi),  0,-std::sin(psi),std::cos(psi));
	M = product_matrices(R2inv,M);
	
	if (std::abs(M.Get(0,1)) > std::abs(M.Get(1,1)))	
		S[10] = D2;	
	else
		S[10] = D1;
	M = product_matrices(S[10],M);

	if (M.Get(1,1) < 0)
		S[11] = D3;
	else
		S[11] = D1;
	M = product_matrices(S[11],M);
	double theta = -100;
	
	/// Last rotation angle theta.
	if (CGAL::square(M.Get(0,1)) + CGAL::square(M.Get(1,1)) == 0)
		theta = 0;
	else
		theta = Signe(-1 * M.Get(0,1)) * std::acos(M.Get(1,1) / std::sqrt(CGAL::square(M.Get(0,1)) + CGAL::square(M.Get(1,1))));

	S[12] = Matrix(1,-std::tan(theta/2),0,  0,1,0,  0,0,1);
	S[13] = Matrix(1,0,0,  std::sin(theta),1,0,  0,0,1);
	S[14] = S[12];
	
	Matrix R3(std::cos(theta),-std::sin(theta),0,  std::sin(theta),std::cos(theta),0,  0,0,1);	
	Matrix R3inv = inverse_matrix(R3);
	/*Matrix S16 = */product_matrices(R3inv,M); // MT
		
	Vector u(Dist.x,Dist.y,Dist.z);	
	Matrix m_inter;	

	// Procedure of the bijection.
	for (int i=0;i<15;i++)
	{
		if (( i == 0) || ( i == 1) || ( i == 5) || ( i == 6) || ( i == 10) || ( i == 11))
			m_inter = S[i];
		else
			m_inter = Matrix(S[i].Get(0,0), -S[i].Get(0,1), -S[i].Get(0,2),
							 -S[i].Get(1,0), S[i].Get(1,1), -S[i].Get(1,2),
							 -S[i].Get(2,0), -S[i].Get(2,1), S[i].Get(2,2));
		u = m_inter * u;		
		
		int x=0,y=0,z=0;
		x = ceil(u.x() - 0.5);
		y = ceil(u.y() - 0.5);
		z = ceil(u.z() - 0.5);		

		u = Vector((double)x, (double)y,(double)z);				
	}

	Point_Int New_Coordinates;
	New_Coordinates.x = (int)u.x();
	New_Coordinates.y = (int)u.y();
	New_Coordinates.z = (int)u.z();	
	
	delete []S;
	return New_Coordinates;
}


//#define DBG_Inverse_Frenet_Rotation

/**
 \fn	Point_Int Inverse_Frenet_Rotation(const Point_Int &Frenet, const Vector &T1,
 		const Vector &T2,const Vector &normal)

 \brief	Inverse operation of frenet rotation. This permits to refind the original coordinates.

 \param	Frenet	The frenet.
 \param	T1	  	The first const Vector &.
 \param	T2	  	The second const Vector &.
 \param	normal	The normal.

 \return	.
 */

inline Point_Int Inverse_Frenet_Rotation(const Point_Int &Frenet, const Vector &T1,const Vector &T2,const Vector &normal)
{
#ifdef DBG_Inverse_Frenet_Rotation //TODO-elo-rm-dbg 
	static unsigned int dbg_Inverse_Frenet_Rotation_call_cnt = 0;		
	std::cout << __func__ << "  call #" << ++dbg_Inverse_Frenet_Rotation_call_cnt << std::endl;
#endif

	Matrix R(T1.x(),T2.x(),normal.x(),T1.y(),T2.y(),normal.y(),T1.z(),T2.z(),normal.z());
	Matrix M = R;

	Matrix D1(1,0,0,  0,1,0,  0,0,1);
	Matrix D2(0,1,0,  1,0,0,  0,0,1);
	Matrix D3(-1,0,0,  0,-1,0,  0,0,1);
	Matrix D4(1,0,0,  0,0,1,  0,1,0);
	Matrix D5(1,0,0,  0,-1,0,  0,0,-1);
	Matrix *S = new Matrix[16];	
		
	// Verify in order to find the smallest rotation angle.
	if (std::abs(M.Get(0,2)) > std::abs(M.Get(1,2)))	
		S[0] = D2;
	else	
		S[0] = D1;	
	M = product_matrices(S[0],M);
	
	if (M.Get(1,2) < 0)
		S[1] = D3;	
	else
		S[1] = D1;
	
	M = product_matrices(S[1],M);
	
	/// first rotation angle : phi;
	double phi = -100;

	if (CGAL::square(M.Get(0,2)) + CGAL::square(M.Get(1,2)) == 0)
		phi = 0;
	else
		phi = Signe(-1 * M.Get(0,2)) * std::acos(M.Get(1,2) / std::sqrt(CGAL::square(M.Get(0,2)) + CGAL::square(M.Get(1,2))));		
	
	Matrix R1(std::cos(phi),-std::sin(phi),0,  std::sin(phi),std::cos(phi),0,  0,0,1);

	S[2] = Matrix(1,-std::tan(phi/2),0,  0,1,0,  0,0,1);
	S[3] = Matrix(1,0,0,  std::sin(phi),1,0,  0,0,1);
	S[4] = S[2];

	//Matrix R1inv = inverse_matrix(R1);
	Matrix R1inv(std::cos(phi),std::sin(phi),0,  -std::sin(phi),std::cos(phi),0,  0,0,1);

	M = product_matrices(R1inv,M);
	
	if (std::abs(M.Get(1,2)) > std::abs(M.Get(2,2)))
		S[5] = D4;
	else
		S[5] = D1;

	M = product_matrices(S[5],M);

	if (M.Get(2,2) < 0)
		S[6] = D5;
	else
		S[6] = D1;
	
	M = product_matrices(S[6],M);
	double psi = -100;
	
	/// Second rotation angle psi.
	if (CGAL::square(M.Get(1,2)) + CGAL::square(M.Get(2,2)) == 0)
		psi = 0;
	else
		psi = Signe(-1*M.Get(1,2)) * std::acos(M.Get(2,2) / std::sqrt(CGAL::square(M.Get(1,2)) + CGAL::square(M.Get(2,2))));
		
	Matrix R2(1,0,0,  0,std::cos(psi),-std::sin(psi),  0,std::sin(psi),std::cos(psi));
	S[7] = Matrix(1,0,0,  0,1,-std::tan(psi/2),  0,0,1);
	S[8] = Matrix(1,0,0,  0,1,0,  0,std::sin(psi),1);
	S[9] = S[7];

	//Matrix R2inv = inverse_matrix(R2);
	Matrix R2inv(1,0,0,  0,std::cos(psi),std::sin(psi),  0,-std::sin(psi),std::cos(psi));
	M = product_matrices(R2inv,M);
	
	if (std::abs(M.Get(0,1)) > std::abs(M.Get(1,1)))	
		S[10] = D2;	
	else
		S[10] = D1;
	M = product_matrices(S[10],M);


	if (M.Get(1,1) < 0)
		S[11] = D3;
	else
		S[11] = D1;
	M = product_matrices(S[11],M);
	double theta = -100;
	
	/// Last rotation angle theta.
	if (CGAL::square(M.Get(0,1)) + CGAL::square(M.Get(1,1)) == 0)
		theta = 0;
	else
		theta = Signe(-1 * M.Get(0,1)) * std::acos(M.Get(1,1) / std::sqrt(CGAL::square(M.Get(0,1)) + CGAL::square(M.Get(1,1))));
		
	S[12] = Matrix(1,-std::tan(theta/2),0,  0,1,0,  0,0,1);
	S[13] = Matrix(1,0,0,  std::sin(theta),1,0,  0,0,1);
	S[14] = S[12];
	
	Matrix R3(std::cos(theta),-std::sin(theta),0,  std::sin(theta),std::cos(theta),0,  0,0,1);	
	Matrix R3inv = inverse_matrix(R3);
	/*Matrix S16 = */product_matrices(R3inv,M); // MT
	Vector u(Frenet.x, Frenet.y, Frenet.z);
	
	Matrix m_inter;
		
	for (int i=14;i>-1;i--)
	{
		m_inter = S[i];

#ifdef DBG_Inverse_Frenet_Rotation //TODO-elo-rm-dbg
		if( dbg_Inverse_Frenet_Rotation_call_cnt == 1596 )
		{
			std::cout << __func__ << "  mark #1" << "  i=" << i << "  u=" << u << "  m_inter=" << m_inter.Get(0,0) << " " << m_inter.Get(1,0) << " " << m_inter.Get(2,0) << " " << m_inter.Get(0,1) << " " << m_inter.Get(1,1) << " " << m_inter.Get(2,1) << " " << m_inter.Get(0,2) << " " << m_inter.Get(1,2) << " " << m_inter.Get(2,2) << std::endl;
		}
#endif

		u =  -1 * m_inter * u;		

		int x = 0, y = 0, z = 0;
		x = -ceil(u.x() - 0.5);
		y = -ceil(u.y() - 0.5);
		z = -ceil(u.z() - 0.5);

#ifdef DBG_Inverse_Frenet_Rotation //TODO-elo-rm-dbg
		if( dbg_Inverse_Frenet_Rotation_call_cnt == 1596 )
		{
			std::cout << __func__ << "  mark #2" << "  i=" << i << "  u=" << u << "  x=" << x << "  y=" << y << "  z=" << z << std::endl;
		}
#endif

		u = Vector((double)x, (double)y, (double)z);		

#ifdef DBG_Inverse_Frenet_Rotation //TODO-elo-rm-dbg
		if( dbg_Inverse_Frenet_Rotation_call_cnt == 1596 )
		{
			std::cout << __func__ << "  mark #3" << "  i=" << i << "  u=" << u << std::endl;
		}
#endif
	}
	
	Point_Int Dist;
	Dist.x = (int)u.x();
	Dist.y = (int)u.y();
	Dist.z = (int)u.z();

#ifdef DBG_Inverse_Frenet_Rotation //TODO-elo-rm-dbg
	if( dbg_Inverse_Frenet_Rotation_call_cnt == 1596 )
	{
		std::cout << __func__ << "  mark #4" << "  Dist=" << Dist.x << " " << Dist.y << " " << Dist.z << std::endl;
	}
#endif
	
	delete []S;
	return Dist;
}
 
/**
 \fn	bool Is_Geometric_Metric_Violated(const Halfedge_handle &h,const int &type,
 		const unsigned int &valence,const float & Threshold)

 \brief	Query if this object is geometric violated. Here, we define a geometric_metric to preserve maximum the intermediate meshes.


 \param	h		 	The.
 \param	type	 	The type.
 \param	valence  	The valence.
 \param	Threshold	The threshold.

 \return	true if geometric violated, false if not.
 */

bool Is_Geometric_Metric_Violated(const Halfedge_handle &h,const int &type,const unsigned int &valence,const float & Threshold)
{
	bool check = false;
	double volume = 0; // volume;;
	double perimeter = 0; //perimeter;;;
	double area = 0; //area;;
	unsigned int i;

	Halfedge_handle g = h;
	Point3d* Points = new Point3d[valence+1];
	g = g->next(); // g points front vertex

	// Points[0] = coordinates of the front vertex
	// Points[i] = coordinates of the ith vertex in the counterclockwise way.
	Points[0] = g->vertex()->point();
	for (i = 1; i < (valence+1); i++)
	{
		Points[i] = g->opposite()->vertex()->point();
		g = g->prev_on_vertex();
	}

	// caculate perimeter
	Vector* Vectors = new Vector[valence];
	Vectors[0] = Points[1] - Points[valence];
	perimeter = std::sqrt(Vectors[0] * Vectors[0]);

	for (i = 1; i< valence; i++)
	{
		Vectors[i] = Points[(i+1)] - Points[i];
		double per = std::sqrt(Vectors[i] * Vectors[i]);
		perimeter = perimeter + per;
	}

	// calculate volume and area;;;
	if (valence == 3) //[0 1 2 3]
	{
		volume = CGAL::volume(Points[0],Points[1],Points[3],Points[2]);
		volume = abs(volume);

		area = Area_Facet_Triangle(Points[1],Points[2],Points[3]);
		area = abs(area);
	}

	else if (valence == 4)
	{
		if ((type == 5) || (type == 8))// [0 1 2 4], [0 2 3 4]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[4],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[2],Points[4],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[2],Points[3],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;
		}

		else if ((type == 6) || (type == 7))// [0 1 2 3], [0 1 3 4]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[3],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[1],Points[4],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[3]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[1],Points[3],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;
		}
	}

	else if (valence == 5)
	{
		if ((type == 9) || (type == 12)) // [0 1 2 5] , [0 2 4 5], [0 2 3 4]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[5],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[2],Points[5],Points[4]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[2],Points[4],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[2],Points[4],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[2],Points[3],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;


		}

		else if (type == 10) // [0 1 4 5], [0 1 2 4], [0 2 3 4]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[5],Points[4]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[1],Points[4],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[2],Points[4],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[4],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[2],Points[3],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

		}

		else if (type == 11)// [0 1 2 3], [0 1 3 5], [0 3 4 5]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[3],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[1],Points[5],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[3],Points[5],Points[4]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[3]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[1],Points[3],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[3],Points[4],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

		}
	}

	else if (valence == 6)
	{
		if ((type == 13) || (type == 16)) //[0 1 2 6], [0 2 3 4] , [0 4 5 6], [0 2 4 6]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[6],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[2],Points[4],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[4],Points[6],Points[5]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[2],Points[6],Points[4]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[6]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[2],Points[3],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[4],Points[5],Points[6]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[2],Points[6],Points[4]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;
		}

		else if ((type == 14) || (type == 15))// [ 0 1 2 3] [ 0 3 4 5] [ 0 1 5 6] [ 0 1 3 5]
		{
			double vol = CGAL::volume(Points[0],Points[1],Points[3],Points[2]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[3],Points[5],Points[4]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[1],Points[6],Points[5]);
			vol = abs(vol);
			volume = volume + vol;

			vol = CGAL::volume(Points[0],Points[1],Points[5],Points[3]);
			vol = abs(vol);
			volume = volume + vol;

			double area_intermediate = Area_Facet_Triangle(Points[1],Points[2],Points[3]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[3],Points[4],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[1],Points[5],Points[6]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;

			area_intermediate = Area_Facet_Triangle(Points[1],Points[3],Points[5]);
			area_intermediate = abs(area_intermediate);
			area = area + area_intermediate;
		}
	}

	double volume_3_root = 0.0;

	double volume_measured = 0.0;

	if (volume == 0)
		volume_3_root = volume;
	else
		volume_3_root = exp(1.0/3.0*(double)log(volume));


	// volume_mea
	volume_measured = volume_3_root / (perimeter / (double)valence);

	// if Volume calculated is greater than threshold value, we should give up removal of this vertex
	if (volume_measured > Threshold)
		check = true;

	// tp free memory.
	delete[] Points;
	delete[] Vectors;

	return check;
}


 
/**
 \fn	bool Is_Normal_Flipping_Occured(const Halfedge_handle &h,const unsigned &valence)

 \brief	Query if removal of the front vertex can cause a normal flipping problem.
 
 \param	h	   	The.
 \param	valence	The valence.

 \return	true if normal flipping occured, false if not.
 */

bool Is_Normal_Flipping_Occured(const Halfedge_handle &h,const unsigned &valence)
{
	int type = Find_Type(h, valence);
	bool check = false;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Halfedge_handle g = h;
	g = g->next();

    Vector V_normal = h->next()->vertex()->normal();

	double length = std::sqrt(V_normal*V_normal);
	if (length != 0)
		V_normal = V_normal / length;

	Point3d *Points = new Point3d[valence];
  Vector  *Normal = new Vector[valence - 2];


	for (unsigned int i=0;i<valence;i++)
	{
		Points[i] = g->opposite()->vertex()->point();
		g = g->prev_on_vertex();
	}

	for (unsigned int j = 0;j<(valence - 2);j++)
	{
	    Normal[j] = CGAL::NULL_VECTOR;
	}

    // if valence = 3 , there is no normal flipping;
	if (valence == 3)
	{
	    check = false;
	}

	// quadrangle
	else if ((type == 5) || (type == 8)) // 0 1 3 , 1 2 3
	{
		Normal[0] = Triangle_Normal(Points[0],Points[1],Points[3]);
		Normal[1] = Triangle_Normal(Points[1],Points[2],Points[3]);
	}

	else if (( type == 6) || (type == 7)) // 0 1 2 , 0 2 3
	{
		Normal[0] = Triangle_Normal(Points[0],Points[1],Points[2]);
		Normal[1] = Triangle_Normal(Points[0],Points[2],Points[3]);
	}

	// pentagone
	else if ((type == 9) || (type == 12)) // 0 1 4 , 1 2 3 , 1 3 4
	{
		Normal[0] = Triangle_Normal(Points[0],Points[1],Points[4]);
		Normal[1] = Triangle_Normal(Points[1],Points[2],Points[3]);
		Normal[2] = Triangle_Normal(Points[1],Points[3],Points[4]);
	}


	else if (type == 10) // 0 1 3 , 1 2 3 , 0 3 4
	{
		Normal[0] = Triangle_Normal(Points[0],Points[1],Points[3]);
		Normal[1] = Triangle_Normal(Points[1],Points[2],Points[3]);
		Normal[2] = Triangle_Normal(Points[0],Points[3],Points[4]);
	}

	else if (type == 11) // 0 1 2 , 0 2 4 , 2 3 4
	{
		Normal[0] = Triangle_Normal(Points[0],Points[1],Points[2]);
		Normal[1] = Triangle_Normal(Points[0],Points[2],Points[4]);
		Normal[2] = Triangle_Normal(Points[2],Points[3],Points[4]);
	}

	// Hexagone

	else if ((type == 13) || (type == 16)) // 0 1 5,  1 2 3 , 3 4 5, 1 3 5
	{
    Normal[0] = Triangle_Normal(Points[0],Points[1],Points[5]);
		Normal[1] = Triangle_Normal(Points[1],Points[2],Points[3]);
		Normal[2] = Triangle_Normal(Points[3],Points[4],Points[5]);
    Normal[3] = Triangle_Normal(Points[1],Points[3],Points[5]);
	}

	else if ((type == 14) || (type == 15)) // 0 1 2 , 2 3 4 , 4 5 0, 0 2 4
	{
	  Normal[0] = Triangle_Normal(Points[0],Points[1],Points[2]);
		Normal[1] = Triangle_Normal(Points[2],Points[3],Points[4]);
		Normal[2] = Triangle_Normal(Points[4],Points[5],Points[0]);
    Normal[3] = Triangle_Normal(Points[0],Points[2],Points[4]);

	}

	for (unsigned int i=0;i<(valence-2);i++)
	{
        double length_normal = std::sqrt(Normal[i]*Normal[i]);
        if (length_normal != 0)
            Normal[i] = Normal[i] / length_normal;
	}

    for (unsigned int i = 0 ; i < (valence - 2); i++)
    {
        double cosine = V_normal * Normal[i];
        double cosine_rad = std::acos(cosine);

        if (cosine_rad >= PI/2)
            check = true;
    }
	delete []Points;
    delete []Normal;

	return check;
}


 
/**
 \fn	bool Is_Border_Vertex(const Halfedge_handle & h)

 \brief	Query if 'h' is border vertex.
 
 \param	h	The.

 \return	true if border vertex, false if not.
 */

bool Is_Border_Vertex(const Halfedge_handle & h)
{
	bool check = false;
	Halfedge_around_vertex_circulator hvc = h->vertex_begin();
	Halfedge_around_vertex_circulator hvc_end = hvc;
	
	CGAL_For_all(hvc,hvc_end)
	{
		if (hvc->is_border_edge())
		{
			check = true;
			break;
		}
	}
	
	return check;	
}

 
/**
 \fn	Vector Normal_Patch(const Halfedge_handle& const_h, const unsigned int &valence)

 \brief	calculates a normal vector of a patch caused by a removal of a front vertex.

 
 \param	const_h	The constant h.
 \param	valence	The valence.

 \return	.
 */

Vector Normal_Patch(const Halfedge_handle& const_h, const unsigned int &valence)
{
	Halfedge_handle h = const_h;
	int type = Find_Type(h, valence);
	
	double area[5]={0,0,0,0,0};
	
	Vector *normals = new Vector[5];
	Vector normal = CGAL::NULL_VECTOR;
	for (int i = 0;i<5;i++)
	{
		normals[i] = CGAL::NULL_VECTOR;
	}	
	
	// Triangle
	if ((type == 1) || (type == 2) || (type == 4))
	{
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);
	}
	else if (type == 3)
	{
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);
	}

	// quadrangle
	else if ((type == 5) || (type == 8))
	{
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);
		h = h->prev()->opposite();
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);
	}

	else if (( type == 6) || (type == 7))
	{
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);

		h = h->next()->opposite();
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);
	}

	// pentagone
	else if ((type == 9) || (type == 12))
	{
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);

		h = h->prev()->opposite();
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);

		h = h->next()->opposite();
		normals[3] = Triangle_Normal(h);
		area[3] = Area_Facet_Triangle(h);
	}


	else if (type == 10)
	{
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);

		h = h->next()->opposite();
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);

		h = h->prev()->opposite();
		normals[3] = Triangle_Normal(h);
		area[3] = Area_Facet_Triangle(h);

	}

	else if (type == 11)
	{
		Halfedge_handle g = h;
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);

		h = g->prev()->opposite();
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);

		h = g->next()->opposite();
		normals[3] = Triangle_Normal(h);
		area[3] = Area_Facet_Triangle(h);
	}

	// Hexagone
	else if ((type == 13) || (type == 16))
	{
		Halfedge_handle g;
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);

		h = h->prev()->opposite();
		g = h;
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);

		h = g->prev()->opposite();
		normals[3] = Triangle_Normal(h);
		area[3] = Area_Facet_Triangle(h);

		h = g->next()->opposite();
		normals[4] = Triangle_Normal(h);
		area[4] = Area_Facet_Triangle(h);
	}

	else if ((type == 14) || (type == 15))
	{
		Halfedge_handle g;
		normals[1] = Triangle_Normal(h);
		area[1] = Area_Facet_Triangle(h);

		h = h->next()->opposite();
		g = h;
		normals[2] = Triangle_Normal(h);
		area[2] = Area_Facet_Triangle(h);

		h = g->prev()->opposite();
		normals[3] = Triangle_Normal(h);
		area[3] = Area_Facet_Triangle(h);

		h = g->next()->opposite();
		normals[4] = Triangle_Normal(h);
		area[4] = Area_Facet_Triangle(h);
	}

	for (unsigned int i=0;i<(valence-2);i++)
		area[0] = area[0] + area[i+1];

	for (unsigned int i=0;i<(valence-2);i++)
		area[i+1] = area[i+1]/area[0];
	
	for (unsigned int i=0;i<(valence-2);i++)
		normal = normal + area[i+1] * normals[i+1];
	
	if (area[0] == 0.0)
		normal = CGAL::NULL_VECTOR;
	
	double length = std::sqrt(normal*normal);
	if (length != 0)
		normal = normal / length;
	
	delete []normals;
	return normal;
}
 
/**
 \fn	Point3d Barycenter_Patch_After_Removal(const Halfedge_handle & h,const int &valence)

 \brief	gives the position of the barycenter of the patch for decimation conquest.

 
 \param	h	   	The.
 \param	valence	The valence.

 \return	.
 */

Point3d Barycenter_Patch_After_Removal(const Halfedge_handle & h,const int &valence)
{
	Halfedge_handle g = h;
	double x = 0, y = 0, z = 0;    
	for (int i=0;i<valence;i++)
	{
		Point3d pt = g->vertex()->point();
		x += pt.x();
		y += pt.y();
		z += pt.z();
		
		g = g->next();
	}	
	
	x /= valence;
	y /= valence;
	z /= valence;
	
	return Point3d(x,y,z);
}


 
/**
 \fn	Point3d Barycenter_Patch_Before_Removal(const Halfedge_handle & h)

 \brief	gives the position of the barycenter of the patch for regulation conquest.

 
 \param	h	The.

 \return	.
 */

Point3d Barycenter_Patch_Before_Removal(const Halfedge_handle & h)
{
	Halfedge_handle g = h;
	int valence = (int)g->next()->vertex_degree();	

	double x = 0., y = 0., z = 0.;

	Halfedge_around_vertex_circulator vh_it = g->next()->vertex_begin();
    Halfedge_around_vertex_circulator vh_it_end = vh_it;

	CGAL_For_all(vh_it, vh_it_end)
	{
		Point3d pt = vh_it->opposite()->vertex()->point();
		x += pt.x();
		y += pt.y();
		z += pt.z();		
	}	
	
	x = x / (double)valence;
	y = y / (double)valence;
	z = z / (double)valence;
	
	return Point3d(x,y,z);
}

//#define DEBUG_Retriangulation

/**
 \fn	void Retriangulation(Polyhedron &pMesh, const Halfedge_handle & ch,
 		const unsigned &valence, const unsigned & Vertex_number, const int & Component_ID)

 \brief	Retriangulates the hole left by a removal of a vertex.

 
 \param [in,out]	pMesh	The mesh.
 \param	ch				 	The ch.
 \param	valence			 	The valence.
 \param	Vertex_number	 	The vertex number.
 \param	Component_ID	 	Identifier for the component.
 */

void Retriangulation(Polyhedron &pMesh, const Halfedge_handle & ch, const unsigned &valence, const unsigned & Vertex_number, const int & Component_ID)
{
	int type = Find_Type(ch, valence);
	Halfedge_handle h;
	h = ch;
	Halfedge_handle g;
	g = h->next();
	

	// Triangle
	if ((type == 1) || (type == 2) || (type == 4))
	{
		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;
				
		h->facet()->normal() = Triangle_Normal(h);
		
		h->facet()->Component_Number = Component_ID;

		if (g->vertex()->Vertex_Sign == NOSIGN)		
			g->vertex()->Vertex_Sign = PLUS;
	}
	else if (type == 3)
	{
		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;
		
		h->facet()->Component_Number = Component_ID;
		
		h->facet()->normal() = Triangle_Normal(h);
		
		if (g->vertex()->Vertex_Sign == NOSIGN)
			g->vertex()->Vertex_Sign = MINUS;
	}

	// quadrangle
	else if ((type == 5) || (type == 8))
	{
		if (g->vertex()->Vertex_Sign == NOSIGN)
			g->vertex()->Vertex_Sign = PLUS;
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = MINUS;

		h = h->prev();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h, g); 
		
		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		h->facet()->Facet_Flag = TO_BE_REMOVED;		
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;

	}

	else if (( type == 6) || (type == 7))
	{
		if (g->vertex()->Vertex_Sign == NOSIGN)
			g->vertex()->Vertex_Sign = MINUS;
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = PLUS;

		g = g->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);
		
		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		h->facet()->Facet_Flag = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
	}
	// pentagone
	else if ((type == 9) || (type == 12))
	{
		if (g->vertex()->Vertex_Sign == NOSIGN)
			g->vertex()->Vertex_Sign = PLUS;
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = MINUS;
		if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign = PLUS;

		h = h->prev();

#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);
		
		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;		
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		g = h->next()->next();
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);
		
		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
	}

	else if (type == 10)
	{
		if (g->vertex()->Vertex_Sign == NOSIGN)
			g->vertex()->Vertex_Sign = PLUS;
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = MINUS;
		if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign = PLUS;

		g = h;
		h = h->prev()->prev();

#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);	

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		g = h->next();
		h = h->prev();		

#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);
		
		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;

	}

	else if (type == 11)
	{
		if (g->vertex()->Vertex_Sign == NOSIGN)
			g->vertex()->Vertex_Sign = MINUS;
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = PLUS;
		if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign = MINUS;

		g = g->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
		
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		g = h->next()->next();
		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;
		
		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
	}

	// Hexagone
	else if ((type == 13) || (type == 16))
	{
		if (g->vertex()->Vertex_Sign ==NOSIGN)
			g->vertex()->Vertex_Sign = PLUS;
		if (g->next()->vertex()->Vertex_Sign ==NOSIGN)
			g->next()->vertex()->Vertex_Sign = MINUS;
		if (g->next()->next()->vertex()->Vertex_Sign ==NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign = PLUS;
		if (g->next()->next()->next()->vertex()->Vertex_Sign ==NOSIGN)
			g->next()->next()->next()->vertex()->Vertex_Sign = MINUS;

		h = h->prev();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;
		
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());		
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
		
		
		g = h->next()->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
		
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		g = h->next()->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;
		
		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
	}

	else if ((type == 14) || (type == 15))
	{
		if (g->vertex()->Vertex_Sign ==NOSIGN)
			g->vertex()->Vertex_Sign = MINUS;
		if (g->next()->vertex()->Vertex_Sign ==NOSIGN)
			g->next()->vertex()->Vertex_Sign = PLUS;
		if (g->next()->next()->vertex()->Vertex_Sign ==NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign = MINUS;
		if (g->next()->next()->next()->vertex()->Vertex_Sign ==NOSIGN)
			g->next()->next()->next()->vertex()->Vertex_Sign = PLUS;


		g = g->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
		
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next()->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;

		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;

		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		g = h->next()->next();		
#ifdef DEBUG_Retriangulation
		std::cout << "split face" << std::endl; //TODO-elo-rm-dbg
#endif
		h = pMesh.split_facet(h,g);

		h->facet()->Component_Number = Component_ID;
		h->opposite()->facet()->Component_Number = Component_ID;
		
		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());
		
		h->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->facet()->Patch_Index = Vertex_number;
		h->opposite()->facet()->Facet_Flag = TO_BE_REMOVED;
		//h->opposite()->facet()->Patch_Index = Vertex_number;
	}
}

/**
 \fn	bool Is_Border_Manifold_Property_Violated(const Halfedge_handle & g,
 		const Halfedge_handle & First_border_edge, const int & Type, const int & Valence)

 \brief	Query if this object is border manifold property violated.

 
 \param	g				 	The.
 \param	First_border_edge	The first border edge.
 \param	Type			 	The type.
 \param	Valence			 	The valence.

 \return	true if border manifold property violated, false if not.
 */

bool Is_Border_Manifold_Property_Violated(const Halfedge_handle & g, const Halfedge_handle & First_border_edge, const int & Type, const int & Valence)
{	
	bool res = false;

	// Added by Adrien Maglo.
	// Test if we are not in the case of a mesh pimple.
	
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
		   res = true;
		   break;
	   }
   }
	
	
	/*if(Valence == 3)
	{
		int Ref = Border_edge->next()->vertex()->Vertex_Number;
		
		Halfedge_around_vertex_circulator Hvc = Border_edge->opposite()->vertex()->vertex_begin();
		Halfedge_around_vertex_circulator Hvc_end = Hvc;
		
		CGAL_For_all(Hvc, Hvc_end)
		{
			if(Hvc->opposite()->vertex()->Vertex_Number == Ref)
				res = true;
		}
	}

	else if(Valence == 4)
	{
		if( ((Number_jump == 0) && ((Type == 5) || (Type == 8))) ||
			((Number_jump == 1) && ((Type == 6) || (Type == 7))) ||
			((Number_jump == 2) && ((Type == 5) || (Type == 8))) )
		{
			int Ref = Border_edge->opposite()->prev()->opposite()->next()->vertex()->Vertex_Number;

			Halfedge_around_vertex_circulator Hvc = Border_edge->opposite()->vertex()->vertex_begin();
			Halfedge_around_vertex_circulator Hvc_end = Hvc;
			
			CGAL_For_all(Hvc, Hvc_end)
			{
				if(Hvc->opposite()->vertex()->Vertex_Number == Ref)
					res = true;
			}
		}
		else if ( ((Number_jump == 0) && ((Type == 6) || (Type == 7))) ||
				  ((Number_jump == 1) && ((Type == 5) || (Type == 8))) ||
				  ((Number_jump == 2) && ((Type == 6) || (Type == 7))) )	
		{
			int Ref = Border_edge->next()->vertex()->Vertex_Number;

			Halfedge_around_vertex_circulator Hvc = Border_edge->opposite()->next()->vertex()->vertex_begin();
			Halfedge_around_vertex_circulator Hvc_end = Hvc;
			
			CGAL_For_all(Hvc, Hvc_end)
			{
				if(Hvc->opposite()->vertex()->Vertex_Number == Ref)
					res = true;
			}
		}
	}*/
	return res;
}
 
/**
 \fn	bool Is_Manifold_Property_Violated(const Halfedge_handle & h, const int &type,
 		const int &valence)

 \brief	Query if removal of this vertex would violate the manifold_property or not.

 
 \param	h	   	The.
 \param	type   	The type.
 \param	valence	The valence.

 \return	true if manifold property violated, false if not.
 */

bool Is_Manifold_Property_Violated(const Halfedge_handle & h, const int &type,const int &valence)
{
	bool check = false;
	Halfedge_handle g = h;
	int* Points_index = new int[valence];

	// if valence is 3, no new edge is inserted, so always safe to remove.
	if (valence == 3)
		return false;

	else
	{
		// Points_index[  ] contains all boundary vertices' indices (ordered in counterclockwise)
		Points_index[0] = g->vertex()->Vertex_Number;
		g = g->next(); // g points center vertex;

		for (int i=1; i<valence; i++)
		{
			g = g->prev_on_vertex();// around the vertex in the counterclockwise way.
			Points_index[i] = g->opposite()->vertex()->Vertex_Number;
		}

		// quadrangle
		if (valence == 4)
		{
			if ((type == 5) || (type == 8))
			{
				g = h->opposite();
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[1])
						check = true;
				}
			}

			else if (( type == 6) || (type == 7))
			{
				g = h;
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[2])
						check = true;;
				}

			}
		}

		//pendtagone : 2 edges to verify
		if (valence == 5)
		{
			if ((type == 9) || (type == 12))
			{
				g = h->opposite();
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[1])
						check = true;
				}

				g = h->next()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[3])
						check = true;
				}
			}

			else if (type == 10)
			{
				g = h;
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[3])
						check = true;
				}

				g = h->next()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[3])
						check = true;
				}
			}

			else if (type == 11)
			{
				g = h;
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[2])
						check = true;
				}

				g = h->opposite();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[2])
						check = true;
				}

			}
		}

		//hexagone

		if (valence == 6)
		{
			if ((type == 13) || (type == 16))
			{
				g = h->opposite();
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[1])
						check = true;
				}

				g = h->opposite();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[3])
						check = true;
				}

				g = h->next()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[3])
						check = true;
				}
			}

			else if ((type == 14) || (type == 15))
			{
				g = h;
				Halfedge_around_vertex_circulator Hvc = g->vertex_begin();
				Halfedge_around_vertex_circulator Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[2])
						check = true;
				}

				g = h;
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[4])
						check = true;
				}

				g = h->prev()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number == Points_index[2])
						check = true;
				}
			}
		}
	}

	delete []Points_index;
	return check;
}


/*	To check if two 3-tuples of integers are same or not.
	output : true if thw points are the same, false if not. */
/*
bool Is_Same_Coordinate(Point_Int fv, Point_Int sv)
{
	bool check = true;
	
	if (fv.x != sv.x)
		check = false;
	if (fv.y != sv.y)
		check = false;
	if (fv.z != sv.z)
		check = false;
		
	return check;
}
*/


// To match first vertices flags between two meshes. This function is needed for under_quantization().
//void Attibute_Seed_Gate_Flag(Polyhedron &Original_mesh, Polyhedron &New_mesh)
//{
//	Vertex_iterator pVertex = NULL;
//	Vertex_iterator pVertex2 = New_mesh.vertices_begin();
//
//	for (pVertex = Original_mesh.vertices_begin(); pVertex != Original_mesh.vertices_end(); pVertex++)
//	{
//		pVertex2->Seed_Edge = pVertex->Seed_Edge;
//		pVertex2++;
//	}
//}

/**
 \fn	void LAB_To_LCH(const float * LAB, float * LCH)

 \brief	Lab to lch.


 \param	LAB			   	The lab.
 \param [in,out]	LCH	If non-null, the lch.
 */

void LAB_To_LCH(const float * LAB, float * LCH)
{
	// L
	LCH[0] = LAB[0];
	
	// C
	LCH[1] = sqrt(LAB[1]*LAB[1] + LAB[2] * LAB[2]);
	
	// H 
	LCH[2] = atan2(LAB[2], LAB[1]) * 180 / PI;
	
	if (LCH[2] < 0.)
		LCH[2] += 360.;
	if (LCH[2] >= 360.)
		LCH[2] -= 360.;

}

/**
 \fn	double Color_Distance_ICE94(const float * LCH1, const float * LCH2)

 \brief	Color distance ice 94.

 \param	LCH1	The first lch.
 \param	LCH2	The second lch.

 \return	.
 */

double Color_Distance_ICE94(const float * LCH1, const float * LCH2)
{
	double L1 = LCH1[0], L2 = LCH2[0];
	double C1 = LCH1[1], C2 = LCH2[1];
	double H1 = LCH1[2], H2 = LCH2[2];

	double dL = L2 - L1;
	double dC = C2 - C1;
	double dH = H2 - H1;

	//double sL = 1;
	double sC = 1 + 0.045 * C1; //C1
	double sH = 1 + 0.015 * C1; // C1
	
	//double kL = 1.;
	//double kC = 1.;
	//double kH = 1.;	

	return sqrt(dL*dL + dC*dC/sC/sC + dH*dH/sH/sH);
}

/**
 \fn	double Color_Distance_CMC21(const float * LCH1, const float * LCH2)

 \brief	Color distance cmc 21.


 \param	LCH1	The first lch.
 \param	LCH2	The second lch.

 \return	.
 */

double Color_Distance_CMC21(const float * LCH1, const float * LCH2)
{
	double L1 = LCH1[0], L2 = LCH2[0];
	double C1 = LCH1[1], C2 = LCH2[1];
	double H1 = LCH1[2], H2 = LCH2[2];

	double dL = L2 - L1;
	double dC = C2 - C1;
	double dH = H2 - H1;


	double F = sqrt(pow(C1,double(4.))/(pow(C1,double(4.)) + 1900));
	double T;
	if ((H1 <= 345) && (H1 >= 164))
	{
		double rad = (H1 + 168) * PI / (double)180;
		T = abs(0.2*cos(rad)) + 0.56;		
	}
	else
	{
		double rad = (H1 + 35) * PI / (double)180;
		T = abs(0.4*cos(rad)) + 0.36;
	}
		
	double sL, sC, sH;
	
	if (L1 < 16)
		sL = 0.511;
	else
		sL = 0.040975 * L1 / (1 + 0.01765 * L1);

	sC = 0.0638*C1 / (1 + 0.0131*C1) + 0.638;
	
	sH = sC*(F*T+1.-F);
	

	double kL = 2.;
	double kC = 1.;
	//double kH = 1.;	

	double vL = dL/kL/sL;
	double vC = dC/kC/sC;
	double vH = dH/sH;


	return sqrt(vL*vL + vC*vC + vH*vH);
}

/**
 \fn	double Color_Distance_CMC11(const float * LCH1, const float * LCH2)

 \brief	Color distance cmc 11.


 \param	LCH1	The first lch.
 \param	LCH2	The second lch.

 \return	.
 */

double Color_Distance_CMC11(const float * LCH1, const float * LCH2)
{
	double L1 = LCH1[0], L2 = LCH2[0];
	double C1 = LCH1[1], C2 = LCH2[1];
	double H1 = LCH1[2], H2 = LCH2[2];

	double dL = L2 - L1;
	double dC = C2 - C1;
	double dH = H2 - H1;


	double F = sqrt(pow(C1,double(4.))/(pow(C1,double(4.)) + 1900));
	double T;
	if ((H1 <= 345) && (H1 >= 164))
	{
		double rad = (H1 + 168) * PI / (double)180;
		T = abs(0.2*cos(rad)) + 0.56;		
	}
	else
	{
		double rad = (H1 + 35) * PI / (double)180;
		T = abs(0.4*cos(rad)) + 0.36;
	}
		
	double sL, sC, sH;
	
	if (L1 < 16)
		sL = 0.511;
	else
		sL = 0.040975 * L1 / (1 + 0.01765 * L1);

	sC = 0.0638*C1/(1 + 0.0131*C1) + 0.638;
	
	sH = sC*(F*T + 1. - F);
	

	double kL = 1.;
	double kC = 1.;
	//double kH = 1.;	

	double vL = dL/kL/sL;
	double vC = dC/kC/sC;
	double vH = dH/sH;


	return sqrt(vL*vL + vC*vC + vH*vH);
}


//

/**
 \fn	Color_Unit Get_Average_Vertex_Color_Before_Removal(const Halfedge_handle & h,
 		const int & valence)

 \brief	Gets an average color of neighboring vertices ( before removal of front vertex of h)


 \param	h	   	The.
 \param	valence	The valence.

 \return	The average vertex color before removal.
 */

Color_Unit Get_Average_Vertex_Color_Before_Removal(const Halfedge_handle & h, const int & valence)
{
	Halfedge_handle g = h;
	
	Color_Unit Average_color;
	Average_color.c0 = 0;
	Average_color.c1 = 0;
	Average_color.c2 = 0;

	for (int i =0 ; i < valence; i++)
	{			
		Color_Unit col = Get_Vertex_Color(g);
		Average_color.c0 += col.c0;
		Average_color.c1 += col.c1;
		Average_color.c2 += col.c2;
		
		g = g->next()->opposite()->next();
	}
	Color_Unit Resulting_color;
	Resulting_color.c0 = floor((float)Average_color.c0 / valence + 0.5);
	Resulting_color.c1 = floor((float)Average_color.c1 / valence + 0.5);
	Resulting_color.c2 = floor((float)Average_color.c2 / valence + 0.5);	
			
	return Resulting_color;	
}


//

/**
 \fn	Color_Unit Get_Average_Vertex_Color_After_Removal(const Halfedge_handle & h,
 		const int & valence)

 \brief	Gets an average color of neighboring vertices ( After removal of front vertex of h)

 
 \param	h	   	The.
 \param	valence	The valence.

 \return	The average vertex color after removal.
 */

Color_Unit Get_Average_Vertex_Color_After_Removal(const Halfedge_handle & h, const int & valence)
{		
	Halfedge_handle g = h;
	
	Color_Unit Average_color;
	Average_color.c0 = 0;
	Average_color.c1 = 0;
	Average_color.c2 = 0;

	for (int i =0 ; i < valence; i++)
	{			
		Color_Unit col = Get_Vertex_Color(g);
		Average_color.c0 += col.c0;
		Average_color.c1 += col.c1;
		Average_color.c2 += col.c2;
		
		g = g->next();
	}	

	Color_Unit Resulting_color;	
	Resulting_color.c0 = floor((float)Average_color.c0 / valence + 0.5);
	Resulting_color.c1 = floor((float)Average_color.c1 / valence + 0.5);
	Resulting_color.c2 = floor((float)Average_color.c2 / valence + 0.5);	
				
	return Resulting_color;	
}


 

/**
 \fn	double Get_Angle_Weight_Youn(const Halfedge_handle & h)

 \brief	Calculate weight of h->vertex for Youn's method.

 

 \param	h	The.

 \return	The angle weight youn.
 */

double Get_Angle_Weight_Youn(const Halfedge_handle & h)
{
	Point3d P = h->next()->vertex()->point();
	Point3d Q = h->vertex()->point();
	Point3d R = h->opposite()->next()->vertex()->point();
	
	Vector V1 = P - Q;
	Vector V2 = R - Q;

	double Denom = std::sqrt(V1*V1) * std::sqrt(V2*V2);

	double Multi = V1*V2;
	
	double Cosine = 0;
	if (Denom != 0.0)
		Cosine = Multi / Denom;
	else
		return 0.0;

	if (Cosine > 1.0)
		Cosine = 1.0;
	if (Cosine < -1.0)
		Cosine = -1.0;

	double Cosinus_rad = acos(Cosine);

	return Cosinus_rad;	
}



//Youn's method for prediction.
Color_Unit Get_Average_Vertex_Color_Youn(const Halfedge_handle & h,const int & valence)
{	
	Halfedge_handle g = h;	
		
	vector<Color_Unit> Neighbor_color;
	vector<double> Weights;		
	
	g = g->next();
	for (int i = 0; i < valence; i++)
	{		
		Neighbor_color.push_back(Get_Vertex_Color(g->opposite())); // color of neighbor (LUV or RGB by Get_Vertex_Color)
		Weights.push_back(Get_Angle_Weight_Youn(g)); // angle formed by front vertex and two neighbor of actual neighbor vertex
		g = g->opposite()->prev();
	}
	double Sum_weights = 0;
	for (int i=0;i<valence;i++)
	{
		Sum_weights += Weights[i];
	}
	
	float Temp_color[3] = {0.0, 0.0, 0.0};	
	for (int i=0;i<valence;i++)
	{
		Temp_color[0] += Weights[i] / Sum_weights * Neighbor_color[i].c0;
		Temp_color[1] += Weights[i] / Sum_weights * Neighbor_color[i].c1;
		Temp_color[2] += Weights[i] / Sum_weights * Neighbor_color[i].c2;
	}
	
	Color_Unit Resulting_color;	
	
	Resulting_color.c0 = (int)floor(Temp_color[0] + 0.5);	
	Resulting_color.c1 = (int)floor(Temp_color[1] + 0.5);		
	Resulting_color.c2 = (int)floor(Temp_color[2] + 0.5);		
	
	return Resulting_color;	
}


 
/**
 \fn	Color_Unit Get_Average_Vertex_Color_Lee(const Halfedge_handle & h, const int & valence)

 \brief	Gets an average vertex color prediction using Lee's method.

 
 \param	h	   	The.
 \param	valence	The valence.

 \return	The average vertex color lee.
 */

Color_Unit Get_Average_Vertex_Color_Lee(const Halfedge_handle & h, const int & valence)
{	
	Halfedge_handle g = h;	
	
	vector<Color_Unit> Neighbor_color;		//contain differents colors among neighbors 
	int Average_c0 = 0, Average_c1 = 0, Average_c2 = 0;
	
	g = g->next();

	for (int i = 0; i < valence; i++)
	{
		Color_Unit col = Get_Vertex_Color(g->opposite());
		
		bool check = false;
		for (unsigned j=0; j < Neighbor_color.size(); j++)
		{
			if (Neighbor_color[j] == col)
			{
				check = true;
				break;
			}
		}
		
		if (!check)		
			Neighbor_color.push_back(col);
		
		Average_c0 += col.c0;
		Average_c1 += col.c1;
		Average_c2 += col.c2;
		
		g = g->opposite()->prev();
	}
				
	Color_Unit Average_color;

	Average_color.c0 = (int)floor((float)(Average_c0 / valence) + 0.5);		
	Average_color.c1 = (int)floor((float)(Average_c1 / valence) + 0.5);
	Average_color.c2 = (int)floor((float)(Average_c2 / valence) + 0.5);	
	
	double Critere_c0 = 50000.0, Critere_c1 = 50000.0, Critere_c2 = 50000.0;
	int Index_c0 = -1, Index_c1 = -1, Index_c2 = -1;	

	// Find for each component among neighbors' color, the component which is the nearst
	// to the average value.
	for (unsigned i=0; i<Neighbor_color.size(); i++)
	{
		int Color_component_0 = Neighbor_color[i].c0;
		int Color_component_1 = Neighbor_color[i].c1;
		int Color_component_2 = Neighbor_color[i].c2;
		
		if (abs(Color_component_0 - Average_color.c0) < Critere_c0)
		{
			Critere_c0 = abs(Color_component_0 - Average_color.c0);
			Index_c0 = i;
		}
		if (abs(Color_component_1 - Average_color.c1) < Critere_c1)
		{
			Critere_c1 = abs(Color_component_1 - Average_color.c1);
			Index_c1 = i;
		}
		if (abs(Color_component_2 - Average_color.c2) < Critere_c2)
		{
			Critere_c2 = abs(Color_component_2 - Average_color.c2);
			Index_c2 = i;
		}			
	}
	
	Color_Unit Resulting_color;
	Resulting_color.c0 = Neighbor_color[Index_c0].c0;
	Resulting_color.c1 = Neighbor_color[Index_c1].c1;
	Resulting_color.c2 = Neighbor_color[Index_c2].c2;
	
	return Resulting_color;		
}



/*
//My method for prediction front vertex's color.
Color_Unit Get_Average_Vertex_Color_Lee(const Halfedge_handle & h, const int & valence)
{	
	Halfedge_handle g = h;		

	vector<Color_Unit> Neighbor_color;		//contain differents colors among neighbors 
	float Average_c0 = 0., Average_c1 = 0., Average_c2 = 0.;
	
	vector<vector<float>> NC_lch;

	
	g = g->next();

	for (int i = 0; i < valence; i++)
	{
		Color_Unit col = Get_Vertex_Color(g->opposite());
		
		vector<float> LCH;
		float LAB[3];
		float Temp_LCH[3];
		LAB[0] = g->opposite()->vertex()->color_float(0);
		LAB[1] = g->opposite()->vertex()->color_float(1);
		LAB[2] = g->opposite()->vertex()->color_float(2);
			
		LAB_To_LCH(LAB, Temp_LCH);
		
		bool check = false;
		for (unsigned j=0; j < Neighbor_color.size(); j++)
		{
			if (Neighbor_color[j] == col)
			{
				check = true;
				break;
			}
		}
		
		if (!check)
		{			
			LCH.push_back(Temp_LCH[0]);
			LCH.push_back(Temp_LCH[1]);
			LCH.push_back(Temp_LCH[2]);

			Neighbor_color.push_back(col);
			NC_lch.push_back(LCH);
		}
		
		Average_c0 += Temp_LCH[0];
		Average_c1 += Temp_LCH[1];
		Average_c2 += Temp_LCH[2];
		
		g = g->opposite()->prev();
	}
	
	float Ref_LCH[3];
	Ref_LCH[0] = Average_c0 / (float)valence;
	Ref_LCH[1] = Average_c1 / (float)valence;
	Ref_LCH[2] = Average_c2 / (float)valence;

	double Dist = 50000.;
	int Corresponding_index = -1;
	
	float Comp_LCH[3];
	// Find for each component among neighbors' color, the component which is the nearst
	// to the average value.
	for (unsigned i=0; i < NC_lch.size(); i++)
	{
		Comp_LCH[0] = NC_lch[i][0];
		Comp_LCH[1] = NC_lch[i][1];
		Comp_LCH[2] = NC_lch[i][2];
		double Temp_dist = Color_Distance_CMC21(Ref_LCH, Comp_LCH);

		if (Temp_dist < Dist)
		{
			Dist = Temp_dist;
			Corresponding_index = i;
		}
	}
	
	return Neighbor_color[Corresponding_index];		
}
*/


 
/**
 \fn	bool Remove_Edges(Polyhedron & pMesh, const Halfedge_handle & h, const int & type)

 \brief	 Remove edges to create a hole.

 
 \param [in,out]	pMesh	The mesh.
 \param	h				 	The.
 \param	type			 	The type.

 \return	true if it succeeds, false if it fails.
 */

bool Remove_Edges(Polyhedron & pMesh, const Halfedge_handle & h, const int & type)
{
	bool check = false;
	Halfedge_handle g = h;

	//triangle

	if ((type == 1) || (type == 2) || (type == 4))
	{
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = PLUS;
	}

	else if (type == 3)
	{
		if (g->next()->vertex()->Vertex_Sign == NOSIGN)
			g->next()->vertex()->Vertex_Sign = MINUS;
	}


	// quadrangle
	else if ((type == 5) || (type == 8))
	{
		//verification
		if (g->prev()->opposite()->facet()->Facet_Flag != FREE)
			check = true;
		if (check == false)
		{
			g = g->prev();
			pMesh.join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = MINUS;
		}
	}

	else if ((type == 6) || (type == 7))
	{
		//verification
		if (g->next()->opposite()->facet()->Facet_Flag != FREE)
			check = true;
		if (check == false)
		{
			g = g->next();
			pMesh.join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = PLUS;
		}
	}

	//pentagone
	else if ((type == 9) || (type == 12))
	{
		g = g->prev()->opposite();
		if (g->facet()->Facet_Flag != FREE)
			check = true;
		g = g->next()->opposite();
		if (g->facet()->Facet_Flag != FREE)
			check = true;

		if (check == false)
		{
			g = h->prev();
			g = pMesh.join_facet(g);
			g = g->next();
			g = pMesh.join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign = PLUS;
		}
	}

	else if (type == 10)
	{
		g = g->next()->opposite();
		if (g->facet()->Facet_Flag != FREE)
			check = true;

		g = g->prev()->opposite();
		if (g->facet()->Facet_Flag != FREE)
			check = true;

		if (check == false)
		{
			g = h->next()->opposite();
			g = pMesh.join_facet(g);
			g = pMesh.join_facet(g);
			

			g = h;
			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign = PLUS;
		}

	}

	else if (type == 11)
	{
		if (g->next()->opposite()->facet()->Facet_Flag != FREE)
			check = true;
		if (g->prev()->opposite()->facet()->Facet_Flag != FREE)
			check = true;

		if (check == false)
		{
			g = g->next();
			g = pMesh.join_facet(g);
			g = g->prev();
			g = pMesh.join_facet(g);

			g = h;

			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign = MINUS;
		}
	}

	else if ((type == 13) || (type == 16))
	{
		g = g->prev()->opposite();
		if (g->facet()->Facet_Flag != FREE)
			check = true;
		if (g->next()->opposite()->facet()->Facet_Flag != FREE)
			check = true;
		if (g->prev()->opposite()->facet()->Facet_Flag != FREE)
			check = true;

		if (check == false)
		{
			g = h->prev()->opposite();
			g = pMesh.join_facet(g);
			g = pMesh.join_facet(g);
			g = pMesh.join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->next()->vertex()->Vertex_Sign = MINUS;
		}
	}

	else if ((type == 14) || (type == 15))
	{
		g = g->next()->opposite();
		if (g->facet()->Facet_Flag != FREE)
			check = true;
		if (g->next()->opposite()->facet()->Facet_Flag != FREE)
			check = true;
		if (g->prev()->opposite()->facet()->Facet_Flag != FREE)
			check = true;

		if (check == false)
		{
			g = h->next()->opposite();
			g = pMesh.join_facet(g);
			g = pMesh.join_facet(g);
			g = pMesh.join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign = PLUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign = MINUS;
			if (g->next()->next()->next()->next()->vertex()->Vertex_Sign == NOSIGN)
				g->next()->next()->next()->next()->vertex()->Vertex_Sign = PLUS;
		}
	}

	return check;
}

 
/**
 \fn	double Sqrt_Color(const Color_Unit & col1, const Color_Unit &col2)

 \brief	calculates square root distance between two colors.

 
 \param	col1	The first col.
 \param	col2	The second col.

 \return	.
 */

double Sqrt_Color(const Color_Unit & col1, const Color_Unit &col2)
{
	int diff0 = col1.c0 - col2.c0;
	int diff1 = col1.c1 - col2.c1;
	int diff2 = col1.c2 - col2.c2;

	const double factor_c0 = 1.0;
	const double factor_c1 = 1.0;
	const double factor_c2 = 1.0;

	double distance = sqrt(factor_c0 * pow((double)diff0,(int)2) 
						 + factor_c1 * pow((double)diff1,(int)2) 
						 + factor_c2 * pow((double)diff2,(int)2));

	return distance;
}









 
/**
 \fn	void RGB_To_YUV(const float r, const float g, const float b, float* yuv)

 \brief	Conversion from RGB to YUV.

 
 \param	r			   	The.
 \param	g			   	The.
 \param	b			   	The.
 \param [in,out]	yuv	If non-null, the yuv.
 */

void RGB_To_YUV(const float r, const float g, const float b, float* yuv)
{
	yuv[0] = 0.299    * r + 0.587   * g + 0.114   * b;
	yuv[1] = -0.14713 * r - 0.28886 * g + 0.436   * b;
	yuv[2] = 0.615    * r - 0.51499 * g - 0.10001 * b;
}

/**
 \fn	void YUV_To_RGB(const float y, const float u, const float v, float* rgb)

\brief	Conversion from YUV to RGB.

 
 \param	y			   	The y coordinate.
 \param	u			   	The u.
 \param	v			   	The v.
 \param [in,out]	rgb	If non-null, the rgb.
 */

void YUV_To_RGB(const float y, const float u, const float v, float* rgb)
{
	rgb[0] = y               + 1.13983 * v;
	rgb[1] = y - 0.39465 * u - 0.5806  * v;
	rgb[2] = y + 2.03211 * u;
}

 

/**
 \fn	void RCT_RGB_To_YCbCr(const float r, const float g, const float b, int* YCbCr)

 \brief	Reversible color transform based on integer-based computation RGB to YCbCr

 
 \param	r				 	The.
 \param	g				 	The.
 \param	b				 	The.
 \param [in,out]	YCbCr	If non-null, the y coordinate cb carriage return.
 */

void RCT_RGB_To_YCbCr(const float r, const float g, const float b, int* YCbCr)
{
	int R = (int)floor(r * RGB_Range + 0.5);
	int G = (int)floor(g * RGB_Range + 0.5);
	int B = (int)floor(b * RGB_Range + 0.5);		

	YCbCr[0] = (int)floor(float(1.0/4.0*(R + 2*G + B) + 0.5));
	YCbCr[1] = B - G;
	YCbCr[2] = R - G;


	//temp : comparison
	//int y = YCbCr[0];
	//int Cb = YCbCr[1];
	//int Cr = YCbCr[2];

	//int G2 = y - (int)floor(float(1.0/4.0*(Cb + Cr) + 0.5));
	//int R2 = Cr + G;
	//int B2 = Cb + G;	
}

/**
 \fn	void RCT_YCbCr_To_RGB(const int y, const int Cb, const int Cr, float * rgb)

 \brief	Reversible color transform based on integer-based computation YCbCR to RGB.
 
 \param	y			   	The y coordinate.
 \param	Cb			   	The cb.
 \param	Cr			   	The carriage return.
 \param [in,out]	rgb	If non-null, the rgb.
 */

void RCT_YCbCr_To_RGB(const int y, const int Cb, const int Cr, float * rgb)
{
	
	int G = y - (int)floor(float(1.0/4.0*(Cb + Cr) + 0.5));
	int R = Cr + G;
	int B = Cb + G;

	rgb[0] = (float)R / RGB_Range;
	rgb[1] = (float)G / RGB_Range;
	rgb[2] = (float)B / RGB_Range;
}



/*
	
*/

/**
 \fn	void ICT_RGB_To_YCbCr(const float r, const float g, const float b, float* YCbCr)

 \brief	Irreversible Color Transform RGB to YCbCr.


 \param	r				 	The.
 \param	g				 	The.
 \param	b				 	The.
 \param [in,out]	YCbCr	If non-null, the y coordinate cb carriage return.
 */

void ICT_RGB_To_YCbCr(const float r, const float g, const float b, float* YCbCr)
{
	YCbCr[0] =   0.299   * r + 0.587   * g + 0.114   * b;
	YCbCr[1] = - 0.16875 * r - 0.33126 * g + 0.5     * b;
	YCbCr[2] =   0.5     * r - 0.41869 * g - 0.08131 * b;
}

/**
 \fn	void ICT_YCbCr_To_RGB(const float y, const float Cb, const float Cr, float * rgb)

 \brief	Irreversible Color Transform yCbCr to RGB.

 \param	y			   	The y coordinate.
 \param	Cb			   	The cb.
 \param	Cr			   	The carriage return.
 \param [in,out]	rgb	If non-null, the rgb.
 */

void ICT_YCbCr_To_RGB(const float y, const float Cb, const float Cr, float * rgb)
{
	rgb[0] = y                + 1.4021 * Cr;
	rgb[1] = y - 0.34413 * Cb - 0.71414 * Cr;
	rgb[2] = y + 1.7718  * Cb;
}


// l = [0, 100], a = [-0.86,0.86], b = [-1.07, 0.94];

/**
 \fn	void RGB_To_LAB(float r, float g, float b, float* lab)

 \brief	Rgb to lab.

 \param	r			   	The.
 \param	g			   	The.
 \param	b			   	The.
 \param [in,out]	lab	If non-null, the lab.
 */

void RGB_To_LAB(float r, float g, float b, float* lab)
{	
	// R 0 -- 1, G 0 -- 1, B 0 -- 1

	//http://www.brucelindbloom.com
	float X, Y, Z;
	//float L;
	float eps = 216.f/24389.f;
	float k = 24389.f/27.f;

	float Xr = 0.964221f;  // reference white D50
	float Yr = 1.0f;
	float Zr = 0.825211f;

	// RGB to XYZ

	// assuming sRGB (D65)
	if (r <= 0.04045)
		r = r / 12.92;
	else
		r = (float)pow((r+0.055)/1.055,2.4);

	if (g <= 0.04045)
		g = g / 12.92;
	else
		g = (float)pow((g+0.055)/1.055,2.4);

	if (b <= 0.04045)
		b = b / 12.92;
	else
		b = (float)pow((b+0.055)/1.055,2.4);


	X =  0.412424*r     + 0.357579*g	+ 0.180464*b;
	Y =  0.212656*r     + 0.715158 *g	+ 0.0721856*b;
	Z =  0.0193324*r    + 0.119193*g	+ 0.950444*b;
	
	// XYZ to LAB

	float xr = X / Xr;
	float yr = Y / Yr;
	float zr = Z / Zr;

	float fx, fy, fz;

	if (xr > eps)		
		fx = (float)pow(xr, (float)(1.f/3.f));
	else
		fx = (k*xr+16.0)/116.0;

	if (yr > eps)
		fy = (float)pow(yr, (float)(1.f/3.f));
	else
		fy = (k*yr+16.0)/116.f;

	if (zr > eps)
		fz = (float)pow(zr, (float)(1.f/3.f));
	else
		fz = (k * zr + 16.0)/116.f;

	lab[0] = (float)(116.0 * fy - 16.0); // L
	lab[1] = (float)(500.0 * (fx - fy)); // A
	lab[2] = (float)(200.0 * (fy - fz)); // B
}

/**
 \fn	void LAB_To_RGB(float L, float A, float B, float* rgb)

 \brief	Lab to rgb.

 \param	L			   	The.
 \param	A			   	a.
 \param	B			   	The.
 \param [in,out]	rgb	If non-null, the rgb.
 */

void LAB_To_RGB(float L, float A, float B, float* rgb)
{
	// LAB TO XYZ

	float eps = 216.f/24389.f;
	float k = 24389.f/27.f;

	float Xr = 0.964221f;  // reference white D50
	float Yr = 1.0f;
	float Zr = 0.825211f;
	//float Xr = 96.4221f;  // reference white D50
	//float Yr = 100.f;
	//float Zr = 82.5211f;

	float fx, fy, fz;
	fy = (L + 16.0) / 116.0;		
	fx = (A / 500.0) + fy;
	fz = fy - B / 200.0;

	float xr, yr, zr;
	float fx_cube = (float)pow(fx, (float)3.0);
	float fy_cube = (float)pow(fy, (float)3.0);
	float fz_cube = (float)pow(fz, (float)3.0);

	if (fx_cube > eps)
		xr = fx_cube;
	else
		xr = (116.0*fx - 16.0) / k;

	if (fy_cube > eps)
		yr = fy_cube;
	else
		yr = (float)(L / k);

	if (fz_cube > eps)
		zr = fz_cube;
	else
		zr = (116.0*fz - 16.0) / k;

	float X, Y, Z;
	X = xr * Xr;
	Y = yr * Yr;
	Z = zr * Zr;


	// XYZ to RGB

	// R 0..1, G 0..1, B 0..1

	float r, g, b;

	// assuming sRGB (D65)

	r =  3.24071*X      + -1.53726*Y	+ -0.498571*Z;
	g =  -0.969258*X	+ 1.87599 *Y	+ 0.0415557*Z;
	b =  0.0556352*X    + -0.203996*Y	+ 1.05707*Z;

	if (r <= 0.0031308)
		rgb[0] = 12.92*r;
	else
		rgb[0] = 1.055*(float)pow(r,float(1.0)/float(2.4))-0.055;

	if (g <= 0.0031308)
		rgb[1] = 12.92*g;
	else
		rgb[1] = 1.055*(float)pow(g,float(1.0)/float(2.4))-0.055;

	if (b <= 0.0031308)
		rgb[2] = 12.92*b;
	else
		rgb[2] = 1.055*(float)pow(b,(float(1.0)/float(2.4)))-0.055;

	//rgb[0] = R;
	//rgb[1] = G;
	//rgb[2] = B;



}
// Gabriel

/**
 \fn	void RGB_To_LUV(float r, float g, float b, float* luv)

 \brief	Rgb to luv.


 \param	r			   	The.
 \param	g			   	The.
 \param	b			   	The.
 \param [in,out]	luv	If non-null, the luv.
 */

void RGB_To_LUV(float r, float g, float b, float* luv)
{
	// R 0..1, G 0..1, B 0..1

	//http://www.brucelindbloom.com

	/*float rf, gf, bf;
	float X_, Y_, Z_, X, Y, Z, fx, fy, fz, xr, yr, zr;*/
	float X, Y, Z, yr;
	float L;
	float eps = 216.f/24389.f;
	float k = 24389.f/27.f;

	float Xr = 0.964221f;  // reference white D50
	float Yr = 1.0f;
	float Zr = 0.825211f;


	// RGB to XYZ

	// assuming sRGB (D65)
	if (r <= 0.04045)
		r = r / 12.92;
	else
		r = (float)pow((r+0.055)/1.055,2.4);

	if (g <= 0.04045)
		g = g/12.92;
	else
		g = (float)pow((g+0.055)/1.055,2.4);

	if (b <= 0.04045)
		b = b/12.92;
	else
		b = (float)pow((b+0.055)/1.055,2.4);


	X =  0.412424*r     + 0.357579*g	+ 0.180464*b;
	Y =  0.212656*r     + 0.715158 *g	+ 0.0721856*b;
	Z =  0.0193324*r    + 0.119193*g	+ 0.950444*b;

	// XYZ to Luv

	float u, v, u_, v_, ur_, vr_;

	u_ = 4*X / (X + 15*Y + 3*Z);
	v_ = 9*Y / (X + 15*Y + 3*Z);

	ur_ = 4*Xr / (Xr + 15*Yr + 3*Zr);
	vr_ = 9*Yr / (Xr + 15*Yr + 3*Zr);

	yr = Y/Yr;

	if ( yr > eps )
		L = 116.f * (float)pow(yr, (float)(1.f/3.f)) - 16.f;
	else
		L = k * yr;

	u = 13.0*L*(u_ - ur_);
	v = 13.0*L*(v_ - vr_);

	luv[0] = (float)L;
	luv[1] = (float)u;
	luv[2] = (float)v;

}




/*
	Conversion from RGB to YUV and from YUV to RGB, implemented by Gabriel
	This code contains some errors, I have not found and corrected these errors...
*/
/*
// Conversion de LUV vers RGB
void LUV_To_RGB(float l, float u, float v, float* rgb)
{
	//http://www.brucelindbloom.com

	float eps = 216.f/24389.f;
	float k = 24389.f/27.f;

	float Xr = 0.964221f;  // reference white D50
	float Yr = 1.0f;
	float Zr = 0.825211f;
	//float Xr = 96.4221f;  // reference white D50
	//float Yr = 100.f;
	//float Zr = 82.5211f;

	// Luv to XYZ

	float X, Y, Z, ur_, vr_;

	ur_ = 4*Xr / (Xr + 15*Yr + 3*Zr);
	vr_ = 9*Yr / (Xr + 15*Yr + 3*Zr);

	if ( l > k*eps )
	{
		Y = (float)pow((float)((l+16.0)/116.f), (float)3.0);
	}
	else
	{
		Y = l/k;
	}

	float Cd = Y*( (39*l)/(v+13*l*vr_) - 5);
	float Ca = ((52*l)/(u+13*l*ur_) - 1)/3.;
	float Cb = -5*Y;

	X = (Cd-Cb)/(Ca+(1./3.));
	Z = X*Ca + Cb;

	// XYZ to RGB

	// R 0..1, G 0..1, B 0..1

	float R, G, B, r, g, b;

	// assuming sRGB (D65)

	r =  3.24071*X      + -1.53726*Y	+ -0.498571*Z;
	g =  -0.969258*X	+ 1.87599 *Y	+ 0.0415557*Z;
	b =  0.0556352*X    + -0.203996*Y	+ 1.05707*Z;

	if (r <= 0.0031308)
		R = 12.92*r;
	else
		R = 1.055*(float)pow(r,float(1.0)/float(2.4))-0.055;

	if (g <= 0.0031308)
		G = 12.92*g;
	else
		G = 1.055*(float)pow(g,float(1.0)/float(2.4))-0.055;

	if (b <= 0.0031308)
		B = 12.92*b;
	else
		B = 1.055*(float)pow(b,(float(1.0)/float(2.4)))-0.055;

	rgb[0] = R;
	rgb[1] = G;
	rgb[2] = B;
}
*/

int Determine_Color_Configuration(const Halfedge_handle &h)
{
	/*Halfedge_handle g = h->next();	
	vector<Color_Unit> Neighbor_color;
	vector<float*> Neighbor_lab;
	vector<float*> Neighbor_rgb;
	Color_Unit VC = Get_Vertex_Color(g);
	size_t valence = g->vertex_degree();
	float Vert_pos[3];
	float Vert_col_lab[3];
	float Vert_col_rgb[3];

	Vert_pos[0] = g->vertex()->point().x();
	Vert_pos[1] = g->vertex()->point().y();
	Vert_pos[2] = g->vertex()->point().z();
	Vert_col_lab[0] = g->vertex()->color_float(0);
	Vert_col_lab[1] = g->vertex()->color_float(1);
	Vert_col_lab[2] = g->vertex()->color_float(2);

	LAB_To_RGB(Vert_col_lab[0], Vert_col_lab[1], Vert_col_lab[2], Vert_col_rgb);
	FILE * Local_mesh = fopen("local_mesh.off","w");
	fprintf(Local_mesh,"COFF\n");
	fprintf(Local_mesh,"%d %d 0\n",valence + 1, valence);
	fprintf(Local_mesh,"%f %f %f %f %f %f\n", Vert_pos[0], Vert_pos[1], Vert_pos[2], Vert_col_rgb[0], Vert_col_rgb[1], Vert_col_rgb[2]);
	
	for (int i = 0; i < (int)valence; i++)
	{
		g = g->opposite();
		Color_Unit col = Get_Vertex_Color(g);		
		Neighbor_color.push_back(col);
		
		Vert_pos[0] = g->vertex()->point().x();
		Vert_pos[1] = g->vertex()->point().y();
		Vert_pos[2] = g->vertex()->point().z();
		Vert_col_lab[0] = g->vertex()->color_float(0);
		Vert_col_lab[1] = g->vertex()->color_float(1);
		Vert_col_lab[2] = g->vertex()->color_float(2);
		Neighbor_lab.push_back(Vert_col_lab);
		LAB_To_RGB(Vert_col_lab[0], Vert_col_lab[1], Vert_col_lab[2], Vert_col_rgb);
		Neighbor_rgb.push_back(Vert_col_rgb);
		fprintf(Local_mesh,"%f %f %f %f %f %f\n", Vert_pos[0], Vert_pos[1], Vert_pos[2], Vert_col_rgb[0], Vert_col_rgb[1], Vert_col_rgb[2]);

		g = g->prev();
	}

	for (int i = 0; i < (int)valence - 1; i++)
		fprintf(Local_mesh,"3 0 %d %d\n", i+1, i+2);

	fprintf(Local_mesh,"3 0 %d 1\n", valence);
	
	fclose(Local_mesh);*/
	return 0;
}

double Color_Distance(const float * LAB1, const float * LAB2)
{
	float LCH1[3];
	float LCH2[3];

	LAB_To_LCH(LAB1, LCH1);
	LAB_To_LCH(LAB2, LCH2);

	double dist;
	//dist = Color_Distance_ICE94(LCH1, LCH2);
	//dist = Color_Distance_CMC21(LCH1, LCH2);
	dist = sqrt(pow((double)LAB2[0]-(double)LAB1[0], 2.0) + pow((double)LAB2[1] - (double)LAB1[1], 2.0) + pow((double)LAB2[2] - (double)LAB1[2], 2.0));
	return dist;
}

double Check_Color_Importance(const Halfedge_handle & ah)//, const unsigned & valence)
{
	return 0;
	/*float LAB1[3];
	float LAB2[3];

	Halfedge_handle h = ah;
	h = h->next();

	LAB1[0] = h->vertex()->color_float(0);
	LAB1[1] = h->vertex()->color_float(1);
	LAB1[2] = h->vertex()->color_float(2);

	LAB2[0] = h->opposite()->vertex()->color_float(0);
	LAB2[1] = h->opposite()->vertex()->color_float(1);
	LAB2[2] = h->opposite()->vertex()->color_float(2);

	double dist = Color_Distance(LAB1, LAB2);

	return dist;*/
}

#endif

#endif
