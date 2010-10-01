///////////////////////////////////////////////////////////////////////////
// Author: Ho Lee
// Year: 2009
// LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#include <mepp_config.h>
#ifdef BUILD_component_Canonical

#ifdef _MSC_VER
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif

#include "Canonical_Component.h"

#include <cmath>


// tags for vertices and facets.
#define FREE -1
#define CONQUERED 0
#define TO_BE_REMOVED 1

// tags for retriangulation
#define PLUS 1
#define MINUS -1
#define NOSIGN 0

const double PI = 3.14159265359;



// Description : Initialize all flags -verticeces and facets- to FREE and give order to vertices
// This function is called within every conquest.
void Init(PolyhedronPtr pMesh)
{
	int i = 0;

	// vertices flags initialization
	Vertex_iterator pVertex = NULL;
	for(pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); i++,pVertex++)
	{
		pVertex->Vertex_Flag_S = FREE;
		pVertex->Vertex_Number_S = i;
		pVertex->Vertex_Sign_S = NOSIGN;
	}

	// facets flag initialization.
	Facet_iterator pFace = NULL;
	for(pFace = pMesh->facets_begin(); pFace != pMesh->facets_end(); pFace++)
	{
		pFace->Facet_Flag_S = FREE;
	}
}


//Description :: To calculate area of triangle
double Area_Facet_Triangle(const Point3d &P,const Point3d &Q, const Point3d &R)
{
	Vector PQ=Q-P;
	Vector PR=R-P;
	Vector QR=R-Q;


	Vector normal	=	CGAL::cross_product(PQ,QR);
	double area=0.5*sqrt(normal*normal);

	return area;
}

//Description : Gives a normal vector of the triangle which contain the halfedge_handle h
Vector Triangle_Normal(const Halfedge_handle &h)
{
	Point3d P = h->vertex()->point();
	Point3d Q = h->next()->vertex()->point();
	Point3d R = h->next()->next()->vertex()->point();

	Vector PQ=Q-P;
	Vector PR=R-P;
	Vector QR=R-Q;

	Vector normal = CGAL::cross_product(PQ,QR);
	double length = std::sqrt(normal*normal);
	if (length != 0.0)
		normal = normal / length;

	return normal;
}
//Description : Gives a normal vector of the triangle formed by three points P Q R in the counterclockwise way.
Vector Triangle_Normal(const Point3d & P,const Point3d & Q,const Point3d &R)
{
	Vector PQ=Q-P;
	Vector PR=R-P;
	Vector QR=R-Q;

	Vector normal = CGAL::cross_product(PQ,QR);
	double length = std::sqrt(normal*normal);
	if (length != 0.0)
		normal = normal / length;

	return normal;
}

//Description : Gives an area of the triangle which contain the halfedge_handle h
double Area_Facet_Triangle(const Halfedge_handle &h)
{

	Point3d P = h->vertex()->point();
	Point3d Q = h->next()->vertex()->point();
	Point3d R = h->next()->next()->vertex()->point();

	Vector PQ=Q-P;
	Vector PR=R-P;
	Vector QR=R-Q;

	Vector normal =	CGAL::cross_product(PQ,QR);
	double area=0.5*sqrt(normal*normal);

	return area;
}

// Description : To find a correspondent type to retriangulate.
int Find_Type(const Halfedge_handle &h,const unsigned int &valence)
{
	int type = 0;
	if (valence == 3)
	{
		if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 1;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 2;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 3;
		else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 4;
	}

	else if (valence == 4)
	{
		if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 5;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 6;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 7;
		else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 8;
	}

	else if (valence == 5)
	{
		if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 9;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 10;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 11;
		else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 12;
	}

	else if (valence == 6)
	{
		if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 13;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 14;
		else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			type = 15;
		else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			type = 16;
	}
	return type;
}

//Description :: Check if removal of this vertex would violate the manifold_property or not.
bool Check_Manifold_Property(Halfedge_handle h, const int &type,const int &valence)
{
	bool check = false;
	Halfedge_handle g = h;
	int* Points_index = new int[valence];

	// if valence is 3, no new edge is inserted, so always safe to remove.
	if(valence == 3)
	{
		return false;
	}

	else
	{
		// Points_index[] contains all boundary vertices' indices (ordered in counterclockwise)

		Points_index[0] = g->vertex()->Vertex_Number_S;
		g = g->next(); // g points center vertex;

		for(int i=1; i<valence; i++)
		{
			g = g->prev_on_vertex();// around the vertex in the counterclockwise way.
			Points_index[i] = g->opposite()->vertex()->Vertex_Number_S;
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[1])
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[2])
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[1])
						check = true;
				}

				g = h->next()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[3])
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[3])
						check = true;
				}

				g = h->next()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[3])
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[2])
						check = true;
				}

				g = h->opposite();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[2])
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[1])
						check = true;
				}

				g = h->opposite();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[3])
						check = true;
				}

				g = h->next()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[3])
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
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[2])
						check = true;
				}

				g = h;
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[4])
						check = true;
				}

				g = h->prev()->opposite()->next();
				Hvc = g->vertex_begin();
				Hvc_end = Hvc;

				CGAL_For_all(Hvc,Hvc_end)
				{
					if (Hvc->opposite()->vertex()->Vertex_Number_S == Points_index[2])
						check = true;
				}
			}
		}
	}
	return check;
}


// Description :: Here, we define a geometric_metric to preserve maximum the intermediate meshes.
bool Check_Geometric_Metric(const Halfedge_handle &h,const int &type,const unsigned int &valence,const float & Threshold)
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
	for(i = 1; i < (valence+1); i++)
	{
		Points[i] = g->opposite()->vertex()->point();
		g = g->prev_on_vertex();
	}

	// caculate perimeter
	Vector* Vectors = new Vector[valence];
	Vectors[0] = Points[1] - Points[valence];
	perimeter = std::sqrt(Vectors[0] * Vectors[0]);

	for(i = 1; i< valence; i++)
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

	else if(valence == 4)
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

		else if((type == 6) || (type == 7))// [0 1 2 3], [0 1 3 4]
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

	else if(valence == 5)
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

	else if(valence == 6)
	{
		if((type == 13) || (type == 16)) //[0 1 2 6], [0 2 3 4] , [0 4 5 6], [0 2 4 6]
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

		else if((type == 14) || (type == 15))// [ 0 1 2 3] [ 0 3 4 5] [ 0 1 5 6] [ 0 1 3 5]
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
	delete Points;
	delete Vectors;

	return check;
}


// Descrpition :: To check if removal of the front vertex can cause a normal flipping phenomenon.
bool Check_Normal_Flipping(const Halfedge_handle &h,const unsigned &valence)
{
	int type = Find_Type(h,valence);
	bool check = false;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Halfedge_handle g = h;
	g = g->next();

    Vector V_normal = h->next()->vertex()->normal();

	double length = std::sqrt(V_normal*V_normal);
	if (length != 0)
		V_normal = V_normal / length;

	Point3d *Points = new Point3d[valence];
    Vector * Normal = new Vector[valence - 2];


	for(unsigned int i=0;i<valence;i++)
	{
		Points[i] = g->opposite()->vertex()->point();
		g = g->prev_on_vertex();
	}

	for(unsigned int j = 0;j<(valence - 2);j++)
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

	for(unsigned int i=0;i<(valence-2);i++)
	{
        double length_normal = std::sqrt(Normal[i]*Normal[i]);
        if (length_normal != 0)
            Normal[i] = Normal[i] / length_normal;
	}

    for(unsigned int i = 0 ; i < (valence - 2); i++)
    {
        double cosine = V_normal * Normal[i];
        double cosine_rad = std::acos(cosine);

        if (cosine_rad >= PI/2)
            check = true;
    }

	return check;
}

bool Check_Border_Vertex(const Halfedge_handle & h)
{
	bool check = false;
	Halfedge_around_vertex_circulator hvc = h->vertex_begin();
	Halfedge_around_vertex_circulator hvc_end = hvc;

	CGAL_For_all(hvc,hvc_end)
	{
		if (hvc->is_border_edge() == true)
			check = true;
	}

	return check;
}


// Description :: To retriangulate the hole left by a removal of a vertex.
void Retriangulation(PolyhedronPtr pMesh,Halfedge_handle h,unsigned int &valence)
{
	int type = Find_Type(h,valence);
	Halfedge_handle g;
	g = h->next();


	// Triangle
	if ((type == 1) || (type == 2) || (type == 4))
	{
		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		//normal = Triangle_Normal(h);

		h->facet()->normal() = Triangle_Normal(h);

		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = PLUS;
	}
	else if (type == 3)
	{
		h->facet()->Facet_Flag_S = TO_BE_REMOVED;

		h->facet()->normal() = Triangle_Normal(h);
		//normal = Triangle_Normal(h);
		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = MINUS;
	}

	// quadrangle
	else if ((type == 5) || (type == 8))
	{
		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = MINUS;

		h = h->prev();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		//normals[1] = Triangle_Normal(h->opposite());
		//area[1] = Area_Facet_Triangle(h->opposite());

		//normals[2] = Triangle_Normal(h);
		//area[2] = Area_Facet_Triangle(h);
		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

	}

	else if (( type == 6) || (type == 7))
	{
		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = PLUS;

		g = g->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		//normals[1] = Triangle_Normal(h->opposite());
		//area[1] = Area_Facet_Triangle(h->opposite());
		//normals[2] = Triangle_Normal(h);
		//area[2] = Area_Facet_Triangle(h);

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;


	}

	// pentagone
	else if ((type == 9) || (type == 12))
	{
		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign_S = PLUS;

		h = h->prev();

		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next()->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;
	}

	else if (type == 10)
	{
		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign_S = PLUS;

		g = h;
		h = h->prev()->prev();

		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next();
		h = h->prev();

		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);

		h = pMesh->split_facet(h,g);

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

	}

	else if (type == 11)
	{
		if (g->vertex()->Vertex_Sign_S == NOSIGN)
			g->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign_S = MINUS;

		g = g->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next()->next();

		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

	}

	// Hexagone

	else if ((type == 13) || (type == 16))
	{
		if (g->vertex()->Vertex_Sign_S ==NOSIGN)
			g->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->vertex()->Vertex_Sign_S ==NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->next()->vertex()->Vertex_Sign_S ==NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->next()->next()->vertex()->Vertex_Sign_S ==NOSIGN)
			g->next()->next()->next()->vertex()->Vertex_Sign_S = MINUS;

		h = h->prev();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;


		g = h->next()->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next()->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;


	}

	else if ((type == 14) || (type == 15))
	{
		if (g->vertex()->Vertex_Sign_S ==NOSIGN)
			g->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->vertex()->Vertex_Sign_S ==NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = PLUS;
		if (g->next()->next()->vertex()->Vertex_Sign_S ==NOSIGN)
			g->next()->next()->vertex()->Vertex_Sign_S = MINUS;
		if (g->next()->next()->next()->vertex()->Vertex_Sign_S ==NOSIGN)
			g->next()->next()->next()->vertex()->Vertex_Sign_S = PLUS;


		g = g->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next()->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;

		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		g = h->next()->next();
		CGAL_precondition( h->next() != g);
		CGAL_precondition( g->next() != h);
		h = pMesh->split_facet(h,g);

		h->facet()->normal() = Triangle_Normal(h);
		h->opposite()->facet()->normal() = Triangle_Normal(h->opposite());

		h->facet()->Facet_Flag_S = TO_BE_REMOVED;
		h->opposite()->facet()->Facet_Flag_S = TO_BE_REMOVED;
	}
}

// Description : Regulation conquest
void Canonical_Component::Regulation(PolyhedronPtr pMesh,const bool Normal_flipping,const bool Use_metric,
											const float metric_threshold,const bool Use_forget_metric,const unsigned int &Forget_value)
{
    Init(pMesh);


	Halfedge_iterator hi = pMesh->halfedges_begin();
	std::queue<Halfedge_handle> Halfedges2;

	hi->vertex()->Vertex_Flag_S = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag_S = CONQUERED;
	Halfedges2.push(hi);

	Halfedge_handle h;

	while(!Halfedges2.empty())
	{
		h = Halfedges2.front();
		Halfedges2.pop();

		unsigned valence;
		valence = h->next()->vertex_degree();

		if ((h->facet()->Facet_Flag_S == CONQUERED) || (h->facet()->Facet_Flag_S == TO_BE_REMOVED))
		{
			continue;
		}

		else if ((h->next()->vertex()->Vertex_Flag_S == FREE) && (valence == 3) && (Check_Border_Vertex(h->next()) == false)) // if valence is 3, remove the front vertex.
		{
			Halfedge_handle g = h;
			int type = 1; // ant type of valence 3

			bool Global_check = true;

			// Check if the manifold property is violated.
			bool Check_Manifold = Check_Manifold_Property(g,type,valence);

			g = h;

			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Check_Metric = true;

			if (Use_metric == true)
			{
				if (Use_forget_metric == true)
				{
					if (pMesh->size_of_vertices() > Forget_value)
					{
						Check_Metric = false;
					}
					else
					{
						Check_Metric = Check_Geometric_Metric(g, type, valence,metric_threshold);//this->Metric_threadhold);
					}
				}

				else
				{
					Check_Metric = Check_Geometric_Metric(g, type, valence,metric_threshold);//this->Metric_threadhold);
				}
			}
			else
			{
				Check_Metric = false;
			}

			Check_Metric = false;

			if ((Check_Manifold == false) && (Check_Metric == false))
				Global_check = false;
			else
				Global_check = true;

			if (Global_check == false)
			{
				Halfedge_handle pass = h;

				g = h->next();
				g->vertex()->Vertex_Flag_S = TO_BE_REMOVED;
				g->facet()->Facet_Flag_S = TO_BE_REMOVED;

				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag_S = CONQUERED;
				g->facet()->Facet_Flag_S = TO_BE_REMOVED;
				if(!g->prev()->is_border_edge())
				{
					Halfedge_handle h1 = g->prev()->opposite();
					h1->facet()->Facet_Flag_S = CONQUERED;
					h1->next()->vertex()->Vertex_Flag_S = CONQUERED;
					if(!h1->next()->is_border_edge())
						Halfedges2.push(h1->next()->opposite());
					if(!h1->prev()->is_border_edge())
						Halfedges2.push(h1->prev()->opposite());

				}
				g = g->prev_on_vertex();
				g->opposite()->vertex()->Vertex_Flag_S = CONQUERED;
				g->facet()->Facet_Flag_S = TO_BE_REMOVED;
				if(!g->prev()->is_border_edge())
				{
					Halfedge_handle h2 = g->prev()->opposite();
					h2->facet()->Facet_Flag_S = CONQUERED;
					h2->next()->vertex()->Vertex_Flag_S = CONQUERED;
					if(!h2->next()->is_border_edge())
						Halfedges2.push(h2->next()->opposite());
					if(!h2->prev()->is_border_edge())
						Halfedges2.push(h2->prev()->opposite());
				}
			}

			else
			{
				h->facet()->Facet_Flag_S = CONQUERED;
				h->next()->vertex()->Vertex_Flag_S = CONQUERED;
				if(!h->next()->is_border_edge())
					Halfedges2.push(h->next()->opposite());
				if(!h->prev()->is_border_edge())
					Halfedges2.push(h->prev()->opposite());
			}
		}
		else  // NULL triangle
		{
			//NULL PATCH


			h->facet()->Facet_Flag_S = CONQUERED;
			h->next()->vertex()->Vertex_Flag_S = CONQUERED;
			if(!h->next()->is_border_edge())
					Halfedges2.push(h->next()->opposite());
				if(!h->prev()->is_border_edge())
					Halfedges2.push(h->prev()->opposite());
		}
	}

	Vertex_iterator pVertex = NULL;
	for(pVertex = pMesh->vertices_begin(); pVertex != pMesh->vertices_end(); )
	{
		Vertex_handle vh = pVertex;
		pVertex++;
		Halfedge_handle del = vh->halfedge();

		if (vh->Vertex_Flag_S == TO_BE_REMOVED)
		{
			pMesh->erase_center_vertex(del);
		}

	}

}

// Description : This function select a set of independent vertices to be removed
// Method used by Alliez in "Progressive Compression for Lossless Transmission of Triangle Meshes"
void Canonical_Component::Decimation_Conquest(PolyhedronPtr pMesh,const bool Normal_flipping,const bool Use_metric,
											const float metric_threshold,const bool Use_forget_metric,const unsigned int &Forget_value)
{
	Init(pMesh);
	Halfedge_iterator hi = pMesh->halfedges_begin();
	std::queue<Halfedge_handle> Halfedges;


	// Two vertices of seed edges are flaged CONQUERED
	hi->vertex()->Vertex_Flag_S = CONQUERED;
	hi->opposite()->vertex()->Vertex_Flag_S = CONQUERED;

	// These vertices are also flaged with sign flags for retriangulation
	hi->vertex()->Vertex_Sign_S = PLUS;
	hi->opposite()->vertex()->Vertex_Sign_S = MINUS;

	//push the first halfedge in the queue.
	Halfedges.push(&(*(hi)));

	Halfedge_handle h;

	while(!Halfedges.empty())
	{
		h = Halfedges.front();
		Halfedges.pop();

		unsigned type = 0; // define type of retriangulation

		unsigned valence = h->next()->vertex_degree(); // valence

		if (h->vertex()->point() == h->opposite()->vertex()->point())
		{
			continue;
		}

		// if its front face is not tagged CONQUERED nor TO_BE_REMOVED
		else if ((h->facet()->Facet_Flag_S == CONQUERED) || (h->facet()->Facet_Flag_S == TO_BE_REMOVED))
		{
			continue;
		}

		// if its front vertex is free and has a valence <= 6
		else if ((h->next()->vertex()->Vertex_Flag_S == FREE) && (valence <= 6) && (Check_Border_Vertex(h->next()) == false))
		{
			type = Find_Type(h,valence);

			// remove the front vertex if its removal does not viloate the manifold property
			Halfedge_handle g = h;

			bool Global_check = true;

			// Check if the manifold property is violated.
			bool Check_Manifold = Check_Manifold_Property(g,type,valence);

			g = h;
            bool Check_NF = false;

            if (Normal_flipping == true)
            {
                Check_NF = Check_Normal_Flipping(g,valence);
            }
            else
            {
                Check_NF = false;
            }

			g = h;
			// calculate error caused by the removal. This metric decides if the vertex can be removed or not.
			bool Check_Metric = false;
			if (Use_metric == true)
			{
				if (Use_forget_metric == true)
				{
					if (pMesh->size_of_vertices() > Forget_value)
					{
						Check_Metric = false;
					}
					else
					{
						Check_Metric = Check_Geometric_Metric(g, type, valence,metric_threshold);//this->Metric_threadhold);
					}
				}

				else
				{
					Check_Metric = Check_Geometric_Metric(g, type, valence,metric_threshold);//this->Metric_threadhold);
				}
			}
			else
			{
				Check_Metric = false;
			}

			if ((Check_Manifold == false) && (Check_Metric == false) && (Check_NF == false))
				Global_check = false;
			else
				Global_check = true;



			if (Global_check == true) // violated manifold property or metric -> NULL patch
			{


				h->facet()->Facet_Flag_S = CONQUERED;
				h->next()->vertex()->Vertex_Flag_S = CONQUERED;

				if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
						h->next()->vertex()->Vertex_Sign_S = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
						h->next()->vertex()->Vertex_Sign_S = PLUS;
				}
				else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
				{
					if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
						h->next()->vertex()->Vertex_Sign_S = MINUS;
				}
				else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
				{
					if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
						h->next()->vertex()->Vertex_Sign_S = PLUS;
				}
				if(h->next()->is_border_edge() == false)
					Halfedges.push(h->next()->opposite());
				if(h->prev()->is_border_edge() == false)
					Halfedges.push(h->prev()->opposite());
			}

			else if(Global_check == false) // All conditions are good. -> Remove the center vertex.
			{
				g = h;
				g = pMesh->erase_center_vertex(g->next()); // remove the front vertex
				g = h;
				g->facet()->Facet_Flag_S = CONQUERED;

				for(unsigned int j=0; j<(valence-1);j++)
				{
					g = g->next();
					g->vertex()->Vertex_Flag_S = CONQUERED;

					if(g->is_border_edge() == false)
					{
						g->opposite()->vertex()->Vertex_Flag_S = CONQUERED;
						Halfedges.push(g->opposite());
					}
				}

				Halfedge_handle pass = h;
				Retriangulation(pMesh, pass, valence);
			}
		}

		// if its front face is NULL PATCH
		else if ((h->next()->vertex()->Vertex_Flag_S == CONQUERED) || ((h->next()->vertex()->Vertex_Flag_S == FREE) && (valence > 6)))
		{
			// its front face is tagged CONQUERED
			h->facet()->Facet_Flag_S = CONQUERED;
			h->next()->vertex()->Vertex_Flag_S = CONQUERED;///////////////////////////

			if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
					h->next()->vertex()->Vertex_Sign_S = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
					h->next()->vertex()->Vertex_Sign_S = PLUS;
			}
			else if ((h->vertex()->Vertex_Sign_S == PLUS) && (h->opposite()->vertex()->Vertex_Sign_S == PLUS))
			{
				if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
					h->next()->vertex()->Vertex_Sign_S = MINUS;
			}
			else if ((h->vertex()->Vertex_Sign_S == MINUS) && (h->opposite()->vertex()->Vertex_Sign_S == MINUS))
			{
				if (h->next()->vertex()->Vertex_Sign_S == NOSIGN)
					h->next()->vertex()->Vertex_Sign_S = PLUS;
			}

			// two other output gates are pused into the fifo queue.
			if(h->next()->is_border_edge() == false)
				Halfedges.push(h->next()->opposite());
			if(h->prev()->is_border_edge() == false)
				Halfedges.push(h->prev()->opposite());
		}
	}
}


// Description : Remove edges to create a hole.
bool Remove_Edges(PolyhedronPtr pMesh,const Halfedge_handle &h,const int &type)
{
	bool check = false;
	Halfedge_handle g = h;

	//triangle

	if ((type == 1) || (type == 2) || (type == 4))
	{
		if(g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = PLUS;
	}

	else if (type == 3)
	{
		if(g->next()->vertex()->Vertex_Sign_S == NOSIGN)
			g->next()->vertex()->Vertex_Sign_S = MINUS;
	}


	// quadrangle
	else if ((type == 5) || (type == 8))
	{
		//verification
		if(g->prev()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;
		if (check == false)
		{
			g = g->prev();
			pMesh->join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = MINUS;
		}
	}

	else if ((type == 6) || (type == 7))
	{
		//verification
		if(g->next()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;
		if (check == false)
		{
			g = g->next();
			pMesh->join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = PLUS;
		}
	}

	//pentagone
	else if ((type == 9) || (type == 12))
	{
		g = g->prev()->opposite();
		if (g->facet()->Facet_Flag_S != FREE)
			check = true;
		g = g->next()->opposite();
		if (g->facet()->Facet_Flag_S != FREE)
			check = true;

		if (check == false)
		{
			g = h->prev();
			g = pMesh->join_facet(g);
			g = g->next();
			g = pMesh->join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign_S = PLUS;
		}
	}

	else if (type == 10)
	{
		g = g->next()->opposite();
		if (g->facet()->Facet_Flag_S != FREE)
			check = true;
		g = g->prev()->opposite();
		if (g->facet()->Facet_Flag_S != FREE)
			check = true;

		if (check == false)
		{
			g = h->next()->opposite();
			g = pMesh->join_facet(g);
			g = pMesh->join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign_S = PLUS;
		}

	}

	else if (type == 11)
	{
		if(g->next()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;
		if(g->prev()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;

		if (check == false)
		{
			g = g->next();
			g = pMesh->join_facet(g);
			g = g->prev();
			g = pMesh->join_facet(g);

			g = h;

			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign_S = MINUS;
		}
	}

	else if ((type == 13) || (type == 16))
	{
		g = g->prev()->opposite();
		if (g->facet()->Facet_Flag_S != FREE)
			check = true;
		if (g->next()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;
		if (g->prev()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;

		if (check == false)
		{
			g = h->prev()->opposite();
			g = pMesh->join_facet(g);
			g = pMesh->join_facet(g);
			g = pMesh->join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->next()->vertex()->Vertex_Sign_S = MINUS;
		}
	}

	else if ((type == 14) || (type == 15))
	{
		g = g->next()->opposite();
		if(g->facet()->Facet_Flag_S != FREE)
			check = true;
		if (g->next()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;
		if (g->prev()->opposite()->facet()->Facet_Flag_S != FREE)
			check = true;

		if(check == false)
		{
			g = h->next()->opposite();
			g = pMesh->join_facet(g);
			g = pMesh->join_facet(g);
			g = pMesh->join_facet(g);

			g = h;
			if (g->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->vertex()->Vertex_Sign_S = PLUS;
			if (g->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->vertex()->Vertex_Sign_S = MINUS;
			if (g->next()->next()->next()->next()->vertex()->Vertex_Sign_S == NOSIGN)
				g->next()->next()->next()->next()->vertex()->Vertex_Sign_S = PLUS;
		}
	}

	return check;
}

Canonical_Component::Canonical_Component(Viewer* v, PolyhedronPtr p):mepp_component(v, p)
{
	// MEPP 2
	componentName = "Canonical_Component";
	init = 1;
}

#endif
