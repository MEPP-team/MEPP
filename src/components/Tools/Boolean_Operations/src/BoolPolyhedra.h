#ifndef BOOLPOLYHEDRA_H
#define BOOLPOLYHEDRA_H

#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include "Boolean_Operations_Definitions.h"
#include "CPolyhedron_from_polygon_builder_3.h"
#include "Boolean_Operations_triangulation.h"

#ifdef BOOLEAN_OPERATIONS_DEBUG
#include "Time_measure.h"
#endif // BOOLEAN_OPERATIONS_DEBUG

typedef Polyhedron::HalfedgeDS																	HDS;
typedef CGAL::Simple_cartesian<num_type>														AABB_Kernel;

class Enriched_Triangle : public AABB_Kernel::Triangle_3
{
public:

	typedef AABB_Kernel::Point_3				Point_3;
	typedef AABB_Kernel::FT						FT;

	Enriched_Triangle(Facet_handle _f)
		: AABB_Kernel::Triangle_3(	to_K(_f->facet_begin()->vertex()->point() + (_f->facet_begin()->vertex()->point() - _f->facet_begin()->next()->vertex()->point()) / 1000),
									to_K(_f->facet_begin()->next()->vertex()->point() + (_f->facet_begin()->next()->vertex()->point() - _f->facet_begin()->next()->next()->vertex()->point()) / 1000),
									to_K(_f->facet_begin()->next()->next()->vertex()->point() + (_f->facet_begin()->next()->next()->vertex()->point() - _f->facet_begin()->vertex()->point()) / 1000)
									), f(_f) {}
	
	Facet_handle facet() {return f;}
	inline Point_3 to_K(Enriched_kernel::Point_3 p) {return Point_3((FT)p.x(), (FT)p.y(), (FT)p.z());}

private:
	Facet_handle f;
};

typedef Enriched_Triangle															Triangle;
typedef CGAL::AABB_triangle_primitive<AABB_Kernel,std::list<Triangle>::iterator>	AABB_Primitive;
typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive>								AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits>												AABB_Tree;

class BoolPolyhedra {

private:
	struct Triangle_Cut {								//contenu d'un triangle a decouper
		bool Facet_from_A;								//permet de savoir a quel polyedre appartient ce triangle
		Vector_exact norm_dir;							//direction de la normale
		std::vector<std::vector<InterId>> CutList;		//liste des decoupes a faire (aretes contraintes lors de la triangulation)
		std::vector<InterId> PtList;					//liste des decoupes a faire (points contraints lors de la triangulation)
		std::map<HalfedgeId, InterId> RefInter;			//liste des intersections (classées par halfedges)

		Triangle_Cut() {}
		Triangle_Cut(Vector_exact V, bool ffA) : norm_dir(V), Facet_from_A(ffA) {} 
		~Triangle_Cut() {}
	};
	
public:
	BoolPolyhedra(PolyhedronPtr pMA, PolyhedronPtr pMB, PolyhedronPtr pMout, Bool_Op BOOP) : m_BOOP(BOOP)
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		std::ofstream ofstrMA("input_A_boolsum.off"); ofstrMA << *pMA;
		std::ofstream ofstrMB("input_B_boolsum.off"); ofstrMB << *pMB;
		Time_measure Timer_total, Timer;	
		N_IFA = 0;
		N_IFB = 0;
		N_FFA = 0;
		N_FFB = 0;
		N_LFFA = 0;
		N_LFFB = 0;
		
		d_i_recup_couple = 0;
		d_i_tri_tri = 0;
		d_i_tt_init = 0;
		d_i_tt_verif_particular_case = 0;
		d_i_tt_inter_seg_tri = 0;
		d_i_tt_add2sol = 0;
		d_Cut_init = 0;
		d_Cut_add_pts = 0;
		d_Cut_add_segments = 0;
		d_Cut_triangulate = 0;
		d_Cut_add_tri = 0;
		d_Cut_validation = 0;

		Timer_total.Start();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		Init(pMA, pMB);

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Init = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		FindCouples();

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_FindCouples = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		if(m_Couples.size() != 0)
		{
			ComputeIntersections();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_ComputeIntersections = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			CutIntersectedFacets();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_CutIntersectedFacets = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			PropagateFacets();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_PropagateFacets = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			ppbuilder.clear_unused_vertices();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_clear_unused_vertices = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			pMout->delegate(ppbuilder);

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_delegate = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			std::ofstream ofstrMS("output_boolsum.off"); ofstrMS << *pMout;

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_write_off = Timer.GetDiff();
			duration_total = Timer_total.GetDiff();
			ofstrtime.open("output_time.txt");
			WriteData(pMout);
			ColorType();
#endif // BOOLEAN_OPERATIONS_DEBUG

		}
	}
	~BoolPolyhedra() {}
	
private:

	void Init(PolyhedronPtr pMA, PolyhedronPtr pMB)
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		Time_measure Timer;
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG
		
		//copier les 2 polyedres source dans pA et pB
		//m_pA = PolyhedronPtr(new Polyhedron());
		//m_pB = PolyhedronPtr(new Polyhedron());
		//*m_pA = *pMA;
		//*m_pB = *pMB;
		m_pA = pMA;
		m_pB = pMB;
		
#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Copy_Polyhedra = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		//triangulation des polyedres afin de n'avoir que des intersections de type triangle-triangle a traiter
		//cette etape est de toute facon necessaire au bon fonctionnement des AABB Trees.
		if(!m_pA->is_pure_triangle()) m_pA->triangulate();
		if(!m_pB->is_pure_triangle()) m_pB->triangulate();

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_triangulate = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		//initialisation des tags sur les vertex :
		for(Vertex_iterator pVertex = m_pA->vertices_begin();pVertex != m_pA->vertices_end();++pVertex)
		{
			pVertex->Label = 0xFFFFFFFF;
		}

		for(Vertex_iterator pVertex = m_pB->vertices_begin();pVertex != m_pB->vertices_end();++pVertex)
		{
			pVertex->Label = 0xFFFFFFFF;
		}

		//initialisation des tags sur les facettes et les halfedges
		for(Facet_iterator pFacet = m_pA->facets_begin();pFacet != m_pA->facets_end();++pFacet)
		{
			pFacet->Label = 0xFFFFFFFF;
			pFacet->IsExt = false;
			pFacet->IsOK = false;
		}

		for(Facet_iterator pFacet = m_pB->facets_begin();pFacet != m_pB->facets_end();++pFacet)
		{
			pFacet->Label = 0xFFFFFFFF;
			pFacet->IsExt = false;
			pFacet->IsOK = false;
		}

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Init_tags = Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

	}


	void FindCouples()
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		Time_measure Timer;
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		//on cherche les facettes de B intersectées par toutes les facettes de A
		//on procede différemment selon la taille des polyedres : l'arbre sera
		//construit sur le polyedre qui aura le moins de facettes
		Facet_iterator pFacet =	NULL;
		std::list<AABB_Tree::Primitive_id> primitives;
		std::list<Triangle> triangles;

		HalfedgeId i = 0;
		FacetId j = 0;

		if(m_pA->size_of_facets() < m_pB->size_of_facets())
		{
			//construction AABB Tree sur le plus petit Polyedre (A)
			for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++) triangles.push_back(Triangle(pFacet));
			tree.rebuild(triangles.begin(),triangles.end());

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_Build_AABBTree = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			//test de chaque facette de B avec A
			for (pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
			{
				tree.all_intersected_primitives(Triangle(pFacet), std::back_inserter(primitives));
				if(primitives.size() !=0)
				{
					Facet_Handle.push_back(pFacet);
					pFacet->Label = j++;
					pFacet->facet_begin()->Label = i++;
					pFacet->facet_begin()->next()->Label = i++;
					pFacet->facet_begin()->next()->next()->Label = i++;
					Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(pFacet->facet_begin()), false));
					do {
						if(primitives.back()->facet()->Label == 0xFFFFFFFF)
						{
							Facet_Handle.push_back(primitives.back()->facet());
							primitives.back()->facet()->Label = j++;
							primitives.back()->facet()->facet_begin()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->next()->Label = i++;
							Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), true));
						}
						m_Couples.push_back(pair<FacetId, FacetId>(primitives.back()->facet()->Label, pFacet->Label));
						primitives.pop_back();
					}
					while(primitives.size() != 0);
				}
			}

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_Test_AABBTree = Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

		}
		else
		{
			//construction AABB Tree sur le plus petit Polyedre (B)
			for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++) triangles.push_back(Triangle(pFacet));
			tree.rebuild(triangles.begin(),triangles.end());

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_Build_AABBTree = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			//test de chaque facette de A avec B
			for (pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
			{
				tree.all_intersected_primitives(Triangle(pFacet), std::back_inserter(primitives));
				if(primitives.size() !=0)
				{
					Facet_Handle.push_back(pFacet);
					pFacet->Label = j++;
					pFacet->facet_begin()->Label = i++;
					pFacet->facet_begin()->next()->Label = i++;
					pFacet->facet_begin()->next()->next()->Label = i++;
					Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(pFacet->facet_begin()), true));
					do {
						if(primitives.back()->facet()->Label == 0xFFFFFFFF)
						{
							Facet_Handle.push_back(primitives.back()->facet());
							primitives.back()->facet()->Label = j++;
							primitives.back()->facet()->facet_begin()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->next()->Label = i++;
							Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), false));
						}
						m_Couples.push_back(pair<FacetId, FacetId>(pFacet->Label, primitives.back()->facet()->Label));
						primitives.pop_back();
					}
					while(primitives.size() != 0);
				}
			}

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_Test_AABBTree = Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

		}
	}

	void ComputeIntersections()
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		Time_measure Timer;
#endif // BOOLEAN_OPERATIONS_DEBUG

		FacetId Ind_A, Ind_B;
		while(!m_Couples.empty())
		{

#ifdef BOOLEAN_OPERATIONS_DEBUG
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			Ind_A = m_Couples.back().first;
			Ind_B = m_Couples.back().second;
			m_Couples.pop_back();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			d_i_recup_couple += Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			InterTriangleTriangle(Ind_A, Ind_B);

#ifdef BOOLEAN_OPERATIONS_DEBUG
			d_i_tri_tri += Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

		}
	}

	void CutIntersectedFacets()
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		Time_measure Timer;
#endif // BOOLEAN_OPERATIONS_DEBUG

		Triangle_Cut TriCut;
		Halfedge_handle he;

		for(FacetId Facet = 0 ; Facet != Inter_tri.size() ; ++Facet)
		{
			if(!Inter_tri[Facet].CutList.empty())
			{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

				TriCut = Inter_tri[Facet];
				he = Facet_Handle[Facet]->facet_begin();
				Triangulation<Exact_Kernel> T(he, TriCut.norm_dir);
				bool IsExt[3];

#ifdef BOOLEAN_OPERATIONS_DEBUG
		d_Cut_init += Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

				for(int i = 0;i!=TriCut.PtList.size();++i)
				{
					T.add_new_pt(InterPts[TriCut.PtList[i]], TriCut.PtList[i]);
				}

#ifdef BOOLEAN_OPERATIONS_DEBUG
		d_Cut_add_pts += Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

				for(int i = 0;i!=TriCut.CutList.size();++i)
				{
					T.add_segment(InterPts[TriCut.CutList[i][0]], InterPts[TriCut.CutList[i][1]], TriCut.CutList[i][0], TriCut.CutList[i][1]);
				}

#ifdef BOOLEAN_OPERATIONS_DEBUG
		d_Cut_add_segments += Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

				vector<vector<unsigned long>> Tri_set = T.get_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false, IsExt);

#ifdef BOOLEAN_OPERATIONS_DEBUG
		d_Cut_triangulate += Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

				ppbuilder.add_triangle(Tri_set, he);

#ifdef BOOLEAN_OPERATIONS_DEBUG
		d_Cut_add_tri += Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

				Facet_Handle[Facet]->IsOK = true;
				if(IsExt[0]) he->opposite()->facet()->IsExt = true;
				if(IsExt[1]) he->next()->opposite()->facet()->IsExt = true;
				if(IsExt[2]) he->next()->next()->opposite()->facet()->IsExt = true;

#ifdef BOOLEAN_OPERATIONS_DEBUG
		d_Cut_validation += Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

			}
		}
	}

	void PropagateFacets()
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		Time_measure Timer;
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		Facet_handle pFacet = NULL, f = NULL, nf = NULL;
		stack<Facet_handle> tmpTriangles;

		for (pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
		{
			if(pFacet->IsOK) tmpTriangles.push(pFacet);
		}

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Prop_InitA = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		while(!tmpTriangles.empty())
		{
			f = tmpTriangles.top();
			tmpTriangles.pop();
			nf = f->facet_begin()->opposite()->facet();
			if(!nf->IsOK)
			{
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt == true) add_facet_to_solution(nf, true);
			}
			nf = f->facet_begin()->next()->opposite()->facet();
			if(!nf->IsOK)
			{
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt == true) add_facet_to_solution(nf, true);
			}
			nf = f->facet_begin()->next()->next()->opposite()->facet();
			if(!nf->IsOK)
			{
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt == true) add_facet_to_solution(nf, true);
			}
		}

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Prop_A = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		for (pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
		{
			if(pFacet->IsOK) tmpTriangles.push(pFacet);
		}

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Prop_InitB = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		while(!tmpTriangles.empty())
		{
			f = tmpTriangles.top();
			tmpTriangles.pop();
			nf = f->facet_begin()->opposite()->facet();
			if(!nf->IsOK)
			{
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt == true) add_facet_to_solution(nf, false);
			}
			nf = f->facet_begin()->next()->opposite()->facet();
			if(!nf->IsOK)
			{
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt == true) add_facet_to_solution(nf, false);
			}
			nf = f->facet_begin()->next()->next()->opposite()->facet();
			if(!nf->IsOK)
			{
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt == true) add_facet_to_solution(nf, false);
			}
		}
		
#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Prop_B = Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

	}
	
	void InterTriangleTriangle(FacetId A, FacetId B)
	{
#ifdef BOOLEAN_OPERATIONS_DEBUG
		Time_measure Timer;
#endif // BOOLEAN_OPERATIONS_DEBUG

		//on a une intersection entre les facettes A et B
		//calcul des normales
		Vector_exact nA, nB;
		nA = Inter_tri[A].norm_dir;
		nB = Inter_tri[B].norm_dir;

		if(CGAL::cross_product(nA, nB) != CGAL::NULL_VECTOR)
		{
			//intersection de type segment ou point

#ifdef BOOLEAN_OPERATIONS_DEBUG
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			//calcul des facets:
			Facet_handle fA, fB, fA2, fB2;
			fA = Facet_Handle[A];
			fB = Facet_Handle[B];

			//calcul des halfedges
			Halfedge_handle heA[3], heB[3];
			heA[0] = fA->facet_begin();
			heA[1] = heA[0]->next();
			heA[2] = heA[1]->next();
			heB[0] = fB->facet_begin();
			heB[1] = heB[0]->next();
			heB[2] = heB[1]->next();

			//calcul des points
			Point3d_exact ptA[3], ptB[3];
			ptA[0] = point_to_exact(heA[0]->vertex()->point());
			ptA[1] = point_to_exact(heA[1]->vertex()->point());
			ptA[2] = point_to_exact(heA[2]->vertex()->point());
			ptB[0] = point_to_exact(heB[0]->vertex()->point());
			ptB[1] = point_to_exact(heB[1]->vertex()->point());
			ptB[2] = point_to_exact(heB[2]->vertex()->point());

			//calcul des positions des points par rapport aux plans des triangles
			//positif : au dessus du plan
			//negatif : en dessous du plan
			//nul : sur le plan
			num_type posA[3], posB[3];
			posA[0] = nB * (ptA[0] - ptB[0]);
			posA[1] = nB * (ptA[1] - ptB[0]);
			posA[2] = nB * (ptA[2] - ptB[0]);
			posB[0] = nA * (ptB[0] - ptA[0]);
			posB[1] = nA * (ptB[1] - ptA[0]);
			posB[2] = nA * (ptB[2] - ptA[0]);

#ifdef BOOLEAN_OPERATIONS_DEBUG
			d_i_tt_init += Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			//calcul du code correspondant à la position du triangle par rapport au plan de l'autre
			//code sur 6 bits (2 bits pour chacun des 3 points)
			//10 -> point au dessus ; 01 -> point en dessous ; 00 -> point sur le plan (si au dessus ni en dessous)
			unsigned short posAbin, posBbin;
			posAbin =	  ( (posA[0] > 0)? 32 : 0 )
						+ ( (posA[0] < 0)? 16 : 0 )
						+ ( (posA[1] > 0)? 8 : 0 )
						+ ( (posA[1] < 0)? 4 : 0 )
						+ ( (posA[2] > 0)? 2 : 0 )
						+ ( (posA[2] < 0)? 1 : 0 );

			posBbin =	  ( (posB[0] > 0)? 32 : 0 )
						+ ( (posB[0] < 0)? 16 : 0 )
						+ ( (posB[1] > 0)? 8 : 0 )
						+ ( (posB[1] < 0)? 4 : 0 )
						+ ( (posB[2] > 0)? 2 : 0 )
						+ ( (posB[2] < 0)? 1 : 0 );

			//on determine si on a une intersection de type arrete//plan et on memorise l'edge correspondant
			unsigned short edgeA = 3, edgeB = 3;
			if(     posAbin == 1  || posAbin == 2 ) edgeA = 1; //points 0 et 1 sur le plan
			else if(posAbin == 16 || posAbin == 32) edgeA = 2; //points 1 et 2 sur le plan
			else if(posAbin == 4  || posAbin == 8 ) edgeA = 0; //points 2 et 0 sur le plan
			if(     posBbin == 1  || posBbin == 2 ) edgeB = 1; //points 0 et 1 sur le plan
			else if(posBbin == 16 || posBbin == 32) edgeB = 2; //points 1 et 2 sur le plan
			else if(posBbin == 4  || posBbin == 8 ) edgeB = 0; //points 2 et 0 sur le plan

			//si on a detecté une telle intersection, on verifie si il y a bel et bien un chevauchement des 2 polyedres
			//en cas de simple contact "arrete-plan" ou "arrete-arrete", on ne calcule pas l'intersection
			Vector_exact nA2, nB2;
			num_type p;
			bool dnc = false;
			bool invert_direction = false;

			if(edgeA != 3 && edgeB == 3)
			{
				fA2 = heA[edgeA]->opposite()->facet();
				nA2 = Inter_tri[fA2->Label].norm_dir;
				p = CGAL::cross_product(nA, nB) * CGAL::cross_product(nA2, nB);
				if(p < 0) return;
				if(p == 0)
				{
					// intersection limite de cas coplanaire de 2 triangles
					switch(m_BOOP)
					{
					case UNION:
						if(posA[(edgeA+1)%3] * (nA2 * nB) > 0) dnc = true;
						break;
					case INTER:
						if(posA[(edgeA+1)%3] > 0) dnc = true;
						break;
					case MINUS:
						if(posA[(edgeA+1)%3] * (nA2 * nB) < 0) dnc = true;
						break;
					}
				}
			}
			else if(edgeA == 3 && edgeB != 3)
			{
				fB2 = heB[edgeB]->opposite()->facet();
				nB2 = Inter_tri[fB2->Label].norm_dir;
				p = CGAL::cross_product(nA, nB) * CGAL::cross_product(nA, nB2);
				if(p < 0) return; //contact edge/plan 
				if(p == 0)
				{
					switch(m_BOOP)
					{
					case UNION:
						if(posB[(edgeB+1)%3] < 0) dnc = true;
						break;
					case INTER:
						if(posB[(edgeB+1)%3] * (nB2 * nA) < 0) dnc = true;
						break;
					case MINUS:
						if(posB[(edgeB+1)%3] > 0) dnc = true;
						break;
					}
				}
			}
			else if(edgeA != 3 && edgeB != 3)
			{
				bool Intersection = false;
				Vector_exact nAcnB2, nA2cnB;
				num_type nAnB2, nA2nB, nA2nB2;
				num_type posA2_A, posB_A, posB2_A, posB_B2, posA_B, posB2_B, posB_A2, posB2_A2, posA2_B, posA2_B2;
				Point3d_exact ptA2, ptB2;

				fA2 = heA[edgeA]->opposite()->facet();
				fB2 = heB[edgeB]->opposite()->facet();
				nA2 = Inter_tri[fA2->Label].norm_dir;
				nB2 = Inter_tri[fB2->Label].norm_dir;

				nAcnB2 = CGAL::cross_product(nA, nB2);
				nA2cnB = CGAL::cross_product(nA2, nB);

				nAnB2 = nA * nB2;
				nA2nB = nA2 * nB;
				nA2nB2 = nA2 * nB2;

				ptA2 = point_to_exact(heA[edgeA]->opposite()->next()->vertex()->point());
				ptB2 = point_to_exact(heB[edgeB]->opposite()->next()->vertex()->point());

				posA_B = posA[(edgeA+1)%3];
				posB_A = posB[(edgeB+1)%3];
				posB_A2 = nA2 * (ptB[(edgeB+1)%3] - ptA[edgeA]);
				posB_B2 = nB2 * (ptB[(edgeB+1)%3] - ptA[edgeA]);
				posA2_A = nA * (ptA2 - ptA[edgeA]);
				posA2_B = nB * (ptA2 - ptA[edgeA]);
				posA2_B2 = nB2 * (ptA2 - ptA[edgeA]);
				posB2_A = nA * (ptB2 - ptA[edgeA]);
				posB2_A2 = nA2 * (ptB2 - ptA[edgeA]);
				posB2_B = nB * (ptB2 - ptA[edgeA]);

				if(nAcnB2 == CGAL::NULL_VECTOR && nA2cnB == CGAL::NULL_VECTOR
					&& nAnB2 * nA2nB > 0) return;

				if(posB_A * posB2_A > 0)
				{
					if(posB_B2 > 0) Intersection = true;
				}
				else if(posB_A * posB2_A < 0)
				{
					if(posA_B < 0) Intersection = true;
				}
				else
				{
					if(posA_B * posB2_B < 0) //non coplanaire
					{
						if(posB_B2 > 0) Intersection = true;
					}
					else //coplanaire
					{
						if(nAnB2 < 0)
						{
							if(m_BOOP == UNION) Intersection = true;
						}
						else
						{
							if(m_BOOP == MINUS) Intersection = true;
						}
					}
				}

				if(posB_A2 * posB2_A2 > 0)
				{
					if(posB_B2 > 0) Intersection = !Intersection;
				}
				else if(posB_A2 * posB2_A2 < 0)
				{
					if(posA2_B < 0) Intersection = !Intersection;
				}
				else if(posB2_A2 == 0)
				{
					if(posA2_B * posB2_B < 0) //non coplanaire
					{
						if(posB_B2 > 0) Intersection = !Intersection;
					}
					else //coplanaire
					{
						if(nA2nB2 < 0)
						{
							if(m_BOOP == UNION) Intersection = !Intersection;
						}
						else
						{
							if(m_BOOP == MINUS) Intersection = !Intersection;
						}
					}
				}
				else if(posB_A2 == 0)
				{
					if(posA2_B2 * posB_B2 < 0) //non coplanaire
					{
						if(posB_B2 > 0) Intersection = !Intersection;
					}
					else //coplanaire
					{
						if(nA2nB < 0)
						{
							if(m_BOOP == UNION) Intersection = !Intersection;
						}
						else
						{
							if(m_BOOP == MINUS) Intersection = !Intersection;
						}
					}
				}

				if(!Intersection) dnc = true;

				//inversion ?

				if(posB_A * posA2_A > 0 && ((posA2_B != 0)?posA2_B:posB2_B) * posA_B > 0) invert_direction = true;
			}

			if(posAbin == 5 || posAbin == 10 || posAbin == 17 || posAbin == 34 || posAbin == 20 || posAbin == 40
			|| posBbin == 5 || posBbin == 10 || posBbin == 17 || posBbin == 34 || posBbin == 20 || posBbin == 40) return;

#ifdef BOOLEAN_OPERATIONS_DEBUG
			d_i_tt_verif_particular_case += Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			bool res;
			bool ptA_on_triB = false;
			bool ptB_on_triA = false;
			InterId ptI;
			vector<InterId> ptInter;

			switch(posBbin)
			{
			//intersections courantes
			case 26:
			case 37:
				res = InterTriangleSegment(heB[0], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = InterTriangleSegment(heB[1], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			case 25:
			case 38:
				res = InterTriangleSegment(heB[1], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = InterTriangleSegment(heB[2], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			case 22:
			case 41:
				res = InterTriangleSegment(heB[2], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = InterTriangleSegment(heB[0], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			//cas particulier d'un sommet de B sur le plan (intersection de type arete-plan)
			case 6:
			case 9:
				res = InterTriangleSegment(heB[2], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = IsInTriangle(heB[0], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			case 18:
			case 33:
				res = InterTriangleSegment(heB[0], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = IsInTriangle(heB[1], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			case 24:
			case 36:
				res = InterTriangleSegment(heB[1], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = IsInTriangle(heB[2], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			//intersection arete-plan
			case 1:
			case 2:
				res = IsInTriangle(heB[0], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = IsInTriangle(heB[1], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			case 16:
			case 32:
				res = IsInTriangle(heB[1], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = IsInTriangle(heB[2], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			case 4:
			case 8:
				res = IsInTriangle(heB[2], fA, ptI);
				if(res) ptInter.push_back(ptI);
				res = IsInTriangle(heB[0], fA, ptI);
				if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
				break;
			default:
				return;
			//rien a calculer
			//les autres nombres possible sont 5, 10, 17, 34, 20, 40 (intersection point sur plan)
			//et 21, 42 (aucune intersection : theoriquement impossible)
			}
			if(ptInter.size() != 2)
			{
				switch(posAbin)
				{
				//intersections courantes
				case 26:
				case 37:
					res = InterTriangleSegment(heA[0], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = InterTriangleSegment(heA[1], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				case 25:
				case 38:
					res = InterTriangleSegment(heA[1], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = InterTriangleSegment(heA[2], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				case 22:
				case 41:
					res = InterTriangleSegment(heA[2], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = InterTriangleSegment(heA[0], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				//cas particulier d'un sommet de B sur le plan (intersection de type arete-plan)
				case 6:
				case 9:
					res = InterTriangleSegment(heA[2], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = IsInTriangle(heA[0], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				case 18:
				case 33:
					res = InterTriangleSegment(heA[0], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = IsInTriangle(heA[1], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				case 24:
				case 36:
					res = InterTriangleSegment(heA[1], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = IsInTriangle(heA[2], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				//intersection arete-plan
				case 1:
				case 2:
					res = IsInTriangle(heA[0], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = IsInTriangle(heA[1], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				case 16:
				case 32:
					res = IsInTriangle(heA[1], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = IsInTriangle(heA[2], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				case 4:
				case 8:
					res = IsInTriangle(heA[2], fB, ptI);
					if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					if(ptInter.size() != 2)
					{
						res = IsInTriangle(heA[0], fB, ptI);
						if(res && (ptInter.size() == 0 || ptInter[0] != ptI)) ptInter.push_back(ptI);
					}
					break;
				default:
					return;
				//rien a calculer
				//les autres nombres possible sont 5, 10, 17, 34, 20, 40 (intersection point sur plan)
				//et 21, 42 (aucune intersection : theoriquement impossible)
				}
			}

#ifdef BOOLEAN_OPERATIONS_DEBUG
			d_i_tt_inter_seg_tri += Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			if(ptInter.size() == 2)
			{
				if(dnc)
				{
					return;
				}
				vector<InterId> ptInterInv;
				ptInterInv.push_back(ptInter[1]);
				ptInterInv.push_back(ptInter[0]);
				
				if(CGAL::cross_product(nA, nB) * (InterPts[ptInter[1]] - InterPts[ptInter[0]]) * ((invert_direction == true)?-1:1) > 0)
				{
					switch(m_BOOP)
					{
					case UNION:
						Inter_tri[fA->Label].CutList.push_back(ptInter);
						if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInter);
						Inter_tri[fB->Label].CutList.push_back(ptInterInv);
						if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
						break;
					case INTER:
						Inter_tri[fA->Label].CutList.push_back(ptInterInv);
						if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
						Inter_tri[fB->Label].CutList.push_back(ptInter);
						if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInter);
						break;
					case MINUS:
						Inter_tri[fA->Label].CutList.push_back(ptInter);
						if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInter);
						Inter_tri[fB->Label].CutList.push_back(ptInter);
						if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInter);
						break;
					}
				}
				else
				{
					switch(m_BOOP)
					{
					case UNION:
						Inter_tri[fA->Label].CutList.push_back(ptInterInv);
						if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
						Inter_tri[fB->Label].CutList.push_back(ptInter);
						if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInter);
						break;
					case INTER:
						Inter_tri[fA->Label].CutList.push_back(ptInter);
						if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInter);
						Inter_tri[fB->Label].CutList.push_back(ptInterInv);
						if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
						break;
					case MINUS:
						Inter_tri[fA->Label].CutList.push_back(ptInterInv);
						if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
						Inter_tri[fB->Label].CutList.push_back(ptInterInv);
						if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
						break;
					}
				}
			}

#ifdef BOOLEAN_OPERATIONS_DEBUG
			d_i_tt_add2sol += Timer.GetDiff();
#endif // BOOLEAN_OPERATIONS_DEBUG

		}
	}

	bool InterTriangleSegment(Halfedge_handle& he, Facet_handle& f, InterId& I)
	{
		if(Inter_tri[f->Label].RefInter.count(he->Label) != 0)
		{
			I = Inter_tri[f->Label].RefInter[he->Label];
		}
		else
		{
			unsigned short res = 0;
			Vector_exact e1, e2, dir, p, s, q;
			num_type u, v, tmp;

			Point3d_exact s1 = point_to_exact(he->opposite()->vertex()->point());
			Point3d_exact s2 = point_to_exact(he->vertex()->point());
			Point3d_exact v0 = point_to_exact(f->facet_begin()->vertex()->point());
			Point3d_exact v1 = point_to_exact(f->facet_begin()->next()->vertex()->point());
			Point3d_exact v2 = point_to_exact(f->facet_begin()->next()->next()->vertex()->point());

			e1 = v1 - v0;
			e2 = v2 - v0;
			dir = s2 - s1;
			p = CGAL::cross_product(dir, e2);
			tmp = (num_type)1/(p*e1);
			s = s1 - v0;
			u = tmp * s * p;
			if(u < 0 || u > 1) return false;
			q = CGAL::cross_product(s, e1);
			v = tmp * dir * q;
			if(v < 0 || v > 1) return false;
			if(u + v > 1) return false;


			I = InterPts.size();
			InterPts.push_back(s1+(tmp*e2*q)*dir);
			ppbuilder.add_vertex(point_to_double(InterPts[I]), I);

			if(u == 0) res += 1;	//intersection sur he(0)
			if(v == 0) res += 2;	//intersection sur he(1)
			if(u+v == 1) res += 4;	//intersection sur he(2)

			switch(res)
			{
			case 0:
				{
					Inter_tri[f->Label].RefInter[he->Label] = I;
					Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
				}
				break;
			case 1:
				{
					Inter_tri[he->facet()->Label].PtList.push_back(I);
					Inter_tri[he->opposite()->facet()->Label].PtList.push_back(I);
					Inter_tri[f->Label].PtList.push_back(I);
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].PtList.push_back(I);
					Inter_tri[f->Label].RefInter[he->Label] = I;
					Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[he->Label] = I;
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
				}
				break;
			case 2:
				{
					Inter_tri[he->facet()->Label].PtList.push_back(I);
					Inter_tri[he->opposite()->facet()->Label].PtList.push_back(I);
					Inter_tri[f->Label].PtList.push_back(I);
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].PtList.push_back(I);
					Inter_tri[f->Label].RefInter[he->Label] = I;
					Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[he->Label] = I;
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
				}
				break;
			case 3:
				{
					Inter_tri[he->facet()->Label].PtList.push_back(I);
					Inter_tri[he->opposite()->facet()->Label].PtList.push_back(I);
					f->facet_begin()->vertex()->Label = I;
					//marquer toutes les facettes contenant le point f->facet_begin()->vertex()
					Halfedge_around_vertex_circulator	H_circ = f->facet_begin()->vertex_begin(),
														H_end = f->facet_begin()->vertex_begin();
					do {
						Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
						Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						H_circ++;
					} while(H_circ != H_end);
				}
				break;
			case 4:
				{
					Inter_tri[he->facet()->Label].PtList.push_back(I);
					Inter_tri[he->opposite()->facet()->Label].PtList.push_back(I);
					Inter_tri[f->Label].PtList.push_back(I);
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].PtList.push_back(I);
					Inter_tri[f->Label].RefInter[he->Label] = I;
					Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[he->Label] = I;
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
				}
				break;
			case 5:
				{
					Inter_tri[he->facet()->Label].PtList.push_back(I);
					Inter_tri[he->opposite()->facet()->Label].PtList.push_back(I);
					f->facet_begin()->next()->next()->vertex()->Label = I;
					//marquer toutes les facettes contenant le point f->facet_begin()->next()->next()->vertex()
					Halfedge_around_vertex_circulator	H_circ = f->facet_begin()->next()->next()->vertex_begin(),
														H_end = f->facet_begin()->next()->next()->vertex_begin();
					do {
						Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
						Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						H_circ++;
					} while(H_circ != H_end);
				}
				break;
			case 6:
				{
					Inter_tri[he->facet()->Label].PtList.push_back(I);
					Inter_tri[he->opposite()->facet()->Label].PtList.push_back(I);
					f->facet_begin()->next()->vertex()->Label = I;
					//marquer toutes les facettes contenant le point f->facet_begin()->next()->vertex()
					Halfedge_around_vertex_circulator	H_circ = f->facet_begin()->next()->vertex_begin(),
														H_end = f->facet_begin()->next()->vertex_begin();
					do {
						Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
						Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						H_circ++;
					} while(H_circ != H_end);
				}
				break;
			}
		}

		return true;
	}

	bool IsInTriangle(Halfedge_handle& he, Facet_handle& f, InterId& I)
	{
		unsigned short res = 0;
		Point3d_exact p = point_to_exact(he->vertex()->point());
		Point3d_exact v0 = point_to_exact(f->facet_begin()->vertex()->point());
		Point3d_exact v1 = point_to_exact(f->facet_begin()->next()->vertex()->point());
		Point3d_exact v2 = point_to_exact(f->facet_begin()->next()->next()->vertex()->point());

		Vector_exact N = Inter_tri[f->Label].norm_dir;
		num_type u, v, w;

		u = N * CGAL::cross_product(v0 - v2, p - v2);
		if(u < 0) return false;
		v = N * CGAL::cross_product(v1 - v0, p - v0);
		if(v < 0) return false;
		w = N * CGAL::cross_product(v2 - v1, p - v1);
		if(w < 0) return false;

		I = he->vertex()->Label;
		if(I == 0xFFFFFFFF)
		{
			I = InterPts.size();
			he->vertex()->Label = I;
			InterPts.push_back(p);
			ppbuilder.add_vertex(he->vertex()->point(), I);
		}

		if(u == 0) res += 1;	//intersection sur he(0)
		if(v == 0) res += 2;	//intersection sur he(1)
		if(w == 0) res += 4;	//intersection sur he(2)

		switch(res)
		{
		case 0:		//intersection pointe d'edge dans le triangle
			{
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		case 1:
			{
				Inter_tri[f->Label].PtList.push_back(I);
				Inter_tri[f->facet_begin()->opposite()->facet()->Label].PtList.push_back(I);
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		case 2:
			{
				Inter_tri[f->Label].PtList.push_back(I);
				Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].PtList.push_back(I);
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		case 3:
			{
				f->facet_begin()->vertex()->Label = I;
				Halfedge_around_vertex_circulator	H_circ = he->vertex_begin(),
													H_end = he->vertex_begin();
				do {
					Halfedge_around_vertex_circulator	F_circ = f->facet_begin()->vertex_begin(),
														F_end = f->facet_begin()->vertex_begin();
					do {
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
						F_circ++;
					} while(F_circ != F_end);
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		case 4:
			{
				Inter_tri[f->Label].PtList.push_back(I);
				Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].PtList.push_back(I);
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		case 5:
			{
				f->facet_begin()->next()->next()->vertex()->Label = I;
				Halfedge_around_vertex_circulator 	H_circ = he->vertex_begin(),
													H_end = he->vertex_begin();
				do {
					Halfedge_around_vertex_circulator 	F_circ = f->facet_begin()->next()->next()->vertex_begin(),
														F_end = f->facet_begin()->next()->next()->vertex_begin();
					do {
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
						F_circ++;
					} while(F_circ != F_end);
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		case 6:
			{
				f->facet_begin()->next()->vertex()->Label = I;
				Halfedge_around_vertex_circulator	H_circ = he->vertex_begin(),
													H_end = he->vertex_begin();
				do {
					Halfedge_around_vertex_circulator	F_circ = f->facet_begin()->next()->vertex_begin(),
														F_end = f->facet_begin()->next()->vertex_begin();
					do {
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
						F_circ++;
					} while(F_circ != F_end);
					H_circ++;
				} while(H_circ != H_end);
			}
			break;
		}
		return true;
	}

	void add_facet_to_solution(Facet_handle pFacet, bool facet_from_A)
	{
		if(pFacet->Label < Inter_tri.size())
		{
			Triangle_Cut TriCut = Inter_tri[pFacet->Label];
			Halfedge_handle he = pFacet->facet_begin();
			Triangulation<Exact_Kernel> T(he, TriCut.norm_dir);
			for(int i = 0;i!=TriCut.PtList.size();++i)
			{
				T.add_new_pt(InterPts[TriCut.PtList[i]], TriCut.PtList[i]);
			}
			vector<vector<unsigned long>> Tri_set = T.get_all_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false);
			ppbuilder.add_triangle(Tri_set, he);
		}
		else
		{
			if(m_BOOP == MINUS && !facet_from_A) ppbuilder.add_triangle(pFacet, true);
			else ppbuilder.add_triangle(pFacet, false);
		}
		pFacet->facet_begin()->opposite()->facet()->IsExt = true;
		pFacet->facet_begin()->next()->opposite()->facet()->IsExt = true;
		pFacet->facet_begin()->next()->next()->opposite()->facet()->IsExt = true;
	} 
	
#ifdef BOOLEAN_OPERATIONS_DEBUG
	void ColorType()
	{
		Facet_iterator pFacet =	NULL;
		Point3d tests[6];
		for (pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++)
		{
			if(pFacet->Label < Inter_tri.size()) pFacet->color(1.0, 0.0 ,0.0);
			else if(pFacet->IsExt) pFacet->color(0.0, 1.0 ,0.0);
			else pFacet->color(0.0, 0.0 ,1.0);
		}
		for (pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++)
		{
			if(pFacet->Label < Inter_tri.size()) pFacet->color(1.0, 0.0 ,0.0);
			else if(pFacet->IsExt) pFacet->color(0.0, 1.0 ,0.0);
			else pFacet->color(0.0, 0.0 ,1.0);
		}
	}

	void WriteData(PolyhedronPtr pMout)
	{
		for(Facet_iterator pFacet = m_pA->facets_begin();pFacet != m_pA->facets_end();++pFacet)
		{
			if(pFacet->Label < Inter_tri.size()) N_IFA++;
			else if(pFacet->IsExt) N_FFA++;
			else N_LFFA++;
		}
		for(Facet_iterator pFacet = m_pB->facets_begin();pFacet != m_pB->facets_end();++pFacet)
		{
			if(pFacet->Label < Inter_tri.size()) N_IFB++;
			else if(pFacet->IsExt) N_FFB++;
			else N_LFFB++;
		}
		N_CF = pMout->size_of_facets() - N_FFA - N_FFB;

		ofstrtime << "Temps de calcul :"																	<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "algorithme"																			<< std::endl;
		ofstrtime << "   Initialisation :               "			<< tr(duration_Init)					<< std::endl;
		ofstrtime << " + Recherche des intersections :  "			<< tr(duration_FindCouples)				<< std::endl;
		ofstrtime << " + Calcul des Intersections :     "			<< tr(duration_ComputeIntersections)	<< std::endl;
		ofstrtime << " + Decoupe des facettes :         "			<< tr(duration_CutIntersectedFacets)	<< std::endl;
		ofstrtime << " + Propagation de la solution :   "			<< tr(duration_PropagateFacets)			<< std::endl;
		ofstrtime << " + Suppression vertex innutiles : "			<< tr(duration_clear_unused_vertices)	<< std::endl;
		ofstrtime << " + Creation du polyedre           "			<< tr(duration_delegate)				<< std::endl;
		ofstrtime << " + Ecriture du fichier            "			<< tr(duration_write_off)				<< std::endl;
		ofstrtime << "-----------------------------------------"											<< std::endl;
		ofstrtime << " Total :                          "			<< tr(duration_total)					<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "  Initialisation"																		<< std::endl;
		ofstrtime << "     Copie des polyedres :          "			<< tr(duration_Copy_Polyhedra)			<< std::endl;
		ofstrtime << "   + Triangulation :                "			<< tr(duration_triangulate)				<< std::endl;
		ofstrtime << "   + Initialisation des tags :      "			<< tr(duration_Init_tags)				<< std::endl;
		ofstrtime << "  -----------------------------------------"											<< std::endl;
		ofstrtime << "   Total :                          "			<< tr(duration_Init)					<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "  Recherche des intersections"														<< std::endl;
		ofstrtime << "     Construction de l'arbre :      "			<< tr(duration_Build_AABBTree)			<< std::endl;
		ofstrtime << "   + Test d'intersection :          "			<< tr(duration_Test_AABBTree)			<< std::endl;
		ofstrtime << "  -----------------------------------------"											<< std::endl;
		ofstrtime << "   Total :                          "			<< tr(duration_FindCouples)				<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "  Calcul des Intersections"															<< std::endl;
		ofstrtime << "     extraire couple a traiter :    "			<< tr(d_i_recup_couple)					<< std::endl;
		ofstrtime << "   + calcul d'inter tri-tri :       "			<< tr(d_i_tri_tri)						<< std::endl;
		ofstrtime << "  -----------------------------------------"											<< std::endl;
		ofstrtime << "   Total :                          "			<< tr(duration_ComputeIntersections)	<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "    calcul d'inter tri-tri"															<< std::endl;
		ofstrtime << "       Initialisation :               "		<< tr(d_i_tt_init)						<< std::endl;
		ofstrtime << "     + verif cas particuliers :       "		<< tr(d_i_tt_verif_particular_case)		<< std::endl;
		ofstrtime << "     + inter triangle-segment :       "		<< tr(d_i_tt_inter_seg_tri)				<< std::endl;
		ofstrtime << "     + memorisation de l'inter :      "		<< tr(d_i_tt_add2sol)					<< std::endl;
		ofstrtime << "    -----------------------------------------"										<< std::endl;
		ofstrtime << "     Total :                          "		<< tr(d_i_tri_tri)						<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "  Decoupe des facettes"																<< std::endl;
		ofstrtime << "     Initialisation :               "			<< tr(d_Cut_init)						<< std::endl;
		ofstrtime << "   + ajout des points :             "			<< tr(d_Cut_add_pts)					<< std::endl;
		ofstrtime << "   + ajout des segments :           "			<< tr(d_Cut_add_segments)				<< std::endl;
		ofstrtime << "   + selection des bons triangles : "			<< tr(d_Cut_triangulate)				<< std::endl;
		ofstrtime << "   + ajout a la solution :          "			<< tr(d_Cut_add_tri)					<< std::endl;
		ofstrtime << "   + validation des voisins :       "			<< tr(d_Cut_validation)					<< std::endl;
		ofstrtime << "  -----------------------------------------"											<< std::endl;
		ofstrtime << "   Total :                          "			<< tr(duration_CutIntersectedFacets)	<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "  Propagation de la solution"															<< std::endl;
		ofstrtime << "     Initialisation Propagation A : "			<< tr(duration_Prop_InitA)				<< std::endl;
		ofstrtime << "   + Propagation A :                "			<< tr(duration_Prop_A)					<< std::endl;
		ofstrtime << "   + Initialisation Propagation B : "			<< tr(duration_Prop_InitB)				<< std::endl;
		ofstrtime << "   + Propagation B :                "			<< tr(duration_Prop_B)					<< std::endl;
		ofstrtime << "  -----------------------------------------"											<< std::endl;
		ofstrtime << "   Total :                          "			<< tr(duration_PropagateFacets)			<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "Details :"																			<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "Polyedre A :"																			<< std::endl;
		ofstrtime << "Nombre de facettes :              "			<< m_pA->size_of_facets()				<< std::endl;
		ofstrtime << "Nombre de facettes intersectees : "			<< N_IFA								<< std::endl;
		ofstrtime << "Nombre de facettes perdues :      "			<< N_LFFA								<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "Polyedre B :"																			<< std::endl;
		ofstrtime << "Nombre de facettes :              "			<< m_pB->size_of_facets()				<< std::endl;
		ofstrtime << "Nombre de facettes intersectees : "			<< N_IFB								<< std::endl;
		ofstrtime << "Nombre de facettes perdues :      "			<< N_LFFB								<< std::endl;
		ofstrtime																							<< std::endl;
		ofstrtime << "Polyedre Solution :"																	<< std::endl;
		ofstrtime << "Nombre de facettes :              "			<< pMout->size_of_facets()				<< std::endl;
		ofstrtime << "Nombre de facettes de A :         "			<< N_FFA								<< std::endl;
		ofstrtime << "Nombre de facettes de B :         "			<< N_FFB								<< std::endl;
		ofstrtime << "Nombre de facettes Creees:        "			<< N_CF									<< std::endl;
	}
#endif // BOOLEAN_OPERATIONS_DEBUG
	
	//attributs

	//elements initiaux
	Bool_Op m_BOOP;										//operation booleenne en cours
	PolyhedronPtr m_pA, m_pB;							//pointeurs vers polyhedres source
	CPolyhedron_from_polygon_builder_3<HDS> ppbuilder;	//pour la construction du polyedre solution

	//intersections a traiter
	vector<pair<FacetId, FacetId>> m_Couples;			//liste des couples de triangles (dans l'ordre A - B)
	vector<Point3d_exact> InterPts;						//liste des points d'intersection calculés
	vector<Triangle_Cut> Inter_tri;						//informations de decoupage pour les triangles intersectés
	vector<Facet_handle> Facet_Handle;					//index des pointeurs vers facette

	AABB_Tree tree;

	//statistiques
#ifdef BOOLEAN_OPERATIONS_DEBUG
	std::ofstream ofstrtime;

	unsigned int N_IFA;						//nombre de facettes intersectées du polyedre A
	unsigned int N_IFB;						//nombre de facettes intersectées du polyedre B
	unsigned int N_FFA;						//nombre de facettes provenant du polyedre A
	unsigned int N_FFB;						//nombre de facettes provenant du polyedre B
	unsigned int N_LFFA;					//nombre de facettes perdues du polyedre A
	unsigned int N_LFFB;					//nombre de facettes perdues du polyedre B
	unsigned int N_CF;						//nombre de facettes créees

	double duration_Init;					//temps total de la phase d'initialisation
	double duration_Copy_Polyhedra;			//temps de copie des polyedres initiaux
	double duration_triangulate;			//temps de triangulation
	double duration_Init_tags;				//temps d'initialisation des tags

	double duration_FindCouples;			//temps total de recherche d'intersection
	double duration_Build_AABBTree;			//temps de construction de l'arbre
	double duration_Test_AABBTree;			//temps de test d'intersections

	double duration_ComputeIntersections;	//temps total pour le calcul des intersections
	double d_i_recup_couple;				//temps de recuperation du couple a traiter dans la liste
	double d_i_tri_tri;						//temps de calcul de l'intersection triangle triangle
	double d_i_tt_init;						//temps d'initialisation de l'intersection triangle triangle 
	double d_i_tt_verif_particular_case;	//temps de verification de certains cas particuliers
	double d_i_tt_inter_seg_tri;			//temps de calcul des intersections triangle segment
	double d_i_tt_add2sol;					//temps d'ajout de memorisation de la decoupe a faire

	double duration_CutIntersectedFacets;	//temps total de decoupe des facettes
	double d_Cut_init;						//temps d'initialisation de la triangulation
	double d_Cut_add_pts;					//temps d'ajout de points
	double d_Cut_add_segments;				//temps d'ajout de segments
	double d_Cut_triangulate;				//temps de calcul de la triangulation
	double d_Cut_add_tri;					//temps d'ajout des triangles a la solution
	double d_Cut_validation;				//temps de validation des voisins

	double duration_PropagateFacets;		//temps total de propagation
	double duration_Prop_InitA;				//temps d'initialisation pour la propagation sur A
	double duration_Prop_A;					//temps de propagation sur A
	double duration_Prop_InitB;				//temps d'initialisation pour la propagation sur B
	double duration_Prop_B;					//temps de propagation sur B

	double duration_clear_unused_vertices;  //temps de suppression des vertex inutiles

	double duration_delegate;				//temps de creation du polyedre

	double duration_write_off;				//temps d'ecriture du fichier off

	double duration_total;					//temps total d'execution
#endif // BOOLEAN_OPERATIONS_DEBUG
};

#endif
