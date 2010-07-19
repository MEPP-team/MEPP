#ifndef BOOLEAN_OPERATIONS_TRIANGULATION_H
#define BOOLEAN_OPERATIONS_TRIANGULATION_H

#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

template <class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt>>
class Enriched_vertex_base : public Vb
{
public:
	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef Enriched_vertex_base<Gt,Vb2> Other;
	};

private:
	unsigned long m_Label;

public:
	void set_Label(unsigned long Label) {m_Label = Label;}
	unsigned long get_Label() {return m_Label;}
};

template <class Gt, class Fb = CGAL::Constrained_triangulation_face_base_2<Gt>>
class Enriched_face_base : public Fb
{
	typedef Fb Base;

public:
	typedef typename Fb::Triangulation_data_structure::Vertex_handle          Vertex_handle;
	typedef typename Fb::Triangulation_data_structure::Face_handle            Face_handle;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Fb::template Rebind_TDS<TDS2>::Other	Fb2;
		typedef Enriched_face_base<Gt,Fb2>						Other;
	};

private:
	bool m_OK;
	bool m_Ext;

public:
	Enriched_face_base() : Base() {}
	Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2) : Base(v0,v1,v2) {}
	Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
		Face_handle n0, Face_handle n1, Face_handle n2) : Base(v0,v1,v2,n0,n1,n2) {}
	Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2,
					bool c0, bool c1, bool c2 ) : Base(v0,v1,v2,n0,n1,n2) {}

	void set_Ext(bool Ext) {m_Ext = Ext;}
	bool get_Ext() {return m_Ext;}
	void set_OK(bool OK) {m_OK = OK;}
	bool get_OK() {return m_OK;}
};
	
template <class K>
class Triangulation
{
	typedef typename K::Point_3														Point_3;
	typedef typename Enriched_vertex_base<K>										Tri_vb;
	typedef typename Enriched_face_base<K>											Tri_fb;
	typedef typename CGAL::Triangulation_data_structure_2<Tri_vb,Tri_fb>			Tri_DS;
	typedef typename CGAL::No_intersection_tag										Itag;
	typedef typename CGAL::Constrained_Delaunay_triangulation_2<K, Tri_DS, Itag>	Constrained_Delaunay_tri;
	typedef typename Constrained_Delaunay_tri::Vertex_handle						Vertex_handle_tri;
	typedef typename Constrained_Delaunay_tri::Face_handle							Face_handle_tri;
	typedef typename Constrained_Delaunay_tri::Point								Point_tri;
	typedef typename Constrained_Delaunay_tri::Face_iterator						Face_iterator_tri;

public:
	Triangulation(Halfedge_handle he, Vector_exact norm_dir)
	{
		double x = to_double(norm_dir.x());
		double y = to_double(norm_dir.y());
		double z = to_double(norm_dir.z());
		double absx = std::abs(x);
		double absy = std::abs(y);
		double absz = std::abs(z);

		if (absx >= absy && absx >= absz) max_coordinate = (x>0)?0:3;
		else if (absy >= absx && absy >= absz) max_coordinate = (y>0)?1:4;
		else if (absz >= absx && absz >= absy) max_coordinate = (z>0)?2:5;

		v1 = add_new_pt(point_to_exact(he->vertex()->point()), he->vertex()->Label);
		v2 = add_new_pt(point_to_exact(he->next()->vertex()->point()), he->next()->vertex()->Label);
		v3 = add_new_pt(point_to_exact(he->next()->next()->vertex()->point()), he->next()->next()->vertex()->Label);

		if(v2->get_Label() == 0xFFFFFFFF) v2->set_Label(0xFFFFFFFE);
		if(v3->get_Label() == 0xFFFFFFFF) v3->set_Label(0xFFFFFFFD);
	}

	Point_tri get_minvar_point_2(Point_3 p)
	{
		switch(max_coordinate)
		{
		case 0:
			return Point_tri(p.y(),p.z());
			break;
		case 1:
			return Point_tri(p.z(),p.x());
			break;
		case 2:
			return Point_tri(p.x(),p.y());
			break;
		case 3:
			return Point_tri(p.z(),p.y());
			break;
		case 4:
			return Point_tri(p.x(),p.z());
			break;
		case 5:
			return Point_tri(p.y(),p.x());
			break;
		default:
			return Point_tri(p.y(),p.z());
		}
	}

	Vertex_handle_tri add_new_pt(Point_3 p, unsigned long Label)
	{
		if(Label != 0xFFFFFFFF)
			for(unsigned int i = 0;i != pts_point.size();++i)
				if(Label == pts_point[i])
					return pts_vertex[i];
		Vertex_handle_tri v;
		v = ct.insert(get_minvar_point_2(p));
		v->set_Label(Label);
		pts_point.push_back(Label);
		pts_vertex.push_back(v);
		return v;
	}

	void add_segment(Point_3 p1, Point_3 p2, unsigned long Label1, unsigned long Label2)
	{
		c1 = add_new_pt(p1, Label1);
		c2 = add_new_pt(p2, Label2);
		ct.insert_constraint(c1, c2);
	}

	vector<vector<unsigned long>> get_triangles(bool inv_triangles, bool *IsExt)
	{
		IsExt[0] = false;
		IsExt[1] = false;
		IsExt[2] = false;
		vector<vector<unsigned long>> tris;
		//recherche de la premiere face :
		for(Face_iterator_tri fi = ct.faces_begin();fi != ct.faces_end();fi++)
			fi->set_OK(false);

		Face_handle_tri f, f2 = c1->face();
		int i;
		do {
			f = f2;
			f->has_vertex(c1,i);
			f2 = f->neighbor(f->ccw(i));
		} while( ! ( f->has_vertex(c2) && f2->has_vertex(c2) ) );
		if(f->has_vertex(ct.infinite_vertex()))
		{
			f = f2;
			f->set_Ext(false);
		}
		else
		{
			f->set_Ext(true);
		}
		//la facette f est la facette de depart
		std::stack<Face_handle_tri> sfh;
		f->set_OK(true);
		sfh.push(f);
		
		while(!sfh.empty())
		{
			f = sfh.top();
			sfh.pop();

			if(f->get_Ext())
			{
				vector<unsigned long> tri;
				int i;
				tri.push_back(f->vertex(0)->get_Label());
				
				if(f->has_vertex(v1,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[0] = true;
				if(f->has_vertex(v2,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[1] = true;
				if(f->has_vertex(v3,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[2] = true;
				
				if(inv_triangles)
				{
					tri.push_back(f->vertex(2)->get_Label());
					tri.push_back(f->vertex(1)->get_Label());
				}
				else
				{
					tri.push_back(f->vertex(1)->get_Label());
					tri.push_back(f->vertex(2)->get_Label());
				}
				tris.push_back(tri);
			}
			for(i = 0;i!=3;i++)
			{
				if(!(f->neighbor(i)->get_OK() || f->neighbor(i)->has_vertex(ct.infinite_vertex())))
				{
					f->neighbor(i)->set_OK(true);
					f->neighbor(i)->set_Ext((f->is_constrained(i))?!f->get_Ext():f->get_Ext());
					sfh.push(f->neighbor(i));
				}
			}
		}
		return tris;
	}
	vector<vector<unsigned long>> get_all_triangles(bool inv_triangles)
	{
		vector<vector<unsigned long>> tris;
		for(Face_iterator_tri f = ct.faces_begin();f != ct.faces_end();f++)
		{
			vector<unsigned long> tri;
			tri.push_back(f->vertex(0)->get_Label());
			if(inv_triangles)
			{
				tri.push_back(f->vertex(2)->get_Label());
				tri.push_back(f->vertex(1)->get_Label());
			}
			else
			{
				tri.push_back(f->vertex(1)->get_Label());
				tri.push_back(f->vertex(2)->get_Label());
			}
			tris.push_back(tri);
		}
		return tris;
	}
	
private:
	Constrained_Delaunay_tri ct;
	vector<InterId> pts_point;
	vector<Vertex_handle_tri> pts_vertex;
	Vertex_handle_tri v1, v2, v3, c1, c2;
	int max_coordinate;
};

#endif // BOOLEAN_OPERATIONS_TRIANGULATION_H
