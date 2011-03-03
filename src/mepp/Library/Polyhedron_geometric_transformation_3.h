/***************************************************************************
Polyhedron_geometric_transformation_3.h - Transformations affecting only the
polyhedron's geometry
----------------------------------------------------------------------------
begin                : March 2009
copyright            : (C) 2008 Hichem Barki, LIRIS M2DisCo
email                : hichem.barki@liris.cnrs.fr
***************************************************************************/

#ifndef _POLYHEDRON_GEOMETRIC_TRANSFORMATION_3_H
#define _POLYHEDRON_GEOMETRIC_TRANSFORMATION_3_H

#include <cmath>
#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/number_utils.h>

namespace Polyhedron_geometric_transformation_3
{
	/*******************************************************************************************************
	*                                  Polyhedron_3 translation                                            *
	*******************************************************************************************************/
	template <typename Polyhedron_3>
	inline void translate(Polyhedron_3& P, typename const Polyhedron_3::Point_3& source, typename const Polyhedron_3::Point_3& target)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits::Vector_3												Vector_3;
		translate<Polyhedron_3>(P, Vector_3(source, target));
	}
	
	template <typename Polyhedron_3>
	inline void translate(Polyhedron_3& P, typename const Polyhedron_3::Traits::Vector_3& transl)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits::Aff_transformation_3									Aff_transf;
		if (transl != CGAL::NULL_VECTOR)
		{
			Aff_transf translation(CGAL::TRANSLATION, transl);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation);
		}
	}
	/*******************************************************************************************************
	*                                  Polyhedron_3 scaling                                                *
	*******************************************************************************************************/
	template <typename Polyhedron_3>
	inline void scale(Polyhedron_3& P, typename const Polyhedron_3::Traits::FT& factor, typename const Polyhedron_3::Point_3& center = CGAL::ORIGIN)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits::Aff_transformation_3									Aff_transf;
		typedef typename Polyhedron_3::Traits::Vector_3												Vector_3;
		Aff_transf scaling(CGAL::SCALING, factor);
		if (center == CGAL::ORIGIN)
		{
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), scaling);
		}
		else
		{
			Aff_transf translation(CGAL::TRANSLATION, Vector_3(center,CGAL::ORIGIN));
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), scaling);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation.inverse());
		}
	}
	/*******************************************************************************************************
	*                                  Polyhedron_3 rotation                                               *
	*******************************************************************************************************/
	template <typename Number_type>
	inline void set_identity_matrix(Number_type matrix[4][4])
	{
		for (int r = 0; r != 4; ++r)
		{
			for (int c = 0; c != 4; ++c)
			{
				matrix_4X4[r][c] = Number_type(r == c);
			}
		}
	}
	template <typename Number_type>
	inline void get_transpose_matrix(const Number_type matrix[4][4], Number_type transpose[4][4])
	{
		for (int r = 0; r != 4; ++r)
			for (int c = 0; c != 4; ++c)
				transpose[r][c] = matrix[c][r];
	}
	template <typename Kernel>
	inline void get_rotation_matrix(typename const Kernel::Vector_3& u, typename const Kernel::FT& radAngle, typename Kernel::FT rot_matrix[4][4])
	{
		typedef typename Kernel::FT																	Number_type;
		typedef typename Kernel::Vector_3															Vector_3;

		Number_type length = Number_type(sqrt(CGAL::to_double(u.squared_length())));
		Number_type cosT2 = Number_type(cos(CGAL::to_double(radAngle / 2.)));
		Number_type sinT2 = Number_type(sin(CGAL::to_double(radAngle / 2.)));
		Number_type a = sinT2 * u.x() / length;
		Number_type b = sinT2 * u.y() / length;
		Number_type c = sinT2 * u.z() / length;
		rot_matrix[0][0] = 1 - 2*b*b - 2*c*c;
		rot_matrix[0][1] = 2*a*b - 2 * cosT2 * c;
		rot_matrix[0][2] = 2*a*c + 2*cosT2*b;
		rot_matrix[0][3] = 0;
		rot_matrix[1][0] = 2*a*b + 2*cosT2*c;
		rot_matrix[1][1] = 1 - 2*a*a - 2*c*c;
		rot_matrix[1][2] = 2*b*c - 2*cosT2*a;
		rot_matrix[1][3] = 0;
		rot_matrix[2][0] = 2*a*c - 2*cosT2*b;
		rot_matrix[2][1] = 2*b*c + 2*cosT2*a;
		rot_matrix[2][2] = 1 - 2*a*a - 2*b*b;
		rot_matrix[2][3] = 0;
		rot_matrix[3][0] = 0;
		rot_matrix[3][1] = 0;
		rot_matrix[3][2] = 0;
		rot_matrix[3][3] = 1;
	}
	
	template <typename Polyhedron_3>
	inline void rotate(Polyhedron_3& P, typename const Polyhedron_3::Traits::Vector_3& rot_axis,
		typename const Polyhedron_3::Traits::FT& radAngle, typename const Polyhedron_3::Traits::Point_3 center = CGAL::ORIGIN)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits::FT													Number_type;
		typedef typename Polyhedron_3::Traits::Vector_3												Vector_3;
		typedef typename Polyhedron_3::Traits::Aff_transformation_3									Aff_transf;

		Number_type rot_matrix[4][4];

		get_rotation_matrix<Polyhedron_3::Traits>(rot_axis, radAngle, rot_matrix);
		Aff_transf rotation(rot_matrix[0][0],rot_matrix[0][1],rot_matrix[0][2],
							rot_matrix[1][0],rot_matrix[1][1],rot_matrix[1][2],
							rot_matrix[2][0],rot_matrix[2][1],rot_matrix[2][2]);
		if (center == CGAL::ORIGIN)
		{
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), rotation);
		}
		else
		{			
			Aff_transf translation(CGAL::TRANSLATION, Vector_3(center,CGAL::ORIGIN));
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), rotation);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation.inverse());
		}
	}
	/*******************************************************************************************************
	*                                  Polyhedron_3 reflexion                                              *
	*******************************************************************************************************/
	template <typename Polyhedron_3>
	inline void reflect(Polyhedron_3& P, typename const Polyhedron_3::Point_3& center = CGAL::ORIGIN)
	{
		CGAL_assertion(&P != NULL);

		typedef typename Polyhedron_3::Traits::Aff_transformation_3									Aff_transf;
		typedef typename Polyhedron_3::Traits::Vector_3												Vector_3;
		Aff_transf reflection(-1, 0, 0, 0, -1, 0, 0, 0, -1, 1);
		
		if (center == CGAL::ORIGIN)
		{
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), reflection);
		}
		else
		{
			Aff_transf translation(CGAL::TRANSLATION, Vector_3(center,CGAL::ORIGIN));
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), reflection);
			std::transform(P.points_begin(), P.points_end(), P.points_begin(), translation.inverse());
		}
	}

} // end of namespace Polyhedron_geometric_transformation_3
#endif // end of _POLYHEDRON_GEOMETRIC_TRANSFORMATION_3_H