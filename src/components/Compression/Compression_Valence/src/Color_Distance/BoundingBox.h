/***************************************************************************
                                BoundingBox.h
                             -------------------
    update               : 2003-01-16
    copyright            : (C) 2002-2003 by Michaël Roy
    email                : michaelroy@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef _BOUNDINGBOX_
#define _BOUNDINGBOX_

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "VectorT.h"
#include <float.h>
#include <iostream>
#include <vector>

//--
//
// BoundingBox
//
//--
// Bounding box of points in nD.
// Points must have double precision float values
template<int Size>
class BoundingBox
{
	//--
	//
	// Type definition
	//
	//--
	protected :

		typedef VectorT<double,Size> Point;
		
	//--
	//
	// Member data
	//
	//--
	protected :

		Point minpt;
		Point maxpt;
		Point total;
		int point_number;

	//--
	//
	// Member functions
	//
	//--
	public :

		// Default construtor
		inline BoundingBox()
		: minpt(DBL_MAX), maxpt(DBL_MIN), total(0.0), point_number(0) {
		}

		// Copy construtor
		inline BoundingBox(const BoundingBox& bb)
		: minpt(bb.minpt), maxpt(bb.maxpt), total(bb.total), point_number(bb.point_number) {
		}

		//--
		//
		// Data management
		//
		//--
		
		// Reset data
		inline BoundingBox& Reset() {
			minpt = DBL_MAX;
			minpt = DBL_MIN;
			total = 0.0;
			point_number = 0;
			return *this;
		}

		// Add a point
		inline BoundingBox<Size>& AddPoint( const Point& p ) {
			for ( int i=0; i<Size; i++ ) {
				if ( p[i] < minpt[i] ) minpt[i] = p[i];
				if ( p[i] > maxpt[i] ) maxpt[i] = p[i];
			}
			total += p;
			point_number++;
			return *this;
		}

		// Add an array of points
		inline BoundingBox<Size>& AddPoints( const std::vector<Point>& v ) {
			for ( int i=0; i<(int)v.size(); i++ ) {
				AddPoint( v[i] );
			}
			return *this;
		}
		
		//--
		//
		// Data accessors
		//
		//--

		// Minimum bounds (constant)
		inline const Point& Min() const {
			return minpt;
		}

	     	// Maximum bounds (constant)
		inline const Point& Max() const {
			return maxpt;
		}

		// Size of the bounding box (constant)
		inline Point Length() const {
			return maxpt - minpt;
		}

		// Diagonal of the bounding box (constant)
		inline double Diagonal() const {
			return Length().Length();
		}

		// Center of the bounding box (constant)
		inline Point Center() const {
			return total / (double)point_number;
		}

		//--
		//
		// Operators
		//
		//--
		inline BoundingBox<Size>& operator=(const BoundingBox<Size>& bb) {
			minpt = bb.minpt;
			maxpt = bb.maxpt;
			total = bb.total;
			point_number = bb.point_number;
			return *this;
		}

		inline BoundingBox<Size>& operator+=(const BoundingBox<Size>& bb) {
			for ( int i=0; i<Size; i++ ) {
				if ( bb.minpt[i] < minpt[i] ) minpt[i] = bb.minpt[i];
				if ( bb.maxpt[i] > maxpt[i] ) maxpt[i] = bb.maxpt[i];
			}
			total += bb.total;
			point_number += bb.point_number;
			return *this;
		}
	
		inline BoundingBox<Size>& operator+=(const Point& p) {
			AddPoint(p);
			return *this;
		}

		inline BoundingBox<Size>& operator+=(const std::vector<Point>& v) {
			AddPoints(v);
			return *this;
		}
};

//--
//
//Other operators
//
//--
template<int Size>
inline BoundingBox<Size> operator+(const BoundingBox<Size>& bb1, const BoundingBox<Size>& bb2)
{
	BoundingBox<Size> result(bb1);
	return result += bb2;
}

template<int Size>
inline std::ostream& operator<<(std::ostream& out, const BoundingBox<Size>& bb)
{
	// Output minimum and maximum points of the bounding box
	return out<<"Min: "<<bb.Min()<<" - Max: "<<bb.Max()<<" - Center: "<<bb.Center();
}

//--
//
// Type definitions
//
//--
typedef BoundingBox<1> BoundingBox1d;
typedef BoundingBox<2> BoundingBox2d;
typedef BoundingBox<3> BoundingBox3d;
typedef BoundingBox<4> BoundingBox4d;

#endif

#endif // _BOUNDINGBOX_

