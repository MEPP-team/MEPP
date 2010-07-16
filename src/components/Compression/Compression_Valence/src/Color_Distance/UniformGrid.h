/***************************************************************************
                               UniformGrid.h
                             -------------------
    update               : 2002-10-06
    copyright            : (C) 2002 by Micha? ROY
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

#ifndef _UNIFORMGRID_
#define _UNIFORMGRID_

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

//#include "MeshDev.h"
#include "Neighborhood.h"
#include "BoundingBox.h"
#include "ColorMesh.h"
#include "Plane.h"
#include "VectorT.h"


////////////////////////////////////////////////////
// Cell structure
////////////////////////////////////////////////////
struct Cell3D
{
	int v;	// Vertex number
	int f;
	Cell3D* next; // Pointer to next cell
};

////////////////////////////////////////////////////
// CUniformGrid Class
////////////////////////////////////////////////////
class UniformGrid
{

	//////////////////////////////////////////////
	// Public member functions
	public :

		// Constructor
		UniformGrid( Mesh_roy* m, BoundingBox3d* bbox, double d );

		// Destructor
		~UniformGrid();
		
		// Register faces
		void SetFaces();
		
		// Compute distance from given point to nearest neighbors
		Neighborhood* NearestNeighbors( const Vector3d& point );
		
		int FacesTestedNumber() { return _FacesTested; }
		
	//////////////////////////////////////////////
	// Protected member function
	private :

		void DistancePoint2Face( const Vector3d& p, int f );

		// Add point to cell
		void AddOnePoint( int n, int x, int y, int z);
		
		// Add face to cell
		void AddOneFace( int n, int x, int y, int z );

		// Cells size (same size for each dimension)
		double m_rSize;

		// Cells array
		Cell3D**** m_pCell;

		// Min cell coordinates
		Vector3d m_pMin;

		// Cells number in the 3 dimensions
		Vector3i m_pCellNum;

		// Faces tested number
		int _FacesTested;
		
		// Nearest Neighbors
		Neighborhood*	neighbors;
		
		std::vector<Vector3d>& mv;
		std::vector<Vector3i>& mf;
		std::vector<Vector3d>& mfn;
		std::vector<double>    mp; // Mesh_roy face planes

		double Area2D( const Vector2d& a, const Vector2d& b, const Vector2d& c );
		bool InTriangle( const Vector2d& a, const Vector2d& b, const Vector2d& c, const Vector2d& p );
		double DistancePoint2Plane(const Vector3d &v, const Vector3d& n, const double& h);
		int Clamp( int x, int max );
};

//========================================================
//
// Area2D
//
//========================================================
// Area of 2D triangle
inline double UniformGrid::Area2D(const Vector2d& a, const Vector2d& b, const Vector2d& c)
{
	return( 	(b[0] - a[0]) * (c[1] - a[1]) -
				(c[0] - a[0]) * (b[1] - a[1]) );
}


//========================================================
//
// DistancePoint2Plane
//
//========================================================
// Distance point to plane
inline double UniformGrid::DistancePoint2Plane(const Vector3d &v, const Vector3d& n, const double& h)
{
	return ( v | n ) + h;
}

//========================================================
//
// Clamp
//
//========================================================
// Clamp the input to the specified range
inline int UniformGrid::Clamp( int x, int max )
{ 
	return( (x < 0) ? 0 : (x >= max) ? (max - 1) : x );
}


inline void UniformGrid::DistancePoint2Face( const Vector3d &p, int f )
{
	int i, j, k=0;
	double d, l=0, m, n;
	Vector2d aa, bb, cc, pp;
	Vector3d u, v;


	// Save Current Face Vertices Indices
	const int& a = mf[f][0];
	const int& b = mf[f][1];
	const int& c = mf[f][2];


	if ( !(p - mv[a]).Length() ) return;
	if ( !(p - mv[b]).Length() ) return;
	if ( !(p - mv[c]).Length() ) return;

	_FacesTested++;

	/////////////////////////////////////////////
	// Distance Point To Plane
	/////////////////////////////////////////////
	d = DistancePoint2Plane( p, mfn[f], mp[f] );
	// If Distance < Error
	if ( fabs(d) < neighbors->Distance() )
	{
		// Find largest component
		for ( i=0; i<3; i++ )
		{ 
			m = fabs( mfn[f][i] ); // Current Component
			if ( m > l )				// Biggest component
			{
				l = m;		// Save value
				k = i;		// Save component indice
			}
		}
		// Projected Point on plane
		u = p - mfn[f] * d;	
		// project out coordinate "k"
		j = 0;
		for ( i=0; i<3; i++ ) if ( i != k )
		{
			aa[j] = mv[ a ][i];
			bb[j] = mv[ b ][i];
			cc[j] = mv[ c ][i];
			pp[j] = u[i];
			j++;
		}
		// compute areas
		l = Area2D( pp, aa, bb );
		m = Area2D( pp, bb, cc );
		n = Area2D( pp, cc, aa );
		// Test if projected point is in face
		if ( ((l > 0) && (m > 0) && (n > 0)) || ((l < 0) && (m < 0) && (n < 0)) )
		{
		
			v = mv[b] - mv[a];
			v.Normalize();
			
			l = (v^(u - mv[a])).Length() / (v^(mv[c] - mv[a])).Length();

			if ( l > 1 ) l = 1;
			if ( l < 0 ) l = 0;

			
			v = Lerp(mv[a], mv[c], l);
			
			m = (u - v).Length() / (Lerp(mv[b], mv[c], l) - v).Length();
				
			d = fabs(d);
			if ( !((l < 0) || (l > 1) || (m < 0) || (m > 1)) )
			{
				if ( d == neighbors->Distance() )
					neighbors->AddFace( u, f, l, m );
				else neighbors->NewFace( d, u, f, l, m );
				// if distance = 0 -> quit
				if ( !d ) return;
			}
		}
	}

	/////////////////////////////////////////////
	// Distance Point To Edge
	/////////////////////////////////////////////
	u = p - mv[ a ];
	v = mv[ b ] - mv[ a ];
	l = v.Length();
	d = (v^u).Length()/l;
	if ( d < neighbors->Distance() )
	{
		
		m = (u|v)/l;
		if ( (m > 0.0f) && (m < l) )
		{
			m /= l;
			v = Lerp(mv[a], mv[b], m);
			if ( d == neighbors->Distance() )
				neighbors->AddEdge( v, f, 0, m );
			else neighbors->NewEdge( d, v, f, 0, m );
			// if distance = 0 -> quit
			if ( !d ) return;
		}
	}
	
	u = p - mv[ b ];
	v = mv[ c ] - mv[ b ];
	l = v.Length();
	d = (v^u).Length()/l;
	if ( d < neighbors->Distance() )
	{
		m = (u|v)/l;

		if ( (m > 0.0f) && (m < l) )
		{
			m /= l;
			v = Lerp(mv[b], mv[c], m);
			if ( d == neighbors->Distance() )
				neighbors->AddEdge( v, f, 1, m );
			else neighbors->NewEdge( d, v, f, 1, m );
			// if distance = 0 -> quit
			if ( !d ) return;
		}
	}

	u = p - mv[ c ];
	v = mv[ a ] - mv[ c ];
	l = v.Length();
	d = (v^u).Length()/l;
	if ( d < neighbors->Distance() )
	{
		m = (u|v)/l;
		if ( (m > 0.0f) && (m < l) )
		{
			m /= l;
			v = Lerp(mv[c], mv[a], m);
			if ( d == neighbors->Distance() )
				neighbors->AddEdge( v, f, 2, m );
			else neighbors->NewEdge( d, v, f, 2, m );
			// if distance = 0 -> quit
			if ( !d ) return;
		}
	}
}

#endif

#endif
