/***************************************************************************
                               Neighborhood.h
                             -------------------
    update               : 2002-09-03
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

#ifndef _NEIGHBORHOOD_
#define _NEIGHBORHOOD_

#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

//#include "MeshDev.h"
#include "VectorT.h"

///////////////////////////////////////////////////
// Neighbor structure
struct Neighbor
{
	Vector3d	c;		// Coordinates
	int			v;		// Vertex if similar
	int			f;		// Face containing the point
	int			e;		// Edge containing the point
	double		r1;		// Ratio 1
	double		r2;		// Ratio 2
	Neighbor*	next;	// Next neighbor
};


////////////////////////////////////////////////////
//
// Neighborhood Class
//
////////////////////////////////////////////////////
class Neighborhood
{

	//////////////////////////////////////////////
	// Public member functions
	public :

		// Constructor
		Neighborhood();

		// Destructor
		~Neighborhood();		
		
		// Register a neighbor in a new list
		void NewVertex( double dist, const Vector3d& coord, int vertex );
		void NewFace( double dist, const Vector3d& coord, int face, double ratio1, double ratio2 );
		void NewEdge( double dist, const Vector3d& coord, int face, int edge, double ratio );
		
		// Add a neighbor in the list
		void AddVertex( const Vector3d& coord, int vertex );
		void AddFace( const Vector3d& coord, int face, double ratio1, double ratio2 );
		void AddEdge( const Vector3d& coord, int face, int edge, double ratio );

		// Return neighbors distance
		double Distance() const;
		
		// Return neighbors
		Neighbor* Neighbors();
		
		// Delete every registered neighbors
		void Reset();
			
	//////////////////////////////////////////////
	// Protected member function
	private :
	
		// The neighbors
		Neighbor*	neighbors;
		
		// Distance to given point
		double		distance;
};

#endif

#endif
