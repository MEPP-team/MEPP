#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

/***************************************************************************
                              Neighborhood.cxx
                             -------------------
    update               : 2002-02-07
    copyright            : (C) 2002 by Michael ROY
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
#include "Neighborhood.h"
#include <float.h>


//////////////////////////////////////////////////////////////////////////////
//
// Neighborhood
//
//////////////////////////////////////////////////////////////////////////////
// Default constructor
Neighborhood::Neighborhood()
: neighbors( 0 ), distance( DBL_MAX )
{
}

//////////////////////////////////////////////////////////////////////////////
//
// ~Neighborhood
//
//////////////////////////////////////////////////////////////////////////////
// Destructor
Neighborhood::~Neighborhood()
{
	Reset();
}

//////////////////////////////////////////////////////////////////////////////
//
// Reset
//
//////////////////////////////////////////////////////////////////////////////
// Delete every registered neighbors
void Neighborhood::Reset()
{
	if ( neighbors )
	{
		Neighbor* temp;
		while( neighbors )
		{
			temp = neighbors;
			neighbors = neighbors->next;
			delete temp;
		}
		neighbors = 0;
	}
	distance = DBL_MAX;
}
	



//////////////////////////////////////////////////////////////////////////////
//
// New
//
//////////////////////////////////////////////////////////////////////////////
// Delete previous registeres neighbors
// Create a new list with the given neighbor
// Initialize neighbors distance
void Neighborhood::NewVertex( double dist, const Vector3d& coord, int vertex )
{
	Reset();
	neighbors = new Neighbor;
	neighbors->c = coord;
	neighbors->v = vertex;
	neighbors->next = 0;
	distance = dist;
}

void Neighborhood::NewFace( double dist, const Vector3d& coord, int face, double ratio1, double ratio2 )
{
	Reset();
	neighbors = new Neighbor;
	neighbors->c = coord;
	neighbors->v = -1;
	neighbors->f = face;
	neighbors->e = -1;
	neighbors->r1 = ratio1;
	neighbors->r2 = ratio2;
	neighbors->next = 0;
	distance = dist;
}

void Neighborhood::NewEdge( double dist, const Vector3d& coord, int face, int edge, double ratio )
{
	Reset();
	neighbors = new Neighbor;
	neighbors->c = coord;
	neighbors->v = -1;
	neighbors->f = face;
	neighbors->e = edge;
	neighbors->r1 = ratio;
	neighbors->next = 0;
	distance = dist;
}

//////////////////////////////////////////////////////////////////////////////
//
// Add
//
//////////////////////////////////////////////////////////////////////////////
// Add a new neighbor in the list
// If this neighbors is not already registered
void Neighborhood::AddVertex( const Vector3d& coord, int vertex )
{
	Neighbor* temp = neighbors;
	if ( temp->c == coord ) return; // Already registered ?
	while( temp->next )
	{
		temp = temp->next;
		if ( temp->c == coord ) return; // Already registered ?
	}
	temp->next = new Neighbor; // Add neighbor to the end of list
	temp = temp->next;
	temp->c = coord;
	temp->v = vertex;
	temp->next = 0;
}

void Neighborhood::AddFace( const Vector3d& coord, int face, double ratio1, double ratio2 )
{
	Neighbor* temp = neighbors;
	if ( temp->c == coord ) return; // Already registered ?
	while( temp->next )
	{
		temp = temp->next;
		if ( temp->c == coord ) return; // Already registered ?
	}
	temp->next = new Neighbor; // Add neighbor to the end of list
	temp = temp->next;
	temp->c = coord;
	temp->v = -1;
	temp->f = face;
	temp->e = -1;
	temp->r1 = ratio1;
	temp->r2 = ratio2;
	temp->next = 0;
}

void Neighborhood::AddEdge( const Vector3d& coord, int face, int edge, double ratio )
{
	Neighbor* temp = neighbors;
	if ( temp->c == coord ) return; // Already registered ?
	while( temp->next )
	{
		temp = temp->next;
		if ( temp->c == coord ) return; // Already registered ?
	}
	temp->next = new Neighbor; // Add neighbor to the end of list
	temp = temp->next;
	temp->c = coord;
	temp->v = -1;
	temp->f = face;
	temp->e = edge;
	temp->r1 = ratio;
	temp->next = 0;
}
	
//////////////////////////////////////////////////////////////////////////////
//
// distance
//
//////////////////////////////////////////////////////////////////////////////
// Return neighbors distance
double Neighborhood::Distance() const
{
	return( distance );
}

//////////////////////////////////////////////////////////////////////////////
//
// neighbors
//
//////////////////////////////////////////////////////////////////////////////
// Return neighbors list
Neighbor* Neighborhood::Neighbors()
{
	return( neighbors );
}
#endif
