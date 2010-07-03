#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

/***************************************************************************
                              UniformGrid.cxx
                             -------------------
    update               : 2003-02-16
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

#include "UniformGrid.h"


static Cell3D* pCell;



//////////////////////////////////////////////////////////////////////
//
// UniformGrid
//
//////////////////////////////////////////////////////////////////////
// Default constructor
// pC   -> pointer to vertices
// bbox -> vertices bounding box
// size -> cells size
UniformGrid::UniformGrid( Mesh_roy* m, BoundingBox3d* bbox, double dim )
: mv(m->Vertices()), mf(m->Faces()), mfn(m->FaceNormals())
{
	int i, j, k, l;

	// Assume that face normals are computed and normalized
 	assert( m->FaceNormalNumber() == m->FaceNumber() );
 	mp.resize( m->FaceNumber() );
	for ( i=0; i<m->FaceNumber(); i++ )
	{
		mp[i] = - mfn[i] | mv[mf[i][0]];
	}

	m_pMin = bbox->Min();

	m_rSize = dim;
/*
	m_pCellNum[0] = (int)( bbox->Size()[0] / m_rSize ) + 1;
	m_pCellNum[1] = (int)( bbox->Size()[1] / m_rSize ) + 1;
	m_pCellNum[2] = (int)( bbox->Size()[2] / m_rSize ) + 1;
*/
	m_pCellNum = (int)( bbox->Diagonal() / m_rSize + 1.0 );
	
	m_pCell = new Cell3D*** [m_pCellNum[0]];
	for ( i=0; i<m_pCellNum[0]; i++ )
	{
		m_pCell[i] = new Cell3D** [m_pCellNum[1]];
		for ( j=0; j<m_pCellNum[1]; j++ )
		{
			m_pCell[i][j] = new Cell3D* [m_pCellNum[2]];
			for ( k=0; k<m_pCellNum[2]; k++ ) m_pCell[i][j][k] = 0;
		}
	}

	for ( i = 0; i < (int)mv.size(); i++ )
	{
		//////////////////////////////////////////////////////
		// Compute cell position
		j = (int)( (mv[i][0] - m_pMin[0]) / m_rSize );
		k = (int)( (mv[i][1] - m_pMin[1]) / m_rSize );
		l = (int)( (mv[i][2] - m_pMin[2]) / m_rSize );
		/////////////////////////////////////////////////////
		// Register point
		AddOnePoint( i, j, k, l );
	}
	
	SetFaces();

	_FacesTested = 0;
	neighbors = 0;

}

//////////////////////////////////////////////////////////////////////
//
// ~UniformGrid
//
//////////////////////////////////////////////////////////////////////
// Default destructor
UniformGrid::~UniformGrid()
{
	if ( m_pCell )
	{
		register int i, j, k;
		// Free memory allocated for all cells
		for ( i=0; i<m_pCellNum[0]; i++ )
		{
			for ( j=0; j<m_pCellNum[1]; j++ )
			{
				for ( k=0; k<m_pCellNum[2]; k++ )
				{
					while( m_pCell[i][j][k] )
					{
						pCell = m_pCell[i][j][k];
						m_pCell[i][j][k] = pCell->next;
						delete pCell;
					}
				}
				delete [] m_pCell[i][j];
			}
			delete [] m_pCell[i];
		}
		delete [] m_pCell;
	}
	
	delete neighbors;
}

//////////////////////////////////////////////////////////////////////
//
// AddOnePoint
//
//////////////////////////////////////////////////////////////////////
// Add given point in cell
void UniformGrid::AddOnePoint(int n, int x, int y, int z)
{
	//////////////////////////////////////////////////////
	// Is there an existing registered point
	if ( !m_pCell[x][y][z] )
	{
		m_pCell[x][y][z] = new Cell3D;
		m_pCell[x][y][z]->v = n;
		m_pCell[x][y][z]->f = -1;
		m_pCell[x][y][z]->next = 0;
	}
	else
	{
		pCell = m_pCell[x][y][z];
		while( pCell->next ) pCell = pCell->next;
		pCell->next = new Cell3D;
		pCell->next->v = n;
		pCell->next->f = -1;
		pCell->next->next = 0;
	}
}

#define TestRegisterFace() \
if ( s == 1 ) { if ( DistancePoint2Plane( p, mfn[i], mp[i] ) >= 0 ) \
	{ AddOneFace( i, xx, yy, zz ); continue; } } \
else { if ( DistancePoint2Plane( p, mfn[i], mp[i] ) < 0 ) \
	{ AddOneFace( i, xx, yy, zz ); continue; } }

void UniformGrid::SetFaces()
{
	register int x1, x2, y1, y2, z1, z2, xx, yy, zz;
	register int a ,b, c, i;
	register char s;
	Vector3d p;

	/////////////////////////////////////
	// Begin to work with vertices
	for ( i=0; i<(int)mf.size(); i++ )
	{
		/////////////////////////////////////////
		// Set vertices index
		a = mf[i][0];
		b = mf[i][1];
		c = mf[i][2];
		//////////////////////////////////////////////////////
		// Compute cell position
		x1 = x2 = (int)( (mv[a][0] - m_pMin[0]) / m_rSize );
		y1 = y2 = (int)( (mv[a][1] - m_pMin[1]) / m_rSize );
		z1 = z2 = (int)( (mv[a][2] - m_pMin[2]) / m_rSize );
		//////////////////////////////////////////////////////
		// Compute cell position
		xx = (int)( (mv[b][0] - m_pMin[0]) / m_rSize );
		yy = (int)( (mv[b][1] - m_pMin[1]) / m_rSize );
		zz = (int)( (mv[b][2] - m_pMin[2]) / m_rSize );
		// Check for x
		if ( xx < x1 ) x1 = xx;
		else 
		if ( xx > x2 ) x2 = xx;
		// Check for y
		if ( yy < y1 ) y1 = yy;
		else
		if ( yy > y2 ) y2 = yy;
		// Check for z
		if ( zz < z1 ) z1 = zz;
		else
		if ( zz > z2 ) z2 = zz;
		//////////////////////////////////////////////////////
		// Compute cell position
		xx = (int)( (mv[c][0] - m_pMin[0]) / m_rSize );
		yy = (int)( (mv[c][1] - m_pMin[1]) / m_rSize );
		zz = (int)( (mv[c][2] - m_pMin[2]) / m_rSize );
		// Check for x
		if ( xx < x1 ) x1 = xx;
		else if ( xx > x2 ) x2 = xx;
		// Check for y
		if ( yy < y1 ) y1 = yy;
		else if ( yy > y2 ) y2 = yy;
		// Check for z
		if ( zz < z1 ) z1 = zz;
		else if ( zz > z2 ) z2 = zz;
		/////////////////////////////////////////////////////:
		// Compute intersection Plane-Cube
		for ( xx=x1; xx<=x2; xx++ )
		{
			for ( yy=y1; yy<=y2; yy++ )
			{
				for ( zz=z1; zz<=z2; zz++ )
				{
					////////////////////////////////////////
					p[0] = m_pMin[0] + xx*m_rSize;
					p[1] = m_pMin[1] + yy*m_rSize;
					p[2] = m_pMin[2] + zz*m_rSize;
					if ( DistancePoint2Plane( p, mfn[i], mp[i] ) < 0 ) s = 1;
					else s = 2;
					//////////////////////////////////
					p[0] += m_rSize;
					TestRegisterFace()
					//////////////////////////////////
					p[2] += m_rSize;
					TestRegisterFace()
					//////////////////////////////////
					p[0] -= m_rSize;
					TestRegisterFace()
					//////////////////////////////////
					p[1] += m_rSize;
					p[2] -= m_rSize;
					TestRegisterFace()
					//////////////////////////////////
					p[0] += m_rSize;
					TestRegisterFace()
					//////////////////////////////////
					p[2] += m_rSize;
					TestRegisterFace()
					//////////////////////////////////
					p[1] -= m_rSize;
					TestRegisterFace()
				}
			}
		}
	}

}

void UniformGrid::AddOneFace(int n, int x, int y, int z)
{
	//////////////////////////////////////////////////////
	// Is there an existing registered point
	if ( !m_pCell[x][y][z] )
	{
		m_pCell[x][y][z] = new Cell3D;
		m_pCell[x][y][z]->v = -1;
		m_pCell[x][y][z]->f = n;
		m_pCell[x][y][z]->next = 0;
	}
	else
	{
		pCell = m_pCell[x][y][z];
		while( (pCell->f != -1) && (pCell->next) )
			pCell = pCell->next;
		if ( pCell->f == -1 ) pCell->f = n;
		else
		{
			pCell->next = new Cell3D;
			pCell->next->v = -1;
			pCell->next->f = n;
			pCell->next->next = 0;
		}
	}
}

//////////////////////////////////////////////////////////////////////
//
// DistanceToClosestPoint
//
//////////////////////////////////////////////////////////////////////
// Compute distance from given point to closest point
Neighborhood* UniformGrid::NearestNeighbors(const Vector3d &point)
{
	static double d;
	static int i, j, k, n;
	static int xx, yy, zz;
	static int ia, ib, ja, jb, ka, kb;

	//////////////////////////////////////////////////////
	// Compute cell position
	xx = (int)( (point[0] - m_pMin[0]) / m_rSize );
	yy = (int)( (point[1] - m_pMin[1]) / m_rSize );
	zz = (int)( (point[2] - m_pMin[2]) / m_rSize );

	n = 0;
	
	delete neighbors;
	neighbors = 0;
	neighbors = new Neighborhood;
	
	//////////////////////////////////////////////////////////
	// Check for point
	//////////////////////////////////////////////////////////
	do
	{
		n++;
		
		ia = Clamp( xx-n, m_pCellNum[0] );
		ib = Clamp( xx+n, m_pCellNum[0] );
		ja = Clamp( yy-n, m_pCellNum[1] );
		jb = Clamp( yy+n, m_pCellNum[1] );
		ka = Clamp( zz-n, m_pCellNum[2] );
		kb = Clamp( zz+n, m_pCellNum[2] );

		for ( i=ia; i<=ib; i++ )
		{
			for ( j=ja; j<=jb; j++ )
			{
				for ( k=ka; k<=kb; k++ )
				{
					if ( !m_pCell[i][j][k] ) continue;
					pCell = m_pCell[i][j][k];
					do
					{
						/////////////////////////////////////////////
						// compute euclidian distance
						// between current sampled point &
						// current reference mesh point
						if ( pCell->v != -1 )
						{
							d = ( point - mv[ pCell->v ] ).Length();
							if ( d <= neighbors->Distance() )
							{
								if ( d == neighbors->Distance() )
								{
									neighbors->AddVertex( mv[ pCell->v ], pCell->v );
								}
								else
								{
									neighbors->NewVertex( d, mv[ pCell->v ], pCell->v );
									if ( !d ) return( neighbors );
								}
							}
						}

						/////////////////////////////////////////////
						// compute euclidian distance
						// between current sampled point &
						// current reference mesh point
						if ( pCell->f != -1 )
						{
							DistancePoint2Face( point, pCell->f );
							if ( !neighbors->Distance() ) return( neighbors );
						}

						pCell = pCell->next;
	
					}
					while( pCell );
				}
			}
		}
	}
	while( neighbors->Distance() > (n * m_rSize) );
	
	return( neighbors );
}
#endif
