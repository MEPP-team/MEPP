#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

/***************************************************************************
                                 Sample.cpp
                             -------------------
    update               : 2002-01-22
    copyright            : (C) 2002 by Michael ROY
    email                : m.roy@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Sample.h"
#include <iostream>
using namespace std;
//////////////////////////////////////////////////////////////////////////////
//
// Sample
//
//////////////////////////////////////////////////////////////////////////////
// Default constructor
Sample::Sample()
: height( 0 ), lw( 0 ), dev( 0 )
{
}

//////////////////////////////////////////////////////////////////////////////
//
// Sample
//
//////////////////////////////////////////////////////////////////////////////
// Destructor
Sample::~Sample()
{
	delete [] lw;
	for ( int i = 0; i < height; i++ )
		delete [] dev[i];
	delete [] dev;
}

//////////////////////////////////////////////////////////////////////////////
//
// SetHeight
//
//////////////////////////////////////////////////////////////////////////////
// Define scan line number
// Allocate memory for Active Edge List
void Sample::SetHeight( int Height )
{
	// Debug test
	#ifndef NDEBUG
		if ( !Height )
			cerr<<"Sample Info: null height"<<endl;
	#endif
	// Size
	height = Height;
	// Allocate line width array
	lw = new int [height];
}

void Sample::InitDev()
{
	// Allocate deviation array
	dev = new double* [height];
	for ( register int i=0; i<height; i++ )
	{
		dev[i] = new double [lw[i]];
	}
}
#endif
