/***************************************************************************
                                  Sample.h
                             -------------------
    update               : 2002-02-22
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

#ifndef _SAMPLE_
#define _SAMPLE_

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

//#include "MeshDev.h"

//===============================================================
//
// Sample class
//
//===============================================================
// Register all informations for the sampling of a face
class Sample 
{
	public :
	
		// Default constructor
		Sample();
		// Destructor
		~Sample();
		int  operator []( int i ) const { return lw[i]; }
		int& operator []( int i ) { return lw[i]; }
		double  operator ()( int i, int j ) const { return dev[j][i]; }
		double& operator ()( int i, int j ) { return dev[j][i]; }
		void SetHeight( int Height );
		int Height() const { return height; }
		void InitDev();

		
	private :

		int		height;
		int*	    lw;
		double**	dev;
};

#endif

#endif
