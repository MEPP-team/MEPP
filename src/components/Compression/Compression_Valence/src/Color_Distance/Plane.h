/***************************************************************************
                                   Plane_roy.h
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

#ifndef _PLANE_
#define _PLANE_

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "VectorT.h"

struct Plane_roy
{
	Vector3d n;
	double   h;

	Plane_roy() : n(0,0,0), h(0) {}
	Plane_roy( const Plane_roy& p ) : n( p.n ), h( p.h ) {}
	Plane_roy( const Vector3d& n, const double& h ) : n( n ), h( h ) {}
};

#endif

#endif // _PLANE_
