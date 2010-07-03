#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

/***************************************************************************
                                PlaneMesh.cxx
                             -------------------
    update               : 2002-09-25
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

#include "PlaneMesh.h"

void PlaneMesh::ClearAll()
{
	ClearFacePlanes();
	Mesh_roy::ClearAll();
}

void PlaneMesh::ComputeFacePlanes()
{
	// Assume that face normals are computed and normalized
 	assert( FaceNormalNumber() == FaceNumber() );
 	face_planes.resize( FaceNumber() );
	for ( int i=0; i<FaceNumber(); i++ )
	{
		FacePlane(i) = - FaceNormal(i) | Vertex(i, 0);
	}
}
#endif
