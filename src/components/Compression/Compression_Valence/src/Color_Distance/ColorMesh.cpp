#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

/***************************************************************************
                                  Mesh_roy.cpp
                             -------------------
    update               : 2003-04-22
    copyright            : (C) 2002-2003 by Micha? Roy
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

#include "ColorMesh.h"
//#include "FileVrml1.h"
//#include "FileVrml2.h"
#include <fstream>
#include <iostream>

//--
//
// ClearAll
//
//--
void Mesh_roy::ClearAll()
{
	// Initialize mesh
	// Clear all data
	ClearVertices();
	ClearFaces();
	ClearColors();
	ClearTextures();
	ClearFaceNormals();
	ClearVertexNormals();
	ClearTextureName();
}

//--
//
// ComputeFaceNormals
//
//--
// Compute unit normal of every faces
void Mesh_roy::ComputeFaceNormals()
{
	// Resize face normal array
	face_normals.resize( FaceNumber() );

	// For every face
	for ( int i=0; i<FaceNumber(); i++ )
	{
		// Compute unit face normal
		FaceNormal(i) = ComputeFaceNormal(i);
	}
}

//--
//
// ComputeVertexNormals
//
//--
// Compute unit normal of every vertex
void Mesh_roy::ComputeVertexNormals()
{
	int i;

	// Assume that face normals are computed
	assert( FaceNormalNumber() == FaceNumber() );

	// Resize and initialize vertex normal array
	vertex_normals.assign( VertexNumber(), Vector3d(0,0,0) );
	
	// For every face
	for ( i=0 ; i<FaceNumber() ; i++ )
	{
		// Add face normal to vertex normal
		VertexNormal(i,0) += FaceNormal(i);
		VertexNormal(i,1) += FaceNormal(i);
		VertexNormal(i,2) += FaceNormal(i);
	}
	
	// For every vertex
	for ( i=0 ; i<VertexNumber() ; i++)
	{
		// Normalize vertex normal
		VertexNormal(i).Normalize();
	}
}

//--
//
// UpperCase
//
//--
// Upper case string of a given one
static std::string UpperCase( const std::string& s )
{
	// Upper case string
	std::string us;
	
	// For every character in the string
	std::string::const_iterator it = s.begin();
	while( it != s.end() )
	{
		// Convert character to upper case
		us += toupper( *it );
		
		// Next character
		++it;
	}
	
	// Return upper case string
	return us;
}

//--
//
// ReadFile
//
//--
// Read mesh from a file
/*
bool Mesh_roy::ReadFile( const std::string& file_name )
{
	FileFormat file_format(UNKNOWN_FILE);

	// Find file extension
	int pos = file_name.find_last_of('.');
	if ( pos == -1 )
	{
		// File has no extension
		return false;
	}
	
	// Format extension string
	std::string extension = UpperCase( file_name.substr( ++pos ) );
	
	// WRL extension
	if ( extension == "WRL" )
	{
		//
		// Check VRML file version
		// 
		std::ifstream file(file_name.c_str());
		if ( file.is_open() == false ) return false;
		std::string word;
		file>>word;
		if ( word != "#VRML" )
		{
			file.close();
			return false;
		}
		file>>word;
		file.close();
		if ( word == "V1.0" )
		{
			file_format = VRML_1_FILE;
		}
		else if ( word == "V2.0" )
		{
			file_format = VRML_2_FILE;
		}
		else
		{
			return false;
		}
	}
	// IV extension
	else if ( extension == "IV" )
	{
		file_format = INVENTOR_FILE;
	}
	// Other extension
	else
	{
		// Unknown file format
		return false;
	}
	
	// Read file
	switch( file_format )
	{
		// OpenInventor file
		case INVENTOR_FILE :
		
		// VRML 1.0 file
		case VRML_1_FILE :
			return ReadVrml1( *this, file_name );

		// VRML 2.0 file
		case VRML_2_FILE :
			return ReadVrml2( *this, file_name );

		// Other file
		default :
			// Unknown file format
			return false;
	}
}
*/

//--
//
// WriteFile
//
//--
// Write mesh to a file
/*
bool Mesh_roy::WriteFile( const std::string& file_name, const FileFormat& file_format ) const
{
	// Write file
	switch( file_format )
	{
		// OpenInventor file
		case INVENTOR_FILE :
			return WriteVrml1( *this, file_name, true );

		// VRML 1.0 file
		case VRML_1_FILE :
			return WriteVrml1( *this, file_name );

		// VRML 2.0 file
		case VRML_2_FILE :
			return WriteVrml2( *this, file_name );

		// Other file
		default :
			// Unknown file format
			return false;
	}
}
*/
//--
//
// Double2Color
//
//--
// Convert double value in [0,1] to color
Vector3d Double2Color( const double& d )
{
	assert( (d>=0.0) && (d<=1.0) );
	if ( d < 0.0  ) return Vector3d(0,0,1);
	if ( d < 0.25 ) return Vector3d(               0,                d * 4.0,                    1 );
	if ( d < 0.50 ) return Vector3d(               0,                      1, 1 - (d - 0.25) * 4.0 );
	if ( d < 0.75 ) return Vector3d( (d - 0.5) * 4.0,                      1,                    0 );
	if ( d < 1.0  ) return Vector3d(               1, 1.0 - (d - 0.75) * 4.0,                    0 );
	return Vector3d(1,0,0);
}
#endif
