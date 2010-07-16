/***************************************************************************
                                    Mesh_roy.h
                             -------------------
    update               : 2003-04-02
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

#ifndef _MESH_ROY_
#define _MESH_ROY_

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "VectorT.h"
#include <assert.h>
#include <string>
#include <vector>

//--
//
// FileFormat
//
//--
enum FileFormat
{
	UNKNOWN_FILE,
	INVENTOR_FILE,
	VRML_1_FILE,
	VRML_2_FILE
};

//--
//
// Utility function
//
//--

// Convert double value in [0,1] to color
extern Vector3d Double2Color( const double& d );

//--
//
// Mesh_roy
//
//--
class Mesh_roy
{
		
	//--
	//
	// Member Data
	//
	//--
	protected :

		// Vertex array
		std::vector<Vector3d> vertices;
		
		// Face array
		std::vector<Vector3i> faces;
		
		// Color array
		std::vector<Vector3d> colors;

		// Texture coordinate array
		std::vector<Vector2d> textures;

		// Face normal array
		std::vector<Vector3d> face_normals;

		// Vertex normal array
		std::vector<Vector3d> vertex_normals;

		// Texture file name
		std::string texture_name;

	//--
	//
	// Member Functions
	//
	//--
	public :

		//--
		//
		// Constructor / Destructor
		//
		//--
		inline Mesh_roy() {
		}

		inline Mesh_roy(const Mesh_roy& m) : vertices(m.vertices), faces(m.faces),
		colors(m.colors), textures(m.textures), face_normals(m.face_normals),
		vertex_normals(m.vertex_normals), texture_name(m.texture_name) {
		}

		inline ~Mesh_roy() {
		}

		//--
		//
		// File input/output
		//
		//--
        
		// Read mesh from a file
		bool ReadFile( const std::string& file_name );
        
		// Write mesh in a file
		bool WriteFile( const std::string& file_name, const FileFormat& file_format=VRML_1_FILE ) const;
        
		//--
		//
		// Vertex Interface
		//
		//--
		
		// Vertex number
		inline int VertexNumber() const {
			return (int)vertices.size();
		}

		// Vertex array
		inline std::vector<Vector3d>& Vertices() {
			return vertices;
		}

		// Vertex array (constant)
		inline const std::vector<Vector3d>& Vertices() const {
			return vertices;
		}

		// Add vertex v in vertex array
		inline void AddVertex( const Vector3d& v ) {
			vertices.push_back( v );
		}

		// Clear vertex array
		inline void ClearVertices() {
			vertices.clear();
		}

		// Vertex #i
		inline Vector3d& Vertex(int i) {
			assert( (i>=0) && (i<VertexNumber()) );
			return vertices[i];
		}

		// Vertex #i (constant)
		inline const Vector3d& Vertex(int i) const {
			assert( (i>=0) && (i<VertexNumber()) );
			return vertices[i];
		}

		// Vertex of face #f with index #v
		// v is in range 0 to 2
		inline Vector3d& Vertex(int f, int v) {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<VertexNumber()) );
			return vertices[faces[f][v]];
		}

		// Vertex of face #f with index #v (constant)
		// v is in range 0 to 2
		inline const Vector3d& Vertex(int f, int v) const {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<VertexNumber()) );
			return vertices[faces[f][v]];
		}

		//--
		//
		// Face Interface
		//
		//--

		// Face number
		inline int FaceNumber() const {
			return (int)faces.size();
		}

		// Face array
		inline std::vector<Vector3i>& Faces() {
			return faces;
		}

		// Face array (constant)
		inline const std::vector<Vector3i>& Faces() const {
			return faces;
		}

		// Add face f to face array
		inline void AddFace( const Vector3i& f ) {
			faces.push_back( f );
		}

		// Clear face array
		inline void ClearFaces() {
			faces.clear();
		}

		// Face #i
		inline Vector3i& Face(int i) {
			assert( (i>=0) && (i<FaceNumber()) );
			return faces[i];
		}

		// Face #i (constant)
		inline const Vector3i& Face(int i) const {
			assert( (i>=0) && (i<FaceNumber()) );
			return faces[i];
		}

		// Index #v of face #f
		// v is in range 0 to 2
		inline int& Face(int f, int v) {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			return faces[f][v];
		}

		// Index #v of face #f (constant)
		// Index #v is in range 0 to 2
		inline const int& Face(int f, int v) const {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			return faces[f][v];
		}

		// Test if face #f has vertex #v
		inline bool FaceHasVertex( int f, int v ) const {
			return (Face(f,0) == v) || (Face(f,1) == v) || (Face(f,2) == v);
		}

		// Test if face is valid
		inline bool IsValidFace( int f ) const {
			assert( (f>=0) && (f<FaceNumber()) );
			return (faces[f][0]>=0) && (faces[f][1]>=0) && (faces[f][2]>=0) && (faces[f][0]!=faces[f][1]) && (faces[f][0]!=faces[f][2]) && (faces[f][1]!=faces[f][2]);
		}
		
		//--
		//
		// Color Interface
		//
		//--

		// Color number
		inline int ColorNumber() const {
			return (int)colors.size();
		}

		// Color array
		inline std::vector<Vector3d>& Colors() {
			return colors;
		}

		// Color array (constant)
		inline const std::vector<Vector3d>& Colors() const {
			return colors;
		}

		// Add color c in color array
		inline void AddColor( const Vector3d& c ) {
			colors.push_back( c );
		}

		// Clear color array
		inline void ClearColors() {
			colors.clear();
		}

		// Color #i
		inline Vector3d& Color(int i) {
			assert( (i>=0) && (i<ColorNumber()) );
			return colors[i];
		}

		// Color #i (constant)
		inline const Vector3d& Color(int i) const {
			assert( (i>=0) && (i<ColorNumber()) );
			return colors[i];
		}

		// Color of face #f with index #v
		// v is in range 0 to 2
		inline Vector3d& Color(int f, int v) {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<ColorNumber()) );
			return colors[faces[f][v]];
		}

		// Color of face #f with index #v (constant)
		// v is in range 0 to 2
		inline const Vector3d& Color(int f, int v) const {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<ColorNumber()) );
			return colors[faces[f][v]];
		}

		//--
		//
		// Texture Interface
		//
		//--

		// Texture file name
		inline std::string& TextureName() {
			return texture_name;
		}

		// Texture file name (constant)
		inline const std::string& TextureName() const {
			return texture_name;
		}

		// Clear texture file name
		inline void ClearTextureName() {
			texture_name = "";
		}

		// Texture coordinate number
		inline int TextureNumber() const {
			return (int)textures.size();
		}

		// Texture coordinate array
		inline std::vector<Vector2d>& Textures() {
			return textures;
		}

		// Texture coordinate array (constant)
		inline const std::vector<Vector2d>& Textures() const {
			return textures;
		}

		// Add texture coordinate t in texture coordinate array
		inline void AddTexture( const Vector2d& t ) {
			textures.push_back( t );
		}

		// Clear texture coordinate array
		inline void ClearTextures() {
			textures.clear();
		}

		// Texture coordinate #i
		inline Vector2d& Texture(int i) {
			assert( (i>=0) && (i<TextureNumber()) );
			return textures[i];
		}

		// Texture coordinate #i (constant)
		inline const Vector2d& Texture(int i) const {
			assert( (i>=0) && (i<TextureNumber()) );
			return textures[i];
		}

		// Texture coordinate of face #f with index #v
		// v is in range 0 to 2
		inline Vector2d& Texture(int f, int v) {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<TextureNumber()) );
			return textures[faces[f][v]];
		}

		// Texture coordinate of face #f with index #v (constant)
		// v is in range 0 to 2
		inline const Vector2d& Texture(int f, int v) const {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<TextureNumber()) );
			return textures[faces[f][v]];
		}

		//--
		//
		// Face Normal Interface
		//
		//--

		// Face normal number
		inline int FaceNormalNumber() const {
			return (int)face_normals.size();
		}

		// Face normal array
		inline std::vector<Vector3d>& FaceNormals() {
			return face_normals;
		}

		// Face normal array (constant)
		inline const std::vector<Vector3d>& FaceNormals() const {
			return face_normals;
		}

		// Add face normal n to face array
		inline void AddFaceNormal( const Vector3d& n ) {
			face_normals.push_back( n );
		}

		// Clear face array
		inline void ClearFaceNormals() {
			face_normals.clear();
		}

		// Face normal #i
		inline Vector3d& FaceNormal(int i) {
			assert( (i>=0) && (i<FaceNormalNumber()) );
			return face_normals[i];
		}

		// Face normal #i (constant)
		inline const Vector3d& FaceNormal(int i) const {
			assert( (i>=0) && (i<(int)FaceNormalNumber()) );
			return face_normals[i];
		}

		//--
		//
		// Vertex Normal Interface
		//
		//--

		// Vertex normal number
		inline int VertexNormalNumber() const {
			return (int)vertex_normals.size();
		}

		// Vertex normal array
		inline std::vector<Vector3d>& VertexNormals() {
			return vertex_normals;
		}

		// Vertex normal array (constant)
		inline const std::vector<Vector3d>& VertexNormals() const {
			return vertex_normals;
		}

		// Add vertex normal n to vertex normal array
		inline void AddVertexNormal( const Vector3d& n ) {
			vertex_normals.push_back( n );
		}

		// Clear vertex normal array
		inline void ClearVertexNormals() {
			vertex_normals.clear();
		}

		// Vertex normal #i
		inline Vector3d& VertexNormal(int i) {
			assert( (i>=0) && (i<VertexNormalNumber()) );
			return vertex_normals[i];
		}

		// Vertex normal #i (constant)
		inline const Vector3d& VertexNormal(int i) const {
			assert( (i>=0) && (i<VertexNormalNumber()) );
			return vertex_normals[i];
		}

		// Vertex normal of face #f with index #v
		// v is in range 0 to 2
		inline Vector3d& VertexNormal(int f, int v) {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<VertexNormalNumber()) );
			return vertex_normals[faces[f][v]];
		}

		// Vertex normal of face #f with index #v (constant)
		// v is in range 0 to 2
		inline const Vector3d& VertexNormal(int f, int v) const {
			assert( (f>=0) && (f<FaceNumber()) );
			assert( (v>=0) && (v<=2) );
			assert( (faces[f][v]>=0) && (faces[f][v]<VertexNormalNumber()) );
			return vertex_normals[faces[f][v]];
		}

		//--
		//
		// Normal computation
		//
		//--

		// Compute face normals
		void ComputeFaceNormals();

		// Compute vertex normals
		// Assume that face normals are computed
		void ComputeVertexNormals();

		// Compute raw normal of face #f
		inline Vector3d ComputeRawFaceNormal( int f ) const {
			return (Vertex(f,1)-Vertex(f,0)) ^ (Vertex(f,2)-Vertex(f,0));
		}

		// Compute raw normal of face {va, vb, vc}
		inline Vector3d ComputeRawFaceNormal( int va, int vb, int vc ) const {
			return (Vertex(vb)-Vertex(va)) ^ (Vertex(vc)-Vertex(va));
		}

		// Compute unit normal of face #f
		inline Vector3d ComputeFaceNormal( int f ) const {
			assert( (f>=0) && (f<(int)FaceNumber()) );
			return ComputeRawFaceNormal(f).Normalize();
		}

		// Compute unit normal of face {va, vb, vc}
		inline Vector3d ComputeFaceNormal( int va, int vb, int vc ) const {
			return ComputeRawFaceNormal(va, vb, vc).Normalize();
		}

		// Compute area of face #i
		inline double ComputeFaceArea( int i ) const {
			return 0.5 * ComputeRawFaceNormal(i).Length();
		}

		// Compute area of face {va, vb ,vc}
		inline double ComputeFaceArea( int va, int vb, int vc ) const {
			return 0.5 * ComputeRawFaceNormal(va, vb, vc).Length();
		}

		//--
		//
		// Miscellaneaous function
		//
		//--

		// Initialize mesh
		// Clear all data
		void ClearAll();
};

#endif

#endif // _MESH_

