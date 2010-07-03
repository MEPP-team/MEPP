/***************************************************************************
                                  VectorT.h
                             -------------------
    update               : 2003-02-07
    copyright            : (C) 2002-2003 by Michaël Roy
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

#ifndef _VECTORT_
#define _VECTORT_

#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

#include <assert.h>
#include <math.h>
#include <algorithm>
#include <iostream>

#include <string.h>

template<typename Type, int Size>
class VectorT
{
	public :

		static inline int Dimension() { return Size; }

		inline VectorT() {
		}

		inline VectorT(const VectorT<Type,Size>& v) {
			memcpy(values, v.values, Size*sizeof(Type));
		}

		explicit inline VectorT(const Type val[Size]) {
			memcpy(values, val, Size*sizeof(Type));
		}

		explicit inline VectorT(const Type& v) {
			for ( int i=0; i<Size; i++ ) values[i] = v;
		}

		explicit inline VectorT(const Type& v0,const Type& v1) {
			assert(Size >= 2);
			values[0] = v0; values[1] = v1;
		}

		explicit inline VectorT(const Type& v0,const Type& v1,const Type& v2) {
			assert(Size >= 3);
			values[0]=v0; values[1]=v1; values[2]=v2;
		}

		explicit inline VectorT(const Type& v0,const Type& v1,const Type& v2, const Type& v3) {
			assert(Size >= 4);
			values[0]=v0; values[1]=v1; values[2]=v2; values[3]=v3;
		}

		inline operator Type*() { return values; }

		inline operator const Type*() const { return values; }

		inline Type& operator[](int i) {
			assert((i>=0) && (i<Size));
			return values[i]; 
		}

		inline const Type& operator[](int i) const {
			assert((i>=0) && (i<Size));
			return values[i]; 
		}

		inline VectorT<Type,Size>& operator=(const VectorT<Type,Size>& v) {
			memcpy(values, v.values, Size*sizeof(Type)); 
			return *this;
		}
		
		inline VectorT<Type,Size>& operator=(const Type& s) {
			for ( int i=0; i<Size; i++ ) values[i] = s; 
			return *this;
		}

		inline bool operator==(const VectorT<Type,Size>& v) const {
			for (int i=0; i<Size; i++) if (values[i] != v.values[i]) return false;
			return true; 
		}
		
		inline bool operator==(const Type& s) const {
			for (int i=0; i<Size; i++) if (values[i] != s) return false;
			return true; 
		}

		inline bool operator!=(const VectorT<Type,Size>& v) const {
			for (int i=0; i<Size; i++) if (values[i] != v.values[i]) return true;
			return false; 
		}
		
		inline bool operator!=(const Type& s) const {
			for (int i=0; i<Size; i++) if (values[i] != s) return true;
			return false; 
		}

		inline Type operator|(const VectorT<Type,Size>& v) const {
			Type p(values[0]*v.values[0]);
			for (int i=1; i<Size; i++) p += values[i]*v.values[i]; 
			return p; 
		}

		inline const VectorT<Type,Size>& operator*=(const Type& s) {
			for (int i=0; i<Size; i++) values[i] *= s; 
			return *this;
		}

		inline const VectorT<Type,Size>& operator*=(const VectorT<Type,Size>& v) {
			for (int i=0; i<Size; i++) values[i] *= v[i]; 
			return *this;
		}

		inline VectorT<Type,Size> operator*(const Type& s) const {
			return VectorT<Type,Size>(*this) *= s;
		}

		inline VectorT<Type,Size> operator*(const VectorT<Type,Size>& v) const {
			return VectorT<Type,Size>(*this) *= v;
		}

		inline const VectorT<Type,Size>& operator/=(const Type &s) {
			for (int i=0; i<Size; i++) values[i] /= s; 
			return *this;
		}

		inline const VectorT<Type,Size>& operator/=(const VectorT<Type,Size>& v) {
			for (int i=0; i<Size; i++) values[i] /= v[i]; 
			return *this;
		}

		inline VectorT<Type,Size> operator/(const Type &s) const {
			return VectorT<Type,Size>(*this) /= s;
		}

		inline VectorT<Type,Size> operator/(const VectorT<Type,Size>& v) const {
			return VectorT<Type,Size>(*this) /= v;
		}

		inline VectorT<Type,Size>& operator-=(const VectorT<Type,Size>& v) {
			for (int i=0; i<Size; i++) values[i] -= v.values[i]; 
			return *this;
		}

		inline VectorT<Type,Size>& operator-=(const Type& s) {
			for (int i=0; i<Size; i++) values[i] -= s; 
			return *this;
		}

		inline VectorT<Type,Size> operator-(const VectorT<Type,Size>& v) const {
			return VectorT<Type,Size>(*this) -= v;
		}

		inline VectorT<Type,Size> operator-(const Type& s) const {
			return VectorT<Type,Size>(*this) -= s;
		}

		inline VectorT<Type,Size>& operator+=(const VectorT<Type,Size>& v) {
			for (int i=0; i<Size; i++) values[i] += v.values[i]; 
			return *this; 
		}

		inline VectorT<Type,Size>& operator+=(const Type& s) {
			for (int i=0; i<Size; i++) values[i] += s; 
			return *this; 
		}

		inline VectorT<Type,Size> operator+(const Type& s) const {
			return VectorT<Type,Size>(*this) += s; 
		}

		inline VectorT<Type,Size> operator+(const VectorT<Type,Size>& v) const {
			return VectorT<Type,Size>(*this) += v;
		}

		inline VectorT<Type,Size> operator-(void) const {
			VectorT<Type,Size> v;
			for (int i=0; i<Size; i++) v.values[i] = -values[i]; 
			return v; 
		}

		inline VectorT<Type,Size> operator^(const VectorT<Type,Size> &v) const;

		inline double Length() const {
			return sqrt((double)SquareLength());
		}

		inline Type SquareLength() const {
			return *this | *this;
		}

		inline VectorT<Type,Size>& Normalize() {
			double n = Length();
			if (n != 0.0) *this *= 1.0/n;
			return *this;
		}

		inline Type Max() const {
			Type m(values[0]);
			for (int i=1; i<Size; i++) if (values[i]>m) m=values[i];
			return m; 
		}

		inline Type Min() const {
			Type m(values[0]); 
			for (int i=1; i<Size; i++) if (values[i]<m) m=values[i];
			return m; 
		}

		inline Type Mean() const {
			Type m(values[0]); 
			for (int i=1; i<Size; i++) m+=values[i];
			return m/Type(Size); 
		}

		inline VectorT<Type,Size> Min(const VectorT<Type,Size>& v) const {
			VectorT<Type,Size> res;
			for (int i=0; i<Size; i++)	res[i] = std::min(values[i], v[i]);
			return res;
		}

		inline VectorT<Type,Size> Max(const VectorT<Type,Size>& v) const {
			VectorT<Type,Size> res;
			for (int i=0; i<Size; i++)	res[i] = std::max(values[i], v[i]);
			return res;
		}

		inline VectorT<Type,Size>& Minimize(const VectorT<Type,Size>& v) {
			for (int i=0; i<Size; i++)	if (v[i]<values[i]) values[i] = v[i];
			return *this;
		}

		inline VectorT<Type,Size>& Maximize(const VectorT<Type,Size>& v) {
			for (int i=0; i<Size; i++)	if (v[i]>values[i]) values[i] = v[i];
			return *this;
		}

		template<typename func>
		inline VectorT<Type,Size> Apply(const func& f) const {
			VectorT<Type,Size> result;
			for (int i=0; i<Size; i++) result[i] = f(values[i]);
			return result; 
		}

		static inline VectorT<Type,Size> Vectorize(const Type& s) {
			VectorT<Type,Size> result;
			for (int i=0; i<Size; i++) result[i] = s;
			return result; 
		}

		inline VectorT<Type,Size>& Invert() {
			for (int i=0; i<Size; i++) values[i] = 1.0 / values[i];
			return *this;
		}
		
		inline VectorT<Type,Size> Lerp(const VectorT<Type,Size>& v, const double& t) const {
			VectorT<Type,Size> result(*this);
			for (int i=0; i<Size; i++) result[i] += t * (v[i] - values[i]);
			return result;
		}
			
		inline double Cotan(const VectorT<Type,Size>& v) const
		{
			// Cotangent between two vector in nD
			// Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
			// Mark Meyer, Mathieu Desbrun, Peter Schroder, Alan H. Barr
			// VisMath '02, Berlin (Germany)
			double dot_prod = *this | v;
			double denom = SquareLength()*v.SquareLength() - dot_prod*dot_prod;
			// denom can be 0 if *this is colinear to v
			assert( denom > 0 );
			return dot_prod / sqrt(denom);
		}
		
		// Angle From Cotan
		inline double AngleFromCotan(const VectorT<Type,Size>& v) const {
			double udotv = *this | v;
			double denom = SquareLength()*v.SquareLength() - udotv*udotv;
			// denom can be 0 if *this is colinear to v
			assert( denom != 0 );
			return (fabs (atan2 (sqrt(denom), udotv)));
		}

		inline bool IsColinear(const VectorT<Type,Size>& v) const
		{
			double dot_prod = *this | v;
			return (SquareLength()*v.SquareLength() - dot_prod*dot_prod) <= 0;
	//		return ((*this)^v).SquareLength() == 0;
		}
		
		// Obtuse angle
		inline bool ObtuseAngle(const VectorT<Type,Size>& v) const {
			return (*this|v) < 0.0;
		}


	protected :

		Type values[Size];
};

//
// Specializations
//
template<>
inline VectorT<float,3>
VectorT<float,3>::operator^(const VectorT<float,3>& v) const 
{
	return VectorT<float,3>(
		values[1]*v[2]-values[2]*v[1],
		values[2]*v[0]-values[0]*v[2],
		values[0]*v[1]-values[1]*v[0] );
}
  
template<>
inline VectorT<double,3>
VectorT<double,3>::operator^(const VectorT<double,3>& v) const
{
	return VectorT<double,3>(
		values[1]*v[2]-values[2]*v[1],
		values[2]*v[0]-values[0]*v[2],
		values[0]*v[1]-values[1]*v[0]);
}

template<typename Type,int Size>
inline std::istream& operator>>(std::istream& is, VectorT<Type,Size>& vec)
{
	for (int i=0; i<Size; i++) is >> vec[i];
	return is;
}

template<typename Type,int Size>
inline std::ostream& operator<<(std::ostream& os, const VectorT<Type,Size>& vec)
{
	for (int i=0; i<Size-1; i++) os<<vec[i]<<" ";
	os<<vec[Size-1];
	return os;
}

template<typename Type,int Size>
inline VectorT<Type,Size> operator*(Type s, const VectorT<Type,Size>& v )
{
	return v * s;
}

template<typename Type,int Size>
inline VectorT<Type,Size> operator/(Type s, const VectorT<Type,Size>& v )
{
	return v / s;
}

template<typename Type,int Size>
inline VectorT<Type,Size> operator+(Type s, const VectorT<Type,Size>& v )
{
	return v + s;
}

template<typename Type,int Size>
inline VectorT<Type,Size> operator-(Type s, const VectorT<Type,Size>& v )
{
	return v - s;
}

template<typename Type, int Size>
inline Type Dot(const VectorT<Type,Size>& v1, const VectorT<Type,Size>& v2)
{
	return v1 | v2; 
}

template<typename Type, int Size>
inline VectorT<Type,Size> Cross(const VectorT<Type,Size>& v1, const VectorT<Type,Size>& v2)
{
	return v1 ^ v2;
}

template<typename Type, int Size>
inline VectorT<Type,Size> Lerp(const VectorT<Type,Size>& v1, const VectorT<Type,Size>& v2, const double& t)
{
	return v1.Lerp(v2, t);
}

template<typename Type, int Size>
inline double Cotan(const VectorT<Type,Size>& v1, const VectorT<Type,Size>& v2)
{
	return v1.Cotan(v2);
}

template<typename Type, int Size>
inline double Cotan(const VectorT<Type,Size>& vo, const VectorT<Type,Size>& va, const VectorT<Type,Size>& vb)
{
	return Cotan(va-vo, vb-vo);
}

// Angle From Cotan
template<typename Type, int Size>
inline double AngleFromCotan(const VectorT<Type,Size>& u, const VectorT<Type,Size>& v)
{
	return u.AngleFromCotan(v);
}

// Angle From Cotan
template<typename Type, int Size>
inline double AngleFromCotan(const VectorT<Type,Size>& vo, const VectorT<Type,Size>& va, const VectorT<Type,Size>& vb)
{
	return AngleFromCotan( va-vo, vb-vo );
}

// Obtuse angle
template<typename Type, int Size>
inline bool ObtuseAngle(const VectorT<Type,Size>& v1, const VectorT<Type,Size>& v2)
{
	return v1.ObtuseAngle(v2);
}

// ObtuseAngle
template<typename Type, int Size>
inline bool ObtuseAngle(const VectorT<Type,Size>& vo, const VectorT<Type,Size>& va, const VectorT<Type,Size>& vb)
{
	return ObtuseAngle( va-vo, vb-vo );
}

template<typename Type, int Size>
inline bool IsColinear(const VectorT<Type,Size>& v1, const VectorT<Type,Size>& v2)
{
	return v1.IsColinear(v2);
}

template<typename Type, int Size>
inline bool IsColinear(const VectorT<Type,Size>& vo, const VectorT<Type,Size>& va, const VectorT<Type,Size>& vb)
{
	return IsColinear(va-vo, vb-vo);
}

template<typename Type, int Size>
inline double Length(const VectorT<Type,Size>& v)
{
	return v.Length();
}

template<typename Type, int Size>
inline VectorT<Type,Size>& Normalize(VectorT<Type,Size>& v)
{
	return v.Normalize();
}

template<typename Type, int Size>
inline VectorT<Type,Size>& Invert(VectorT<Type,Size>& v)
{
	return v.Invert();
}


typedef VectorT<signed char,1>        Vector1c;
typedef VectorT<unsigned char,1>      Vector1uc;
typedef VectorT<signed short int,1>   Vector1s;
typedef VectorT<unsigned short int,1> Vector1us;
typedef VectorT<signed int,1>         Vector1i;
typedef VectorT<unsigned int,1>       Vector1ui;
typedef VectorT<float,1>              Vector1f;
typedef VectorT<double,1>             Vector1d;

typedef VectorT<signed char,2>        Vector2c;
typedef VectorT<unsigned char,2>      Vector2uc;
typedef VectorT<signed short int,2>   Vector2s;
typedef VectorT<unsigned short int,2> Vector2us;
typedef VectorT<signed int,2>         Vector2i;
typedef VectorT<unsigned int,2>       Vector2ui;
typedef VectorT<float,2>              Vector2f;
typedef VectorT<double,2>             Vector2d;

typedef VectorT<signed char,3>        Vector3c;
typedef VectorT<unsigned char,3>      Vector3uc;
typedef VectorT<signed short int,3>   Vector3s;
typedef VectorT<unsigned short int,3> Vector3us;
typedef VectorT<signed int,3>         Vector3i;
typedef VectorT<unsigned int,3>       Vector3ui;
typedef VectorT<float,3>              Vector3f;
typedef VectorT<double,3>             Vector3d;

typedef VectorT<signed char,4>        Vector4c;
typedef VectorT<unsigned char,4>      Vector4uc;
typedef VectorT<signed short int,4>   Vector4s;
typedef VectorT<unsigned short int,4> Vector4us;
typedef VectorT<signed int,4>         Vector4i;
typedef VectorT<unsigned int,4>       Vector4ui;
typedef VectorT<float,4>              Vector4f;
typedef VectorT<double,4>             Vector4d;

#endif

#endif // _VECTORT_

