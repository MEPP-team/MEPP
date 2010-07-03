#include "../../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Compression_Valence

/***************************************************************************
                                Deviation.cxx
                             -------------------
    update               : 2003-02-16
    copyright            : (C) 2002-2003 by Micha? ROY
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

#include "Deviation.h"
#include <float.h>
//======================================================
//
// Edge2Vertex
//
//======================================================
// Return vertex number for a given edge
// Assume conterclockwise order
static const int Edge2Vertex[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };



Deviation::Deviation()
: ma(0), mb(0), mavn(0), mbvn(0), mafn(0), mbfn(0), bb(0), ug(0), samples(0)
{
}

Deviation::~Deviation()
{
	if ( samples ) delete samples;
	if ( ug ) delete ug;
	if ( bb ) delete bb;
}


bool Deviation::Initialization( Mesh_roy* a, Mesh_roy* b, double SampleStep, double GridSize )
{
	if ( !a || !b || !a->VertexNumber() || !a->FaceNumber() || !b->VertexNumber() || !b->FaceNumber() )
	{
		cerr << "Input meshes are invalid" << endl;
		return false;
	}
	if ( samples )
	{
		delete samples;
		samples = 0;
	}
	if ( bb )
	{
		delete bb;
		bb = 0;
	}
	if ( ug )
	{
		delete ug;
		ug = 0;
	}

	ma = a;
	mb = b;

	mavn = ma->VertexNumber();
	mbvn = mb->VertexNumber();

	mafn = ma->FaceNumber();
	mbfn = mb->FaceNumber();

	// Compute face normals for second mesh
	// Needed by Uniform Grid
	mb->ComputeFaceNormals();

	bb = new BoundingBox3d;
	bb->AddPoints( ma->Vertices() );

	if ( SampleStep )
	{
		step = bb->Diagonal() * SampleStep * 0.01;
		samples = new Sample [mafn];
	}
	else
	{
		step = 0;
	}
	snum = 0;

	bb->AddPoints( mb->Vertices() );
	ug = new UniformGrid( mb, bb, bb->Diagonal() * GridSize * 0.01 );

	delete bb;
	bb = 0;

	return true;
}

bool Deviation::Compute( DeviationType type )
{
	if ( !ma || !mb )
	{
		cerr<<"No meshes"<<endl;
		return false;
	}

	switch( type )
	{
		//
		// Geometric
		//
		case GEOMETRIC_DEVIATION :
			if ( step ) 
				return GeometricDeviationSample();
			else
				return GeometricDeviation();
			break;

		//
		// Normal
		//
		case NORMAL_DEVIATION :
			if ( step ) 
				return MeshDeviationSample( ma->VertexNormals(), mb->VertexNormals() );
			else
				return MeshDeviation( ma->VertexNormals(), mb->VertexNormals() );
			break;

		//
		// Color
		//
		case COLOR_DEVIATION :
			if ( step ) 
				return MeshDeviationSample( ma->Colors(), mb->Colors() );
			else
				return MeshDeviation( ma->Colors(), mb->Colors() );
			break;

		//
		// Texture
		//
		case TEXTURE_DEVIATION :
			if ( step ) 
				return MeshDeviationSample( ma->Textures(), mb->Textures() );
			else
				return MeshDeviation( ma->Textures(), mb->Textures() );
			break;

		//
		// Others
		//
		default :
			cerr<<"Wrong deviation type"<<endl;
			return false;
	}
}

bool Deviation::GeometricDeviation()
{
	dev.resize( mavn );

	for ( register int i = 0; i < mavn; i++ )
	{
		dev[i] = ug->NearestNeighbors(ma->Vertex(i))->Distance();
	}
		
	delete ug;
	ug = 0;

	Statistics();
	Deviation2Material();
	
	return true;
}

template<class Type>
bool Deviation::MeshDeviation( std::vector<Type>& aa, std::vector<Type>& ab )
{
	if ( ((int)aa.size() != mavn) || ((int)ab.size() != mbvn) )
	{
		cerr<<"Wrong attribute number"<<endl;
		return false;
	}
	Neighbor* Nearest;
	double d;

	dev.resize( mavn );
	dev = DBL_MAX;

	for ( register int i = 0; i < mavn; i++ )
	{
		Nearest = ug->NearestNeighbors(ma->Vertex(i))->Neighbors();
		while( Nearest )
		{
			if ( Nearest->v != -1 )
			{
				d = (aa[i]-ab[Nearest->v]).Length();
			}
			else if ( Nearest->e != -1 )
			{			

				const int& a = mb->Face(Nearest->f, Edge2Vertex[Nearest->e][0]);
				const int& b = mb->Face(Nearest->f, Edge2Vertex[Nearest->e][1]);

				d = (aa[i]- Lerp(ab[a], ab[b], Nearest->r1)).Length();
			}
			else
			{			
				const int& a = mb->Face(Nearest->f, 0);
				const int& b = mb->Face(Nearest->f, 1);
				const int& c = mb->Face(Nearest->f, 2);

				d = (aa[i] - Lerp( Lerp(ab[a], ab[c], Nearest->r1), Lerp(ab[b], ab[c], Nearest->r1), Nearest->r2)).Length();
			}
			if ( d < dev[i] )
			{
				dev[i] = d;
			}
			Nearest = Nearest->next;
		}
	}
		
	delete ug;
	ug = 0;

	Statistics();
	Deviation2Material();
	
	return true;
}

bool Deviation::GeometricDeviationSample()
{
	int i, j, k;
	Vector3d p; // Sample point

	for ( k=0; k<mbfn; k++ )
	{
		SampleFace( k );
		// for each scan line
		for ( j=0; j<samples[k].Height(); j++ )
		{
			// Start point
			p = ma->Vertex(k, 0) + samplev * (double)j;
			// For each point in current scan line
			for ( i=0; i<samples[k][j]; i++ )
			{
				samples[k](i, j) = ug->NearestNeighbors( p )->Distance();
				snum++;
				// Next sampled point
				p += sampleu;
			}
		}
	}

	delete ug;
	ug = 0;

	StatisticsSample();

	delete samples;
	samples = 0;

	return true ;
}

template<class Type>
bool Deviation::MeshDeviationSample( vector<Type>& aa, std::vector<Type>& ab )
{
	if ( ((int)aa.size() != mavn) || ((int)ab.size() != mbvn) )
	{
		cerr<<"Wrong attribute number"<<endl;
		return false;
	}

	int i, j, k;
	Vector3d p, ref; // Sample point
	Type mstep, nstep1, nstep2, an, bn, n;
	Neighbor* Nearest;
	double d;


	for ( k=0; k<mbfn; k++ )
	{
		SampleFace( k );
		// Reference point
		ref = ma->Vertex(k, 0);
		nstep1 = (ab[ ma->Face(k, 2) ] - ab[ ma->Face(k, 0) ]) / (double)samples[k].Height();
		nstep2 = (ab[ ma->Face(k, 2) ] - ab[ ma->Face(k, 1) ]) / (double)samples[k].Height();
		an = ab[ ma->Face(k, 0) ];
		bn = ab[ ma->Face(k, 1) ];
		// for each scan line
		for ( j=0; j<samples[k].Height(); j++ )
		{
			// Start point
			p = ma->Vertex(k, 0) + samplev * (double)j;
			mstep = (bn - an) / (double)samples[k][j];
			n = an;
			// For each point in current scan line
			for ( i=0; i<samples[k][j]; i++ )
			{
				Nearest = ug->NearestNeighbors( p )->Neighbors();

				samples[k](i, j) = DBL_MAX;

				while( Nearest )
				{
					if ( Nearest->v != -1 )
					{
						d = (aa[i]-ab[Nearest->v]).Length();
					}
					else if ( Nearest->e != -1 )
					{

						const int& a = mb->Face(Nearest->f, Edge2Vertex[Nearest->e][0]);
						const int& b = mb->Face(Nearest->f, Edge2Vertex[Nearest->e][1]);

						d = (aa[i]- Lerp(ab[a], ab[b], Nearest->r1)).Length();

					}
					else
					{
						const int& a = mb->Face(Nearest->f, 0);
						const int& b = mb->Face(Nearest->f, 1);
						const int& c = mb->Face(Nearest->f, 2);

						d = (aa[i] - Lerp( Lerp(ab[a], ab[c], Nearest->r1), Lerp(ab[b], ab[c], Nearest->r1), Nearest->r2)).Length();
					}
					if ( d < samples[k](i, j) )
					{
						samples[k](i, j) = d;
					}
					Nearest = Nearest->next;
				}
				snum++;
				// Next sampled point
				p += sampleu;
				n += mstep;
			}
			an += nstep1;
			bn += nstep2;
		}
	}
	
	delete ug;
	ug = 0;

	StatisticsSample();

	delete samples;
	samples = 0;
	
	return true;
}

bool Deviation::Statistics()
{
	register int i;
	// No deviation
	if ( !dev.size() )
	{
		cerr<<"Cannot compute statistics: no deviation value"<<endl;
		return false;
	}
	// Minimum
	mindev = dev.min();
	// Maximum
	maxdev = dev.max();
	// Mean
	meandev = dev.sum() / dev.size();
	// Variance
	vardev = 0;
	rmsdev = 0;
	for ( i = 0; i < (int)dev.size(); i++ )
	{
		rmsdev += dev[i] * dev[i];
		vardev += Sqr( dev[i] - meandev );
	}
	vardev /= (double)(dev.size() - 1 );
	rmsdev /= (double)dev.size();
	rmsdev = sqrt(rmsdev);
/*
	// Histogram
	valarray<int> hist(0, 1000);
	for ( i=0; i < (int)dev.size(); i++ )
	{
		hist[ (int)( dev[i] / maxdev * 999.0 ) ]++;
	}
	// Output deviation value histogram to "Histogram.dat" file
	ofstream f( "Histogram.dat" );
	for ( i=0; i<1000; i++ )
	{
		f<<(double)i/999.0*maxdev<<" "<<hist[i]<<endl;
	}
	f.close();
	// Mediane
	int m = hist[0];
	i = 0;
	while( m < (int)dev.size() / 2 )
	{
		i++;
		m += hist[i];
	}
	meddev = (double)i/199.0*maxdev;
*/
	return true;
}

void Deviation::StatisticsSample()
{
	double d;
	int i, j;
	int k;
	vardev = 0;
	maxdev = 0;
	mindev = DBL_MAX;
	meandev = 0;
	rmsdev = 0;
    
	// Histogram
	valarray<int> hist(0, 201);
	for ( k=0; k<mafn; k++ )
	{
	    for ( j=0; j<samples[k].Height(); j++ )
	    {
		for ( i=0; i<samples[k][j]; i++ )
		{
		    d = samples[k](i, j);
		    meandev += d;
		    if ( d > maxdev ) maxdev = d;
		    if ( d < mindev ) mindev = d;
		}
	    }
	}
	meandev /= snum;
	
	for ( k=0; k<mafn; k++ )
	{
		for ( j=0; j<samples[k].Height(); j++ )
		{
			for ( i=0; i<samples[k][j]; i++ )
			{
				d = samples[k](i, j);
				vardev += Sqr(d - meandev);
				rmsdev += d*d;
		//	hist[ (int)( d / maxdev * 200.0 ) ]++;
			}
		}
	}
	vardev /= (double)(snum - 1);
	rmsdev /= (double)snum;
	rmsdev = sqrt(rmsdev);
	
/*
	ofstream f( "Histogram.dat" );
	for ( i=0; i<201; i++ )
		f << ((double)i/200.0) * maxdev << " " << hist[i] << endl;
	f.close();
*/
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Deviation2Color
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Define color for input deviation
// Deviation must be in [0,1]
// Returned color is RGB color values in [0,1]
Vector3d Deviation::Deviation2Color( const double& d )
{
	if     ( d < 0    ) return Vector3d(0,0,0);
	else if ( d < 0.25 ) return Vector3d( 0, d * 4.0, 1 );
	else if ( d < 0.50 ) return Vector3d( 0, 1, 1 - (d - 0.25) * 4.0 );
	else if ( d < 0.75 ) return Vector3d( (d - 0.5) * 4.0, 1, 0 );
	else if ( d < 1    ) return Vector3d( 1, 1.0 - (d - 0.75) * 4.0, 0 );
	return Vector3d(1,0,0);
}

void Deviation::Deviation2Material()
{
	// Per vertex Color binding
	ma->Colors().resize( mavn );
	if ( dev_bound <= 0 )
	{
		dev_bound = maxdev;
	}
	if ( dev_bound == 0 )
	{
		// Avoid division by 0
		// if maximum deviation is null
		// (no deviation)
 		dev_bound = 1;
	}
	// Set deviation values as material
	for ( int i=0; i<mavn; i++ )
	{
		// Normalize deviation values
		ma->Color(i) = Deviation2Color( dev[i] / dev_bound );
	}
}

//===============================================================
//
// SampleFace
//
//===============================================================
// Sample a given face
void Deviation::SampleFace( int face )
{
	// Shortcut for vertex indices
	const int& a = ma->Face(face, 0);
	const int& b = ma->Face(face, 1);
	const int& c = ma->Face(face, 2);

	double x, y;
	int i;

	sampleu = ma->Vertex(b) - ma->Vertex(a);
	x = sampleu.Length();
 	sampleu *= 1.0 / x;
	sampleu *= step;

	samplev = ma->Vertex(c) - ma->Vertex(a);
	y = samplev.Length();
 	samplev *= 1.0/ y;
	samplev *= step;

	x /= step;
	y /= step;
	samples[face].SetHeight( (int)ceil(y) );

	y = x / y;

	for ( i=0; i<samples[face].Height(); i++ )
	{
		samples[face][i] = (int)ceil(x);
		x -= y;
	}

	samples[face].InitDev();
}
#endif
