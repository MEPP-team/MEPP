// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/Surface_mesh_simplification/include/CGAL/Cartesian/MatrixC33.h $
// $Id: MatrixC33.h 32079 2006-06-26 16:30:40Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_CARTESIAN_MATRIXC33_H
#define CGAL_CARTESIAN_MATRIXC33_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <CGAL/Vector_3.h>
#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>
#include <CGAL/Null_matrix.h>

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/optional/optional.hpp>
#endif

CGAL_BEGIN_NAMESPACE

template <class R_>

/**
 \class	MatrixC33

 \brief	Matrix c 33. 
*/

class MatrixC33 
{
public:
  
  typedef R_ R ;
  
  typedef typename R::RT       RT ;
  typedef typename R::Vector_3 Vector_3;

  MatrixC33 ( Null_matrix )
   :
    mR0(NULL_VECTOR)
   ,mR1(NULL_VECTOR)
   ,mR2(NULL_VECTOR)
  {}
  MatrixC33()
  {
  }
  MatrixC33 ( RT const& r0x, RT const& r0y, RT const& r0z
            , RT const& r1x, RT const& r1y, RT const& r1z
            , RT const& r2x, RT const& r2y, RT const& r2z
            )
    :
     mR0(r0x,r0y,r0z)
    ,mR1(r1x,r1y,r1z)
    ,mR2(r2x,r2y,r2z)
  {}
  
  MatrixC33 ( Vector_3 const& r0, Vector_3 const& r1, Vector_3 const& r2 )
    :
     mR0(r0)
    ,mR1(r1)
    ,mR2(r2)
  {}
  
  Vector_3 const& r0() const { return mR0; }
  Vector_3 const& r1() const { return mR1; }
  Vector_3 const& r2() const { return mR2; }
  
  Vector_3& r0() { return mR0; }
  Vector_3& r1() { return mR1; }
  Vector_3& r2() { return mR2; }

  Vector_3 const& operator[] ( int row ) const { return row == 0 ? mR0 : ( row == 1 ? mR1 : mR2 ) ; }  
  Vector_3&       operator[] ( int row )       { return row == 0 ? mR0 : ( row == 1 ? mR1 : mR2 ) ; }  
  
  MatrixC33& operator+= ( MatrixC33 const& m )
  {
    mR0 = mR0 + m.r0() ;
    mR1 = mR1 + m.r1() ;
    mR2 = mR2 + m.r2() ;
    return *this ;
  }
  /*
  const MatrixC33 & operator=(MatrixC33 const & m) const  
  {
	  MatrixC33 n(m.mR0,m.mR1,m.mR2);	  
	  return n;
  }
*/
   void round_elements(void)   
  {
	int rounded[3][3];
	rounded[0][0] = floor(this->mR0.x() + 0.5);
	rounded[1][0] = floor(this->mR0.y() + 0.5);
	rounded[2][0] = floor(this->mR0.z() + 0.5);
	rounded[0][1] = floor(this->mR1.x() + 0.5);
	rounded[1][1] = floor(this->mR1.y() + 0.5);
	rounded[2][1] = floor(this->mR1.z() + 0.5);
	rounded[0][2] = floor(this->mR2.x() + 0.5);
	rounded[1][2] = floor(this->mR2.y() + 0.5);
	rounded[2][2] = floor(this->mR2.z() + 0.5);
	
	Vector n1(rounded[0][0],rounded[1][0],rounded[2][0]);
	Vector n2(rounded[0][1],rounded[1][1],rounded[2][1]);
	Vector n3(rounded[0][2],rounded[1][2],rounded[2][2]);
	
	this->mR0 = n1;
	this->mR1 = n2;
	this->mR2 = n3;	
  }
  RT Get(int i,int j)
  {
	  Vector_3 vec;
	  if (i == 0)
		  vec = mR0;
	  else if (i == 1)
		  vec = mR1;
	  else if (i==2)
		  vec = mR2;

	  if (j == 0)
		  return vec.x();
	  else if (j == 1)
		  return vec.y();
	  else
		  return vec.z();
  }   

  /**
   \fn	MatrixC33& MatrixC33::operator-= ( MatrixC33 const& m )
  
   \brief	Subtraction assignment operator for MatrixC33.
  
   \param	m	The.
  
   \return	The result of the operation.
   */

  MatrixC33& operator-= ( MatrixC33 const& m )
  {
    mR0 = mR0 - m.r0() ;
    mR1 = mR1 - m.r1() ;
    mR2 = mR2 - m.r2() ;
    return *this ;
  }

  /**
   \fn	MatrixC33& MatrixC33::operator*= ( RT const& c )
  
   \brief	Multiplication assignment operator for MatrixC33.
  

   \param	c	The.
  
   \return	The result of the operation.
   */

  MatrixC33& operator*= ( RT const& c )
  {
    mR0 = mR0 * c ;
    mR1 = mR1 * c ;
    mR2 = mR2 * c ;
    return *this ;
  }

  /**
   \fn	MatrixC33& MatrixC33::operator/= ( RT const& c )
  
   \brief	Division assignment operator for MatrixC33.
  

   \param	c	The.
  
   \return	The result of the operation.
   */

  MatrixC33& operator/= ( RT const& c )	
  {
    mR0 = mR0 / c ;
    mR1 = mR1 / c ;
    mR2 = mR2 / c ;
    return *this ;
  }

  /**
   \fn	friend MatrixC33 MatrixC33::operator+ ( MatrixC33 const& a, MatrixC33 const& b )
  
   \brief	Addition operator for MatrixC33.

  
   \param	a	The first value.
   \param	b	A value to add to it.
  
   \return	The result of the operation.
   */

  friend MatrixC33 operator+ ( MatrixC33 const& a, MatrixC33 const& b )
  {
    return MatrixC33(a.r0()+b.r0()
                    ,a.r1()+b.r1()
                    ,a.r2()+b.r2()
                    );
  }

  /**
   \fn	friend MatrixC33 MatrixC33::operator- ( MatrixC33 const& a, MatrixC33 const& b )
  
   \brief	Subtraction operator for MatrixC33.

   \param	a	The first value.
   \param	b	A value to subtract from it.
  
   \return	The result of the operation.
   */

  friend MatrixC33 operator- ( MatrixC33 const& a, MatrixC33 const& b )
  {
    return MatrixC33(a.r0()-b.r0()
                    ,a.r1()-b.r1()
                    ,a.r2()-b.r2()
                    );
  }

  /**
   \fn	friend MatrixC33 MatrixC33::operator* ( MatrixC33 const& m, RT const& c )
  
   \brief	Multiplication operator for MatrixC33.
  
  
   \param	m	The first value to multiply.
   \param	c	The second value to multiply.
  
   \return	The result of the operation.
   */

  friend MatrixC33 operator* ( MatrixC33 const& m, RT const& c )
  {
    return MatrixC33(m.r0()*c,m.r1()*c,m.r2()*c);
  }
  friend MatrixC33 operator* ( RT const& c, MatrixC33 const& m )
  {
    return MatrixC33(m.r0()*c,m.r1()*c,m.r2()*c);
  }
  
  friend MatrixC33 operator/ ( MatrixC33 const& m, RT const& c )
  {
    return MatrixC33(m.r0()/c,m.r1()/c,m.r2()/c);
  }
  
  friend Vector_3 operator* ( MatrixC33 const& m, Vector_3 const& v )
  {
    return Vector_3(m.r0()*v,m.r1()*v,m.r2()*v);
  }
  friend Vector_3 operator* ( Vector_3 const& v, MatrixC33 const& m )
  {
    return Vector_3(v*m.r0(),v*m.r1(),v*m.r2());
  }

  /**
   \fn	RT MatrixC33::determinant() const
  
   \brief	Gets the determinant of the Matrix.
  
   \return	.
   */

  RT determinant() const
  {
	  return CGAL::determinant/*det3x3_by_formula*/(r0().x(),r0().y(),r0().z()	// MT: changement CGAL 3.5
                            ,r1().x(),r1().y(),r1().z()
                            ,r2().x(),r2().y(),r2().z()
                            );
  }

  /**
   \fn	MatrixC33& MatrixC33::transpose()
  
   \brief	Gets the transpose of Matrix.

   \return	.
   */

  MatrixC33& transpose()
  {
    mR0 = Vector_3(r0().x(),r1().x(),r2().x());
    mR1 = Vector_3(r0().y(),r1().y(),r2().y()); 
    mR2 = Vector_3(r0().z(),r1().z(),r2().z());
    return *this ;
  }
    
private:

  Vector_3 mR0 ;
  Vector_3 mR1 ;
  Vector_3 mR2 ;
} ;

template<class R>
inline
MatrixC33<R> direct_product ( Vector_3<R> const& u, Vector_3<R> const& v )
{
  return MatrixC33<R>( v * u.x()
                     , v * u.y()
                     , v * u.z()
                     ) ;
}

template<class R>
MatrixC33<R> transposed_matrix ( MatrixC33<R> const& m )
{
  MatrixC33<R> copy = m ;
  copy.Transpose();
  return copy ;
}

template<class R>
MatrixC33<R> cofactors_matrix ( MatrixC33<R> const& m )
{
  typedef typename R::RT RT ;
  
  RT c00 =  det2x2_by_formula(m.r1().y(),m.r1().z(),m.r2().y(),m.r2().z());
  RT c01 = -det2x2_by_formula(m.r1().x(),m.r1().z(),m.r2().x(),m.r2().z());
  RT c02 =  det2x2_by_formula(m.r1().x(),m.r1().y(),m.r2().x(),m.r2().y());
  
  RT c10 = -det2x2_by_formula(m.r0().y(),m.r0().z(),m.r2().y(),m.r2().z());
  RT c11 =  det2x2_by_formula(m.r0().x(),m.r0().z(),m.r2().x(),m.r2().z());
  RT c12 = -det2x2_by_formula(m.r0().x(),m.r0().y(),m.r2().x(),m.r2().y());
  
  RT c20 =  det2x2_by_formula(m.r0().y(),m.r0().z(),m.r1().y(),m.r1().z());
  RT c21 = -det2x2_by_formula(m.r0().x(),m.r0().z(),m.r1().x(),m.r1().z());
  RT c22 =  det2x2_by_formula(m.r0().x(),m.r0().y(),m.r1().x(),m.r1().y());
  
  return MatrixC33<R>(c00,c01,c02
                     ,c10,c11,c12
                     ,c20,c21,c22
                     );
}

template<class R>
MatrixC33<R> adjoint_matrix ( MatrixC33<R> const& m )
{
  return cofactors_matrix(m).transpose()  ;
}
/*
template<class R>
boost::optional< MatrixC33<R> > inverse_matrix ( MatrixC33<R> const& m )
{
  typedef typename R::RT RT ;
  
  typedef MatrixC33<R> Matrix ;
  
  typedef boost::optional<Matrix> result_type ;
   
  result_type rInverse ;
  
  RT det = m.determinant();
  
  if ( ! CGAL_NTS is_zero(det) )
  {
    RT c00 = (m.r1().y()*m.r2().z() - m.r1().z()*m.r2().y()) / det; 
    RT c01 = (m.r2().y()*m.r0().z() - m.r0().y()*m.r2().z()) / det;
    RT c02 = (m.r0().y()*m.r1().z() - m.r1().y()*m.r0().z()) / det; 
  
    RT c10 = (m.r1().z()*m.r2().x() - m.r1().x()*m.r2().z()) / det; 
    RT c11 = (m.r0().x()*m.r2().z() - m.r2().x()*m.r0().z()) / det; 
    RT c12 = (m.r1().x()*m.r0().z() - m.r0().x()*m.r1().z()) / det; 
  
    RT c20 = (m.r1().x()*m.r2().y() - m.r2().x()*m.r1().y()) / det; 
    RT c21 = (m.r2().x()*m.r0().y() - m.r0().x()*m.r2().y()) / det; 
    RT c22 = (m.r0().x()*m.r1().y() - m.r0().y()*m.r1().x()) / det; 
    
    rInverse = result_type( Matrix(c00,c01,c02
                                  ,c10,c11,c12
                                  ,c20,c21,c22
                                  )
                          ) ;
  }
  
  return rInverse ;
}
*/
template<class R>
MatrixC33<R>  inverse_matrix( MatrixC33<R> const& m )
{
  typedef typename R::RT RT ;
  
  typedef MatrixC33<R> Matrix ;
  
  //typedef boost::optional<Matrix> result_type ;
   
  Matrix rInverse(0,0,0,0,0,0,0,0,0);
  
  RT det = m.determinant();
  
  if ( ! CGAL_NTS is_zero(det) )
  {
    RT c00 = (m.r1().y()*m.r2().z() - m.r1().z()*m.r2().y()) / det; 
    RT c01 = (m.r2().y()*m.r0().z() - m.r0().y()*m.r2().z()) / det;
    RT c02 = (m.r0().y()*m.r1().z() - m.r1().y()*m.r0().z()) / det; 
  
    RT c10 = (m.r1().z()*m.r2().x() - m.r1().x()*m.r2().z()) / det; 
    RT c11 = (m.r0().x()*m.r2().z() - m.r2().x()*m.r0().z()) / det; 
    RT c12 = (m.r1().x()*m.r0().z() - m.r0().x()*m.r1().z()) / det; 
  
    RT c20 = (m.r1().x()*m.r2().y() - m.r2().x()*m.r1().y()) / det; 
    RT c21 = (m.r2().x()*m.r0().y() - m.r0().x()*m.r2().y()) / det; 
    RT c22 = (m.r0().x()*m.r1().y() - m.r0().y()*m.r1().x()) / det; 
    
    rInverse = Matrix(c00,c01,c02,c10,c11,c12,c20,c21,c22);                          
  }
  
  return rInverse;
}


template<class R>
MatrixC33<R> product_matrices(MatrixC33<R> const& m,MatrixC33<R> const& n)
{
    typedef typename R::RT RT ;
  
    typedef MatrixC33<R> Matrix ;
  
    //typedef boost::optional<Matrix> result_type ;
   
    Matrix result(0,0,0,0,0,0,0,0,0);
	RT mT[3][3];
	mT[0][0] = m.r0().x();mT[0][1] = m.r0().y();mT[0][2] = m.r0().z();
	mT[1][0] = m.r1().x();mT[1][1] = m.r1().y();mT[1][2] = m.r1().z();
	mT[2][0] = m.r2().x();mT[2][1] = m.r2().y();mT[2][2] = m.r2().z();

	RT nT[3][3];
	nT[0][0] = n.r0().x();nT[0][1] = n.r0().y();nT[0][2] = n.r0().z();
	nT[1][0] = n.r1().x();nT[1][1] = n.r1().y();nT[1][2] = n.r1().z();
	nT[2][0] = n.r2().x();nT[2][1] = n.r2().y();nT[2][2] = n.r2().z();

	RT C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

	for (int i=0; i< 3 ; i++)
	{
		for (int j=0; j<3; j++)
		{
			for (int k=0;k<3;k++)
			{
				C[i][j] += mT[i][k]*nT[k][j];
			}
		}
	}
	
    
	
	return Matrix(C[0][0],C[0][1],C[0][2],C[1][0],C[1][1],C[1][2],C[2][0],C[2][1],C[2][2]);                          
}

CGAL_END_NAMESPACE

#endif

#endif // CGAL_CARTESIAN_MATRIXC33_H //
// EOF //
 
 
