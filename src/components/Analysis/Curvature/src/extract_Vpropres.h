#ifndef _EXTRACT_VPROPRES
#define _EXTRACT_VPROPRES


/*!
 * \file extract_Vpropres.h
 * \brief Eigenvalues and eigenvectors of a real symmetric square matrix by Jacobi's method 
 * \author Jean-Pierre Moreau
 */
 
#include <mepp_config.h>
#ifdef BUILD_component_Curvature

/**
 \fn	int ValPro(int N,double **A,double ER,int IM,double**U,double **AD);

 \brief	Eigenvalues and eigenvectors of a double symmetric matrix by JACOBI'S METHOD
 
 \param	N	NMAX of matrix A 
 \param A 	double symmetric matrix
 \param	ER	maximum number of iterations
 \param	IM	maximum number of iterations
 \param [out]	U 	eigenvectors given in lines
 \param [out]	AD	eigenvalues stored in ascending order on main diagonal 

 */
int ValPro(int N,double **A,double ER,int IM,double**U,double **AD);


/**
 \fn	void EigSrt(double **D,double **V,int N);

 \brief	Given the eigenvalues D and eigenvectors V as output from VALPRO,  
   this routine sorts the eigenvalues into ascending order, and       
  rearranges the lines of V accordingly. The method used is straight 
   insertion.

 \param [in,out]	D	eigenvalues stored in ascending order on main diagonal 
 \param [in,out]	V	eigenvectors given in lines
 \param	N			 	NMAX 
 */
void EigSrt(double **D,double **V,int N);

#endif

#endif
