#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include<math.h>
#include"extract_Vpropres.h"

#define ABS_GUY  fabs
#define SQRT_GUY sqrt
#define MACH_EPS 2e-16


#define NB_VP_PRIS 10

/*********************************************************************
*      Eigenvalues and eigenvectors of a double symmetric matrix     *
*                      by JACOBI'S METHOD                            *
* ------------------------------------------------------------------ *
* INPUTS:                                                            *
*           A(N,N) : double symmetric matrix                         *
*           N      : NMAX of matrix A                                *
*           ER     : desired précision                               *
*           IM     : maximum number of iterations                    *
* OUTPUTS:                                                           *
*           AD(N,N): eigenvalues stored in ascending order on main   *
*                    diagonal                                        *
*           U(N,N) : eigenvectors given in lines                     *
*                                                                    *
*                          C++ version from FORTRAN by J-P Moreau    *
* ------------------------------------------------------------------ *
* NOTA: index 0 not used here.                                       *
*********************************************************************/
int ValPro(int N,double **A,double ER,int IM,double**U,double **AD)  {
	double XNA,XND,ST,C,S,P,PR,Q,B,TE;
	int I,IT,J,K,NI,NJ=0;
	NI=-1;
	for (I=1; I<N+1; I++)
		for (J=1; J<N+1; J++) {
		  U[I][J]=0.0;
		  AD[I][J]=0.0;
		}
	IT=0;
	XNA=0.0;
	XND=0.0;
	ST=0.0;
	for (I=1; I<N+1; I++) {
		U[I][I]=1.0;
		XND=XND+A[I][I]*A[I][I];
		if (I==N) goto s10;
		for (J=I+1; J<N+1; J++) {
			ST=ST+(A[I][J]*A[I][J]);
			U[I][J]=0.0;
			U[J][I]=0.0;
		}
	}
 s10: XNA=XND+2.0*ST;
      do {
        ST=0.0;
        for (I=1; I<N; I++) {
          for (J=I+1; J<N+1; J++) {
            TE=ABS_GUY(A[I][J]);
            if (TE > ST) {
              ST=TE;
              NI=I;
              NJ=J;
            }
          }
        }
		if (NI==-1)
			return -1;
        B=A[NI][NI]-A[NJ][NJ];
        if (ABS_GUY(B) > MACH_EPS) goto s15;
        C=1.0/SQRT_GUY(2.0);

        if (A[NI][NJ]>0) S=C;
        else if (A[NI][NJ]<0) S=-C;
        else S=0;
        // S=C*DSIGN(0.1D+01,A[NI][NJ]) fortran
        goto s20;
   s15: Q=ABS_GUY(B);
        if (B<0) P=-2.0*A[NI][NJ];
        else if (B>0) P=2.0*A[NI][NJ];
        else P=0.0;
        // P=2.0*A[NI][NJ]*DSIGN(0.1D+01,B)  fortran
        ST=SQRT_GUY(P*P+Q*Q);
        C=SQRT_GUY((1.0+Q/ST)/2.0);
        S=P/(2.0*ST*C);
   s20: for (K=1; K<N+1; K++) {
          ST=U[K][NI];
          U[K][NI]=C*ST+S*U[K][NJ];
          U[K][NJ]=C*U[K][NJ]-S*ST;
        }
        if (NI==1) goto s30;
        for (K=1; K<NI; K++) {
          ST=A[K][NI];
          A[K][NI]=C*ST+S*A[K][NJ];
          A[K][NJ]=C*A[K][NJ]-S*ST;
        }
   s30: if (NJ==NI+1) goto s40;
        for (K=NI+1; K<NJ; K++) {
          ST=A[NI][K];
          A[NI][K]=C*ST+S*A[K][NJ];
          A[K][NJ]=C*A[K][NJ]-S*ST;
        }
   s40: if (NJ==N) goto s50;
        for (K=NJ+1; K<N+1; K++) {
          ST=A[NI][K];
          A[NI][K]=C*ST+S*A[NJ][K];
          A[NJ][K]=C*A[NJ][K]-S*ST;
        }
   s50: XND=XND+2.0*A[NI][NJ]*A[NI][NJ];
        ST=A[NI][NI];
        A[NI][NI]=C*C*ST+2.0*S*C*A[NI][NJ]+S*S*A[NJ][NJ];
        A[NJ][NJ]=C*C*A[NJ][NJ]+S*S*ST-2.0*S*C*A[NI][NJ];
        A[NI][NJ]=0.0;
        IT=IT+1;
        PR=ABS_GUY(1.0-XNA/XND);
        if (PR<ER) goto s60;
      } while (IT <= IM);

  s60:
      for (I=1; I<N+1; I++)
        for (J=1; J<N+1; J++)
          AD[I][J]=A[I][J];
    return 0;
    } // valpro

  /*********************************************************************
  * Given the eigenvalues D and eigenvectors V as output from VALPRO,  *
  * this routine sorts the eigenvalues into ascending order, and       *
  * rearranges the lines of V accordingly. The method used is straight *
  * insertion.                                                         *
  *********************************************************************/
    void EigSrt(double **D,double **V,int N) {
	  double P;
	  int I,J,K;
	  for (I=1; I<N; I++) {
        K=I;
	    P=D[I][I];
	    for (J=I+1; J<N+1; J++)
		  if (fabs(D[J][J]) >= fabs(P))  {
	        K=J;
	        P=D[J][J];
		  }

		if (K!=I)  {
	      D[K][K]=D[I][I];
	      D[I][I]=P;
          for (J=1; J<N+1; J++)  {
	        P=V[J][I];
	        V[J][I]=V[J][K];
	        V[J][K]=P;
		  }
		}
	  }
	}
/*
void MATREAD(double A[NMAX][NMAX],int N) {
	int I,J; double temp;
    for (I=1; I<N+1; I++)
		for (J=1; J<N+1; J++) {
			fscanf(fp1,"%lg",&temp);
			A[I][J] = temp;
        }
}
void MATPRINT(char *titre,double A[NMAX][NMAX],int N) {
	int I, J;
    printf("\n\n  %s  \n\n",titre);
    for (I=1; I<N+1; I++)
		for (J=1; J<N+1; J++)  {
			printf("%13.6f",A[I][J]);
			if (J==N) printf("\n");
		}
}
int WriteHead (FILE *fp, char * string){
	if (string == NULL) return (-2);
	if (fprintf (fp,"\n%s\n%s\n%s\n\n", Separator, string, Separator) <= 0)
		return (-1);
	return 0;
}
int WriteEnd (FILE *fp){
	if (fprintf (fp,"\n\n%s\n", Separator) <= 0) return (-1);
	 return 0;
}  */
#endif
