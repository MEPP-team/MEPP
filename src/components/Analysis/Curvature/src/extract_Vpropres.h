#ifndef _EXTRACT_VPROPRES
#define _EXTRACT_VPROPRES

#include "../../../../mepp/mepp_config.h"
#ifdef BUILD_component_Curvature

int ValPro(int N,double **A,double ER,int IM,double**U,double **AD);
void EigSrt(double **D,double **V,int N);

#endif

#endif
