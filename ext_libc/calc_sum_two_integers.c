#ifndef __SUMINT_H
#define __SUMINT_H

#include <stdio.h>
#include <math.h>


#ifdef PARALLEL
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif

static int c_sumtwointegers(double* ref, double* other, int s1, int s2 )

{

    int r = 0;
    r = s1 + s2;
    printf("Hola PEPEPEPEPEPE %i \n",r);
    printf("\n");
    return s1+s2;

}
#endif