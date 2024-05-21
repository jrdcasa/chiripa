#ifndef __DISCRIMINANT_H
#define __DISCRIMINANT_H

#include <stdio.h>
#include <math.h>


static void c_discriminant(double ref[],
                           double conf[],
                           double n[],
                           int Ni, int Nj, int dim,
                           float Rvdw_i[], float Rvdw_j[], double rmaxmin[2],
                           int max_index[2], int min_index[2])
{

    double rijx, rijy, rijz, dij;
    double b, s1, s2;
    double rmax, rmin;
    int imax, jmax, imin, jmin;
    double discriminant;
    double tmp;
    rijx = rijy = rijz = b = 0.0;
    rmin = LONG_MAX;
    rmax = LONG_MIN;
    imax = jmax = -1;
    imin = jmin = -1;

    for (int i=0; i<Ni; i++) {
        for (int j=0; j<Nj; j++) {
            // Calculate the distance of the two atoms
            rijx = conf[j*dim+0]-ref[i*dim+0];
            rijy = conf[j*dim+1]-ref[i*dim+1];
            rijz = conf[j*dim+2]-ref[i*dim+2];
            dij = sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));
            // Dot product of n and rij
            b = rijx*n[0] + rijy*n[1] + rijz*n[2];
            //if(i==1 && j==0) printf("%f\n", b);
            tmp = Rvdw_i[i] + Rvdw_j[j];
            discriminant = b*b - dij*dij + tmp*tmp;
            if (discriminant >= 0.0)
            {
                s1 = -b + sqrt(discriminant);
                s2 = -b - sqrt(discriminant);
                if (s1 > rmax) {
                    //printf("%f, %f, %i, %i\n", s1, rmax, imax, jmax);
                    rmax = s1; imax = i; jmax = j;
                    //printf("** %f, %f, %i, %i\n", s1, rmax, imax, jmax);
                }
                else if (s2 < rmin) {
                    rmin = s2; imin = i; jmin = j;
                }
            }

        }
    }

    rmaxmin[0] = rmax;
    rmaxmin[1] = rmin;
    max_index[0] = imax;
    max_index[1] = jmax;
    min_index[0] = imin;
    min_index[1] = jmin;

}

static void c_discriminant_okuwaki(double ref[],
                                   double conf[],
                                   double n[],
                                   int Ni, int Nj, int dim,
                                   double dij_base[],
                                   double rmaxmin[2],
                                   int max_index[2], int min_index[2])
{

    double rijx, rijy, rijz, dij;
    double b, s1, s2;
    double rmax, rmin;
    int imax, jmax, imin, jmin;
    double discriminant;
    double tmp;
    rijx = rijy = rijz = b = 0.0;
    rmin = LONG_MAX;
    rmax = LONG_MIN;
    imax = jmax = -1;
    imin = jmin = -1;

    for (int i=0; i<Ni; i++) {
        for (int j=0; j<Nj; j++) {
            // Calculate the distance of the two atoms
            rijx = conf[j*dim+0]-ref[i*dim+0];
            rijy = conf[j*dim+1]-ref[i*dim+1];
            rijz = conf[j*dim+2]-ref[i*dim+2];
            dij = sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));
            // Dot product of n and rij
            b = rijx*n[0] + rijy*n[1] + rijz*n[2];
            //if(i==1 && j==0) printf("%f\n", b);
//            tmp = Rvdw_i[i] + Rvdw_j[j];
            tmp = dij_base[i*Nj+j];

            //printf("%d %d %d %f\n", i, j, i*Nj+j, tmp);
            discriminant = b*b - dij*dij + tmp*tmp;
            if (discriminant >= 0.0)
            {
                s1 = -b + sqrt(discriminant);
                s2 = -b - sqrt(discriminant);
                if (s1 > rmax) {
                    //printf("%f, %f, %i, %i\n", s1, rmax, imax, jmax);
                    rmax = s1; imax = i; jmax = j;
                    //printf("** %f, %f, %i, %i\n", s1, rmax, imax, jmax);

                }
                else if (s2 < rmin) {
                    rmin = s2; imin = i; jmin = j;
                }
            }
            //if(i==17 && j==1) printf("%f %f %f %f\n", dij_base[i*Nj+j], dij, s1, s2);

        }
    }

    rmaxmin[0] = rmax;
    rmaxmin[1] = rmin;
    max_index[0] = imax;
    max_index[1] = jmax;
    min_index[0] = imin;
    min_index[1] = jmin;

}



#endif