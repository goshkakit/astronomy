#include <stdio.h>
#include <stdlib.h>
#include "hermite.h"
//---------------------------------------------------------------------------------//
//C Interpolation by a polynomial using NVAL out of NMAX given data points
//C
//C       Input : ITYP  : +-1: use Lagrange polynomial of degree NVAL-1
//C                         2: use Hermite formula:
//C                            Polynomial of degree 2*NVAL-1
//C
//C               NVAL  : number of points to use for interpolation
//C               NMAX  : number of given points in list
//C               X(I)  : arguments of given values (I=1,...,NMAX)
//C               Y(I)  : functional values Y=f(X)
//C               Z(I)  : derivatives       Z=f(X) (ITYP=2)
//C               XP    : interpolation argument
//C
//C       Output: YP    : interpolated value at XP
//C               ZP    : first derivative of YP (ITYP=-1 only)
//C               IRCODE: return code
//C                       0 = ok
//C                       1 = interpolation not in center interval
//C                       2 = error, argument out of range
//C
//C  The function selects the NVAL values to be used for interpolation such
//C  that the interpolated data point is located in the center interval.
//C  (Works best for NVAL = even number, of course).
//C
//c  An error code is generated if the data point is too close to the
//c  edges or outside of the given arguments.
//C
//C  Author:   W. Gurtner
//C            Astronomical Institute
//C            University of Bern
//C
//C  Created:  June 2002
//C
//C  Modified: 25-JUL-2005 : Check boundaries
//---------------------------------------------------------------------------------//
void hermite (int ityp, double *x, double *y, double *z, int nmax, int nval, double xp, double *yp, double *zp, int *ircode)
{
  int i = 0;
  int i0, j, k, n;
  int debug= 0;
  double pj, sk, vi, ui;

/**
  for (i=0;i<nmax;i++)
    {
      printf("i %d x %f y %f\n",i,x[i],y[i]);
    }
**/
  *ircode= 0;
  *yp= 0.e0;
  *zp= 0.e0;
  n=nval-1;
  /**printf("nval= %d nmax= %d xp= %f\n",nval,nmax,xp);**/

/* Look for given values to be used */ 
  if (xp < x[0] || xp > x[nmax-1])
    {
      *ircode= 2;
      if (debug) printf ("Her1 %d %f %f\n",i,x[i],xp);
      return;
    }

/* Look for given value immediately preceeding interpolation argument */
  for (i=0; i<nmax; i++)
    {
      if (x[i] >= xp)
        {
          *ircode= 0;
          break;
        }
    }

/*  Start index in vectors x,y,z */
  i0=i-(n+1)/2;
  if (i0 < 0)
    {
      if (debug) printf ("Her2 %d %d %d %d %d\n",i,n,i0,i0+n,nmax);
      i0= 0; /* or 1? */
      *ircode=1;
    }
  if (i0+n > nmax)
    {
      if (debug) printf ("Her2 %d %d %d %d %d\n",i,n,i0,i0+n,nmax);
      i0= nmax- n;
      *ircode=1;
    }

/* Lagrange formula for polynomial interpolation */ 
  if (abs(ityp) == 1)
    {
 
      for (i=0; i<n+1; i++) /* i=0 or 1?? */
        {
          pj=1.e0;
          for (j=0; j<n+1; j++)
            {
              if (j != i) pj*=(xp-x[j+i0])/(x[i+i0]-x[j+i0]);
            }
          *yp+=y[i+i0]*pj;
        }
      /**printf("yp= %f\n",*yp);**/

/* COMPUTE DERIVATIVE OF THE LAGRANGE POLYNOMIAL */
  if (ityp == -1)
    {
      for (i=0; i<n+1; i++) /* i=0 or 1?? */
        {
          sk=0.e0;
          for (k=0; k<n+1; k++) /* i=0 or 1?? */
            {
              if (k != i) 
                {
                  pj=1.e0;
                  for (j=0; j<n+1; j++) /* i=0 or 1?? */
                    {
                      if (j != i && j != k)
                        {
                          pj*=(xp-x[j+i0])/(x[i+i0]-x[j+i0]);
                        }
                    } 
                  sk+=pj/(x[i+i0]-x[k+i0]);
                }
            }
          *zp+=y[i+i0]*sk;
        }
    }


/* HERMITE INTERPOLATION (ADDITIONAL USE OF THE DERIVATIVES) */ 
    } else
      if (ityp == 2)
        { 
			 /*printf("Vel:\n");*/
          for (i=0; i<n+1; i++)
            {
              sk=0.e0;
              for (k=0; k<n+1; k++)
                {
                  if (k != i) sk+= 1.e0/(x[i+i0]-x[k+i0]);
                }
              vi= 1.e0-2.e0*(xp-x[i+i0])*sk;
              ui= xp-x[i+i0];

              pj=1.e0;
              for (j=0; j<n+1; j++)
                {
                  if (j != i) pj*=(xp-x[j+i0])/(x[i+i0]-x[j+i0]);
                }

              *yp+=(y[i+i0]*vi+z[i+i0]*ui)*pj*pj;
				  /*printf("vel end: %d %d %f %f %f %f\n",i,i0,ui,pj,z[i+i0],*yp);*/
            }
        }

      *ircode=0;
 }
