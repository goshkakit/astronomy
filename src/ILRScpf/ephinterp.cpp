#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ephi.h"
#include "ephinterp.h"
#include "hermite.h"

int debug = 0;
/*-----------------------------------------------------------------------
     void cinterp_pvac (jdint, seceph, outvec, backvec,
     	aberout,aberback, corr12, corr23, ircode)

  Interpolate position, velocity, aberration, and relativistic corrections
  from cpf file arrays stored in memory.

       jdint        Julian date (integer day)
       seceph       Julian day  (seconds)
       outvec       outv(i,1) ,...,outv(i,6): x,y,z,vx,vy,vz   outbound vector
       backvec      backv(i,1),...,outv(i,6): x,y,z,vx,vy,vz   inbound vector
       aberout      outbound aberration correction
       aberback     inbound aberration correction
       corr12       outbound relativistic correction
       corr23       inbound relativistic correction
       ircode       error code: 0=OK, 1,2=error (time out of range)

  Author: R. L. Ricklefs / UT CSR
  Date: 1 May, 2003
  History:
-----------------------------------------------------------------------*/
void cinterp_pvac (double jdint, double seceph,
	      double outvec[6], double backvec[6],
	      double aberout[3], double aberback[3],
	      double *corr12, double *corr23, int fvel, int *ircode, eh &ephihdr, ei &ephi )
{
  double itime, z, vz;		/* z is dummy for derivative field */
  int j;
  //void hermite ();

  itime = ((seceph / 86400.e0 - ephi.jdf[0]) + (jdint - ephi.jdi[0])) * 86400.e0;
/**
  debug= 1;
  printf("itime %f seceph %f jfd[0] %f jding %f jdi[0] %f\n",
		  itime,seceph,ephi.jdf[0],jdint,ephi.jdi[0]);
**/

  hermite (1, ephi.jds, ephi.c12, &z, ephi.nv, 10, itime, corr12, &vz, ircode);
  if (*ircode != 0)
    {
      if (debug)
	{
	  printf ("readi jds(0): %f jds(ephi.nv-1): %f %d\n",
		  ephi.jds[0], ephi.jds[ephi.nv - 1], ephi.nv);
	  printf ("itime= %f\n", itime);
	  printf ("hermite error: c12\n");
	}
      return;
    }
  hermite (1, ephi.jds, ephi.c23, &z, ephi.nv, 10, itime, corr23, &vz, ircode);
  if (*ircode != 0)
    {
      if (debug)
	printf ("hermite error: c23\n");
      return;
    }
  for (j = 0; j < 3; j++)
    {
      hermite (1, ephi.jds, &ephi.aberoutv[j][0], &z, ephi.nv, 10, itime, &aberout[j], &vz, ircode); //MAXENT
      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: aberout %d\n", j);
	  return;
	}
      hermite (1, ephi.jds, &ephi.aberbackv[j][0], &z, ephi.nv, 10, itime, &aberback[j], &vz, ircode);
      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: aberback %d\n", j);
	  return;
	}
    }
/* 3 or 6 elements? Handle velocity if present */
  for (j = 0; j < 6; j++)
    {
      /* Out */
      /* Position */
      if (j < 3)
	{
	  if (fabs (ephi.outv[j + 3][0]) > 1.e-10)
	    {
	      /* Use velocities if they exist (hermite) */
	      hermite( 2, ephi.jds, &ephi.outv[j][0], &ephi.outv[j + 3][0], ephi.nv, 10, itime, &outvec[j], &vz, ircode );
	    }
	  else
	    {
	      /* Otherwise, ignore them (fvel = -1) or create them (fvel = 1) */
	      hermite( -fvel, ephi.jds, &ephi.outv[j][0], &z, ephi.nv, 10, itime, &outvec[j], &outvec[j + 3], ircode);
	    }
	}
      /* Velocity */
      /* Interpolate if they exist */
      else if (fabs (ephi.outv[j][0]) > 1.e-10)
	{
	  hermite (1, ephi.jds, &ephi.outv[j][0], &z, ephi.nv, 10, itime,  &outvec[j], &vz, ircode);
	}

      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: outvec %d\n", j);
	  return;
	}

      /* Back */
      /* Position */
      if (j < 3)
	{
	  if (fabs (ephi.backv[j + 3][0]) > 1.e-10)
	    {
	      /* Use velocities if they exist (hermite) */
	      hermite (2, ephi.jds, &ephi.backv[j][0], &ephi.backv[j + 3][0], ephi.nv, 10, itime, &backvec[j], &vz, ircode);
	    }
	  else
	    {
	      /* Otherwise, ignore them (fvel = -1) or create them (fvel = 1) */
	      hermite (-fvel, ephi.jds, &ephi.backv[j][0], &z, ephi.nv,  10, itime, &backvec[j], &backvec[j + 3], ircode);
	    }
	}
      /* Velocity */
      /* Interpolate if they exist */
      else if (fabs (ephi.backv[j][0]) > 1.e-10)
	{
	  hermite (1, ephi.jds, &ephi.backv[j][0], &z, ephi.nv, 10, itime,  &backvec[j], &vz, ircode);
	}
      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: backvec %d\n", j);
	  return;
	}
      if (debug)
	printf ("readi-> %d %f %f\n", j, outvec[j], backvec[j]);
    }
}

/*-----------------------------------------------------------------------
    void interp_off (jdint, seceph, offsetout,offsetback,
         ircode)

  Interpolate offset from center of main body from cpf file arrays
  stored in memory.

       jdint        Julian date (integer day)
       seceph       Julian day  (seconds)
       offsetout    outbound (or bounce) offset in m
       offsetback   inbound offset in m
       ircode       error code: 0=OK, 1,2=error (time out of range)

  Author: R. L. Ricklefs / UT CSR
  Date: 1 May, 2004
  History:
-----------------------------------------------------------------------*/
void cinterp_off (double jdint, double seceph,  double offsetout[6], double offsetback[6], int *ircode, eh &ephihdr, ei &ephi)
{
  double itime, z, vz;		/* z is dummy for derivative field */
  int debug = 0;
  int j;
  //void hermite ();

  itime =
    ((seceph / 86400.e0 - ephi.offsetjdf[0]) +
     (jdint - ephi.offsetjdi[0])) * 86400.e0;
  for (j = 0; j < 3; j++)
    {
      hermite (1, ephi.offsetjds, &ephi.offsetoutv[j][0], &z,   ephi.nvoff, 10, itime, &offsetout[j], &vz, ircode);
      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: offsetout %d\n", j);
	  return;
	}
      hermite (1, ephi.offsetjds, &ephi.offsetbackv[j][0], &z, ephi.nvoff, 10, itime, &offsetback[j], &vz, ircode);
      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: offsetback %d\n", j);
	  return;
	}
    }
}

void cinterp_rot (double jdint, double seceph, int *arottype, double *agast, double arotang[3], int *ircode, eh &ephihdr, ei &ephi)
{
  double itime, z, vz;		/* z is dummy for derivative field */
  int debug = 0;
  int j;
  //void hermite ();

  itime =
    ((seceph / 86400.e0 - ephi.rotjdf[0]) +
     (jdint - ephi.rotjdi[0])) * 86400.e0;
  *arottype = ephihdr.rot_type;

  hermite (1, ephi.rotjds, ephi.gast, &z, ephi.nvrot, 10, itime, agast, &vz, ircode);
  if (*ircode != 0)
    {
      if (debug)
	printf ("hermite error: gast\n");
      return;
    }

  for (j = 0; j < 3; j++)
    {
      hermite (1, ephi.rotjds, &ephi.rotangv[j][0], &z,  ephi.nvrot, 10, itime, &arotang[j], &vz, ircode);
      if (*ircode != 0)
	{
	  if (debug)
	    printf ("hermite error: rotang %d\n", j);
	  return;
	}
    }
  if (debug)
    printf ("rot: %12.6f %12.6f %12.6f\n",
	    arotang[0], arotang[1], arotang[2]);
}

/*-----------------------------------------------------------------------
      void interp_rot (jdint, seceph, arottype, agast, arotang,
                            ircode)

  Interpolate rotation angles (euler or north pole RA, Dec, and W)
  from cpf file arrays stored in memory.

       jdint        Julian date (integer day)
       seceph       Julian day  (seconds)
       arotype      Rotation type (int): 0=n/a; 1=Lunar euler; 2=ra/dec/W
       agast        Greenwich apparent sidereal time
       arotang      x, y, z, or ra, dec, w (3 elements, d.p. in deg)
       ircode       error code: 0=OK, 1,2=error (time out of range)

  Author: R. L. Ricklefs / UT CSR
  Date: 1 May, 2004
  History:
-----------------------------------------------------------------------*/
void cinterp_eop (double jdint, double seceph, double x, double y, double dutc, int *ircode, eh &ephihdr, ei &ephi)
{
  double itime, z, vz;		/* z is dummy for derivative field */
  int debug = 0;
  //void hermite ();

  itime = ((seceph / 86400.e0 - ephi.jdf[0]) + (jdint - ephi.jdi[0])) * 86400.e0;

  hermite( 1, ephi.eopjds, &ephi.eopxy[0][0], &z, ephi.nveop, 10, itime, &x, &vz, ircode);
  if (*ircode != 0)
    {
      if (debug)
	printf ("hermite error: eop x\n");
      return;
    }
  hermite( 1, ephi.eopjds, &ephi.eopxy[1][0], &z, ephi.nveop, 10, itime, &y, &vz, ircode);
  if (*ircode != 0)
    {
      if (debug)
	printf ("hermite error: eop y\n");
      return;
    }
 // hermite( 1, ephi.eopjds, &ephi.eopdutc, &z, ephi.nveop, 10, itime, &dutc, ircode);
 // if (*ircode != 0)
 //   {
 //     if (debug)
	//printf ("hermite error: eop dutc\n");
 //     return;
 //   }
}
