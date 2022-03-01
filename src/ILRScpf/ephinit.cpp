#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ephinit.h"

/*-----------------------------------------------------------------------
  int cephinit (ephname)

  Open cpf ephemeris file and read into memory.
  Convenience function using char name as input; calls ephinitu.

  Input
       ephname - Name of cpf file (character string)

  Output
       ierr    - error code: 0= OK, -1= error

  Author: R. L. Ricklefs / UT CSR
  Date: 1 Apr, 2005
  History:
  12/06/2005 - Corrected handing of corr12 and corr23 relativistic correction
               to match the FORTRAN. rlr.

-----------------------------------------------------------------------*/
int cephinit (char *str, eh &ephihdr, ei &ephi )
{
  FILE *str_in;

  if ((str_in = fopen (str, "r")) == NULL)
    {
      return (-1);
    }

  cephinitu (str_in, ephihdr, ephi );
  fclose (str_in);
  return (0);
}

/*-----------------------------------------------------------------------
 int cephinitu (unit)

  Read cpf ephemeris file into memory.

  Input
       unit    - fortran unit number

  Output
       ierr    - error code: 0= OK, -1= error
       contents of /ephihdr/ (header info)
       contents of /ephi/ (record info)

  Author: R. L. Ricklefs / UT CSR
  Date: 1 May 2003
  History: 
     02.02.2006 S.Riepl changed function return type to the same as
           given in this documentation above and changed scanf format string
           from i5 to %5d. Cleaned up some scanf statements producing
           warnings when compiled with -Wall option.

  Bug: function always returns 0, regardless of file reading error status!

-----------------------------------------------------------------------*/
int cephinitu (FILE * str_in, eh &ephihdr, ei &ephi )
{
  int i, j;
  int idir;
  int mjd;
  double sod;
  double dgast = 0.e0;
  double pi;
  char ephstr[256];

  /* initialize all the arrays, in case there is nothing to read */
  for (i = 0; i < MAXENT; i++)
    {
      for (j = 0; j < 6; j++)
	{
	  ephi.outv[j][i] = 0.e0;
	  ephi.backv[j][i] = 0.e0;
	}
      for (j = 0; j < 3; j++)
	{
	  ephi.aberoutv[j][i] = 0.e0;
	  ephi.aberbackv[j][i] = 0.e0;
	  ephi.offsetoutv[j][i] = 0.e0;
	  ephi.offsetbackv[j][i] = 0.e0;
	  ephi.rotangv[j][i] = 0.e0;
	}
      ephi.c12[i] = 0.e0;
      ephi.c23[i] = 0.e0;
      ephi.leapflag[i] = 0;
      ephi.oscrelcorr[i] = 0.e0;
      ephi.gast[i] = 0.e0;
      ephi.eopxy[0][i] = 0.e0;
      ephi.eopxy[1][i] = 0.e0;
      ephi.eopdutc[i] = 0.e0;
    }

  ephihdr.atrk_ro_0 = 0;
  ephihdr.xtrk_ro_0 = 0;
  ephihdr.rtrk_ro_0 = 0;
  ephihdr.atrk_ro_1 = 0;
  ephihdr.xtrk_ro_1 = 0;
  ephihdr.rtrk_ro_1 = 0;
  ephihdr.xtrk_ro_2 = 0;
  ephihdr.rtrk_ro_2 = 0;
  ephihdr.prf = 0.e0;
  ephihdr.xpond_xmit_delay = 0.e0;
  ephihdr.xpond_utc_off = 0.e0;
  ephihdr.xpond_osc_drift = 0.e0;
  ephihdr.cofm_corr = 0;

  ephi.nv = -1;
  ephi.nvoff = -1;
  ephi.nvrot = -1;
  ephi.nveop = -1;

  pi = 4.e0 * atan (1.e0);

  while (fgets (ephstr, 256, str_in) != NULL)
    {
      if (strncmp (ephstr, "h1", 2) == 0 || strncmp (ephstr, "H1", 2) == 0)
	{
	  sscanf (&ephstr[7], "%2d  %3c %4d %2d %2d %2d  %4d %10c %10c",
		  &ephihdr.version, ephihdr.eph_source, &ephihdr.epyear,
		  &ephihdr.epmon, &ephihdr.epday, &ephihdr.ephour,
		  &ephihdr.eph_seq, ephihdr.tar_name, ephihdr.eph_notes);
	}
      else if (strncmp (ephstr, "h2", 2) == 0
	       || strncmp (ephstr, "H2", 2) == 0)
	{
	  sscanf (&ephstr[3],
		  "%8d %4d %8d %4d %2d %2d %2d %2d %2d %4d %2d %2d %2d %2d %2d %5d %1d %1d %2d %1d %1d",
		  &ephihdr.cospar_id, &ephihdr.sic, &ephihdr.norad,
		  &ephihdr.sttpyear, &ephihdr.sttpmonth,
		  &ephihdr.sttpday, &ephihdr.sttphour, &ephihdr.sttpmin,
		  &ephihdr.sttpsec, &ephihdr.endpyear, &ephihdr.endpmonth,
		  &ephihdr.endpday, &ephihdr.endphour, &ephihdr.endpmin,
		  &ephihdr.endpsec, &ephihdr.ephsep, &ephihdr.compat,
		  &ephihdr.tar_type, &ephihdr.ref_frame, &ephihdr.rot_type,
		  &ephihdr.cofm_app);
	}
      else if (strncmp (ephstr, "h3", 2) == 0
	       || strncmp (ephstr, "H3", 2) == 0)
	{
	  sscanf (&ephstr[3],
            "%5d %5d %5d %5d %5d %5d %5d %5d %5d",
		  &ephihdr.atrk_ro_0, &ephihdr.xtrk_ro_0, &ephihdr.rtrk_ro_0,
		  &ephihdr.atrk_ro_1, &ephihdr.xtrk_ro_1, &ephihdr.rtrk_ro_1,
		  &ephihdr.atrk_ro_2, &ephihdr.xtrk_ro_2, &ephihdr.rtrk_ro_2);
	}
      else if (strncmp (ephstr, "h4", 2) == 0
	       || strncmp (ephstr, "H4", 2) == 0)
	{
	  sscanf (&ephstr[3],
		  "%lf %lf %lf %lf",
		  &ephihdr.prf, &ephihdr.xpond_xmit_delay,
		  &ephihdr.xpond_utc_off, &ephihdr.xpond_osc_drift);
	}
      else if (strncmp (ephstr, "h5", 2) == 0
	       || strncmp (ephstr, "H5", 2) == 0)
	{
	  /*sscanf (&ephstr[3], "%7.4f", &ephihdr.cofm_corr); */
	  sscanf (&ephstr[3], "%lf", &ephihdr.cofm_corr);
	}
      else if (strncmp (ephstr, "h9", 2) == 0
	       || strncmp (ephstr, "H9", 2) == 0)
	{
	}
      else if (strncmp (ephstr, "10", 2) == 0)
	{
	  sscanf (&ephstr[3], "%d", &idir);
	  if (idir == 0 || idir == 1)
	    {
	      (ephi.nv)++;
				  /**printf("nv= %d\n",ephi.nv);**/
	      if (ephi.nv > MAXENT)
		{
		  printf ("ephemeris read error: too many records\n");
		  break;
		}
	      sscanf (&ephstr[3], "%d %d %lf %d %lf %lf %lf", &ephihdr.ltcorr,	/* = idir */
		      &mjd, &sod, &ephi.leapflag[ephi.nv],
		      &ephi.outv[0][ephi.nv], &ephi.outv[1][ephi.nv],
		      &ephi.outv[2][ephi.nv]);
/**
	  printf(
	    "%ld %f %d %lf %lf %lf",
	    mjd,sod,ephi.leapflag[ephi.nv],ephi.outv[0][ephi.nv],
	    ephi.outv[1][ephi.nv],ephi.outv[2][ephi.nv]);
**/
	      ephi.jdi[ephi.nv] = mjd + 2400000.0e0;
	      ephi.jdf[ephi.nv] = sod / 86400.0e0 + 0.5e0;
	      if (ephi.jdf[ephi.nv] > 1.0e0)
		{
		  ephi.jdi[ephi.nv] += 1.e0;
		  ephi.jdf[ephi.nv] -= 1.e0;
		}
	      /* In case there is not back vector (e.g., satellites)... */
	      ephi.backv[0][ephi.nv] = -ephi.outv[0][ephi.nv];
	      ephi.backv[1][ephi.nv] = -ephi.outv[1][ephi.nv];
	      ephi.backv[2][ephi.nv] = -ephi.outv[2][ephi.nv];
	    }
	  else if (idir == 2)
	    {
	      /* For now, ignore the return time and leapflag. */
	      /* Assume index nv is the same. */
	      sscanf (&ephstr[3],
		      "%d %d %lf %d %lf %lf %lf", &idir,
		      &mjd, &sod, &ephi.leapflag[ephi.nv],
		      &ephi.backv[0][ephi.nv], &ephi.backv[1][ephi.nv],
		      &ephi.backv[2][ephi.nv]);
	      if (ephi.backv[0][ephi.nv] * ephi.outv[0][ephi.nv] > 0.e0)
		{
		  ephi.backv[0][ephi.nv] = -ephi.backv[0][ephi.nv];
		  ephi.backv[1][ephi.nv] = -ephi.backv[1][ephi.nv];
		  ephi.backv[2][ephi.nv] = -ephi.backv[2][ephi.nv];
		}
	    }
	}
      else if (strncmp (ephstr, "20", 2) == 0)
	{
	  sscanf (&ephstr[3], "%d", &idir);
	  if (idir == 0 || idir == 1)
	    {
	      sscanf (&ephstr[3],
		      "%d %lf %lf %lf", &idir,
		      &ephi.outv[3][ephi.nv], &ephi.outv[4][ephi.nv],
		      &ephi.outv[5][ephi.nv]);
	      /* In case there is not back vector (e.g., satellites)... */
	      ephi.backv[3][ephi.nv] = -ephi.outv[3][ephi.nv];
	      ephi.backv[4][ephi.nv] = -ephi.outv[4][ephi.nv];
	      ephi.backv[5][ephi.nv] = -ephi.outv[5][ephi.nv];
	    }
	  else if (idir == 2)
	    {
	      sscanf (&ephstr[3],
		      "%d %lf %lf %lf", &idir,
		      &ephi.backv[0][ephi.nv], &ephi.backv[1][ephi.nv],
		      &ephi.backv[2][ephi.nv]);
	    }
	}
      else if (strncmp (ephstr, "30", 2) == 0)
	{
	  sscanf (&ephstr[3], "%d", &idir);
	  if (idir == 0 || idir == 1)
	    {
	      sscanf (&ephstr[3],
		      "%d %lf %lf %lf %lf", &idir,
		      &ephi.aberoutv[0][ephi.nv], &ephi.aberoutv[1][ephi.nv],
		      &ephi.aberoutv[2][ephi.nv], &ephi.c12[ephi.nv]);
	      /* This is a scalar: Be sure it is positive... */
	      ephi.c12[ephi.nv] = fabs (ephi.c12[ephi.nv]);
	      /* In case there is no return information... */
	      ephi.aberbackv[0][ephi.nv] = -ephi.aberoutv[0][ephi.nv];
	      ephi.aberbackv[1][ephi.nv] = -ephi.aberoutv[1][ephi.nv];
	      ephi.aberbackv[2][ephi.nv] = -ephi.aberoutv[2][ephi.nv];
	      ephi.c23[ephi.nv] = -ephi.c12[ephi.nv];	/* scalar */
	    }
	  else if (idir == 2)
	    {
	      sscanf (&ephstr[3],
		      "%d %lf %lf %lf %lf", &idir,
		      &ephi.aberbackv[0][ephi.nv],
		      &ephi.aberbackv[1][ephi.nv],
		      &ephi.aberbackv[2][ephi.nv], &ephi.c23[ephi.nv]);
	      /* These are vectors: outbound and back vectors should be of
	       * opposite sign. */
	      if (ephi.backv[0][ephi.nv] * ephi.outv[0][ephi.nv] > 0.e0)
		{
		  ephi.aberbackv[0][ephi.nv] = -ephi.aberbackv[0][ephi.nv];
		  ephi.aberbackv[1][ephi.nv] = -ephi.aberbackv[1][ephi.nv];
		  ephi.aberbackv[2][ephi.nv] = -ephi.aberbackv[2][ephi.nv];
		}
	      /* This is a scalar: Be sure it is positive... */
	      ephi.c23[ephi.nv] = ephi.c23[ephi.nv];
	    }
	  /* WEAK TEST. REDO THIS! */
	  if (ephi.c23[ephi.nv] * ephi.c12[ephi.nv] > 0.e0)
	    ephi.c23[ephi.nv] = -ephi.c12[ephi.nv];
	}
      else if (strncmp (ephstr, "40", 2) == 0)
	{
	  sscanf (&ephstr[3], "%lf", &ephi.oscrelcorr[ephi.nv]);
	}
      else if (strncmp (ephstr, "50", 2) == 0)
	{
	  sscanf (&ephstr[3], "%d", &idir);
	  if (idir == 0 || idir == 1)
	    {
	      (ephi.nvoff)++;
	      if (ephi.nvoff > MAXENT)
		{
		  printf ("ephemeris read error: too many records\n");
		  break;
		}
	      sscanf (&ephstr[3],
		      "%d %d %lf %lf %lf %lf", &idir,
		      &mjd, &sod, &ephi.offsetoutv[0][ephi.nvoff],
		      &ephi.offsetoutv[1][ephi.nvoff],
		      &ephi.offsetoutv[2][ephi.nvoff]);
/**
	  printf(
	    "%ld %f %d %lf %lf %lf",
	    mjd,sod,ephi.leapflag[ephi.nv],ephi.offsetoutv[0][ephi.nv],
	    ephi.offsetoutv[1][ephi.nv],ephi.offsetoutv[2][ephi.nv]);
**/
	      ephi.offsetjdi[ephi.nvoff] = mjd + 2400000.0e0;
	      ephi.offsetjdf[ephi.nvoff] = sod / 86400.0e0 + 0.5e0;
	      if (ephi.offsetjdf[ephi.nvoff] > 1.0e0)
		{
		  ephi.offsetjdi[ephi.nvoff] += 1.e0;
		  ephi.offsetjdf[ephi.nvoff] -= 1.e0;
		}
	      /* In case there is not back vector */
	      ephi.offsetbackv[0][ephi.nv] = -ephi.offsetoutv[0][ephi.nv];
	      ephi.offsetbackv[1][ephi.nv] = -ephi.offsetoutv[1][ephi.nv];
	      ephi.offsetbackv[2][ephi.nv] = -ephi.offsetoutv[2][ephi.nv];
	    }
	  else if (idir == 2)
	    {
	      sscanf (&ephstr[3],
		      "%d %d %lf %lf %lf %lf", &idir,
		      &mjd, &sod, &ephi.offsetbackv[0][ephi.nv],
		      &ephi.offsetbackv[1][ephi.nv],
		      &ephi.offsetbackv[2][ephi.nv]);
	      if (ephi.offsetbackv[0][ephi.nvoff] *
		  ephi.offsetoutv[0][ephi.nvoff] > 0.e0)
		{
		  ephi.offsetbackv[0][ephi.nvoff] =
		    -ephi.offsetbackv[0][ephi.nvoff];
		  ephi.offsetbackv[1][ephi.nvoff] =
		    -ephi.offsetbackv[1][ephi.nvoff];
		  ephi.offsetbackv[2][ephi.nvoff] =
		    -ephi.offsetbackv[2][ephi.nvoff];
		}
	    }
	}
      else if (strncmp (ephstr, "60", 2) == 0)
	{
	  (ephi.nvrot)++;
	  if (ephi.nvrot > MAXENT)
	    {
	      printf ("ephemeris read error: too many records\n");
	      break;
	    }
	  sscanf (&ephstr[3],
		  "%d %lf %lf %lf %lf %lf",
		  &mjd, &sod, &ephi.rotangv[0][ephi.nvrot],
		  &ephi.rotangv[1][ephi.nvrot], &ephi.rotangv[2][ephi.nvrot],
		  &ephi.gast[ephi.nvrot]);
/**
	  printf(
	    "%ld %lf %lf %lf %lf %lf\n",
	    mjd,sod,ephi.rotangv[0][ephi.nvrot],ephi.rotangv[1][ephi.nvrot],ephi.rotangv[2][ephi.nvrot],ephi.gast[ephi.nvrot]);
**/
	  ephi.rotjdi[ephi.nvrot] = mjd + 2400000.0e0;
	  ephi.rotjdf[ephi.nvrot] = sod / 86400.0e0 + 0.5e0;
	  if (ephi.rotjdf[ephi.nvrot] > 1.0e0)
	    {
	      ephi.rotjdi[ephi.nvrot] += 1.e0;
	      ephi.rotjdf[ephi.nvrot] -= 1.e0;
	    }
/*       Sidereal time wraps around from 24 to 0 hours... */
	  ephi.gast[ephi.nvrot] += dgast;

	  if (ephi.nvrot > 0
	      && ephi.gast[ephi.nvrot] < ephi.gast[ephi.nvrot - 1])
	    {
	      dgast += 24.e0;
	      ephi.gast[ephi.nvrot] += 24.e0;
	    }
	  for (i = 0; i < 3; i++)
	    {
	      ephi.rotangv[i][ephi.nvrot] *= pi / 180.0e0;
	    }

	}
      else if (strncmp (ephstr, "70", 2) == 0)
	{
	  (ephi.nveop)++;
	  if (ephi.nveop > MAXENT)
	    {
	      printf ("ephemeris read error: too many records\n");
	      break;
	    }
	  sscanf (&ephstr[3],
		  "%d %lf %lf %lf %lf",
		  &mjd, &sod, &ephi.eopxy[0][ephi.nveop],
		  &ephi.eopxy[1][ephi.nveop], &ephi.eopdutc[ephi.nveop]);
/**
	  printf(
	    "%ld %f %d %lf %lf %lf",
	    mjd,sod,ephi.leapflag[ephi.nv],ephi.outv[0][ephi.nv],
	    ephi.outv[1][ephi.nv],ephi.outv[2][ephi.nv]);
**/
	  ephi.eopjdi[ephi.nveop] = mjd + 2400000.0e0;
	  ephi.eopjdf[ephi.nveop] = sod / 86400.0e0 + 0.5e0;
	  if (ephi.eopjdf[ephi.nveop] > 1.0e0)
	    {
	      ephi.eopjdi[ephi.nveop] += 1.e0;
	      ephi.eopjdf[ephi.nveop] -= 1.e0;
	    }
	}
      else if (strncmp (ephstr, "99", 2) == 0)
	{
	}
      /* Comment */
      else if (strncmp (ephstr, "00", 2) == 0)
	{
	}
    }

/*     create a small usable time tag */
  (ephi.nv)++;
  for (i = 0; i < ephi.nv; i++)
    {
      ephi.jds[i] = ((ephi.jdf[i] - ephi.jdf[0]) +
		     (ephi.jdi[i] - ephi.jdi[0])) * 86400.e0;
    }
  (ephi.nvoff)++;
  for (i = 0; i < ephi.nvoff; i++)
    {
      ephi.offsetjds[i] = ((ephi.offsetjdf[i] - ephi.offsetjdf[0]) +
			   (ephi.offsetjdi[i] -
			    ephi.offsetjdi[0])) * 86400.e0;
    }
  (ephi.nvrot)++;
  for (i = 0; i < ephi.nvrot; i++)
    {
      ephi.rotjds[i] = ((ephi.rotjdf[i] - ephi.rotjdf[0]) +
			(ephi.rotjdi[i] - ephi.rotjdi[0])) * 86400.e0;
    }
  (ephi.nveop)++;
  for (i = 0; i < ephi.nveop; i++)
    {
      ephi.eopjds[i] = ((ephi.eopjdf[i] - ephi.eopjdf[0]) +
			(ephi.eopjdi[i] - ephi.eopjdi[0])) * 86400.e0;
    }
  return (0);
}
