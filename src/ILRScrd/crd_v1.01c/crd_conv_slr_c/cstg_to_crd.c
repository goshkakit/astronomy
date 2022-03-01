#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define FNLEN   256
#define FMT_VERSION 0

#include "crd.h"
#include "cstg.h"

char *version = "1.01b 02/02/2012";
char *site = "UTX";
FILE *str_in, *str_out_np, *str_out_qlk;
fpos_t startpos;
struct cstg_hdr hdr;
struct cstg_sed sed;
struct cstg_np np;
struct rh1 h1;
struct rh2 h2;
struct rh3 h3;
struct rh4 h4;
struct rc0 c0;
struct rc1 c1;
struct rc2 c2;
struct rc3 c3;
struct rc4 c4;
struct rd10 d10;
struct rd11 d11;
struct rd12 d12;
struct rd20 d20;
struct rd21 d21;
struct rd30 d30;
struct rd40 d40;
struct rd50 d50;
struct rd60 d60;
struct rd00 d00;

void get_sat_ids ();
void get_station_id ();
void move_hdr_to_crd ();
void move_np_to_crd ();
void move_sed_to_crd ();
void read_cstg_hdr ();
void read_cstg_np ();
void read_cstg_sampled_eng ();
void setup_files ();
void sodtohms ();
void write_headers ();
void write_data ();
void write_end_of_data_block ();
void write_end_of_file ();

/*-------------------------------------------------------------------------
 * Program: cstg_to_crd
 *
 * Purpose:
 * Converts old ILRS (CSTG) quicklook and normalpoint files to the 
 * ILRS CRD format
 *
 * Calling sequence:
 *   cstg_to_crd -i cstg_np_qlk_filename -o crd_filename_base
 *	where crd_filename_base has no extension. The extension
 *	'.npt' is added for the normalpoint file and
 *	'.qlk' is addded for the quicklook (sampled engineering) file.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   June 29, 2007 - Initial version
 *   March 2, 2011 - Add comment to converted file to distinguish from native
 *                   CRD creation. NOTE: User should change "site" to reflect 
 *                   their site or organization. rlr.
 *   Feb 2, 2010 - Because of round-off problems, change 1.e7 to
 *                 (long double)1.e7 in second of days calulations.
 *                 1.01b. rlr.
 *
**-----------------------------------------------------------------------*/
int
main (argc, argv)
     int argc;
     char *argv[];
{
  int first_record = 0;
  int data_type = 0;
  long data_release = 0;
  long filpos;
  int n_np = 0, n_sed = 0;
  int reset_def = 1;
  int status;
  time_t lclock;
  char str[256];
  double sec_of_day;
  double ssec_of_day;
  double esec_of_day;
  struct tm gmt;

/*  Greetings */
  printf("Convert CSTG Normalpoints and Quicklook to CRD Format: version %s\n",version);

/*  get file names and open files  */
  setup_files (argc, argv);

/* get processing time */
  time (&lclock);
  gmt = *gmtime (&lclock);

/* read and process the file... */
  while ((status = fgets (str, 256, str_in)) != NULL)
    {
      /* Section separator */
/* =========================== Normal points  ============================= */
      if (strncmp (str, "99999", 5) == 0)
	{
	  /* Normal points */
	  data_type = 1;

	  /* Header */
	  if ((status = fgets (str, 256, str_in)) != NULL)
	    {
	      if (strlen (str) <= 1)
		continue;	/* Blank lines... */

	      /**printf ("header: [%s]\n", str); **/
	      read_cstg_hdr (str, &hdr);
	      first_record = 1;

	      /* Get the start and end times of file (and release version) */
	      fgetpos (str_in, &startpos);
	      while ((status = fgets (str, 256, str_in)) != NULL
		     && strncmp (str, "88888", 5) != 0
		     && strncmp (str, "99999", 5) != 0)
		{

		  if (strlen (str) <= 1)
		    continue;	/* Blank lines... */
		  getfield (str, 0, 12, &sec_of_day);
		  getifield (str, 47, 1, &data_release);
		  sec_of_day /= 1.e7;
		  if (first_record)
		    {
		      ssec_of_day = sec_of_day;
		      first_record = 0;
		    }
		  else
		    {
		      esec_of_day = sec_of_day;
		    }
		}
	      /**printf ("s sec %f e sec %f\n", ssec_of_day,
	         esec_of_day); **/
	      fsetpos (str_in, &startpos);

	      /* Write the header(s) */
	      move_hdr_to_crd (gmt, data_type, ssec_of_day, esec_of_day,
			       data_release);
	      write_headers (str_out_np);

	      /* Process the normalpoitn records */
	      reset_def = 1;
	      while ((status = fgets (str, 256, str_in)) != NULL
		     && strncmp (str, "88888", 5) != 0
		     && strncmp (str, "99999", 5) != 0)
		{
		  if (strlen (str) <= 1)
		    continue;	/* Blank lines... */
		  read_cstg_np (str, &np);
		  move_np_to_crd ();
		  write_data (str_out_np, reset_def);
		  reset_def = 0;
		}

	      /* May need to "unread" that last record... */
	      if (strncmp (str, "88888", 5) == 0
               || strncmp (str, "99999", 5) == 0)
		fseek (str_in, (long) (-6), SEEK_CUR);

	      /* end the session... */
	      write_end_of_data_block (str_out_np);
	      n_np++;
	    }
	}
/* =========================== Eng Quicklook ============================= */
      else if (strncmp (str, "88888", 5) == 0)
	{
	  /* Sampled Eng = quicklook */
	  data_type = 2;

	  /* Header */
          filpos = ftell (str_in);

	  if ((status = fgets (str, 256, str_in)) != NULL)
	    {
	      read_cstg_hdr (str, &hdr);
	      first_record = 1;

	      /* read to end of file to get final date and time */
	      while ((status = fgets (str, 256, str_in)) != NULL)
		{
		  getfield (str, 0, 12, &sec_of_day);
		  sec_of_day /= 1.e7;
		  if (first_record)
		    {
		      ssec_of_day = sec_of_day;
		      first_record = 0;
		    }
		  else
		    {
		      esec_of_day = sec_of_day;
		    }
		}
              fseek (str_in, filpos, SEEK_SET);
	      fgets (str, 256, str_in);
	    }
	  /* Write headers */
	  move_hdr_to_crd (gmt, data_type, ssec_of_day, esec_of_day, 0);
	  write_headers (str_out_qlk);

	  /*Process all the sampled eng data records */
	  reset_def = 1;
	  while ((status = fgets (str, 256, str_in)) != NULL)
	    {
	      read_cstg_sampled_eng (str, &sed);
	      move_sed_to_crd (ssec_of_day);
	      write_data (str_out_qlk, reset_def);
	      reset_def = 0;
	    }

	  /* and end the session */
	  write_end_of_data_block (str_out_qlk);
	  n_sed++;
	}
    }

  /* clean up and go home... */
  if (n_np > 0)
    write_end_of_file (str_out_np);
  if (n_sed > 0)
    write_end_of_file (str_out_qlk);
  fclose (str_in);
  fclose (str_out_np);
  fclose (str_out_qlk);
}

/*-------------------------------------------------------------------------
**
**	move_hdr_to_crd - Copy info from the input header records 
**			  and elsewhere to the CRD header records.
**
**-----------------------------------------------------------------------*/
void
move_hdr_to_crd (struct tm gmt, int data_type, double ssec_of_day,
		 double esec_of_day, int data_release)
{
  char mon_name[4];
  int mon, day;
  int i;

/* Comment */
  sprintf (d00.comment,"Converted from the CSTG format by %s using cstg_to_crd V%s", site, version);

/* Format header */
  strcpy (h1.crd_literal,"CRD");
  h1.format_version = FMT_VERSION;
  h1.prod_year = gmt.tm_year + 1900;
  h1.prod_mon = gmt.tm_mon + 1;
  h1.prod_day = gmt.tm_mday;
  h1.prod_hour = gmt.tm_hour;

/* Station header */
  get_station_id (hdr.cdp_pad_id, &h2.stn_name);
  h2.cdp_pad_id = hdr.cdp_pad_id;
  h2.cdp_sys_num = hdr.cdp_sys_num;
  h2.cdp_occ_num = hdr.cdp_occ_num;
  h2.stn_timescale = hdr.stn_timescale;

/* Target header */
  h3.ilrs_id = hdr.ilrs_id;
  get_sat_ids (0, &h3.ilrs_id, &h3.norad, &h3.sic, &h3.target_name);
  h3.SC_timescale = 0;
  if (hdr.np_window_ind == 2)	/* llr */
    h3.target_type = 2;
  else
    h3.target_type = 1;

  h4.data_type = data_type;
  doytogr (hdr.year - 1900, hdr.doy, &mon, &day, mon_name);
  h4.start_year = h4.end_year = hdr.year;
  h4.start_mon = h4.end_mon = mon;
  h4.start_day = h4.end_day = day;
  sodtohms ((long) ssec_of_day, &h4.start_hour, &h4.start_min, &h4.start_sec);
  sodtohms ((long) esec_of_day, &h4.end_hour, &h4.end_min, &h4.end_sec);
  h4.data_release = data_release;	/* from np detail recds; none for sed */
  h4.refraction_app_ind = 0;
  h4.CofM_app_ind = 0;
  h4.xcv_amp_app_ind = 0;
  if (data_type == 1)
    {
      h4.stn_sysdelay_app_ind = 1;
    }
  else
    {
      h4.stn_sysdelay_app_ind = 0;
    }
  h4.SC_sysdelay_app_ind = 0;
  h4.range_type_ind = 2;
  h4.data_qual_alert_ind= 0;

/* Config header */
  c0.detail_type = 0;
  c0.xmit_wavelength = hdr.xmit_wavelength;
  strcpy (c0.config_ids[0], "std");
  for (i = 1; i < 10; i++)
    {
      c0.config_ids[i][0] = '\0';
    }

/* Normalpoint recd */
/* Calibratin Record */
  d40.type_of_data = 0;
  strcpy (d40.sysconfig_id, c0.config_ids[0]);
  d40.num_points_recorded = -1;
  d40.num_points_used = -1;
  d40.one_way_target_dist = -1;
  d40.cal_sys_delay = hdr.cal_sys_delay;
  d40.cal_delay_shift = hdr.cal_delay_shift;
  d40.cal_rms = hdr.cal_rms;
  d40.cal_skew = -1;
  d40.cal_kurtosis = -1;
  d40.cal_PmM = -1;
  d40.cal_shift_type_ind = hdr.cal_type_ind > 4 ? 3 : 2;
  d40.cal_type_ind = hdr.cal_type_ind;
  if (d40.cal_type_ind > 4)
    d40.cal_type_ind -= 5;
  d40.cal_type_ind += 2;
  if (d40.cal_type_ind == 6)
    d40.cal_type_ind = 0;
  d40.detector_channel= 0;

/* Session (pass) statistical Record */
  strcpy (d50.sysconfig_id, c0.config_ids[0]);
  d50.sess_rms = hdr.sess_rms;
  d50.sess_skew = -1;
  d50.sess_kurtosis = -1;
  d50.sess_PmM = -1;
  d50.data_qual_ind = hdr.data_qual_ind;

/* Compatibility record */
  strcpy (d60.sysconfig_id, c0.config_ids[0]);
  d60.sys_change_ind = hdr.sys_change_ind;
  d60.sys_config_ind = hdr.sys_config_ind;
}

/*-------------------------------------------------------------------------
**
**	move_np_to_crd - Copy info from the input normalpoint records 
**			 and elsewhere to the CRD data records.
**
**-----------------------------------------------------------------------*/
void
move_np_to_crd ()
{
/* Normalpoint Record */
  d11.sec_of_day = np.sec_of_day;
  d11.time_of_flight = np.time_of_flight;
  strcpy (d11.sysconfig_id, c0.config_ids[0]);
  d11.epoch_event = 2;
  d11.num_ranges = np.num_ranges;

  if (hdr.np_window_ind == 1)
    d11.np_window_length = 5;
  else if (hdr.np_window_ind == 2)	/* llr */
    {			   /** is doc correct? Is it really 10? **/
      if (np.llr_np_window_ind == 1)
	d11.np_window_length = 300;
      else if (np.llr_np_window_ind == 2)
	d11.np_window_length = 600;
      else if (np.llr_np_window_ind == 3)
	d11.np_window_length = 900;
      else if (np.llr_np_window_ind == 4)
	d11.np_window_length = 1200;
      else if (np.llr_np_window_ind == 5)
	d11.np_window_length = 1500;
      else if (np.llr_np_window_ind == 6)
	d11.np_window_length = 1800;
      else if (np.llr_np_window_ind == 7)
	d11.np_window_length = 2100;
      else if (np.llr_np_window_ind == 8)
	d11.np_window_length = 2400;
      else if (np.llr_np_window_ind == 9)
	d11.np_window_length = 3000;
    }
  else if (hdr.np_window_ind == 3)
    d11.np_window_length = 15;
  else if (hdr.np_window_ind == 4)
    d11.np_window_length = 20;
  else if (hdr.np_window_ind == 5)
    d11.np_window_length = 30;
  else if (hdr.np_window_ind == 6)
    d11.np_window_length = 60;
  else if (hdr.np_window_ind == 7)
    d11.np_window_length = 120;
  else if (hdr.np_window_ind == 8)
    d11.np_window_length = 180;
  else if (hdr.np_window_ind == 9)
    d11.np_window_length = 300;

  if (hdr.np_window_ind == 2)	/* llr */
    d11.time_of_flight += np.scale_or_tof_sec;
  else if (np.scale_or_tof_sec != 0)
    d11.num_ranges *= exp10 ((double) np.scale_or_tof_sec);

  d11.bin_rms = np.bin_rms;
  d11.bin_skew = -1;
  d11.bin_kurtosis = -1;
  d11.bin_PmM = -1;
  d11.return_rate = np.snr;
  d11.detector_channel= 0;

/* Meteorology Record */
  d20.sec_of_day = np.sec_of_day;
  d20.pressure = np.pressure;
  d20.temperature = np.temperature;
  d20.humidity = np.humidity;
  d20.value_origin = 1;

/* Calibration Record */
  d40.sec_of_day = np.sec_of_day;
}

/*-------------------------------------------------------------------------
**
**	move_sed_to_crd -  Copy info from the input sampled enineering 
**			   records and elsewhere to the CRD data records.
**
**-----------------------------------------------------------------------*/
void
move_sed_to_crd (double ssec_of_day)
{
/* Range Record */
  d10.sec_of_day = sed.sec_of_day;
  d10.time_of_flight = sed.time_of_flight;
  strcpy (d10.sysconfig_id, c0.config_ids[0]);
  d10.epoch_event = 2;
  d10.filter_flag = 0;
  d10.detector_channel = 0;
  d10.stop_number = 0;
  d10.xcv_amp = sed.xcv_amp;

/* Meteorology Record */
  d20.sec_of_day = sed.sec_of_day;
  d20.pressure = sed.pressure;
  d20.temperature = sed.temperature;
  d20.humidity = sed.humidity;
  d20.value_origin = 1;

/* POinting Angle Record */
  d30.sec_of_day = sed.sec_of_day;
  d30.azimuth = sed.azimuth;
  d30.elevation = sed.elevation;
  d30.direction_ind = 0;
  d30.angle_origin_ind = sed.angle_origin_ind;
  d30.refraction_corr_ind = 1;	/* Assumption! */

/* Calibration Record */
  d40.sec_of_day = ssec_of_day;
}

/*-------------------------------------------------------------------------
**
**	sodtohms - Convert seconds of day to hour/minute/second 
**
**-----------------------------------------------------------------------*/
void
sodtohms (long sod, int *h, int *m, int *s)
{
  *h = sod / 3600;
  *m = (sod - *h * 3600) / 60;
  *s = sod - *h * 3600 - *m * 60;
}

/*-------------------------------------------------------------------------
**
**	write_headers - Write CRD header records
**
**-----------------------------------------------------------------------*/
void
write_headers (FILE * str_out_crd)
{
  write_00 (str_out_crd, d00);
  write_h1 (str_out_crd, h1);
  write_h2 (str_out_crd, h2);
  write_h3 (str_out_crd, h3);
  write_h4 (str_out_crd, h4);
  write_c0 (str_out_crd, c0);
  write_60 (str_out_crd, d60);
}

/*-------------------------------------------------------------------------
**
**	write_data - Write CRD Data records
**
**-----------------------------------------------------------------------*/
void
write_data (FILE * str_out_crd, int reset_def)
{
  static double opressure = -999;
  static double otemperature = -999;
  static double ohumidity = -999;
  static double orefraction_corr = -999;
  static double otarget_CofM_corr = -999;
  static double ond_value = -999;
  static double otime_bias = -999;
  static double oazimuth = -999;
  static double oelevation = -999;
  static double ocal_delay_shift = -999;

  if (reset_def)
    {
      opressure = otemperature = ohumidity = -999;
      orefraction_corr = otarget_CofM_corr = ond_value = otime_bias = -999;
      oazimuth = oelevation = ocal_delay_shift = -999;
    }

  if (h4.data_type != 1)
    write_10 (str_out_crd, d10);
  else
    write_11 (str_out_crd, d11);

  /* For existing real data! Normally '5' */
  if (fabs (d20.pressure - opressure) > 0.009 
   || fabs (d20.temperature - otemperature) > 0.09
   || fabs (d20.humidity - ohumidity) > 0.9)
    {
      write_20 (str_out_crd, d20);
      opressure = d20.pressure;
      otemperature = d20.temperature;
      ohumidity = d20.humidity;
    }

  if ((fabs (d30.azimuth - oazimuth) > 0.1 ||
       fabs (d30.elevation - oelevation) > 0.1) && h4.data_type != 1)
    {
      write_30 (str_out_crd, d30);
      oazimuth = d30.azimuth;
      oelevation = d30.elevation;
    }

  if (fabs (d40.cal_delay_shift - ocal_delay_shift) > 1)
    {
      write_40 (str_out_crd, d40);
      ocal_delay_shift = d40.cal_delay_shift;
    }
}

/*-------------------------------------------------------------------------
**
**      write_end_of_data_block - Write CRD End-of-data-block record
**
**-----------------------------------------------------------------------*/
void
write_end_of_data_block (FILE * str_out_crd)
{
  write_50 (str_out_crd, d50);
  write_h8 (str_out_crd);
}

/*-------------------------------------------------------------------------
**
**      write_end_of_file - Write CRD End-of-file record
**
**-----------------------------------------------------------------------*/
void
write_end_of_file (FILE * str_out_crd)
{
  write_h9 (str_out_crd);
}

/*-------------------------------------------------------------------------
**
**	setup_files		- Open input and output file names
**
**-----------------------------------------------------------------------*/
void
setup_files (int argc, char *argv[])
{
  char qld_in[FNLEN] = { "" };
  char crd_out1[FNLEN] = { "" };
  char crd_out2[FNLEN] = { "" };

  int found = 0, i;

/*  Get file names  */
  for (i = 1; i < argc; i += 2)
    {
      if (argv[i][0] == '-')
	{
	  if (argv[i][1] == 'i')
	    {
	      strcpy (qld_in, argv[i + 1]);
	      printf ("qld_in = [%s]\n", qld_in);
	      found += 1;
	    }
	  else if (argv[i][1] == 'o')
	    {
	      if (found == 1)
		{
		  strcpy (crd_out1, argv[i + 1]);
		  strcat (crd_out1, ".npt");
		  strcpy (crd_out2, argv[i + 1]);
		  strcat (crd_out2, ".qlk");
		}
	      printf ("crd_out1 = [%s]\n", crd_out1);
	      printf ("crd_out2 = [%s]\n", crd_out2);
	      found += 2;
	    }
	}
    }
  if (found != 3)
    {
      printf ("Usage: cstg_to_crd -i cstg_qlk_np_filename -o crd_filename_base\n");
      exit (1);
    }

/*  open output QLD file  */
  if (strlen (qld_in) > 0)
    {
      if ((str_in = fopen (qld_in, "r")) == NULL)
	{
	  printf ("Could not open file %s\n", qld_in);
	  exit (1);
	}
    }

/*  open output CRD files  */
  if ((str_out_np = fopen (crd_out1, "w")) == NULL)
    {
      printf ("Could not open file %s\n", crd_out1);
      exit (1);
    }
  if ((str_out_qlk = fopen (crd_out2, "w")) == NULL)
    {
      printf ("Could not open file %s\n", crd_out2);
      exit (1);
    }
}

/*-------------------------------------------------------------------------
**
**	get_sat_ids - Given a laser target's ilrs ID, norad ID, sic, or
**		      name, get the other 3 IDs.
**
**-----------------------------------------------------------------------*/
void
get_sat_ids (int mode, int *ilrs_id, int *norad_id, int *sic, char *target)
{
  FILE *sat_id_in;
  char *sat_id_file = "./targets.dat";
  char str[256], ttarget[11];
  int status, tilrs_id, tsic, tnorad_id = 0;
  int i, l;

  if (mode != 3)
    {
      for (i = 0; i < 10; i++)
	{
	  target[i] = ' ';
	}
      target[10] = '\0';
    }

  if ((sat_id_in = fopen (sat_id_file, "r")) == NULL)
    {
      printf ("Could not open file %s\n", sat_id_file);
      exit (1);
    }
  while ((status = fgets (str, 256, sat_id_in)) != NULL)
    {
      sscanf (str, "%s %*s %d %d %d", ttarget, &tsic, &tilrs_id,
	      &tnorad_id);

      if (mode == 0 && tilrs_id == *ilrs_id)
	{
	  *norad_id = tnorad_id;
	  *sic = tsic;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  if (l < 10)
	    target[l] = ' ';
	  target[10] = '\0';
	  break;
	}
      else if (mode == 1 && tsic == *sic)
	{
	  *ilrs_id = tilrs_id;
	  *norad_id = tnorad_id;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  if (l < 10)
	    target[l] = ' ';
	  target[10] = '\0';
	  break;
	}
      else if (mode == 2 && tnorad_id == *norad_id)
	{
	  *ilrs_id = tilrs_id;
	  *sic = tsic;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  target[10] = '\0';
	  break;
	}
      else if (mode == 3 && strcmp (ttarget, target) == 0)
	{
	  *ilrs_id = tilrs_id;
	  *norad_id = tnorad_id;
	  *sic = tsic;
	  break;
	}
    }
  fclose (sat_id_in);
  /*printf ("got sat\n"); */
}

/*-------------------------------------------------------------------------
**
**	 get_station_id - get the 10 character station id by using
**                        the marker id and sites.dat file.
**
**-----------------------------------------------------------------------*/
void
get_station_id (int marker, char *station)
{
  FILE *sttn_id_in;
  char *sttn_id_file = "./sites.dat";
  char str[256], tstation[11];
  int status, tmarker;
  int i;

  for (i = 0; i < 10; i++)
    {
      station[i] = ' ';
    }
  station[10] = '\0';
  if ((sttn_id_in = fopen (sttn_id_file, "r")) == NULL)
    {
      printf ("Could not open file %s\n", sttn_id_file);
      exit (1);
    }
  while ((status = fgets (str, 256, sttn_id_in)) != NULL)
    {
      sscanf (str, "%d %s", &tmarker, tstation);

      if (tmarker == marker)
	{
	  strncpy (station, tstation, 4);
	  break;
	}
    }
  fclose (sttn_id_in);
  /*printf ("got sttn\n"); */

}
