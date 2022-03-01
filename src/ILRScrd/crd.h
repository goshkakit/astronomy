/*  Consolidated laser Ranging Data format (CRD)
 *       record and variable definitions for FORTRAN
 *       R. Ricklefs UT/CSR July 2007
 *  History:
 *  08/xx/07 - added H3 Target type (v0.27)
 *  11/26/07 - added H4 data quality alert
 *	       and #10 stop number
 *	       and #20 origin of values (v0.27) rlr.
 *  05/07/08 - Expand all configuration and data section character strings
 *             to allow up to 40 characters (plus NULL). 
 *             - Added detector channel to normalpoint (11) and calibration (40)
 *             records. 
 *             - Added field for 'crd' literal to H1.
 *             - Record '21' sky_clarity is not double rather than int.
 *               (v1.00 rlr)
 *
 * ----------------------------------------------------------------------
*/

/**
#define NSTATIONS 5
#define NSPACECRAFT 10
#define NRC10 10000
#define NRC11 1000
#define NRC12 10000
#define NRC20 10000
**/

/* Common name abbreviations:
 *  app = applied
 *  cdp = Crustal Dynamics Project (old NASA name)
 *  CofM = Center of Mass
 *  corr = correction or corrected
 *  est = estimated
 *  ind = indicator
 *  num = number
 *  osc = oscillator
 *  off = offset
 *  PmM = peak minus mean
 *  sic = (Goddard) Satellite ID Code
 *  stn = station
 *  SC = spacecraft
 *  sys = system
 *  utc = Universal Time Coordinated
 *  xcv = receive
 */

/* Ranging data header fields */
  /* H1 - format header */
  struct rh1 {
    char crd_literal[4];
    int format_version;
    int prod_year;
    int prod_mon;
    int prod_day;
    int prod_hour;
  };

  /* H2 - station header */
  struct rh2 {
    char stn_name[11];
    int cdp_pad_id;
    int cdp_sys_num;
    int cdp_occ_num;
    int stn_timescale;
  };

  /* H3 - spacecraft header */
  struct rh3 {
    char target_name[11];
    int ilrs_id;
    int sic;
    int norad;
    int SC_timescale;
    int target_type;
  };

  /* H4 - Session header */
  struct rh4 {
    int data_type;
    int start_year;
    int start_mon;
    int start_day;
    int start_hour;
    int start_min;
    int start_sec;
    int end_year;
    int end_mon;
    int end_day;
    int end_hour;
    int end_min;
    int end_sec;
    int data_release;
    int refraction_app_ind;
    int CofM_app_ind;
    int xcv_amp_app_ind;
    int stn_sysdelay_app_ind;
    int SC_sysdelay_app_ind;
    int range_type_ind;
    int data_qual_alert_ind;
  };

  /* Need indicators that these have been read? */
  /* H8 - End of Session footer */
  /* H9 - End of File footer */

/* Ranging data configuration fields (1 of n) */
  /* C0 - System Configuration Record */
  struct rc0 {
    int detail_type;
    double xmit_wavelength;
/**
    char sysconfig_id[5];
    char laserconfig_id[4];
    char detectorconfig_id[4];
    char timingconfig_id[4];
    char xponderconfig_id[4];
**/
    char config_ids[10][41];
  };

  /* C1 - Laser Configuration Record */
  struct rc1 {
    int detail_type;
    char laser_config_id[41];
    char laser_type[41];
    double prim_wavelength;	/* Primary wavelength of laser */
    double nom_fire_rate;	/* Nominal fire rate of laser */
    double pulse_energy;
    double pulse_width;
    double beam_div;
    int pulses_in_semitrain;	/* for multi-pulse systems */
  };

  /* C2 - Detector Configuration Record */
  struct rc2 {
    int detail_type;
    char detector_config_id[41];
    char detector_type[41];
    double app_wavelength;
    double qe;			/* quantum efficiency (in %) */
    double voltage;
    double dark_count;
    char output_pulse_type[41];
    double output_pulse_width;
    double spectral_filter;
    double spectral_filter_xmission;	/* % transmission of filter */
    double spatial_filter;
    char signal_proc[41];	/* signal processing algorithm or pgm name */
  };

  /* C3 - Timing Configuration Record */
  struct rc3 {
    int detail_type;
    char timing_config_id[41];
    char time_source[41];
    char freq_source[41];
    char timer[41];
    char timer_serial_num[41];
    double epoch_delay_corr;
  };

  /* C4 - Transponder Configuration Record */
  struct rc4 {
    int detail_type;
    char xponder_config_id[41];
    long double est_stn_utc_offset;
    double est_stn_osc_drift;
    long double est_xponder_utc_offset;
    double est_xponder_osc_drift;
    long double xponder_clock_ref_time;
    int stn_off_drift_app_ind;
    int SC_off_drift_app_ind;
    int SC_time_simplified_ind;
  };

/* Ranging data fields */
/* Secofday: need int sec and int psec? */
  /* 10 - Range Record */
  struct rd10 {
    long double sec_of_day;
    long double time_of_flight;
    char sysconfig_id[41];
    int epoch_event;
    int filter_flag;
    int detector_channel;
    int stop_number;
    int xcv_amp;
  };

  /* 11 - Normal Point Record */
  struct rd11 {
    long double sec_of_day;
    long double time_of_flight;
    char sysconfig_id[41];
    int epoch_event;
    double np_window_length;
    int num_ranges;
    double bin_rms;
    double bin_skew;
    double bin_kurtosis;
    double bin_PmM;
    double return_rate;
    int detector_channel;
  };

  /* 12 - Range Supplement Record */
  struct rd12 {
    long double sec_of_day;
    char sysconfig_id[41];
    double refraction_corr;
    double target_CofM_corr;
    double nd_value;
    double time_bias;
  };

  /* 20 - Meteorological Record */
  struct rd20 {
    long double sec_of_day;
    double pressure;
    double temperature;
    double humidity;
    int value_origin;
  };

  /* 21 - Meteorological Supplement Record */
  struct rd21 {
    long double sec_of_day;
    double wind_speed;
    double wind_direction;
    char precip_type[41];
    int visibility;
    double sky_clarity;
    int atmospheric_seeing;
    int cloud_cover;
  };

  /* 30 - Pointing Angles Record */
  struct rd30 {
    long double sec_of_day;
    double azimuth;
    double elevation;
    int direction_ind;
    int angle_origin_ind;
    int refraction_corr_ind;
  };

  /* 40 - Calibration Record */
  struct rd40 {
    long double sec_of_day;
    int type_of_data;
    char sysconfig_id[41];
    int num_points_recorded;
    int num_points_used;
    double one_way_target_dist;
    double cal_sys_delay;
    double cal_delay_shift;
    double cal_rms;
    double cal_skew;
    double cal_kurtosis;
    double cal_PmM;
    int cal_type_ind;
    int cal_shift_type_ind;
    int detector_channel;
  };

  /* 50 - Session Statistics Record */
  struct rd50 {
    char sysconfig_id[41];
    double sess_rms;
    double sess_skew;
    double sess_kurtosis;
    double sess_PmM;
    int data_qual_ind;
  };

  /* 60 - Compatibility Record */
  struct rd60 {
    char sysconfig_id[41];
    int sys_change_ind;
    int sys_config_ind;
  };

  /* 9X - User Defined Record */
  struct rd9x {
    /********** 
        Add userdefined record types and fields here 
    **********/
  };

  /* 00 - Comment Record */
  struct rd00 {
    char comment[81];
  };


  /* Ranging data header/footer records */
  /* H1 - format header */
  int read_h1(char * str, struct rh1 *header);

  /* H2 - station header */
  int read_h2(char * str, struct rh2 *header);

  /* H3 - spacecraft header */
  int read_h3(char * str, struct rh3 *header);

  /* H4 - Session header */
  int read_h4(char * str, struct rh4 *header);

  /* Need indicators that these have been read? */
  /* H8 - End of Session footer */
  void read_h8(char * str);

  /* H9 - End of File footer */
  void read_h9(char * str);

  /* Ranging data configuration records (1 of n) */
  /* C0 - System Configuration Record */
  int read_c0(char * str, struct rc0 *config);

  /* C1 - Laser Configuration Record */
  int read_c1(char * str, struct rc1 *config);

  /* C2 - Detector Configuration Record */
  int read_c2(char * str, struct rc2 *config);

  /* C3 - Timing Configuration Record */
  int read_c3(char * str, struct rc3 *config);

  /* C4 - Transponder Configuration Record */
  int read_c4(char * str, struct rc4 *config);

  /* Ranging data records */
  /* Secofday: need int sec and int psec? */
  /* 10 - Range Record */
  int read_10(char * str, struct rd10 *data_recd);

  /* 11 - Normal Point Record */
  int read_11(char * str, struct rd11 *data_recd);

  /* 12 - Range Supplement Record */
  int read_12(char * str, struct rd12 *data_recd);

  /* 20 - Meteorological Record */
  int read_20(char * str, struct rd20 *data_recd);

  /* 21 - Meteorological Supplement Record */
  int read_21(char * str, struct rd21 *data_recd);

  /* 30 - Pointing Angles Record */
  int read_30(char * str, struct rd30 *data_recd);

  /* 40 - Calibration Record */
  int read_40(char * str, struct rd40 *data_recd);

  /* 50 - Session Statistics Record */
  int read_50(char * str, struct rd50 *data_recd);

  /* 60 - Compatibility Record */
  int read_60(char * str, struct rd60 *data_recd);

  /* 9X - User Defined Records 90-99 */
  int read_9x(char * str, struct rd9x *data_recd);

  /* 00 - Comment Record */
  int read_00(char * str, struct rd00 *data_recd);

//##################################################
  /* Ranging data header/footer records */
  /* H1 - format header */
  void write_h1(FILE * str_out, struct rh1 header);

  /* H2 - station header */
  void write_h2(FILE * str_out, struct rh2 header);

  /* H3 - spacecraft header */
  void write_h3(FILE * str_out, struct rh3 header);

  /* H4 - Session header */
  void write_h4(FILE * str_out, struct rh4 header);

  /* Need indicators that these have been read? */
  /* H8 - End of Session footer */
  void write_h8(FILE * str_out);

  /* H9 - End of File footer */
  void write_h9(FILE * str_out);

  /* Ranging data configuration records (1 of n) */
  /* C0 - System Configuration Record */
  void write_c0(FILE * str_out, struct rc0 config);

  /* C1 - Laser Configuration Record */
  void write_c1(FILE * str_out, struct rc1 config);

  /* C2 - Detector Configuration Record */
  void write_c2(FILE * str_out, struct rc2 config);

  /* C3 - Timing Configuration Record */
  void write_c3(FILE * str_out, struct rc3 config);

  /* C4 - Transponder Configuration Record */
  void write_c4(FILE * str_out, struct rc4 config);

  /* Ranging data records */
  /* 10 - Range Record */
  void write_10(FILE * str_out, struct rd10 data_recd);

  /* 11 - Normal Point Record */
  void write_11(FILE * str_out, struct rd11 data_recd);

  /* 12 - Range Supplement Record */
  void write_12(FILE * str_out, struct rd12 data_recd);

  /* 20 - Meteorological Record */
  void write_20(FILE * str_out, struct rd20 data_recd);

  /* 21 - Meteorological Supplement Record */
  void write_21(FILE * str_out, struct rd21 data_recd);

  /* 30 - Pointing Angles Record */
  void write_30(FILE * str_out, struct rd30 data_recd);

  /* 40 - Calibration Record */
  void write_40(FILE * str_out, struct rd40 data_recd);

  /* 50 - Session Statistics Record */
  void write_50(FILE * str_out, struct rd50 data_recd);

  /* 60 - Compatibility Record */
  void write_60(FILE * str_out, struct rd60 data_recd);

  /* 9X - User Defined Record */
  void write_9x(FILE * str_out, struct rd9x data_recd);

  /* 00 - Comment Record */
  void write_00(FILE * str_out, struct rd00 data_recd);
