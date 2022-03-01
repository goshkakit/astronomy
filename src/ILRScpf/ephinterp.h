#include "ephi.h"

void cinterp_pvac (double jdint, double seceph, double outvec[6], double backvec[6],
	      double aberout[3], double aberback[3],
	      double *corr12, double *corr23, int fvel, int *ircode, eh &ephihdr, ei &ephi);

void cinterp_off (double jdint, double seceph, double offsetout[6], double offsetback[6], int *ircode, eh &ephihdr, ei &ephi);

void cinterp_rot (double jdint, double seceph,  int *arottype, double *agast, double arotang[3], int *ircode, eh &ephihdr, ei &ephi);

// error
void cinterp_eop (double jdint, double seceph, double x, double y, double dutc, int *ircode, eh &ephihdr, ei &ephi);