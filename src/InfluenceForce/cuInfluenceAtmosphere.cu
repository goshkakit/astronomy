//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Earth Atmosphere
// GPU version
//==============================================================================//
#include <math.h>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//==============================================================================//
// некоторын дополнительные функции
//==============================================================================//
__device__ int kernalAnint( double X )
{
	double res;
	int iX = (int)X;
	double delta = X - (double)iX;

	if( delta > 0.5 )
		res = iX + 1;
	else
		res = iX;

	return iX;
};
__device__ double kernalAmax( double X, double Y )
{
	if( X > Y )
		return X;
	else
		return Y;
}
__device__ double kernalAmin( double X, double Y )
{
	if( X < Y )
		return X;
	else
		return Y;
}
__device__ double kernalAmod( double X, double Y )
{
	int s = (int)(X/Y);
	double res = X - ((double)s)*Y;

	return res;
}
//==============================================================================//
// Функция расчет плотности верхних слоев атмосферы согласно ГОСТу.
//==============================================================================//
__device__  double kernalRoa2004_2( double time, double *x, double ajd0, double delt0 )
{
	double f0t[7] = {  75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 200.0 };
	double Re = 6378.136;
	double om = 6.300388008;
	double alpha = 0.0033528037;
	double b1900 = 2415020.31352;
	double  btau = 365.2422;
	double  ro0 = 1.58868E-8;
	double   an0 =  .20580000E+01;
	double   an1 =  .58870000E-02;
	double  an2 = -.40120000E-05;

	double saem[4];
	saem[0] = 100.0;
	saem[1] = 100.0;
	saem[2] = 100.0;
	saem[3] = 3.0;

	double Rmod, sinF, h, F107, F81, aKp, F0;
	double e0, e1, e2, e3, e4, e5, e6, e7, e8;
	double a0, a1, a2, a3, a4, a5, a6;
	double al0, al1, al2, al3, al4;
	double c0, c1, c2, c3, c4;
	double d0, d1, d2, d3, d4;
	double b0, b1, b2, b3, b4;
	double fi1, power, aKs0, aKs1, aKs2, aKs3, aKs4, aKss1, aKss4, ron;
	double ct, sz, xss, yss, zss, rr;
	double cosFi1, sinFi1, app, s_app, c_app, xs, ys, zs, rs, cosFi, cos05, ak0, ak1, d, ad, ak2, ak3, ak4;
	double h0, akst1, akst2, ast0, t;
	int iq, iq6 ,ish_a, ish_b, ish_c, ish_d, ish_e, ish_l;

	t = ajd0+(delt0+time)/86.4;
	Rmod = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
	sinF = x[2]/Rmod;
	h = kernalAmax( Rmod*1.0E+3-Re*(1.0-alpha*sinF*sinF), 0.0);
		
	double roa2004 = 0;

	if (h >= 1500.0) 
	{
		roa2004 = 0.0;
		return roa2004;
	}
	else if (h > 120.0) 
	{
		F107 = saem[1];
		F81  = saem[2];
		aKp  = saem[3];
		iq = kernalAmin( kernalAmax( kernalAnint((F81-75.0)/25.0), 0), 7) +1;
			
		F0 = f0t[iq-1]; //!!!

		if( F81 > 200.0)
		{
			iq = 6;
		}
		else if (F81 > 225.0) 
		{
			F0 = 250.0;
			iq = 7;
		}

		iq6 = (iq-1)*6;

		iq6 = iq6 - 1; //!! коррекция индекса

		ish_a=0;
		ish_b=0;
		ish_c=0;
		ish_d=0;
		ish_e=0;
		ish_l=0;
		if (h > CUtick[iq6+1]) ish_a=7;
		if (h > CUtick[iq6+2]) ish_b=7;
		if (h > CUtick[iq6+3]) ish_c=7;
		if (h > CUtick[iq6+4]) ish_d=7;
		if (h > CUtick[iq6+5]) ish_e=7;
		if (h > CUtick[iq6+6]) ish_l=7;

		// ------------------
		iq = iq - 1; //!! коррекция индекса
		// a - coeffs
		a0= CUaa0[ish_a+iq];
		a1= CUaa1[ish_a+iq];
		a2= CUaa2[ish_a+iq];
		a3= CUaa3[ish_a+iq];
		a4= CUaa4[ish_a+iq];
		a5= CUaa5[ish_a+iq];
		a6= CUaa6[ish_a+iq];
		// al - coeffs
		al0= CUaal0[ish_l+iq];
		al1= CUaal1[ish_l+iq];
		al2= CUaal2[ish_l+iq];
		al3= CUaal3[ish_l+iq];
		al4= CUaal4[ish_l+iq];
		// c - coeffs
		c0= CUcc0[ish_c+iq];
		c1= CUcc1[ish_c+iq];
		c2= CUcc2[ish_c+iq];
		c3= CUcc3[ish_c+iq];
		c4= CUcc4[ish_c+iq];
		// d - coeffs
		d0= CUdd0[ish_d+iq];
		d1= CUdd1[ish_d+iq];
		d2= CUdd2[ish_d+iq];
		d3= CUdd3[ish_d+iq];
		d4= CUdd4[ish_d+iq];
		// b - coeffs
		b0= CUbb0[ish_b+iq];
		b1= CUbb1[ish_b+iq];
		b2= CUbb2[ish_b+iq];
		b3= CUbb3[ish_b+iq];
		b4= CUbb4[ish_b+iq];
		// e - coeffs
		e0= CUee0[ish_e+iq];
		e1= CUee1[ish_e+iq];
		e2= CUee2[ish_e+iq];
		e3= CUee3[ish_e+iq];
		e4= CUee4[ish_e+iq];
		e5= CUee5[iq];
		e6= CUee6[iq];
		e7= CUee7[iq];
		e8= CUee8[iq];
		// ------------------

		fi1 = CUff1[iq];

		power = a0+h*(a1+h*(a2+h*(a3+h*(a4+h*(a5+h*a6)))));
		aKs0 = al0+h*(al1+h*(al2+h*(al3+h*al4)));
		aKs1 = c0+h*(c1+h*(c2+h*(c3+h*c4)));
		aKss1 = an0+h*(an1+h*an2);
		aKs2 = d0+h*(d1+h*(d2+h*(d3+h*d4)));
		aKs3 = b0+h*(b1+h*(b2+h*(b3+h*b4)));
		aKs4 = e0+h*(e1+h*(e2+h*(e3+h*e4)));
		aKss4 = e5+aKp*(e6+aKp*(e7+aKp*e8));
		ron = ro0*exp(power);

		ct  = (t-2415020.0)/36525.0;
		sz  = 628.33195099*ct+1.739935890;
		xss = -cos(sz);
		sinF = -sin(sz);
		yss = 0.91747*sinF;
		zss = 0.397805*sinF;
		rr = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );

		cosFi1 = cos(fi1);
		sinFi1 = sin(fi1);
		app = sz+( kernalAmod(t, 1.0)-0.5)*om;
		s_app = sin(app);
		c_app = cos(app);
		xs = xss*c_app + yss*s_app;
		ys = yss*c_app - xss*s_app;
		zs = zss;
		rs = sqrt(xs*xs+ys*ys+zs*zs);
		cosFi = (x[0]*(xs*cosFi1-ys*sinFi1)+ x[1]*(ys*cosFi1+xs*sinFi1)+ x[2]*zs)/(rs*Rmod);
		cos05 = sqrt((1.0+cosFi)/2.0);
		ak0 = 1.0 + aKs0*(F81-F0)/F0;

		//aK1 = aKs1*(cos05**aKss1); //!!!
		ak1 = aKs1*pow(cos05, aKss1);

		d = kernalAmod(t-b1900, btau);

		ad = -2.53418E-2 + d*(-2.44075E-3+d*(3.08389E-6 +d*( 
				2.90115E-6 + d*(-4.99606E-8+d*(3.36327E-10+d*( 
			-1.0966E-12 + d*(1.73227E-15+d*(-1.06271E-18))))))));
		ak2 = ad*aKs2;
		ak3 = aKs3*(F107-F81)/(F81+abs(F107-F81));
		ak4 = aKs4*aKss4;
		roa2004 = ron*ak0*(1.0 + ak1 + ak2 + ak3 + ak4);
		return roa2004;
	}
	else if ((h > 0.0) && (h < 20.0)) 
	{
		h0=0.0;
		ast0=1.228;
		akst1=-9.0764E-2;
		akst2=-2.0452E-3;
	}
	else if ((h >= 20.0)&&(h < 60.0))
	{
		h0=20.0;
		ast0=9.013E-2;
		akst1=-0.16739;
		akst2=6.2669E-4;
	}
	else if ((h >= 60.0)&&(h < 100.0)) 
	{
		h0=60.0;
		ast0=3.104E-4;
		akst1=-0.137;
		akst2=-7.8653E-4;
	}
	else if ((h >= 100.0)&&(h < 120.0)) 
	{
		h0=100.0;
		ast0=3.66E-7;
		akst1=-0.18553;
		akst2=1.5397E-3;
	}
	power = (h-h0)*(akst1+akst2*(h-h0));
	roa2004 = ast0*exp(power);
	return roa2004;
}

//==============================================================================//
// Процедура расчет возмущения, вызванного сопротивлением атмосферы.
//==============================================================================//
__device__ void kernalAtm_drag( double *x, double t, double *f, double sigma_up, double ajd0, double delt0 )
{
	double v;
	v = sqrt( x[3]*x[3] + x[4]*x[4] + x[5]*x[5] );

	double rc = kernalRoa2004_2( t, x, ajd0, delt0 );
	double coeff =  rc*sigma_up*v*1.0E+6;

	f[0] = -x[3]*coeff;
	f[1] = -x[4]*coeff;
	f[2] = -x[5]*coeff;
}
//==============================================================================//
