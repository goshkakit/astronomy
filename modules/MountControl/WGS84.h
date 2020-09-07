#pragma once

//WGS84 constants
const double pi = 3.14159265358979;										// 32385;
const double RG = 180.0 / pi;											// Radian-degree convertion
const double GR = 1.0 / RG;

const double WEarth = 7.2921158553e-5;									// Earth angular speed [rad/s]
const double RMEarth = 6371.0088;										// Earth mean radius [km]
const double R0Earth = 6378.137;										// Earth equatorial radius [km]
const double R0Earth_m = R0Earth * 1e3;									// Earth equatorial radius [m]
const double RPEarth = 6356.752314245;									// Earth polar radius [km]
const double FlEarth = 1.0 / 298.257223563;								//Earth flattening
const double E2Earth = 1.0 - (RPEarth / R0Earth)*(RPEarth / R0Earth);	// Square eccentricity of Earth
const double GMEarth = 398600.44189;									// Earth gravitational parameter [km^3 / s^2]
const double GMEarth_m = GMEarth * 1e9;									// Earth gravitational parameter [m^3 / s^2]
const double gEarth = 0.00980665;										// Earth gravitational acceleration [km / s^2]
const double gEarth_m = 9.80665;										// Earth gravitational acceleration [m / s^2]
const double gEarth_sm = gEarth * 1e5;									// Earth gravitational acceleration [sm / s^2]

//elev, lat, lon -> xyzITRF
void WGS84_XYZ(double Hw, double Fwg, double Lwg, double & X, double & Y, double & Z);

//xyzITRF -> elev, lat, lon
void XYZ_WGS84(double X, double Y, double Z, double & Hw, double & Fwg, double & Lwg);
void RF_WGS84(double R, double F, double & Hz, double & Fzg);
void SinCos(double a, long double & sina, long double & cosa);
double sign(double Val);
double Sqr(double a);