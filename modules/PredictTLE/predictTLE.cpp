#include <stdio.h>


#include <fstream>
#include <iostream>
#include "json/json.h"
#include <vector>

#include "Norad\coreLib.h"
#include "Norad\cOrbit.h"

#include "common\DataConverter.h"
#include "common\mytypes.h"
#include "common\TLELoader.h"

#include "InfluenceForce\InfluenceForce.h"

const double RG = 180.0 / M_PI;	// Radian-degree convertion
const double GR = 1.0 / RG;


struct Ttelescop
{
	int id;						// Observer number
	double lon, lat, height;	// Telescop  position
	double rate;				// Frame rate
	double FOVwa;				// FOV 05 width on ASC
	double FOVwb;				// FOV 05 width on Dec
	double sig;					// FOR TRACK SIMULATION ONLY :sigma of Ra,Dec measurements errors [ang.sec]
	double max_dur;				// FOR TRACK SIMULATION ONLY :maximum track duration [sec]
	std::string name;
};

struct ViewParams
{
	double az, el;
	double range;
};

Json::Value read_Json(std::string &FileName)
{
	Json::CharReaderBuilder builder;
	Json::Value value;
	JSONCPP_STRING errs;
	builder["collectComments"] = false;

	std::ifstream doc(FileName, std::ifstream::binary);
	bool ok = parseFromStream(builder, doc, &value, &errs);
	if (!ok)
	{
		std::cout << "Failed to parse configuration\n" << errs;
	}

	return value;
}

void read_tel_Json(Ttelescop &tel, Json::Value Jtel, std::string &FileName)
{
	tel.lon = Jtel["site_lon"].asDouble();
	tel.lat = Jtel["site_lat"].asDouble();
	tel.height = Jtel["site_height"].asDouble();
	
	//WGS84_XYZ(height, lat, lon, tel.site.x, tel.site.y, tel.site.z);

	if (Jtel.isMember("id"))	tel.id = Jtel["id"].asInt();
	else tel.id = -1;
	//check_site(lon, lat, height, tel.id);

	if (Jtel.isMember("Frame_rate"))	tel.rate = Jtel["Frame_rate"].asDouble();
	else tel.rate = 5.0;

	if (Jtel.isMember("RMS"))	tel.sig = Jtel["RMS"].asDouble();
	else tel.sig = 0.001389;
	tel.sig = tel.sig * GR;

	if (Jtel.isMember("Max_dur"))	tel.max_dur = Jtel["Max_dur"].asDouble();
	else tel.max_dur = 1e10;

	if (Jtel.isMember("name"))	tel.name = Jtel["name"].asString();
	else tel.name = "Unknown";
}

void read_tels_Json(vector<Ttelescop> &tels, Json::Value Jtels, std::string &FileName)
{
	unsigned i;
	tels.resize(Jtels.size());
	for (i = 0; i < Jtels.size(); ++i) {
		read_tel_Json(tels[i], Jtels[i], FileName);
	}
}

Force::InfluenceForce *IForce = new Force::InfluenceForce();

ViewParams ViewTEMESatFromTelescope( double JD_UTC, double* tel_pos, double *sat_pos )
{
	DataConverter Dconv;

	// MDB
	double jdm = JD_UTC + 0.125;
	double date1 = Dconv.JDtoYYYYMMDD(jdm);
	double time1 = Dconv.SECtoHHMMSS(date1, jdm);

	// установка времени
	double int1, ajd1, delt1;
	IForce->set_time(date1, time1, &ajd1, &delt1, &int1);

	// матрица перевода в земную систему
	double Arot[9];
	IForce->iers_update_matrix(int1, Arot, ajd1, delt1);

	// матрица перехода из земной в нормальную систему
	//double invArot[9];
	//IForce->transpose(Arot, invArot);

	// матрица перехода из ICRF в TEME
	double A_Teme[9];
	IForce->GetTemeMatrix(int1, A_Teme, ajd1, delt1);

	// обратная матрица перехода из TEME в ICRF
	double invA_Teme[9];
	IForce->transpose(A_Teme, invA_Teme);


	// перевод вектора состояния в ICRF 
	double PICRF[3];
	IForce->matVecMul(invA_Teme, sat_pos, PICRF);
	// переход в земную
	double PITRF[3];
	IForce->matVecMul(Arot, PICRF, PITRF);

	// перевод координат телескопа в ICRF
	//double TICRF[3];
	//IForce->matVecMul(invArot, tel_pos, TICRF);

	// вектор телескопа
	S3DCoordinate Rtel;
	Rtel.x = tel_pos[0];
	Rtel.y = tel_pos[1];
	Rtel.z = tel_pos[2];

	// вектор наблюдения
	S3DCoordinate Rn;
	Rn.x = PITRF[0];
	Rn.y = PITRF[1];
	Rn.z = PITRF[2];
	S3DCoordinate Rmeg = Rn - Rtel;

	double cosa = (Rtel*Rmeg) / (Rtel.norm()*Rmeg.norm());
	double zenit = acos(cosa);

	ViewParams vp;

	vp.az = 0;
	vp.el = 90.0 - zenit*RG;
	vp.range = Rmeg.norm();

	return vp;
}

void main()
{
	string DBRoot = "C:\\Users\\Nikolay\\AppData\\Local\\ODSW";
	std::string pinpf_tels = DBRoot + "\\DB\\input\\Telescopes\\tels_simulation.json";
	std::string tle_inpu = DBRoot + "\\DB\\input\\Test\\TLE\\orbit.tle";

	double SAT_elv = 10.0;

	Json::Value inpf, JTels;
	std::vector<Ttelescop> tels;

	inpf = read_Json(pinpf_tels);
	if (inpf.isMember("Telescopes")) {
		JTels = inpf["Telescopes"];
	}
	else std::cout << "Missing parameter 'Telescopes' in " << pinpf_tels << "\n";

	read_tels_Json(tels, JTels, pinpf_tels);
	if (tels.size() == 0) {
		std::cout << "Empty array of telescopes in " << pinpf_tels << "\n";
	};

	std::cout << "Count of telescopes in " << tels.size() << "\n";

	// load TLE
	TLELoader tleLoader;
	tleLoader.LoadData(tle_inpu.c_str(), 3);

	// init rotation
	IForce->Init_CPU();

	DataConverter DC;

	for (int tl = 0; tl < tels.size(); tl++)
	{
		//    Latitude  in degrees (negative south)
		//    Longitude in degrees (negative west)
		//    Altitude  in km
		cSite siteView(tels[tl].lat, tels[tl].lon, tels[tl].height );

		double tel_pos[3];
		double tel_pos_wgs72[3];

		double lon = tels[tl].lon*GR;
		double lat = tels[tl].lat*GR;
		double alt = tels[tl].height;

		double h = 6378.135 + tels[tl].height;

		double theta = lon;
		double c = 1.0 / sqrt(1.0 + F * (F - 2.0) * sqr(sin(lat)));
		double s = sqr(1.0 - F) * c;
		double achcp = (XKMPER_WGS72 * c + alt) * cos(lat);

		tel_pos_wgs72[0] = achcp * cos(theta);         // km
		tel_pos_wgs72[1] = achcp * sin(theta);         // km
		tel_pos_wgs72[2] = (XKMPER_WGS72 * s + alt) * sin(lat);   // km

		tel_pos[0] = h*cos(lat)*cos(lon);
		tel_pos[1] = h*cos(lat)*sin(lon);
		tel_pos[2] = h*sin(lat);

		for (int ti = 0; ti < 3; ti++)
		{
			printf("telesope D72 %f m \n", (tel_pos_wgs72[ti] - tel_pos[ti])*1000.0);
		}

		for (int i = 0; i < tleLoader.NORADList.size(); i++)
		{
			cOrbit *orbit = tleLoader.NORADList[i];
			double JDtle = orbit->Epoch().Date();
			double jds = (int)JDtle;
			double jd_delta = 2.0;
			double jde = (int)JDtle+ jd_delta;

			double dayS = DC.JDtoYYYYMMDD(jds);
			double dayE = DC.JDtoYYYYMMDD(jde);
			printf( "JD start: %f ,  Day start: %f\n", jds, dayS );

			bool SatVisible = false;
			double stepSec = 0.1;
			double stepMin = stepSec / 60.0;

			double mpe_start = (jds - JDtle)*24.0*60.0;
			double mpe_stop = (jde - JDtle)*24.0*60.0;
			double JDstart = 0;

			// Calculate position, velocity
			// mpe = "minutes past epoch
			for (double mpe = mpe_start; mpe <= mpe_stop; mpe += stepMin)
			{
				double currenJd = JDtle + mpe / 24.0 / 60.0;

				// Get the position of the satellite at time "mpe"
				cEciTime eci = orbit->GetPosition(mpe);

				// Now get the "look angle" from the site to the satellite. 
				// Note that the ECI object "eciSDP4" contains a time associated
				// with the coordinates it contains; this is the time at which
				// the look angle is valid.
				cTopo topoLook = siteView.GetLookAngle(eci);

				double sat_pos[3];
				sat_pos[0] = eci.Position().m_x;
				sat_pos[1] = eci.Position().m_y;
				sat_pos[2] = eci.Position().m_z;
				ViewParams vp = ViewTEMESatFromTelescope(currenJd, tel_pos_wgs72, sat_pos);
				

				double current_el = topoLook.ElevationDeg();
				//current_el = vp.el;

				//printf("el %f R %f\n", vp.el, vp.range);

				if (current_el >= SAT_elv && !SatVisible)
				{
					
					SatVisible = true;
					printf("EL: %3.5f [%3.5f] R: %5.5f [%5.5f] ", vp.el, topoLook.ElevationDeg(), vp.range, vp.range - topoLook.RangeKm());
					printf("AZ: %3.3f S: %.7f ", topoLook.AzimuthDeg(), currenJd);
					JDstart = currenJd;
				}
				if (current_el < SAT_elv && SatVisible)
				{
					SatVisible = false;
					double duration = (currenJd - JDstart)*86400.0;
					printf("EL: %3.5f [%3.5f] R: %5.5f [%5.5f] ", vp.el, topoLook.ElevationDeg(), vp.range, vp.range - topoLook.RangeKm());
					printf("AZ: %3.3f E: %.7f D: %.5f\n", topoLook.AzimuthDeg(), currenJd, duration );

				}
			}
		}
	}
}