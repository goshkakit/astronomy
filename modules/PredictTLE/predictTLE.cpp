#include <stdio.h>


#include <fstream>
#include <iostream>
#include <vector>

#include "Norad\coreLib.h"
#include "Norad\cOrbit.h"

#include "common\DataConverter.h"
#include "common\mytypes.h"
#include "common\TLELoader.h"

#include "InfluenceForce\InfluenceForce.h"

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "MountController.h"

using namespace VecMath;

const double RG = 180.0 / M_PI;	// Radian-degree convertion
const double GR = 1.0 / RG;

int ETALON_SAT[104] = {16908,
					36287,
					37384,
					37948,
					41434,
					40938,
					38250,
					40749,
					40748,
					36508,
					19751,
					20026,
					37846,
					37847,
					38857,
					38858,
					40128,
					40129,
					40544,
					40545,
					40889,
					40890,
					41859,
					41175,
					41174,
					41550,
					41549,
					41860,
					41861,
					41862,
					41579,
					37138,
					37372,
					37867,
					37868,
					39155,
					40001,
					40315,
					41554,
					37781,
					39635,
					40269,
					40547,
					41241,
					41384,
					33105,
					41240,
					39227,
					8820,
					22195,
					38077,
					27944,
					37158,
					42738,
					42917,
					42965,
					37755,
					39086,
					41335,
					7646,
					22824,
					39068,
					39452,
					39451,
					39453,
					42829,
					36605,
					31698,
					34661,
					32711,
					32384,
					32260,
					29601,
					29486,
					28874,
					28474,
					28361,
					28190,
					28129,
					27704,
					27663,
					26690,
					26605,
					26407,
					26360,
					25933,
					25030,
					24876,
					24320,
					23953,
					23833,
					23027,
					22877,
					22779,
					22700,
					22657,
					22581,
					22446,
					22275,
					22231,
					22108,
					22014,
					21930,
					21890 };

int GPS_SAT[36] = { 34661,
						32711,
						32384,
						32260,
						29601,
						29486,
						28874,
						28474,
						28361,
						28190,
						28129,
						27704,
						27663,
						26690,
						26605,
						26407,
						26360,
						25933,
						25030,
						24876,
						24320,
						23953,
						23833,
						23027,
						22877,
						22779,
						22700,
						22657,
						22581,
						22446,
						22275,
						22231,
						22108,
						22014,
						21930,
						21890 };




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

bool read_Json(std::string &FileName, rapidjson::Document &document)
{
    std::ifstream stream(FileName, std::ifstream::binary);

    rapidjson::IStreamWrapper wrapper(stream);

    document.ParseStream(wrapper);

    if (document.HasParseError())
    {
        size_t offset = document.GetErrorOffset();
        rapidjson::ParseErrorCode code = document.GetParseError();

        std::cout << "Failed to parse JSON document: code = " << code << ", offset = " << offset << ", path = " << FileName << std::endl;

        return false;
    }

    return true;
}

void read_tel_Json(Ttelescop &tel, const rapidjson::Value &value)
{
    if (value.IsObject()) 
    {
        if (value.HasMember("site_lon"))
        {
            const rapidjson::Value &v_site_lon = value["site_lon"];

            if (v_site_lon.IsNumber())
            {
                tel.lon = v_site_lon.GetDouble();
            }
        }

        if (value.HasMember("site_lat"))
        {
            const rapidjson::Value &v_site_lon = value["site_lat"];

            if (v_site_lon.IsNumber())
            {
                tel.lat = v_site_lon.GetDouble();
            }
        }

        if (value.HasMember("site_height"))
        {
            const rapidjson::Value &v_site_lon = value["site_height"];

            if (v_site_lon.IsNumber())
            {
                tel.height = v_site_lon.GetDouble();
            }
        }

        //WGS84_XYZ(height, lat, lon, tel.site.x, tel.site.y, tel.site.z);

        tel.id = -1;

        if (value.HasMember("id"))
        {
            const rapidjson::Value &v_site_lon = value["id"];

            if (v_site_lon.IsNumber())
            {
                tel.id = v_site_lon.GetInt();
            }
        }

        //check_site(lon, lat, height, tel.id);

        tel.rate = 5.0;

        if (value.HasMember("Frame_rate"))
        {
            const rapidjson::Value &v_site_lon = value["Frame_rate"];

            if (v_site_lon.IsNumber())
            {
                tel.rate = v_site_lon.GetDouble();
            }
        }

        tel.sig = 0.001389;

        if (value.HasMember("RMS"))
        {
            const rapidjson::Value &v_site_lon = value["RMS"];

            if (v_site_lon.IsNumber())
            {
                tel.sig = v_site_lon.GetDouble();
            }
        }

        tel.sig = tel.sig * GR;

        tel.max_dur = 1e10;

        if (value.HasMember("Max_dur"))
        {
            const rapidjson::Value &v_site_lon = value["Max_dur"];

            if (v_site_lon.IsNumber())
            {
                tel.max_dur = v_site_lon.GetDouble();
            }
        }

        tel.name = "Unknown";

        if (value.HasMember("name"))
        {
            const rapidjson::Value &v_site_lon = value["name"];

            if (v_site_lon.IsNumber())
            {
                tel.name = v_site_lon.GetString();
            }
        }
    }
}

void read_tels_Json(vector<Ttelescop> &tels, const rapidjson::Value &Jtels)
{
    if (Jtels.IsArray())
    {
        tels.resize(Jtels.Size());

        for (unsigned idx = 0, end = Jtels.Size(); idx != end; ++idx)
        {
            read_tel_Json(tels[idx], Jtels[idx]);
        }
    }
}

Force::InfluenceForce *IForce = new Force::InfluenceForce();
std::shared_ptr<NewConvertor> convertor;

void ConvertXYZtoRADEC(double* resultPosition, double* inTelescopePosition, double* Ra, double* Dec)
{
	// входные данные задаются в ICRF
	// вектор направления в системе ICRF
	double x = resultPosition[0] - inTelescopePosition[0];
	double y = resultPosition[1] - inTelescopePosition[1];
	double z = resultPosition[2] - inTelescopePosition[2];

	double r = atan2(y, x);
	double d = atan2(z, sqrt(x * x + y * y));

	double pi = 3.1415926535;

	if (r < 0)
		r = 2.0 * pi + r;

	*Ra = r;
	*Dec = d;
}

void ITRFToICRF(double JD_UTC, double* posITRF, double* posICRF)
{
	DataConverter Dconv;
	// MDB
	JD_UTC = JD_UTC + 0.125;
	double dataMDB = Dconv.JDtoYYYYMMDD(JD_UTC);
	double timeMDB = Dconv.SECtoHHMMSS(dataMDB, JD_UTC);

	// установка времени
	double int1, ajd1, delt1;
	IForce->set_time(dataMDB, timeMDB, &ajd1, &delt1, &int1);

	// матрица перевода в земную систему
	double Arot[9];
	IForce->iers_update_matrix(int1, Arot, ajd1, delt1);
	// матрица перехода из земной в нормальную систему
	double invArot[9];
	IForce->transpose(Arot, invArot);
	IForce->matVecMul(invArot, posITRF, posICRF);
}

// получение положения телескопа в момент измерения в системе ICRF
void ConvertTEMEtoICRF(double JD_UTC, double* inPTEME, double* outPICRF)
{
	// установка времени
	DataConverter Dconv;
	// MDB
	JD_UTC = JD_UTC + 0.125;
	double dataMDB = Dconv.JDtoYYYYMMDD(JD_UTC);
	double timeMDB = Dconv.SECtoHHMMSS(dataMDB, JD_UTC);

	double int1, ajd1, delt1;
	IForce->set_time(dataMDB, timeMDB, &ajd1, &delt1, &int1);

	// матрица перехода из ICRF в TEME
	double A_Teme[9];
	IForce->GetTemeMatrix(int1, A_Teme, ajd1, delt1);
	// обратная матрица перехода из TEME в ICRF
	double invA_Teme[9];
	IForce->transpose(A_Teme, invA_Teme);

	// перевод вектора состояния в ICRF 
	IForce->matVecMul(invA_Teme, inPTEME, outPICRF);
}

void getAzElSatellite()
{

}

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
	convertor = std::make_shared<NewConvertor>();

	double tel_lon = 37.4340;
	double tel_lat = 55.8581;
	double tel_h = 0.225;
	double RMEarth = 6371.0088;

	double SAT_MinEl = 30.0;
	double SAT_MaxEl = 70.0;

	double x, y, z;
	convertor->WGS84_XYZ(tel_h, tel_lat, tel_lon, x, y, z);

	// telescope xyz
	double Tel_itrf[3];
	Tel_itrf[0] = x;
	Tel_itrf[1] = y;
	Tel_itrf[2] = z;

	// norad telescope position
	cSite siteView(tel_lat, tel_lon, tel_h);

	double xt, yt, zt;
	xt = (RMEarth + tel_h) * cos(tel_lat* GR) * cos(tel_lon * GR);
	yt = (RMEarth + tel_h) * cos(tel_lat * GR) * sin(tel_lon * GR);
	zt = (RMEarth + tel_h) * sin(tel_lat * GR);

	printf("TEL WGS84: %f %f %f\n", x, y, z);
	printf("TEL Spherical: %f %f %f\n", xt, yt, zt);

	// TLE
	std::string tle_input = "tle/TLE20201003.txt";
	// load TLE
	TLELoader tleLoader;
	//tleLoader.LoadData(tle_input.c_str(), 2, GPS_SAT, 36);
	tleLoader.LoadData(tle_input.c_str(), 2, ETALON_SAT, 104);
	printf("TLE count: %d\n", tleLoader.NORADList.size());

	// init rotation
	IForce->Init_CPU();

	// date time convertor
	DataConverter DC;

	double DateView = 20201008;
	//double TimeView1 = 160000;
	//double TimeView2 = 173000;
	double TimeView1 = 143000;
	double TimeView2 = 160000;
	double D = DC.YYYYMMDDtoJD(DateView) - 0.5;	// юлианская дата начала суток
	double T1 = DC.HHMMSSToSec(TimeView1);
	double T2 = DC.HHMMSSToSec(TimeView2);
	double JD_StartView = D + T1 / 86400.0; // полная юлианская дата по UTC
	double JD_EndView = D + T2 / 86400.0;
	double DT = (JD_EndView - JD_StartView) * 86400.0;
	double stepSec = 120;
	int maxCountPoint = 8;

	int Npos = DT / stepSec;
	std::vector< S3DCoordinate > Plan;
	std::vector< int > PlanIds;
	for (int i = 0; i < Npos; i++)
	{
		Plan.push_back(S3DCoordinate(-1, 0, 0));
		PlanIds.push_back(0);
	}
	
	for (int i = 0; i < tleLoader.NORADList.size(); i++)
	{
		cOrbit* orbit = tleLoader.NORADList[i];
		double JDtle = orbit->Epoch().Date();

		if (orbit->Apogee() > 20000)
		{
			continue;
		}
		
		//double jds = (int)JDtle;
		//double jd_delta = 2.0;
		//double jde = (int)JDtle + jd_delta;
		double jds = JD_StartView;
		double jde = JD_EndView;

		//double dayS = DC.JDtoYYYYMMDD(jds + 0.125);
		//double dayE = DC.JDtoYYYYMMDD(jde);
		//printf("Day time of start: %f : %f\n", dayS, DC.SECtoHHMMSS(dayS, jds + 0.125) );

		bool SatVisible = false;
		double stepMin = stepSec / 60.0;
		double JD_start = 0;

		double mpe_start = (jds - JDtle) * 24.0 * 60.0;
		double mpe_stop = (jde - JDtle) * 24.0 * 60.0;

		// Calculate position, velocity,  mpe = "minutes past epoch
		int timeIndex = 0;
		int countPointAdd = 0;

		for (double mpe = mpe_start; mpe <= mpe_stop; mpe += stepMin)
		{
			double JD_current = JDtle + mpe / 24.0 / 60.0;

			// Get the position of the satellite at time "mpe"
			cEciTime eci = orbit->GetPosition(mpe);

			// Now get the "look angle" from the site to the satellite. Note that the ECI object "eciSDP4" contains a time associated
			// with the coordinates it contains; this is the time at which the look angle is valid.
			cTopo topoLook = siteView.GetLookAngle(eci);
			double tle_el = topoLook.ElevationDeg();
			double tle_az = topoLook.AzimuthDeg();

			// verify view angel
			// перевод в ICRF координат телескопа
			double Tel_icrf[3];
			ITRFToICRF(JD_current, Tel_itrf, Tel_icrf);

			// положение спутника по TLE в момент измерения
			double Ptle[3]; // прогноз положения в километрах
			Ptle[0] = eci.Position().m_x;
			Ptle[1] = eci.Position().m_y;
			Ptle[2] = eci.Position().m_z;
			double PICRF[3];

			ConvertTEMEtoICRF(JD_current, Ptle, PICRF);

			// вычисление просчитанных измерений
			double Ra_rad; // rad
			double Dec_rad; // rad
			ConvertXYZtoRADEC(PICRF, Tel_icrf, &Ra_rad, &Dec_rad);

			Angs RaDec;
			RaDec.ang1 = Ra_rad;
			RaDec.ang2 = Dec_rad;

			on_surface tel_pos;
			tel_pos.height = 225.0;
			tel_pos.latitude = 55.8581;
			tel_pos.longitude = 37.4340;
			tel_pos.temperature = 15;
			tel_pos.pressure = 1000;
			Angs sat_view = convertor->RaDec2AzElev(JD_current, RaDec, tel_pos);

			double sat_az = sat_view.ang1 * RG;
			double sat_el = sat_view.ang2 * RG;

			if (sat_el >= SAT_MinEl && sat_el < SAT_MaxEl)
			{
				if (Plan[timeIndex].x < 0 || countPointAdd < maxCountPoint)
				{
					S3DCoordinate pt;
					pt.x = JD_current;
					pt.y = sat_az;
					pt.z = sat_el;
					Plan[timeIndex] = pt;
					PlanIds[timeIndex] = atoi(orbit->SatId().c_str());
					countPointAdd++;
				}
			}

			if (tle_el >= SAT_MinEl && !SatVisible)
			{
				JD_start = JD_current;
				SatVisible = true;
				double dayS = DC.JDtoYYYYMMDD(JD_start + 0.125);
				printf("%s %.0f:%.03f\t", orbit->SatId().c_str(), dayS, DC.SECtoHHMMSS(dayS, JD_start + 0.125));
				printf("EL: %3.5f [%3.5f] Az: %5.5f [%5.5f] ", tle_el, sat_el,  tle_az, sat_az);
			}
			if (tle_el < SAT_MinEl && SatVisible)
			{
				SatVisible = false;
				double duration = (JD_current - JD_start) * 86400.0 / 60.0;
				printf("EL: %3.5f [%3.5f] Az: %5.5f [%5.5f] %.0f min\n", tle_el, sat_el, tle_az, sat_az, duration);
			}

			timeIndex++;
		}
		if (SatVisible)
		{
			printf("\n");
		}
	}

	printf("PLAN:\n");
	FILE* fre = fopen("plan.txt", "w");
	for (int i = 0; i < Plan.size(); i++)
	{
		if (Plan[i].x > 0)
		{
			double Dtime = DC.JDtoYYYYMMDD(Plan[i].x + 0.125);
			double Ptime = DC.SECtoHHMMSS(Dtime, Plan[i].x + 0.125);

			fprintf(fre, "%d %.0f %.03f\t %.10lf %lf  %lf\n", PlanIds[i], Dtime, Ptime, Plan[i].x, Plan[i].y, Plan[i].z );
			printf("%d %.0f %.03f\t %.10lf  %lf  %lf\n", PlanIds[i], Dtime, Ptime, Plan[i].x, Plan[i].y, Plan[i].z);
		}
		else
		{
			printf("...\n");
		}
	}
	fclose(fre);
}

/*
void main()
{
	string DBRoot = "C:\\Users\\Nikolay\\AppData\\Local\\ODSW";
	std::string pinpf_tels = DBRoot + "\\DB\\input\\Telescopes\\tels_simulation.json";
	std::string tle_inpu = DBRoot + "\\DB\\input\\Test\\TLE\\orbit.tle";

	double SAT_elv = 10.0;

    std::vector<Ttelescop> tels;

    {
        rapidjson::Document d;

        if (read_Json(pinpf_tels, d))
        {
            if (d.IsObject())
            {
                if (d.HasMember("Telescopes"))
                {
                    const rapidjson::Value &v_Telescopes = d["Telescopes"];

                    if (v_Telescopes.IsArray())
                    {
                        read_tels_Json(tels, v_Telescopes);
                    }
                    else
                    {
                        std::cout << "Invalid parameter 'Telescopes' in: " << pinpf_tels << std::endl;
                    }
                }
                else
                {
                    std::cout << "Missing parameter 'Telescopes' in: " << pinpf_tels << std::endl;
                }
            }
            else
            {
                std::cout << "Invalid JSON document: " << pinpf_tels << std::endl;
            }
        }
        else
        {
            std::cout << "Failed to read JSON document: " << pinpf_tels << std::endl;
        }
    }

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
				double current_el = topoLook.ElevationDeg();
				double current_az = topoLook.ElevationDeg();

				double sat_pos[3];
				sat_pos[0] = eci.Position().m_x;
				sat_pos[1] = eci.Position().m_y;
				sat_pos[2] = eci.Position().m_z;
				ViewParams vp = ViewTEMESatFromTelescope(currenJd, tel_pos_wgs72, sat_pos);
				


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

*/