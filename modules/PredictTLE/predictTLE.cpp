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