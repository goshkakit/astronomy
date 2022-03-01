#include  <stdio.h>

#include <fstream>
#include <iostream>
#include <vector>
#include "common/stV_type.h"

#include "VerifyEtalon/CorrespondenceData.h"

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
	double x, y, z;
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

void read_stV_Json(TstV &stV, const rapidjson::Value &JstV)
{
    if (JstV.IsObject())
    {
        if (JstV.HasMember("JD"))
        {
            const rapidjson::Value &v_JD = JstV["JD"];

            if (v_JD.IsNumber())
            {
                stV.jd = v_JD.GetDouble();
            }
        }

        if (JstV.HasMember("id"))
        {
            const rapidjson::Value &v_JD = JstV["id"];

            if (v_JD.IsNumber())
            {
                stV.ID = v_JD.GetInt();
            }
        }

        if (JstV.HasMember("est"))
        {
            const rapidjson::Value &v_arr = JstV["est"];

            if (v_arr.IsArray())
            {
                if (v_arr.Size() >= 6)
                {
                    if (v_arr[0].IsNumber()) stV.x = v_arr[0].GetDouble();
                    if (v_arr[1].IsNumber()) stV.y = v_arr[1].GetDouble();
                    if (v_arr[2].IsNumber()) stV.z = v_arr[2].GetDouble();
                    if (v_arr[3].IsNumber()) stV.Vx = v_arr[3].GetDouble();
                    if (v_arr[4].IsNumber()) stV.Vy = v_arr[4].GetDouble();
                    if (v_arr[5].IsNumber()) stV.Vz = v_arr[5].GetDouble();

                    if (v_arr.Size() > 6)
                    {
                        if (v_arr[6].IsNumber()) stV.BalFactor = v_arr[6].GetDouble();

                        if (v_arr.Size() > 7)
                        {
                            if (v_arr[7].IsNumber()) stV.SMRatio = v_arr[7].GetDouble();
                        }
                    }
                }
            }
        }
    }
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

void main(int argc, char *argv[])
{
	/* 
	std::string pinpf_tels = argv[1];
	// lazer
	char *htsName = argv[2];
	// tle
	char *tleName = argv[3];
	// simple
	char *optfname = argv[4];
	*/

	// 1 tests
	std::string pinpf_tels = "data/etalon/tels.json";
	char *htsName = "data/etalon/lageos2_cpf_131022_7951.sgf";
	char *tleName = "data/etalon/TLE20131020.txt";
	char *optfname = "data/etalon/track_13_10_22_22_23_simple_num199.dat";
	std::string pinpf_stV = "data/etalon/stV_lageos2.json";
	bool useSTV = false;

	// 2
	//std::string pinpf_tels = "data/etalon2/tels.json";
	//char *htsName = "data/etalon2/cryosat2_cpf_200210_5411.esa";
	//char *tleName = "data/etalon2/TLE20200210.txt";
	//char *optfname = "data/etalon2/track_20_02_10_18_33_simple_num00016078.dat";
	//std:string pinpf_stV = "data/etalon2/stV_lageos2.json";
	//bool useSTV = false;

	/*std::string pinpf_tels = "etalon\\Ajisai\\tel_10989.json";
	char *htsName = "etalon\\lageos2_cpf_131022_7951.sgf";
	char *tleName = "etalon\\TLE20131020.txt";
	char *optfname = "etalon\\Ajisai\\track_19_05_10_02_44_simple_num00000113.dat";
	std:string pinpf_stV = "etalon\\Ajisai\\stV_ajisai_track_19_05_09.json";*/

	/*std::string pinpf_tels = "etalon\\Ajisai\\tel_10995.json";
	char *htsName = "etalon\\lageos2_cpf_131022_7951.sgf";
	char *tleName = "etalon\\TLE20131020.txt";
	char *optfname = "etalon\\Ajisai\\track_19_05_12_20_01_simple_num00031991.dat";
	std:string pinpf_stV = "etalon\\Ajisai\\stV_ajisai_track_19_05_09.json";*/

	/*std::string pinpf_tels = "etalon\\Starlette\\tel_10999.json";
	char *htsName = "etalon\\lageos2_cpf_131022_7951.sgf";
	char *tleName = "etalon\\TLE20131020.txt";
	char *optfname = "etalon\\Starlette\\track_19_04_02_19_52_simple_num00044656.dat";
	std:string pinpf_stV = "etalon\\Starlette\\stV_starlette_190402.json";*/

	/*std::string pinpf_tels = "etalon\\Lares\\tel_10068.json";
	char *htsName = "etalon\\lageos2_cpf_131022_7951.sgf";	
	char *tleName = "etalon\\TLE20131020.txt";
	char *optfname = "etalon\\Lares\\track_19_02_01_08_58_simple_num00020619.dat";
	std:string pinpf_stV = "etalon\\Lares\\stV_lares_track_19_02_01_1.json";*/

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

    TstV stV;

	if (useSTV)
	{
        {
            rapidjson::Document d;

            if (read_Json(pinpf_stV, d))
            {
                if (d.IsObject())
                {
                    read_stV_Json(stV, d);
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

		std::cout << stV.jd << "  " << stV.x << "  " << stV.y << "  " << stV.z << "  " << stV.Vx << "  " << stV.Vy << "  " << stV.Vz << "  " << stV.BalFactor << std::endl;
	}

	CorrespondenceData CDD;
	CDD.InitModyle();

	double telpos[3];
	telpos[0] = tels[0].x;
	telpos[1] = tels[0].y;
	telpos[2] = tels[0].z;

	printf("Telecope position: %f %f %f\n", telpos[0], telpos[1], telpos[2]);
	printf("%s\n%s\n%s\n", htsName, tleName, optfname);

	CDD.RunCorrespondenceData(optfname, htsName, tleName, telpos, stV, useSTV);
}