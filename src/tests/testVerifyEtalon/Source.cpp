#include  <stdio.h>

#include <fstream>
#include <iostream>
#include "json/json.h"
#include <vector>

#include "VerifyEtalon/CorrespondenceData.h"

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

	tel.x = Jtel["x"].asDouble();
	tel.y = Jtel["y"].asDouble();
	tel.z = Jtel["z"].asDouble();

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

void read_tels_Json(std::vector<Ttelescop> &tels, Json::Value Jtels, std::string &FileName)
{
	unsigned i;
	tels.resize(Jtels.size());
	for (i = 0; i < Jtels.size(); ++i) {
		read_tel_Json(tels[i], Jtels[i], FileName);
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

	std::string pinpf_tels = "etalon\\tels.json";
	char *htsName = "etalon\\lageos2_cpf_131022_7951.sgf";
	char *tleName = "etalon\\TLE20131020.txt";
	char *optfname = "etalon\\track_13_10_22_22_23_simple_num199.dat";

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

	CorrespondenceData CDD;
	CDD.InitModyle();

	double telpos[3];
	telpos[0] = tels[0].x;
	telpos[1] = tels[0].y;
	telpos[2] = tels[0].z;

	printf("Telecope position: %f %f %f\n", telpos[0], telpos[1], telpos[2]);
	printf("%s\n%s\n%s\n", htsName, tleName, optfname);

	CDD.RunCorrespondenceData(optfname, htsName, tleName, telpos);
}