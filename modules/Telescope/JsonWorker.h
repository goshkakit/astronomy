#pragma once

#include <fstream>
#include <iostream>

#include "json/json.h"
#include "AstroTypes.h"

class JSONWorker
{
private:
	const double RG = 180.0 / M_PI;	// Radian-degree convertion
	const double GR = 1.0 / RG;

public:
	JSONWorker() {}

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

	void read_tel_Json(TelescopObject &tel, Json::Value Jtel, std::string &FileName)
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

	void read_tels_Json(std::vector<TelescopObject> &tels, Json::Value Jtels, std::string &FileName)
	{
		unsigned i;
		tels.resize(Jtels.size());
		for (i = 0; i < Jtels.size(); ++i) {
			read_tel_Json(tels[i], Jtels[i], FileName);
		}
	}

};