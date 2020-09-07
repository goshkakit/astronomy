#pragma once

#include <fstream>
#include <iostream>

#include "AstroTypes.h"

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

class JSONWorker
{
private:
	const double RG = 180.0 / M_PI;	// Radian-degree convertion
	const double GR = 1.0 / RG;

public:
	JSONWorker() {}

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

    void read_tel_Json(TelescopObject &tel, const rapidjson::Value &value)
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

    void read_tels_Json(std::vector<TelescopObject> &tels, const rapidjson::Value &Jtels)
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

};