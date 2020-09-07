#pragma once

// File processing
#include "JsonWorker.h"

// Common
#include "common/TLELoader.h"
#include "common/DataConverter.h"
#include "common/mytypes.h"

class AstroDataManager
{
public:
	AstroDataManager() {}

	JSONWorker jsonWorker;
	TLELoader tleLoader;
	std::vector<TelescopObject> tels;

	void loadTLE(std::string tle_input)
	{
		tleLoader.LoadData(tle_input.c_str(), 2);
	}

	void loadTelescope(std::string pinpf_tels)
    {
        {
            rapidjson::Document d;

            if (jsonWorker.read_Json(pinpf_tels, d))
            {
                if (d.IsObject())
                {
                    if (d.HasMember("Telescopes"))
                    {
                        const rapidjson::Value &v_Telescopes = d["Telescopes"];

                        if (v_Telescopes.IsArray())
                        {
                            jsonWorker.read_tels_Json(tels, v_Telescopes);
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
	}
};