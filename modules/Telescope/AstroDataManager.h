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
		Json::Value inpf, JTels;


		inpf = jsonWorker.read_Json(pinpf_tels);
		if (inpf.isMember("Telescopes")) {
			JTels = inpf["Telescopes"];
		}
		else std::cout << "Missing parameter 'Telescopes' in " << pinpf_tels << "\n";

		jsonWorker.read_tels_Json(tels, JTels, pinpf_tels);
		if (tels.size() == 0) {
			std::cout << "Empty array of telescopes in " << pinpf_tels << "\n";
		};

		std::cout << "Count of telescopes in " << tels.size() << "\n";
	}
};