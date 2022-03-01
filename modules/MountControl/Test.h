#pragma once
#include <iostream>

#include "MountController.h"

int TestSetCurrentDirection() {
	double Telpos[3];
	Telpos[0] = 4604.18;
	Telpos[1] = 3298.4;
	Telpos[2] = 2758.99;

	std::shared_ptr<NewConvertor> convertor = std::make_shared<NewConvertor>();
	MountSpecification mount_spec;
	MountController Controller(convertor, mount_spec, Telpos);

	//Test MountController::SetCurrentDirectionRaDec(...)
	std::cout << "Test MountController::SetCurrentDirectionRaDec(...)" << std::endl;
	{
		double Jd = 2458896.0;
		double Ra = 0.42466583183083417;
		double Dec = 0.87727866072481087;

		Controller.SetCurrentDirectionRaDec(Jd, { Ra, Dec });

		CurrentDirection cur_dir = Controller.GetCurrentDirection();
		std::cout << "Jd = " << Jd << ", Ra = " << Ra << ", Dec = " << Dec
			<< " -> Jd = " << cur_dir.Jd << ", Alph = " << cur_dir.OwnAxes.ang1 << ", Bet = " << cur_dir.OwnAxes.ang2 << std::endl;
		std::cout << "True values for Jd = 2458896.0, Ra = 0.42467, Dec = 0.87727 is Alph = 2.86056, Bet = -1.12296" << std::endl;
	}
	std::cout << std::endl;

	//Test MountController::SetCurrentDirectionAzElev(...)
	std::cout << "Test MountController::SetCurrentDirectionAzElev(...)" << std::endl;
	{
		double Jd = 2458896.0;
		double Az = M_PI / 6;
		double Elev = M_PI / 3;

		Controller.SetCurrentDirectionAzElev(Jd, { Az, Elev });

		CurrentDirection cur_dir = Controller.GetCurrentDirection();
		std::cout << "Jd = " << Jd << ", Az = " << Az << ", Elev = " << Elev
			<< " -> Jd = " << cur_dir.Jd << ", Alph = " << cur_dir.OwnAxes.ang1 << ", Bet = " << cur_dir.OwnAxes.ang2 << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 0.523599, Elev = 1.0472 is Alph = 2.86056, Bet = -1.12296" << std::endl;
	}
	std::cout << endl;

	//Test MountController::SetCurrentDirectionOwnAxes(...)
	std::cout << "Test MountController::SetCurrentDirectionOwnAxes(...)" << std::endl;
	{
		std::cout << "Test 1:" << std::endl;
		Controller.SetCurrentDirectionAzElev(2458896.0, { M_PI / 3 , M_PI / 3 });
		Controller.SetCurrentDirectionOwnAxes(2458896.0, Controller.GetCurrentDirection().OwnAxes);
		CurrentDirection cur_dir = Controller.GetCurrentDirection();
		std::cout << "Jd = " << cur_dir.Jd << ", Az = " << cur_dir.AzElev.ang1 * 180 / M_PI << ", Elev = " << cur_dir.AzElev.ang2 * 180 / M_PI << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 60, Elev = 60" << std::endl;

		std::cout << "Test 2:" << std::endl;
		Controller.SetCurrentDirectionAzElev(2458896.0, { 5 * M_PI / 6 , M_PI / 3 });
		Controller.SetCurrentDirectionOwnAxes(2458896.0, Controller.GetCurrentDirection().OwnAxes);
		cur_dir = Controller.GetCurrentDirection();
		std::cout << "Jd = " << cur_dir.Jd << ", Az = " << cur_dir.AzElev.ang1 * 180 / M_PI << ", Elev = " << cur_dir.AzElev.ang2 * 180 / M_PI << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 150, Elev = 60" << std::endl;

		std::cout << "Test 3:" << std::endl;
		Controller.SetCurrentDirectionAzElev(2458896.0, { 4 * M_PI / 3 , M_PI / 3 });
		Controller.SetCurrentDirectionOwnAxes(2458896.0, Controller.GetCurrentDirection().OwnAxes);
		cur_dir = Controller.GetCurrentDirection();
		std::cout << "Jd = " << cur_dir.Jd << ", Az = " << cur_dir.AzElev.ang1 * 180 / M_PI << ", Elev = " << cur_dir.AzElev.ang2 * 180 / M_PI << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 240, Elev = 60" << std::endl;

		std::cout << "Test 4:" << std::endl;
		Controller.SetCurrentDirectionAzElev(2458896.0, { 5 * M_PI / 3 , M_PI / 3 });
		Controller.SetCurrentDirectionOwnAxes(2458896.0, Controller.GetCurrentDirection().OwnAxes);
		cur_dir = Controller.GetCurrentDirection();
		std::cout << "Jd = " << cur_dir.Jd << ", Az = " << cur_dir.AzElev.ang1 * 180 / M_PI << ", Elev = " << cur_dir.AzElev.ang2 * 180 / M_PI << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 300, Elev = 60" << std::endl;
	}

	return 0;
}