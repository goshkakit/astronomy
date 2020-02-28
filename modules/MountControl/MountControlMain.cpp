#include "Convertor.h"
#include "Mat3x3.h"
#include "Calibration.h"
#include "comPortControl.h"

#include <iostream>

int main() {
	Convertor telescope;
	double Telpos[3];
	Telpos[0] = 4604.18;
	Telpos[1] = 3298.4;
	Telpos[2] = 2758.99;
	telescope.SetTelPosITRF(Telpos);
	telescope.InitMountSpec(8, 1000, 5000, 720, 200);
	//Начальное положение
	telescope.SetAzElevPos(2458896.0, 0.0, 0.0);
	std::pair<double, double> AlphBet = telescope.GetAlphBetPos();
	telescope.SetMotorsPos(0.0, 0.0);

	traject tr;

	//
	printf("Connect\n");
	// Init
	ComPortControl CPC;
	int res = CPC.OpenCOM(L"COM4");
	if (res != 0)
	{
		return 0;
	}
	// Wait
	Sleep(2000);
	printf("Start stepper control!\n");
	// Read data after connect
	std::string result = CPC.ReadResponce();
	std::cout << "Recive: " << result << std::endl;

	CPC.getPosition();

	printf("\nSet 2 absolute position:\n");
	int k = 0;
	while (k < 100)
	{
		double Az = 0.0;
		double Elev = 0.0;
		std::cin >> Az;
		std::cin >> Elev;
		
		int flag = telescope.LimitsAzElev(Az * pi / 180, Elev * pi / 180);
		if (flag == -2) {
			std::cout << "Target position out of range!" << std::endl;
		}
		else if (flag == -1) {
			std::cout << "Current position out of range!" << std::endl;
		}
		else if (flag == 0) {
			tr = telescope.GoToAzElev(2458896.0, Az * pi / 180, Elev * pi / 180);
			CPC.runCommand(tr.endpos.first, -tr.endpos.second);
			printf("\nSet 2 absolute position:\n");
			k++;
		}
		else {
			tr = telescope.GoToAzElev(2458896.0, 0.0, 0.0);
			CPC.runCommand(tr.endpos.first, -tr.endpos.second);
			printf("\nSet 2 absolute position:\n");
			tr = telescope.GoToAzElev(2458896.0, Az * pi / 180, Elev * pi / 180);
			CPC.runCommand(tr.endpos.first, -tr.endpos.second);
			printf("\nSet 2 absolute position:\n");
		}
	}

	std::cout << std::endl;
}