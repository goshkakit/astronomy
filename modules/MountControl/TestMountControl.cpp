#include "Convertor.h"
#include "Mat3x3.h"
//#include "Calibration.h"

#include <iostream>

int Test_() {

	Convertor telescope;
	double Telpos[3];
	Telpos[0] = 4604.18;
	Telpos[1] = 3298.4;
	Telpos[2] = 2758.99;
	telescope.SetTelPosITRF(Telpos);
	std::cout << std::endl;

	//Test Convertor::GetAlphBetPos(...) from Ra, Dec
	std::cout << "Test Convertor::GetAlphBetPos(...) from Ra, Dec" << std::endl;
	{
		double Jd = 2458896.0;
		double Ra = 0.42466583183083417;
		double Dec = 0.87727866072481087;
		telescope.SetRaDecPos({ Jd, Ra, Dec });

		Angs alphbet;
		alphbet = telescope.GetAlphBetPos();
		std::cout << "Jd = " << Jd << ", Ra = " << Ra << ", Dec = " << Dec
			<< " -> Jd = " << alphbet.Jd << ", Alph = " << alphbet.ang1 << ", Bet = " << alphbet.ang2 << std::endl;
		std::cout << "True values for Jd = 2458896.0, Ra = 0.42467, Dec = 0.87727 is Alph = 2.86056, Bet = -1.12296" << std::endl;
	}

	std::cout << std::endl;

	//Test Convertor::GetAlphBetPos(...) from Az, Elev
	std::cout << "Test Convertor::GetAlphBetPos(...) from Az, Elev" << std::endl;
	{
		double Jd = 2458896.0;
		double Az = pi / 6;
		double Elev = pi / 3;
		telescope.SetAzElevPos({ Jd, Az, Elev });

		Angs alphbet;
		alphbet = telescope.GetAlphBetPos();
		std::cout << "Jd = " << Jd << ", Az = " << Az << ", Elev = " << Elev
			<< " -> Jd = " << alphbet.Jd << ", Alph = " << alphbet.ang1 << ", Bet = " << alphbet.ang2 << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 0.523599, Elev = 1.0472 is Alph = 2.86056, Bet = -1.12296" << std::endl;
	}

	std::cout << std::endl;
	
	//Test Convertor::AlphBet2AzElev(...)
	std::cout << "Test Convertor::AlphBet2AzElev(...)" << std::endl;
	{
		std::cout << "Test 1:" << std::endl;
		telescope.SetAzElevPos({ 2458896.0, pi / 3, pi / 3 });
		Angs alphbet = telescope.GetAlphBetPos();
		Angs AzElev = telescope.AlphBet2AzElev(alphbet);
		std::cout << "Jd = " << AzElev.Jd << ", Az = " << AzElev.ang1 * 180 / pi << ", Elev = " << AzElev.ang2 * 180 / pi << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 60, Elev = 60" << std::endl;

		std::cout << "Test 2:" << std::endl;
		telescope.SetAzElevPos({ 2458896.0, 5 * pi / 6, pi / 3 });
		alphbet = telescope.GetAlphBetPos();
		AzElev = telescope.AlphBet2AzElev(alphbet);
		std::cout << "Jd = " << AzElev.Jd << ", Az = " << AzElev.ang1 * 180 / pi << ", Elev = " << AzElev.ang2 * 180 / pi << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 150, Elev = 60" << std::endl;

		std::cout << "Test 3:" << std::endl;
		telescope.SetAzElevPos({ 2458896.0, 4 * pi / 3, pi / 3 });
		alphbet = telescope.GetAlphBetPos();
		AzElev = telescope.AlphBet2AzElev(alphbet);
		std::cout << "Jd = " << AzElev.Jd << ", Az = " << AzElev.ang1 * 180 / pi << ", Elev = " << AzElev.ang2 * 180 / pi << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 240, Elev = 60" << std::endl;

		std::cout << "Test 4:" << std::endl;
		telescope.SetAzElevPos({ 2458896.0, 5 * pi / 3, pi / 3 });
		alphbet = telescope.GetAlphBetPos();
		AzElev = telescope.AlphBet2AzElev(alphbet);
		std::cout << "Jd = " << AzElev.Jd << ", Az = " << AzElev.ang1 * 180 / pi << ", Elev = " << AzElev.ang2 * 180 / pi << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 300, Elev = 60" << std::endl;
	}

	//Test Calibration::CalculateA(...)
	//Calibration test(Telpos, 2458896.0, 0.0, 0.0);
	//std::cout << std::endl;
	{
		/*
		std::cout << "Test Calibration::CalculateA(...)" << std::endl;
		std::vector<std::pair<double, double>> real_points = { { 0, pi },{ pi, pi / 2 },{ 3 * pi / 2, pi } };
		double deltAlph = 10 / 180 * pi;
		double deltBet = 10 / 180 * pi;
		//double deltBet = pi;
		std::vector<std::pair<double, double>> tel_points = { { 0 + deltAlph, pi + deltBet },{ pi + deltAlph, pi / 2 + deltBet },{ 3 * pi / 2 + deltAlph, pi + deltBet } };
		std::vector<double> R(9);
		R = test.CalculateA(tel_points, real_points);
		std::vector<double> Y1 = { real_points[0].first, real_points[0].second, 1 };
		std::vector<double> Y2 = { real_points[1].first, real_points[1].second, 1 };
		std::vector<double> Y3 = { real_points[2].first, real_points[2].second, 1 };
		std::vector<double> X1 = Mat3x3XStolb3x1(R, Y1);
		std::vector<double> X2 = Mat3x3XStolb3x1(R, Y2);
		std::vector<double> X3 = Mat3x3XStolb3x1(R, Y3);

		std::cout << "CalcAlph1 - Alph1 = " << X1[0] - tel_points[0].first << ", CalcBet1 - Bet1 = " << X1[1] - tel_points[0].second << ", CalcAlph2 - Alph2 = " << X2[0] - tel_points[1].first
			<< ", CalcBet2 - Bet2 = " << X2[1] - tel_points[1].second << ", CalcAlph3 - Alph3 = " << X3[0] - tel_points[2].first << ", CalcBet3 - Bet3 = " << X3[1] - tel_points[2].second << std::endl;
		*/
	}

	std::cout << std::endl;

	//Test Convertor::CalcTraject(...)
	{
		std::cout << "Test Convertor::CalcTraject(...)" << std::endl;
		MotPos motor = telescope.GetMotorsPos();

		std::cout << "Test 1:" << std::endl;
		std::cout << "inAlph = 0.0, inBet = 0.0, outAlph = pi/4, outBet = pi/4" << std::endl;
		std::cout << "stepsAlph = 144000, stepsBet = 144000, T = 24.4" << std::endl;
		traject tr = telescope.CalcTraject({ 0.0, 0.0, 0.0 }, { 0.0, pi / 4, pi / 4 });
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T << std::endl;

		std::cout << "Test 2:" << std::endl;
		std::cout << "inAlph = 0.0, inBet = 0.0, outAlph = pi/4, outBet = pi/2" << std::endl;
		std::cout << "stepsAlph = 144000, stepsBet = 288000, T = 38.8" << std::endl;
		tr = telescope.CalcTraject({ 0.0, 0.0, 0.0 }, { 0.0, pi / 4, pi / 2 });
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T << std::endl;

		std::cout << "Test 3:" << std::endl;
		std::cout << "inAlph = 0.0, inBet = 0.0, outAlph = pi/8, outBet = pi/4" << std::endl;
		std::cout << "stepsAlph = 72000, stepsBet = 144000, T = 24.4" << std::endl;
		tr = telescope.CalcTraject({ 0.0, 0.0, 0.0 }, { 0.0, pi / 8, pi / 4 });
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T << std::endl;

		std::cout << "Test 4:" << std::endl;
		std::cout << "inAlph = pi/4, inBet = 0.0, outAlph = 0, outBet = 0" << std::endl;
		std::cout << "stepsAlph = -144000, stepsBet = 0, T = 24.4, VAlph = -10000, VBet = 0, aAlph = -1000, aBet = 0" << std::endl;
		tr = telescope.CalcTraject({ 0.0, pi / 4, 0.0 }, { 0.0, 0.0, 0.0 });
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T
			<< ", VAlph = " << tr.V.first << ", VBet = " << tr.V.second << ", aAlph = " << tr.a.first << ", aBet = " << tr.a.second << std::endl;

		std::cout << "Test 5:" << std::endl;
		std::cout << "inAlph = pi/16, inBet = pi/16, outAlph = 0, outBet = 0" << std::endl;
		std::cout << "stepsAlph = -36000, stepsBet = -36000, T = 12, VAlph = -6000, VBet = -6000, aAlph = -1000, aBet = -1000" << std::endl;
		tr = telescope.CalcTraject({ 0.0, pi / 16, pi / 16 }, { 0.0, 0.0, 0.0 });
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T
			<< ", VAlph = " << tr.V.first << ", VBet = " << tr.V.second << ", aAlph = " << tr.a.first << ", aBet = " << tr.a.second << std::endl;
	}

	std::cout << std::endl;

	return 0;
}