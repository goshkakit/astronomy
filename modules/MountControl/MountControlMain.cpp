#include "Convertor.h"
#include "Mat3x3.h"
#include "Calibration.h"

#include <iostream>

int main() {
	Convertor telescope;
	double Telpos[3];
	Telpos[0] = 4604.18;
	Telpos[1] = 3298.4;
	Telpos[2] = 2758.99;
	telescope.SetTelPosITRF(Telpos);
	std::cout << std::endl;

	Calibration test(Telpos, 2458896.0, 0.0, pi / 2);
	std::cout << std::endl;

	//Test Convertor::GetAlphBetPos(...) from Ra, Dec
	std::cout << "Test Convertor::GetAlphBetPos(...) from Ra, Dec" << std::endl;
	{
		double Jd = 2458896.0;
		double Ra = 0.42466583183083417;
		double Dec = 0.87727866072481087;
		telescope.SetRaDecPos(Jd, Ra, Dec);

		std::pair<double, double> alphbet;
		alphbet = telescope.GetAlphBetPos();
		std::cout << "Jd = " << Jd << ", Ra = " << Ra << ", Dec = " << Dec
			<< " -> Alph = " << alphbet.first << ", Bet = " << alphbet.second << std::endl;
		std::cout << "True values for Jd = 2458896.0, Ra = 0.42467, Dec = 0.87727 is Alph = 0.523599, Bet = 1.0472" << std::endl;
	}

	std::cout << std::endl;

	//Test Convertor::GetAlphBetPos(...) from Az, Elev
	std::cout << "Test Convertor::GetAlphBetPos(...) from Az, Elev" << std::endl;
	{
		double Jd = 2458896.0;
		double Az = pi/6;
		double Elev = pi/3;
		telescope.SetAzElevPos(Jd, Az, Elev);

		std::pair<double, double> alphbet;
		alphbet = telescope.GetAlphBetPos();
		std::cout << "Jd = " << Jd << ", Az = " << Az << ", Elev = " << Elev
			<< " -> Alph = " << alphbet.first << ", Bet = " << alphbet.second << std::endl;
		std::cout << "True values for Jd = 2458896.0, Az = 0.523599, Elev = 1.0472 is Alph = 0.523599, Bet = 1.0472" << std::endl;
	}

	std::cout << std::endl;

	//Test Calibration::CalculateA(...)
	{
		std::cout << "Test Calibration::CalculateA(...)" << std::endl;
		std::vector<std::pair<double, double>> real_points = { { 0.0, pi / 6.0 },{ 2.0 * pi / 3.0, pi / 6.0 + pi / 8.0 },{ 4.0 * pi / 3.0, pi / 6.0 + pi / 4.0 } };
		double deltAlph = 5 * pi / 180;
		double deltBet = 10 * pi / 180;
		std::vector<std::pair<double, double>> tel_points = { { 0.0 + deltAlph, pi / 6.0 + deltBet },{ 2.0 * pi / 3.0 + deltAlph, pi / 6.0 + pi / 8.0 + deltBet },{ 4.0 * pi / 3.0 + deltAlph, pi / 6.0 + pi / 4.0 + deltBet } };
		std::vector<double> A(6);
		A = test.CalculateA(tel_points, real_points);
		double deltAlph_ = A[2] * 180 / pi;
		double deltBet_ = A[5] * 180 / pi;

		std::cout << "deltAlph = 5 deg, deltBet = 10 deg" << std::endl;
		std::cout << "Calc_deltAlph = " << deltAlph_ << " deg, Calc_deltBet = " << deltBet_ << " deg" << std::endl;
	}

	std::cout << std::endl;

	//Test Convertor::CalcTraject(...)
	{
		std::cout << "Test Convertor::CalcTraject(...)" << std::endl;
		std::pair<double, double> motor = telescope.GetMotorsPos();

		std::cout << "Test 1:" << std::endl;
		std::cout << "inAlph = 0.0, inBet = 0.0, outAlph = pi/4, outBet = pi/4" << std::endl;
		std::cout << "stepsAlph = 298000, stepsBet = 298000, T = 39.8" << std::endl;
		traject tr = telescope.CalcTraject(0.0, 0.0, pi / 4, pi / 4);
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T << std::endl;

		std::cout << "Test 2:" << std::endl;
		std::cout << "inAlph = 0.0, inBet = 0.0, outAlph = pi/4, outBet = pi/2" << std::endl;
		std::cout << "stepsAlph = 298000, stepsBet = 596000, T = 69.6" << std::endl;
		tr = telescope.CalcTraject(0.0, 0.0, pi / 4, pi / 2);
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T << std::endl;

		std::cout << "Test 3:" << std::endl;
		std::cout << "inAlph = 0.0, inBet = 0.0, outAlph = pi/8, outBet = pi/4" << std::endl;
		std::cout << "stepsAlph = 149000, stepsBet = 298000, T = 39.8" << std::endl;
		tr = telescope.CalcTraject(0.0, 0.0, pi / 8, pi / 4);
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T << std::endl;

		std::cout << "Test 4:" << std::endl;
		std::cout << "inAlph = pi/4, inBet = 0.0, outAlph = 0, outBet = 0" << std::endl;
		std::cout << "stepsAlph = -298000, stepsBet = 0, T = 39.8, VAlph = -10000, VBet = 0, aAlph = -1000, aBet = 0" << std::endl;
		tr = telescope.CalcTraject(pi / 4, 0.0, 0.0, 0.0);
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T
			<< ", VAlph = " << tr.V.first << ", VBet = " << tr.V.second << ", aAlph = " << tr.a.first << ", aBet = " << tr.a.second << std::endl;

		std::cout << "Test 5:" << std::endl;
		std::cout << "inAlph = pi/16, inBet = pi/16, outAlph = 0, outBet = 0" << std::endl;
		std::cout << "stepsAlph = -74500, stepsBet = -74500, T = 17.2627, VAlph = -8631.34, VBet = -8631.34, aAlph = -1000, aBet = -1000" << std::endl;
		tr = telescope.CalcTraject(pi / 16, pi / 16, 0.0, 0.0);
		std::cout << "Real_stepsAlph = " << tr.endpos.first - tr.startpos.first << ", Real_stepsBet = " << tr.endpos.second - tr.startpos.second << ", Real_T = " << tr.T
			<< ", VAlph = " << tr.V.first << ", VBet = " << tr.V.second << ", aAlph = " << tr.a.first << ", aBet = " << tr.a.second << std::endl;
	}

	std::cout << std::endl;

	return 0;
}