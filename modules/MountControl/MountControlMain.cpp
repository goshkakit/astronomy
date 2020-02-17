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

	//Test telescope.GetAlphBetPos() from Ra, Dec
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

	//Test telescope.GetAlphBetPos() from Az, Elev
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

	//Test Mat3x3
	std::vector<double> M = {2.0, 5.0, 7.0, 6.0, 3.0, 4.0, 5.0, -2.0, -3.0};
	double Det = Det3x3(M);
	std::vector<double> Dop = AlgebDop3x3(M);
	std::vector<double> T = Transpose3x3(Dop);
	std::vector<double> Inv = Inversion3x3(M);

	//Test Calibration::CalculateA(...)
	Calibration test(Telpos, 2458896.0, 0.0, pi/2);
	std::vector<std::pair<double, double>> real_points = { {0.0, pi / 6.0}, {2.0 * pi / 3.0, pi / 6.0  + pi / 8.0}, { 4.0 * pi / 3.0, pi / 6.0 + pi / 4.0 } };
	double deltAlph = 5 * pi / 180;
	double deltBet = 10 * pi / 180;
	std::vector<std::pair<double, double>> tel_points = { { 0.0 + deltAlph, pi / 6.0 + deltBet },{ 2.0 * pi / 3.0 + deltAlph, pi / 6.0 + pi / 8.0 + deltBet },{ 4.0 * pi / 3.0 + deltAlph, pi / 6.0 + pi / 4.0 + deltBet } };
	std::vector<double> A(6);
	A = test.CalculateA(tel_points, real_points);
	double deltAlph_ = A[2] * 180 / pi;
	double deltBet_ = A[5] * 180 / pi;

	return 0;
}