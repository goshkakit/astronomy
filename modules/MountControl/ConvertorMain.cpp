#include "Convertor.h"

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
		std::cout << "True values for Jd = 2458896.0, Ra = 0.42467, Dec = 0.87727 is Alph = 1.0472, Bet = 0.523599" << std::endl;
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
		std::cout << "True values for Jd = 2458896.0, Az = 0.523599, Elev = 1.0472 is Alph = 1.0472, Bet = 0.523599" << std::endl;
	}

	return 0;
}