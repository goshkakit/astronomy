#include "Calibration.h"
#include "Mat3x3.h"

Calibration::Calibration(double* pos, const double& Jd, const double& Alph, const double& Bet) {
	telpos.SetTelPosITRF(pos);
	telpos.SetAlphBetPos(Jd, Alph, Bet);
}

Calibration::Calibration(const double& x, const double& y, const double& z, const double& Jd, const double& Alph, const double& Bet) {
	telpos.SetTelPosITRF(x, y, z);
	telpos.SetAlphBetPos(Jd, Alph, Bet);
}

Calibration::~Calibration() {
}

void Calibration::Calibrate() {

}

std::vector<double> Calibration::CalculateA(const std::vector<std::pair<double, double>>& tel_points, const std::vector<std::pair<double, double>>& real_points) {
	/*
	std::vector<double> A(9), b1(3), b2(3), x1(3), x2(3);
	A[0] = real_points[0].first;
	A[1] = real_points[0].second;
	A[2] = 1.0;
	A[3] = real_points[1].first;
	A[4] = real_points[1].second;
	A[5] = 1.0;
	A[6] = real_points[2].first;
	A[7] = real_points[2].second;
	A[8] = 1;
	
	b1[0] = tel_points[0].first;
	b1[1] = tel_points[1].first;
	b1[2] = tel_points[2].first;

	b2[0] = tel_points[0].second;
	b2[1] = tel_points[1].second;
	b2[2] = tel_points[2].second;

	std::vector<double> InvA(9);
	InvA = Inversion3x3(A);
	x1 = Mat3x3XStolb3x1(InvA, b1);
	x2 = Mat3x3XStolb3x1(InvA, b2);

	std::vector<double> res(6);
	for (int i = 0; i < 3; i++) {
		res[i] = x1[i];
		res[i + 3] = x2[i];
	}
	*/

	std::vector<double> res(6);
	double SumDelt1 = 0, SumDelt2 = 0;
	for (int i = 0; i < tel_points.size(); i++) {
		SumDelt1 += (tel_points[i].first - real_points[i].first);
		SumDelt2 += (tel_points[i].second - real_points[i].second);
	}
	res[0] = 1.0;
	res[1] = 0.0;
	res[2] = SumDelt1 / tel_points.size();
	res[3] = 0.0;
	res[4] = 1.0;
	res[5] = SumDelt2 / tel_points.size();

	return res;
}