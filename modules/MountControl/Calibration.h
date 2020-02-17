#pragma once
#include "Convertor.h"

#include <vector>

class Calibration {
public:
	Calibration(double* pos, const double& Jd, const double& Alph, const double& Bet);
	Calibration(const double& x, const double& y, const double& z, const double& Jd, const double& Alph, const double& Bet);
	~Calibration();

	void Calibrate();

	//������� ���������� ������� ��������������
	std::vector<double> CalculateA(const std::vector<std::pair<double, double>>& tel_points, const std::vector<std::pair<double, double>>& real_points);
private:
	//��������� ��
	Convertor telpos;

	//�����
	double backlash_pp, backlash_pm, backlash_mm, backlash_mp;

	//������� �������������
	std::vector<double> App, Apm, Amm, Amp;
};