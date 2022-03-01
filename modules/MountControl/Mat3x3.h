#pragma once

#include <vector>

//���������� ������������ ��� ������� 3�3
double Det3x3(const std::vector<double>& M);

//���������� ��������������� ���������� 3�3
std::vector<double> AlgebDop3x3(const std::vector<double>& M);

//���������������� ������� 3x3
std::vector<double> Transpose3x3(const std::vector<double>& M);

//���������� �������� ������� 3x3
std::vector<double> Inversion3x3(const std::vector<double>& M);

//��������� ������� 3�3 �� ������� 3�1
std::vector<double> Mat3x3XStolb3x1(const std::vector<double>& M, const std::vector<double>& S);

//��������� ������� 2�2
std::vector<double> Inversion2x2(const std::vector<double>& M);