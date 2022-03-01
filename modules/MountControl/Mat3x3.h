#pragma once

#include <vector>

//Вычисление детерминанта для матрицы 3х3
double Det3x3(const std::vector<double>& M);

//Вычисление алгебраического дополнения 3х3
std::vector<double> AlgebDop3x3(const std::vector<double>& M);

//Транспонирование матрицы 3x3
std::vector<double> Transpose3x3(const std::vector<double>& M);

//Вычисление обратной матрицы 3x3
std::vector<double> Inversion3x3(const std::vector<double>& M);

//Умножение матрицы 3х3 на столбец 3х1
std::vector<double> Mat3x3XStolb3x1(const std::vector<double>& M, const std::vector<double>& S);

//Обращение матрицы 2х2
std::vector<double> Inversion2x2(const std::vector<double>& M);