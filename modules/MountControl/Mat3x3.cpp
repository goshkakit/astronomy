#include "Mat3x3.h"

//Вычисление детерминанта для матрицы 3х3
double Det3x3(const std::vector<double>& M) {
	return (M[0] * (M[4] * M[8] - M[5] * M[7]) - 
			M[1] * (M[3] * M[8] - M[5] * M[6]) +
			M[2] * (M[3] * M[7] - M[4] * M[6]));
}

//Вычисление алгебраического дополнения 3х3
std::vector<double> AlgebDop3x3(const std::vector<double>& M) {
	std::vector<double> Dop(9);
	Dop[0] = M[4] * M[8] - M[5] * M[7];
	Dop[1] = -M[3] * M[8] + M[5] * M[6];
	Dop[2] = M[3] * M[7] - M[4] * M[6];
	Dop[3] = -M[1] * M[8] + M[2] * M[7];
	Dop[4] = M[0] * M[8] - M[2] * M[6];
	Dop[5] = -M[0] * M[7] + M[1] * M[6];
	Dop[6] = M[1] * M[5] - M[2] * M[4];
	Dop[7] = -M[0] * M[5] + M[2] * M[3];
	Dop[8] = M[0] * M[4] - M[1] * M[3];
	return Dop;
}

//Транспонирование матрицы 3x3
std::vector<double> Transpose3x3(const std::vector<double>& M) {
	std::vector<double> T(9);
	T[0] = M[0];
	T[1] = M[3];
	T[2] = M[6];
	T[3] = M[1];
	T[4] = M[4];
	T[5] = M[7];
	T[6] = M[2];
	T[7] = M[5];
	T[8] = M[8];
	return T;
}

//Вычисление обратной матрицы 3x3
std::vector<double> Inversion3x3(const std::vector<double>& M) {
	std::vector<double> Inv(9);
	double Det = Det3x3(M);
	if (Det != 0) {
		std::vector<double> Dop(9);
		Dop = AlgebDop3x3(M);
		std::vector<double> T(9);
		T = Transpose3x3(Dop);
		for (int i = 0; i < T.size(); i++) {
			Inv[i] = T[i] / Det;
		}
	}
	return Inv;
}

//Умножение матрицы 3х3 на столбец 3х1
std::vector<double> Mat3x3XStolb3x1(const std::vector<double>& M, const std::vector<double>& S) {
	std::vector<double> St(3);
	St[0] = M[0] * S[0] + M[1] * S[1] + M[2] * S[2];
	St[1] = M[3] * S[0] + M[4] * S[1] + M[5] * S[2];
	St[2] = M[6] * S[0] + M[7] * S[1] + M[8] * S[2];
	return St;
}

std::vector<double> Inversion2x2(const std::vector<double>& M) {
	std::vector<double> Inv(4);
	double Det = M[0] * M[3] - M[1] * M[2];
	if (Det != 0) {
		Inv[0] = M[3]/Det;
		Inv[1] = -M[2]/Det;
		Inv[2] = -M[1]/Det;
		Inv[3] = M[0]/Det;
	}
	return Inv;
}