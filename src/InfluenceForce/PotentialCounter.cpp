//
//  PotentialCounter.cpp
//  Project
//
//  Created by Георгий on 30.11.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#include "PotentialCounter.h"

//==============================================================================//
// Считывание коэффициентов разложения из строки
//==============================================================================//
void data_row::ReadFromStr(const char* str)
{
    sscanf(str, "%d %d %le %le", &n, &m, &a, &b);
}

//==============================================================================//
// Считывание коэффициентов разложения из файла
// Файл "C:/Users/user/source/repos/potential/Project/data.txt" или "data.txt"
//==============================================================================//
void PotentialCounter::LoadFromFile(string filename){
    string line;
    ifstream in(filename);
    if (!in)
        cout << "Check coefs file location" << endl;
    else
    {
        const long len = 90, strings = 65338;
        const char ch = '\n';
        char mass[len] = {};
        this->data = new data_row[strings];
        for(int r = 0; r<strings; r++)
        {
            memset(mass, 0, len);
            in.getline(mass, len-1,ch);
            data[r].ReadFromStr(mass);
        }
    }
}

//==============================================================================//
// Выбор коэффициента С из массива
//==============================================================================//
double PotentialCounter::Cnm(int n,int m){
    int index = (n+3)*(n-2)/2 + m ;
    return data[index].a;
}

//==============================================================================//
// Выбор коэффициента S из массива
//==============================================================================//
double PotentialCounter::Snm(int n,int m){
    int index = (n+3)*(n-2)/2 + m ;
    return data[index].b;
}

//==============================================================================//
// Нормирующий множитель
//==============================================================================//
long double PotentialCounter::coefCounter(int n, int m)
{
    long double coef = 2.0*(2.0*n+1.0);
    if(m != 0){
        for(int i = n-m+1; i<=n+m; i++){
            coef/=i;
        }
    }
    return sqrt(abs(coef));
}

//==============================================================================//
// Вычисление соответствующего полинома Лежандра
//==============================================================================//
void PotentialCounter::LegendreCounter(int n, int m, double x){
    double fact;
    int j;
    int k;
    legArr = new double[n+1];
    for (j = 0; j < n + 1; j++)
    {
        legArr[j] = 0.0;
    }
    if (m <= n)
    {
        legArr[m] = 1.0;
        fact = 1.0;
        for (k = 0; k < m; k++)
        {
            legArr[m] = - legArr[m] * fact * sqrt (1.0 - x * x);
            fact = fact + 2.0;
        }
    }
    if (m + 1 <= n)
    {
        legArr[m+1] = x * (double) (2*m+1) * legArr[m];
    }
    for (j = m+2; j <= n; j++)
     {
         legArr[j] = ((double)(2 * j     - 1) * x * legArr[j-1]
                     +(double)(  - j - m + 1) *     legArr[j-2] )
                     /(double)(    j - m    );
     }
}

//==============================================================================//
// Вычисление члена ряда
//==============================================================================//
long double PotentialCounter::partSumCounter(double r, double fi, double lmbd, int n, int m){
    this->LegendreCounter(n,m,sin(fi));
    long double partSum = 0;
    partSum = pow((r0/r),(n+1)) * legArr[n] *
              (Cnm(n, m) * cos(m * lmbd) + Snm(n, m) * sin(m * lmbd)) * this->coefCounter(n, m);
    this->cleanLeg();
    return partSum;
}

void PotentialCounter::cleanLeg(){
    delete[] legArr;
    legArr = nullptr;
}

//==============================================================================//
// Вычисление ряда без первых двух гармоник
//==============================================================================//
long double PotentialCounter::potential(double r, double fi, double lmbd){
    long double sum = 0;
    for (int n=3; n<=length; n++) {
        for(int m=0; m<=n; m++){
            sum += this->partSumCounter(r, fi, lmbd, n, m);
        }
    }
    return f_c * M_c * sum/r0;
}

//==============================================================================//
// Вычисление нулевой гармоники
//==============================================================================//
double PotentialCounter::U0(double r){
    return f_c * M_c / r;
}

//==============================================================================//
// Вычисление второй гармоники
//==============================================================================//
double PotentialCounter::U2(double r, double theta, double fi){
    return -f_c * M_c * J2 * r0 * r0 * (3* sin(theta) * sin(theta) - 1) / (2 * r * r * r) +
            this->coefCounter(2, 1) * f_c * M_c * r0 * r0 * (-3 * sin(theta) * sqrt(1 - sin(theta) * sin(theta))) *
                                      (Cnm(2, 1) * cos(fi) + Snm(2, 1) * sin(fi)) / (r * r * r) +
            this->coefCounter(2, 2) * f_c * M_c * r0 * r0 * (3 * (1 - sin(theta) * sin(theta))) *
                                      (Cnm(2, 2) * cos(2*fi) + Snm(2, 2) * sin(2*fi)) / (r * r * r);
}

//==============================================================================//
// Вычисление радиальной компоненты ускорения
//==============================================================================//
double PotentialCounter::accR(double r, double theta, double fi){
    double Up = U0(r + dR) + U2(r + dR, theta, fi) + this->potential(r + dR, theta, fi);
    double Ul = U0(r)      + U2(r, theta, fi)      + this->potential(r,      theta, fi);
    return (Up - Ul) / dR;
}

//==============================================================================//
// Вычисление theta-компоненты ускорения
//==============================================================================//
double PotentialCounter::accTh(double r, double theta, double fi){
    double Up = U2(r, theta + dTh, fi) + this->potential(r, theta + dTh, fi);
    double Ul = U2(r, theta      , fi) + this->potential(r, theta      , fi);
    return (Up - Ul) / (r0*dTh);
}

//==============================================================================//
// Вычисление fi-компоненты ускорения
//==============================================================================//
double PotentialCounter::accFi(double r, double theta, double fi){
    double Up = U2(r, theta, fi + dFi) + this->potential(r, theta, fi + dFi);
    double Ul = U2(r, theta, fi) + this->potential(r, theta, fi);
    return (Up - Ul) / (r0*dFi);
}

//==============================================================================//
// Вычисление первой гармоники радиальной компоненты ускорения
//==============================================================================//
double PotentialCounter::accR0(double r, double theta, double fi) {
    return (U0(r + dR) + U2(r + dR, theta, fi) - U0(r) - U2(r, theta, fi)) / dR;
}

//==============================================================================//
// Вычисление первой гармоники theta-компоненты ускорения
//==============================================================================//
double PotentialCounter::accTh0(double r, double theta, double fi) {
    return (U2(r, theta + dTh, fi) - U2(r, theta, fi)) / (r0 * dTh);
}

//==============================================================================//
// Вычисление первой гармоники fi-компоненты ускорения
//==============================================================================//
double PotentialCounter::accFi0(double r, double theta, double fi) {
    return (U2(r, theta, fi + dFi) - U2(r, theta, fi)) / (r0 * dFi);
}

/*
    Пример запуска:

    PotentialCounter pc = {};
    pc.length = 100;
    pc.LoadFromFile("maps_gr/data.txt");
    //или pc.LoadFromFile("data.txt");

    double U_pc = pc.potential(r, theta, fi);
    double aR = pc.accR(r, theta, fi);
    double aTh = pc.accTh(r, theta, fi);
    ouble aFi = pc.accFi(r, theta, fi);
*/
