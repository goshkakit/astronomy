//
//  Extrapolation.hpp
//  Project
//
//  Created by Георгий on 04.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#ifndef Extrapolation_hpp
#define Extrapolation_hpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <random>
#include <ctime>
#include "omp.h"
#include "PotentialCounter.h"
#include "Triangulation.hpp"

using namespace std;

//==============================================================================//
// Узел сетки
//==============================================================================//
class Vertice
{
public:
    int index;
    double U, r, theta, fi, x, y, z, accR, accTh, accFi;
    void ReadFromStr(const char* str);
};

//==============================================================================//
// Треугольник
//==============================================================================//
class Tr
{
public:
    int index, fatherInd;
    vector<int> V = {-1, -1, -1};
    vector<int> childInd = {-1, -1, -1, -1};
    void ReadFromStr(const char* str);
};

class Extrapolation{
public:
    int stepRad;
    int startRad;
    int maxRad;
    int trIdx;
    int vertAm;
    int trAm;
    void stepSet(int step);
    void radSet(int stR, int maxR);
    void verticeLoader(string filename);
    void triangleLoader(string filename);
    vector<vector<Vertice>> vert_arr;
    vector<Tr> tr_arr;
    double sphDist(double theta, double fi, int idx, int rIdx);
    void loader(int pol, int deg, int maxR, int step, int stR=r0);
    int zeroSearcher(int rIdx, double x, double y, double z);
    int searcher(int fthrIdx, int rIdx, double x, double y, double z);
    int parSearcher(int fthrIdx, int rIdx, double x, double y, double z);
    bool isInTr(int rIdx, int idx, double x, double y, double z);
    vector<double> layerCounter(double theta, double fi, int trIdx, int rIdx);
    vector<double> counter(double r, double theta, double fi, int tr, int rIdx);
    vector<double> extrapolator(double r, double theta, double fi);
};

#endif /* Extrapolation_hpp */
