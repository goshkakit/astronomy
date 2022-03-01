//
//  Triangulation.hpp
//  Project
//
//  Created by Георгий on 09.12.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#ifndef Triangulation_hpp
#define Triangulation_hpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include "omp.h"
#include <vector>
#include <algorithm>
#include "PotentialCounter.h"

using namespace std;

#define pi 3.141592653589793238462643383279502884L

//==============================================================================//
// Узел сетки при триангуляции
//==============================================================================//
class VerticeWithCompared
{
public:
    double x,y,z,r,theta,fi,U,accR,accTh,accFi;
    int index;
    VerticeWithCompared *compared;
    VerticeWithCompared(){};
    void toSph();
    bool operator < (const VerticeWithCompared &rht) const;
};

//==============================================================================//
// Треугольник при триангуляции
//==============================================================================//
class Triangle
{
public: vector<int> childInd = {-1, -1, -1, -1};
        Triangle();
        Triangle(int index, VerticeWithCompared* V1, VerticeWithCompared* V2, VerticeWithCompared* V3, int fatherIndex = -1);
        Triangle(int index, int V1, int V2, int V3, int fatherIndex = -1);
        int index;
        int fatherInd;
        vector<int> V;
        bool operator == (const Triangle &tr) const;
};

class Triangulation
{
public:
    double r;
    void setRad(double r);
    static double distCounter(VerticeWithCompared V1, VerticeWithCompared V2);
    void zero_triangles();
    void zeroMeshCreator(double r);
    int vIndex;
    int globalIndex;
    int localVIndex;
    void mesher(double r, int degree);
    vector<VerticeWithCompared> vert_arr;
    vector<Triangle> tr_arr;
    vector<vector<VerticeWithCompared>> parVertArr;
    VerticeWithCompared verticeCreator(VerticeWithCompared V1, VerticeWithCompared V2);
    bool vertDetector(VerticeWithCompared V);
    void trianglesCreator(int degree);
    void map(int degree, int polynomDegree, int startRad, int maxRad, int step);
};

#endif /* Triangulation_hpp */
