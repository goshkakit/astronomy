//
//  ParallelTriangulation.hpp
//  Project
//
//  Created by Георгий on 22.03.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#ifndef ParallelTriangulation_hpp
#define ParallelTriangulation_hpp

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

class ParallelTriangulation
{
public:
    void parMap(int degree, int polynomDegree, int startRad, int maxRad, int step, int threads);
};

#endif /* ParallelTriangulation_hpp */
