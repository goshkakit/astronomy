//
//  ParallelTriangulation.cpp
//  Project
//
//  Created by Георгий on 22.03.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "ParallelTriangulation.hpp"

//==============================================================================//
// Распараллеленное построение сетки
//==============================================================================//
void ParallelTriangulation::parMap(int degree, int polynomDegree, int startRad, int maxRad, int step, int threads){
    Triangulation tr = {};
    ofstream infileV;
    string deg  = to_string(degree);
    string polD = to_string(polynomDegree);
    string st   = to_string(step);
    string stR = to_string(startRad);
    string mxR = to_string(maxRad);
    infileV.open("M_" + deg + "_" + polD + "_" + st + ".txt", std::ios_base::app);
    cout << "Work has started" << endl;
    clock_t tf;
    tf = clock();
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#pragma omp parallel private(tr) 
    {
        int prNum = omp_get_num_threads();
        int thN = omp_get_thread_num();
        int str = startRad + thN * (maxRad - startRad) / prNum;
        int enr = 0;
        if (thN != (prNum - 1)) {
            enr = startRad + (thN + 1) * (maxRad - startRad) / prNum - step;
        }   
        else {
            enr = startRad + (thN + 1) * (maxRad - startRad) / prNum;
        }
       // cout << "Proc " << thN << " is working on range " << str << "-" << enr << endl;
        tr.map(degree, polynomDegree, str, enr, step);
#pragma omp for ordered
        for(int k=0; k<threads; k++){
#pragma omp ordered
        {
            cout<<"Proc "<<k<<" is writing in file"<<endl;
            for(int i=0; i<tr.parVertArr.size(); i++){
                for(int j=0; j<tr.parVertArr[i].size(); j++){
                    infileV << tr.parVertArr[i][j].index   << " " <<
                        tr.parVertArr[i][j].U              << " " <<
                        tr.parVertArr[i][j].r              << " " <<
                        tr.parVertArr[i][j].theta - pi / 2 << " " <<
                        tr.parVertArr[i][j].fi + pi        << " " <<
                        tr.parVertArr[i][j].x              << " " <<
                        tr.parVertArr[i][j].y              << " " <<
                        tr.parVertArr[i][j].z              << " " << 
                        tr.parVertArr[i][j].accR           << " " <<
                        tr.parVertArr[i][j].accTh          << " " <<
                        tr.parVertArr[i][j].accFi          << "\n";
                }
            }
        }
        }
        
    }
    tf = clock() - tf;
    cout<<"Global map is done in "<<(double)tf/(CLOCKS_PER_SEC*60)<<" min"<<endl;
    infileV.close();
}

/*
    Пример запуска:

    //ParallelTriangulation pTr = {};
    //pTr.parMap(8, 100, r0, r0+180000000, 100000, 12);
*/