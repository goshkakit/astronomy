//
//  Extrapolation.cpp
//  Project
//
//  Created by Георгий on 04.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "InfluenceForce.h"
namespace Force
{
    //==============================================================================//
    // Считывание из файла вершин
    //==============================================================================//
    void Vertice::ReadFromStr(const char* str) {
        sscanf(str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &index, &r, &theta, &fi, &x, &y, &z, &accR, &accTh, &accFi);
        //sscanf(str, "%d %lf %lf %lf %lf %lf %lf %lf", &index, &U, &r, &theta, &fi, &x, &y, &z);
    }

    void Vertice::readFromFile(ifstream& file) {
        file.read(reinterpret_cast<char*>(&index), sizeof(index));
        file.read(reinterpret_cast<char*>(&r), sizeof(r));
        file.read(reinterpret_cast<char*>(&theta), sizeof(theta));
        file.read(reinterpret_cast<char*>(&fi), sizeof(fi));
        file.read(reinterpret_cast<char*>(&x), sizeof(x));
        file.read(reinterpret_cast<char*>(&y), sizeof(y));
        file.read(reinterpret_cast<char*>(&z), sizeof(z));
        file.read(reinterpret_cast<char*>(&accR), sizeof(accR));
        file.read(reinterpret_cast<char*>(&accTh), sizeof(accTh));
        file.read(reinterpret_cast<char*>(&accFi), sizeof(accFi));
        //cout << "Reading: " << index << "|" << &index << "|" << "|" << sizeof(index) << endl;
    }
    
    //==============================================================================//
    // Считывание из файла треугольников
    //==============================================================================//
    void Tr::ReadFromStr(const char* str) {
        sscanf(str, "%d %d %d %d %d %d %d %d %d", &index, &fatherInd, &childInd[0], &childInd[1], &childInd[2], &childInd[3],
            &V[0], &V[1], &V[2]);
    }

    void Tr::readFromFile(ifstream& file) {
        file.read(reinterpret_cast<char*>(&index), sizeof(index));
        file.read(reinterpret_cast<char*>(&fatherInd), sizeof(fatherInd));
        file.read(reinterpret_cast<char*>(&childInd[0]), sizeof(childInd[0]));
        file.read(reinterpret_cast<char*>(&childInd[1]), sizeof(childInd[1]));
        file.read(reinterpret_cast<char*>(&childInd[2]), sizeof(childInd[2]));
        file.read(reinterpret_cast<char*>(&childInd[3]), sizeof(childInd[3]));
        file.read(reinterpret_cast<char*>(&V[0]), sizeof(V[0]));
        file.read(reinterpret_cast<char*>(&V[1]), sizeof(V[1]));
        file.read(reinterpret_cast<char*>(&V[2]), sizeof(V[2]));
    }

    //==============================================================================//
    // Шаг по радиусу
    //==============================================================================//
    void InfluenceForce::stepSet(int step) {
        this->stepRad = step;
    }

    //==============================================================================//
    // Диапазон сетки
    //==============================================================================//
    void InfluenceForce::radSet(int stR, int maxR) {
        this->startRad = stR;
        this->maxRad = maxR;
        this->maxStep = (maxR - stR) / this->stepRad;
    }

    //==============================================================================//
    // Загрузка данных о вершинах
    //==============================================================================//
    void InfluenceForce::verticeLoader(string filename) {
        //int c_r = 0;
        string line;
        string path = "../maps_gr/" + filename;
        ifstream inf;
        inf.open(path, ios::binary | ios::in);
        if (!inf)
            cout << "No such vertices map" << endl;
        else {
            //vert_arr.reserve(int((maxRad - startRad) / stepRad));
            //const char ch = '\n';
            //char mass[500] = {};
            
            for (int rad = 0; rad <= int((maxRad - startRad) / stepRad); rad++) {
                Vertice vertex = {};
                this->vert_arr.emplace_back(vertAm+1, vertex);
                
                for (int i = 0; i <= this->vertAm; i++) {
                    //cout << "Rad " << rad << " Vert " << i << endl;
                    Vertice& newVert = this->vert_arr[rad][i];
                    newVert.readFromFile(inf);
                    //c_r++;
                   // cout << "Cur Read: " << c_r << endl;
                }
                
                /*
                vector<Vertice> v = {};
                this->vert_arr.push_back(v);
                for (int i = 0; i <= this->vertAm; i++) {
                    Vertice newVert;
                    memset(mass, 0, 500);
                    inf.getline(mass, 499, ch);
                    newVert.ReadFromStr(mass);
                    this->vert_arr[rad].push_back(newVert);
                }
                */
            }
        }
    }

    //==============================================================================//
    // Загрузка данных о треугольниках
    //==============================================================================//
    void InfluenceForce::triangleLoader(string filename) {
        //string line;
        ifstream inf;
        inf.open("../maps_gr/" + filename, ios::binary | ios::in);
        if (!inf)
            cout << "No such triangles map" << endl;
        else {
            
            //tr_arr.reserve(trAm);
            for (int i = 0; i <= this->trAm; i++) {
                Tr tr = {};
                tr.readFromFile(inf);
                this->tr_arr.push_back(tr);
            }
            
            /*
            const char ch = '\n';
            char mass[90] = {};
            for (int i = 0; i <= this->trAm; i++) {
                Tr tr = {};
                memset(mass, 0, 90);
                inf.getline(mass, 89, ch);
                tr.ReadFromStr(mass);
                this->tr_arr.push_back(tr);
            }
            */
        }
    }

    //==============================================================================//
    // Сферическое расстояние
    //==============================================================================//
    double InfluenceForce::sphDist(double theta, double fi, int idx, int rIdx) {
        return acos(sin(theta) * sin(this->vert_arr.at(rIdx).at(idx).theta) +
            cos(theta) * cos(this->vert_arr.at(rIdx).at(idx).theta) *
            cos(fi - this->vert_arr.at(rIdx).at(idx).fi));
    }

    //==============================================================================//
    // Общая загрузка данных о сетке
    //==============================================================================//
    void InfluenceForce::loader(int pol, int deg, int maxR, int step, int stR) {
        if (stR > maxR)
            cout << "Wrong rads" << endl;
        else {
            string degree = to_string(deg);
            string st = to_string(step);
            this->stepSet(step);
            this->radSet(stR, maxR);
            
            if (deg == 4) {
                this->vertAm = 2561;
                this->trAm = 6819;
            }
            if (deg == 5) {
                this->vertAm = 10241;
                this->trAm = 27299;
            }
            if (deg == 6) {
                this->vertAm = 40961;
                this->trAm = 109219;
            }
            if (deg == 7) {
                this->vertAm = 163841;
                this->trAm = 436899;
            }
            if (deg == 8) {
                this->vertAm = 655361;
                this->trAm = 1747619;
            }
            if (deg == 9) {
                this->vertAm = 2621441;
                this->trAm = 6990499;
            }
            //if(startRad == r0)
            //    this->verticeLoader("v_A" + degree + "_100_" + st +".txt");
            //else
            this->verticeLoader("v_D" + degree + "_" + to_string(pol) + "_" + st + "_" + to_string(startRad) + "-" + to_string(maxRad) + "F.txt");
            this->triangleLoader("Triangles_" + degree + "F.txt");

            if (this->tr_arr.size() == 0) {
                cout << "No triangles file" << endl;
                exit;
            }
            if (this->vert_arr.size() == 0) {
                cout << "No vertices file" << endl;
                exit;
            }

            //cout << this->tr_arr[25].index << " " << this->tr_arr[25].childInd[0] << " " << this->tr_arr[25].childInd[1] << " " << this->tr_arr[25].childInd[2] << " " << this->tr_arr[25].childInd[3] << endl;
        }
    }

    //==============================================================================//
    // Поиск треугольника на нулевом шаге триангуляции
    //==============================================================================//
    int InfluenceForce::zeroSearcher(int rIdx, double x, double y, double z) {
        for (int i = 0; i < 20; i++) {
            if (isInTr(rIdx, i, x, y, z))
                return i;
        }
        return -1;
    }

    //==============================================================================//
    // Поиск треугольника, которому принадлежит точка
    //==============================================================================//
    int InfluenceForce::searcher(int fthrIdx, int rIdx, double x, double y, double z) {
        Tr& t = this->tr_arr.at(fthrIdx);
        if (t.childInd[0] != -1) {
            for (int j = 0; j < 4; j++) {
                if (isInTr(rIdx, t.childInd[j], x, y, z)) {
                    return searcher(t.childInd[j], rIdx, x, y, z);
                }
            }
        }
        else
            return fthrIdx;
        return fthrIdx;
    }

    int InfluenceForce::parSearcher(int fthrIdx, int rIdx, double x, double y, double z) {
        Tr& t = this->tr_arr.at(fthrIdx);
        if (t.childInd[0] != -1) {
            int ind;
#pragma omp parallel num_threads(4)
            {
                int j = omp_get_thread_num();
                if (isInTr(rIdx, t.childInd[j], x, y, z))
                    ind = j;
            }
//#pragma omp barrier
            return parSearcher(t.childInd[ind], rIdx, x, y, z);
        }
        else
            return fthrIdx;
        return fthrIdx;
    }
    
    //==============================================================================//
    // Проверка принадлежности точки  треугольнику
    //==============================================================================//
    bool InfluenceForce::isInTr(int rIdx, int idx, double x, double y, double z) {
        Tr& t = this->tr_arr[idx];
        Vertice& V0 = this->vert_arr.at(rIdx).at(t.V[0]);
        Vertice& V1 = this->vert_arr.at(rIdx).at(t.V[1]);
        Vertice& V2 = this->vert_arr.at(rIdx).at(t.V[2]);

        double delta = V0.x * (V1.y * V2.z - V2.y * V1.z) - V1.x * (V0.y * V2.z - V2.y * V0.z) + V2.x * (V0.y * V1.z - V1.y * V0.z);
        double deltaX = x * (V1.y * V2.z - V2.y * V1.z) - V1.x * (y * V2.z - V2.y * z) + V2.x * (y * V1.z - V1.y * z);
        double deltaY = V0.x * (y * V2.z - V2.y * z) - x * (V0.y * V2.z - V2.y * V0.z) + V2.x * (V0.y * z - y * V0.z);
        double deltaZ = V0.x * (V1.y * z - y * V1.z) - V1.x * (V0.y * z - y * V0.z) + x * (V0.y * V1.z - V1.y * V0.z);

        double lmbd = delta / (deltaX + deltaY + deltaZ);
        double a = deltaX / (deltaX + deltaY + deltaZ);
        double b = deltaY / (deltaX + deltaY + deltaZ);
        double c = deltaZ / (deltaX + deltaY + deltaZ);

        if (a > 0 && b > 0 && c > 0 && lmbd > 0)
            return true;
        return false;
    }

    //==============================================================================//
    // Интерполяция на слое
    //==============================================================================//
    vector<double> InfluenceForce::layerCounter(double x, double y, double z, double theta, double fi, int trIdx, int rIdx) {
        //double U = 0;
        vector<double> res(3);
        Tr& t = this->tr_arr.at(trIdx);
        double delta = 1e-10;

        double d1 = sphDist(theta, fi, t.V[0], rIdx) + delta;
        double d2 = sphDist(theta, fi, t.V[1], rIdx) + delta;
        double d3 = sphDist(theta, fi, t.V[2], rIdx) + delta;
        
       
        res[0] = (this->vert_arr.at(rIdx).at(t.V[0]).accR / d1 +
            this->vert_arr.at(rIdx).at(t.V[1]).accR / d2 +
            this->vert_arr.at(rIdx).at(t.V[2]).accR / d3) / (1 / d1 + 1 / d2 + 1 / d3);

        res[1] = (this->vert_arr.at(rIdx).at(t.V[0]).accTh / d1 +
            this->vert_arr.at(rIdx).at(t.V[1]).accTh / d2 +
            this->vert_arr.at(rIdx).at(t.V[2]).accTh / d3) / (1 / d1 + 1 / d2 + 1 / d3);

        res[2] = (this->vert_arr.at(rIdx).at(t.V[0]).accFi / d1 +
            this->vert_arr.at(rIdx).at(t.V[1]).accFi / d2 +
            this->vert_arr.at(rIdx).at(t.V[2]).accFi / d3) / (1 / d1 + 1 / d2 + 1 / d3);
        
       
            /*
            * if (this->sTr == trIdx) {
                res[0] = (this->vert_arr.at(rIdx).at(t.V[0]).accR / d1 +
                    this->vert_arr.at(rIdx).at(t.V[1]).accR / d2 +
                    this->vert_arr.at(rIdx).at(t.V[2]).accR / d3) / (1 / d1 + 1 / d2 + 1 / d3);

                res[1] = (this->vert_arr.at(rIdx).at(t.V[0]).accTh / d1 +
                    this->vert_arr.at(rIdx).at(t.V[1]).accTh / d2 +
                    this->vert_arr.at(rIdx).at(t.V[2]).accTh / d3) / (1 / d1 + 1 / d2 + 1 / d3);

                res[2] = (this->vert_arr.at(rIdx).at(t.V[0]).accFi / d1 +
                    this->vert_arr.at(rIdx).at(t.V[1]).accFi / d2 +
                    this->vert_arr.at(rIdx).at(t.V[2]).accFi / d3) / (1 / d1 + 1 / d2 + 1 / d3);
             }
            else {
                Tr& st = this->tr_arr.at(this->sTr);
                vector<int> nV;
                int co = 0;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        
                        if (this->vert_arr.at(rIdx).at(st.V[i]).index != this->vert_arr.at(rIdx).at(t.V[j]).index) {
                            co++;
                        }
                        if (co == 2) {
                            nV.push_back(this->vert_arr.at(rIdx).at(st.V[i]).index);
                            exit;
                       }
                    }
                    co = 0;
                  
                }

                if (nV.size() == 2) {
                    double d1 = sphDist(theta, fi, nV[0], rIdx) + delta;
                    double d2 = sphDist(theta, fi, nV[1], rIdx) + delta;

                    res[0] = (this->vert_arr.at(rIdx).at(nV[0]).accR / d1 +
                        this->vert_arr.at(rIdx).at(nV[1]).accR / d2) / (1 / d1 + 1 / d2);

                    res[1] = (this->vert_arr.at(rIdx).at(nV[0]).accTh / d1 +
                        this->vert_arr.at(rIdx).at(nV[1]).accTh / d2) / (1 / d1 + 1 / d2);

                    res[2] = (this->vert_arr.at(rIdx).at(nV[0]).accFi / d1 +
                        this->vert_arr.at(rIdx).at(nV[1]).accFi / d2) / (1 / d1 + 1 / d2);
                }
                //double d1s = sphDist(theta, fi, st.V[0], rIdx) + delta;
                //double d2s = sphDist(theta, fi, st.V[1], rIdx) + delta;
                //double d3s = sphDist(theta, fi, st.V[2], rIdx) + delta;
            }
            */
        
        /*
        * 
        Tr& tF = this->tr_arr.at(t.fatherInd);
        double d4 = sphDist(theta, fi, tF.V[0], rIdx) + delta;
        double d5 = sphDist(theta, fi, tF.V[1], rIdx) + delta;
        double d6 = sphDist(theta, fi, tF.V[2], rIdx) + delta;

        res[0] = (this->vert_arr.at(rIdx).at(t.V[0]).accR / d1 +
            this->vert_arr.at(rIdx).at(t.V[1]).accR / d2 +
            this->vert_arr.at(rIdx).at(t.V[2]).accR / d3 +
            this->vert_arr.at(rIdx).at(tF.V[0]).accR / d4 +
            this->vert_arr.at(rIdx).at(tF.V[1]).accR / d5 +
            this->vert_arr.at(rIdx).at(tF.V[2]).accR / d6) / (1 / d1 + 1 / d2 + 1 / d3 + 1 / d4 + 1 / d5 + 1 / d6);

        res[1] = (this->vert_arr.at(rIdx).at(t.V[0]).accTh / d1 +
            this->vert_arr.at(rIdx).at(t.V[1]).accTh / d2 +
            this->vert_arr.at(rIdx).at(t.V[2]).accTh / d3 + 
            this->vert_arr.at(rIdx).at(tF.V[0]).accTh / d4 +
            this->vert_arr.at(rIdx).at(tF.V[1]).accTh / d5 +
            this->vert_arr.at(rIdx).at(tF.V[2]).accTh / d6) / (1 / d1 + 1 / d2 + 1 / d3 + 1 / d4 + 1 / d5 + 1 / d6);

        res[2] = (this->vert_arr.at(rIdx).at(t.V[0]).accFi / d1 +
            this->vert_arr.at(rIdx).at(t.V[1]).accFi / d2 +
            this->vert_arr.at(rIdx).at(t.V[2]).accFi / d3 +
            this->vert_arr.at(rIdx).at(tF.V[0]).accFi / d4 +
            this->vert_arr.at(rIdx).at(tF.V[1]).accFi / d5 +
            this->vert_arr.at(rIdx).at(tF.V[2]).accFi / d6) / (1 / d1 + 1 / d2 + 1 / d3 + 1 / d4 + 1 / d5 + 1 / d6);
            */
        
        return res;
    }

    double InfluenceForce::funcFX(double l, double x, double y, Vertice A, Vertice B, Vertice C, Vertice D, Vertice E, Vertice F) {
        return x * x * (2 * B.accR + 2 * A.accR - 4 * D.accR) / (l * l) + 2 * x * y * (2 * E.accR - 2 * F.accR - B.accR + A.accR) / (l * l) +
            y * y * (A.accR / 2.0 + B.accR / 2.0 + 2 * C.accR + D.accR - 2 * E.accR - 2 * F.accR) / (l * l) +
            x * (4 * D.accR - B.accR - 3 * A.accR) / l + y * (-1.5 * A.accR + B.accR / 2.0 - C.accR - 2 * D.accR + 4 * F.accR) / l + A.accR;
    }

    double InfluenceForce::funcFY(double l, double x, double y, Vertice A, Vertice B, Vertice C, Vertice D, Vertice E, Vertice F) {
        return x * x * (2 * B.accTh + 2 * A.accTh - 4 * D.accTh) / (l * l) + 2 * x * y * (2 * E.accTh - 2 * F.accTh - B.accTh + A.accTh) / (l * l) +
            y * y * (A.accTh / 2.0 + B.accTh / 2.0 + 2 * C.accTh + D.accTh - 2 * E.accTh - 2 * F.accTh) / (l * l) +
            x * (4 * D.accTh - B.accTh - 3 * A.accTh) / l + y * (-1.5 * A.accTh + B.accTh / 2.0 - C.accTh - 2 * D.accTh + 4 * F.accTh) / l + A.accTh;
    }

    double InfluenceForce::funcFZ(double l, double x, double y, Vertice A, Vertice B, Vertice C, Vertice D, Vertice E, Vertice F) {
        return x * x * (2 * B.accFi + 2 * A.accFi - 4 * D.accFi) / (l * l) + 2 * x * y * (2 * E.accFi - 2 * F.accFi - B.accFi + A.accFi) / (l * l) +
            y * y * (A.accFi / 2.0 + B.accFi / 2.0 + 2 * C.accFi + D.accFi - 2 * E.accFi - 2 * F.accFi) / (l * l) +
            x * (4 * D.accFi - B.accFi - 3 * A.accFi) / l + y * (-1.5 * A.accFi + B.accFi / 2.0 - C.accFi - 2 * D.accFi + 4 * F.accFi) / l + A.accFi;
    }

    vector<double> InfluenceForce::projectOnLayer(double x, double y, double z, int rIdx) {
        vector<double> res(3);
        double r = this->startRad + rIdx * this->stepRad;
        double theta = acos(z / sqrt(x * x + y * y + z * z));
        double fi;
        if (x == 0 && y == 0)
        	fi = 0;
        else {
        	fi = atan2(y, x);
        	if (fi < 0) {
        		fi = 2 * pi + fi;
        	}
        }

        res[0] = r * sin(theta) * cos(fi);
        res[1] = r * sin(theta) * sin(fi);
        res[2] = r * cos(theta);

        return res;
    }

    vector<double> InfluenceForce::layerCounterF(double x, double y, double z, int trIdx, int rIdx) {
        vector<double> res(3);
        vector<double> xt(3);
        xt = projectOnLayer(x, y, z, rIdx);

        x = xt[0];
        y = xt[1];
        z = xt[2];
        //cout <<"Projected coordinates "<< x << " " << y << " " << z << endl;
        Tr& t = this->tr_arr.at(trIdx);

        Vertice& A = this->vert_arr.at(rIdx).at(t.V[0]);
        Vertice& B = this->vert_arr.at(rIdx).at(t.V[1]);
        Vertice& C = this->vert_arr.at(rIdx).at(t.V[2]);

        double AB = distCounter(A, B);
        double BC = distCounter(B, C);
        double CA = distCounter(C, A);
        double l = (AB + BC + CA) / 3.0;
       // this->accL += l;
       // cout << l << endl;
        double eu_x = B.x - A.x;
        double eu_y = B.y - A.y;
        double eu_z = B.z - A.z;
        if (eu_x == 0)
            eu_x = 1e-10;
        double a_x = C.x - A.x;
        double a_y = C.y - A.y;
        double a_z = C.z - A.z;

        double n_x = a_y * eu_z - a_z * eu_y;
        double n_y = a_z * eu_x - a_x * eu_z;
        double n_z = a_x * eu_y - a_y * eu_x;

        double alph = (n_z - n_x * eu_z / eu_x) / (n_x * eu_y / eu_x - n_y);
        double bet = (-alph * eu_y - eu_z) / eu_x;

        double ev_z = sqrt(3) * AB * CA / (2 * (a_x * bet + a_y * alph + a_z));
        double ev_y = alph * ev_z;
        double ev_x = bet * ev_z;
       // cout << "eu coordinates: " << eu_x << " " << eu_y << " " << eu_z << endl;
       // cout << "ev coordinates: " << ev_x << " " << ev_y << " " << ev_z << endl;
        eu_x = eu_x / AB;
        eu_y = eu_y / AB;
        eu_z = eu_z / AB;

        double ev = sqrt(ev_x * ev_x + ev_y * ev_y + ev_z * ev_z);

        ev_z = ev_z / ev;
        ev_y = ev_y / ev;
        ev_x = ev_x / ev;
        /*
        double ew_x = eu_y * ev_z - eu_z * ev_y;
        double ew_y = eu_z * ev_x - eu_x * ev_z;
        double ew_z = eu_x * ev_y - eu_y * ev_x;
        */
        double xn = (x - A.x) * eu_x + (y - A.y) * eu_y + (z - A.z) * eu_z;
        double yn = (x - A.x) * ev_x + (y - A.y) * ev_y + (z - A.z) * ev_z;

        //cout << "New coordinates " << xn << " " << yn  << endl;

        //double An_x = 0;
        //double An_y = 0;

        double Bn_x = (B.x - A.x) * eu_x + (B.y - A.y) * eu_y + (B.z - A.z) * eu_z;
        //double Bn_y = 0;

        double Cn_x = (C.x - A.x) * eu_x + (C.y - A.y) * eu_y + (C.z - A.z) * eu_z;
        double Cn_y = (C.x - A.x) * ev_x + (C.y - A.y) * ev_y + (C.z - A.z) * ev_z;

       // cout << "New B coordinates " << Bn_x << " " << Bn_y << endl;
      //  cout << "New C coordinates " << Cn_x << " " << Cn_y << endl;
        double in_x = 0;
        if (xn != Cn_x) {
            double line_slope = (yn - Cn_y) / (xn - Cn_x);
            if (line_slope == 0) {
                cout << "Point is not in triangle!" << endl;
                exit;
            }
            else {
                double line_c = yn - line_slope * xn;
                in_x = -line_c / line_slope;
            }
        }
        else {
            in_x = xn;
        }
     //   cout <<"Intersection x cord " << in_x << endl;

        double inC = sqrt((in_x - Cn_x) * (in_x - Cn_x) + Cn_y * Cn_y);
        double inX = sqrt((in_x - xn) * (in_x - xn) + yn * yn);

        double alpha = in_x / Bn_x;
        double beta = inX / inC;
      //  cout << "Alpha " << alpha << " Beta " << beta << endl;

        res[0] = (A.accR * (1 - alpha) + B.accR * alpha) * (1 - beta) + C.accR * beta;
        res[1] = (A.accTh * (1 - alpha) + B.accTh * alpha) * (1 - beta) + C.accTh * beta;
        res[2] = (A.accFi * (1 - alpha) + B.accFi * alpha) * (1 - beta) + C.accFi * beta;

        /*
        Tr& tF = this->tr_arr.at(t.fatherInd);
        Tr& tC = this->tr_arr.at(tF.childInd.at(3));
        
        Vertice& A = this->vert_arr.at(rIdx).at(tF.V[0]);
        Vertice& B = this->vert_arr.at(rIdx).at(tF.V[1]);
        Vertice& C = this->vert_arr.at(rIdx).at(tF.V[2]);

     //   cout << "FTr sides are " << distCounter(A, B) << " "<< distCounter(B, C)<<" "<< distCounter(C, A)<< endl;

        Vertice& D = this->vert_arr.at(rIdx).at(tC.V[0]);
        Vertice& E = this->vert_arr.at(rIdx).at(tC.V[1]);
        Vertice& F = this->vert_arr.at(rIdx).at(tC.V[2]);

      //  cout << "Tr sides are " << distCounter(D, E) << " " << distCounter(E, F) << " " << distCounter(F, D) << endl;
       // cout << "Distances are " << distCounter(A, D) << " " << distCounter(D, B) << " " << distCounter(B, E) << " " << distCounter(E, C)<< " " << distCounter(C, F) << " " << distCounter(F, A) << endl;
        double AB = distCounter(A, B);
        double BC = distCounter(B, C);
        double CA = distCounter(C, A);
        double l = (AB + BC + CA)/3.0;
       // cout << l << " " << l / 2.0 << " " << (distCounter(D, E) + distCounter(E, F) + distCounter(F, D)) / 3.0 << endl;
        double eu_x = B.x - A.x;
        double eu_y = B.y - A.y;
        double eu_z = B.z - A.z;

       // cout << eu_x << " " << eu_y << " " << eu_z << endl;

        double a_x = C.x - A.x;
        double a_y = C.y - A.y;
        double a_z = C.z - A.z;

        double n_x = a_y * eu_z - a_z * eu_y;
        double n_y = a_z * eu_x - a_x * eu_z;
        double n_z = a_x * eu_y - a_y * eu_x;

        double alph = (n_z - n_x * eu_z / eu_x) / (n_x * eu_y / eu_x - n_y);
        double bet = (-alph * eu_y - eu_z) / eu_x;

        double ev_z = sqrt(3) * AB * CA / (2 * (a_x * bet + a_y * alph + a_z));
        double ev_y = alph * ev_z;
        double ev_x = bet * ev_z;

        //cout << eu_x *eu_x+ eu_y *eu_y+ eu_z*eu_z << endl;
        //cout << ev_x * ev_x + ev_y * ev_y + ev_z * ev_z << endl;
        eu_x = eu_x / AB;
        eu_y = eu_y / AB;
        eu_z = eu_z / AB;

        double ev = sqrt(ev_x * ev_x + ev_y * ev_y + ev_z * ev_z);

        ev_z = ev_z / ev;
        ev_y = ev_y / ev;
        ev_x = ev_x / ev;

        double ew_x = eu_y * ev_z - eu_z * ev_y;
        double ew_y = eu_z * ev_x - eu_x * ev_z;
        double ew_z = eu_x * ev_y - eu_y * ev_x;

        //cout << ew_x * ev_x + ew_y * ev_y + ew_z * ev_z << endl;

        double xn = (x-A.x) * eu_x + (y-A.y) * eu_y + (z-A.z) * eu_z;
        double yn = (x-A.x) * ev_x + (y-A.y) * ev_y + (z-A.z) * ev_z;
        double zn = (x-A.x) * ew_x + (y-A.y) * ew_y + (z-A.z) * ew_z;
       
        //cout <<sqrt(xn*xn+yn*yn+zn*zn) << endl;

        res[0] = funcFX(l, xn, yn, A, B, C, D, E, F);
        res[1] = funcFY(l, xn, yn, A, B, C, D, E, F);
        res[2] = funcFZ(l, xn, yn, A, B, C, D, E, F);
        */
        return res;
    }

    //==============================================================================//
    // Интерполяция между слоями
    //==============================================================================//
    vector<double> InfluenceForce::counter(double x, double y, double z, double r, double theta, double fi, int tr, int rIdx) {
        vector<double> res(3), V1(3);
        double delta = 1e-10;
        V1 = layerCounterF(x, y, z, tr, rIdx);
        //V1 = layerCounter(x, y, z, theta, fi, tr, rIdx);
        double w1 = 1 / (abs(r - (startRad + (stepRad * rIdx))) + delta);
        //double U1 = V1[0];
        double accR1 = V1[0];
        double accTh1 = V1[1];
        double accFi1 = V1[2];
        /*
        if (rIdx != 0) {
            vector<double> V0(3);
            V0 = layerCounter(theta, fi, tr, rIdx - 1);
            double accR0 = V0[0];
            double accTh0 = V0[1];
            double accFi0 = V0[2];
            double w0 = 1 / (abs(r - (startRad + (stepRad * (rIdx - 1)))));
            */
            if (rIdx != maxStep) {
                vector<double> V2(3);
                V2 = layerCounterF(x, y, z, tr, rIdx + 1);
                //V2 = layerCounter(x,y,z, theta, fi, tr, rIdx + 1);
                double accR2 = V2[0];
                double accTh2 = V2[1];
                double accFi2 = V2[2];

                double w2 = 1 / (abs(startRad + (stepRad * (rIdx + 1)) - r) + delta);

                double accR = (accR1 * w1 + accR2 * w2) / (w1 + w2);
                double accTh = (accTh1 * w1 + accTh2 * w2) / (w1 + w2);
                double accFi = (accFi1 * w1 + accFi2 * w2) / (w1 + w2);
                res[0] = accR;
                res[1] = accTh;
                res[2] = accFi;
                /*
                if (rIdx + 1 != maxStep) {
                    vector<double> V3(3);
                    V3 = layerCounter(theta, fi, tr, rIdx + 2);
                    double w3 = 1 / (abs(startRad + (stepRad * (rIdx + 2)) - r));
                    double accR3 = V3[0];
                    double accTh3 = V3[1];
                    double accFi3 = V3[2];

                    double accR = (accR1 * w1 + accR2 * w2 + accR0 * w0 + accR3 * w3) / (w1 + w2 + w0 + w3);
                    double accTh = (accTh1 * w1 + accTh2 * w2 + accTh0 * w0 + accTh3 * w3) / (w1 + w2 + w0 + w3);
                    double accFi = (accFi1 * w1 + accFi2 * w2 + accFi0 * w0 + accFi3 * w3) / (w1 + w2 + w0 + w3);
                    res[0] = accR;
                    res[1] = accTh;
                    res[2] = accFi;
                }
                else {
                
                    //double U = (U1 * w1 + U2 * w2) / (w1 + w2);
                    double accR = (accR1 * w1 + accR2 * w2 + accR0 * w0) / (w1 + w2 + w0);
                    double accTh = (accTh1 * w1 + accTh2 * w2 + accTh0 * w0) / (w1 + w2 + w0);
                    double accFi = (accFi1 * w1 + accFi2 * w2 + accFi0 * w0) / (w1 + w2 + w0);
                    res[0] = accR;
                    res[1] = accTh;
                    res[2] = accFi;
                }
                */
            }
            
            else{
                //res.push_back(U1);
               // res[0] = (accR1 * w1 + accR0 * w0) / (w1 + w0);
               // res[1] = (accTh1 * w1 + accTh0 * w0) / (w1 + w0);
               // res[2] = (accFi1 * w1 + accFi0 * w0) / (w1 + w0);
                res[0] = accR1;
                res[1] = accTh1;
                res[2] = accFi1;
            }
            /*
            if (this->cc == 0 || this->cc == 1) {
                return res;
            }
            else {
                vector<double> intF(3), bF(3);
                if (curX[0] == prX[0]) {
                    intF[0] = curF[0];
                }
                else {
                    intF[0] = (curF[0] * (x * 1e-6 - prX[0]) - prF[0] * (x * 1e-6 - curX[0])) / (curX[0] - prX[0]);
                }
                if (curX[1] == prX[1]) {
                    intF[1] = curF[1];
                }
                else {
                    intF[1] = (curF[1] * (y * 1e-6 - prX[1]) - prF[1] * (y * 1e-6 - curX[1])) / (curX[1] - prX[1]);
                }
                if (curX[2] == prX[2]) {
                    intF[2] = curF[2];
                }
                else {
                    intF[2] = (curF[2] * (z * 1e-6 - prX[2]) - prF[2] * (z * 1e-6 - curX[2])) / (curX[2] - prX[2]);
                }
                //cout << intF[0] << " " << intF[1] << " " << intF[2] << endl; 
               // cout << res[0] << " " << res[1] << " " << res[2] << endl;
                
                if (sTr != tr) {
                    alpha = 0.0;
                }
                else {
                    if (alpha <= 0.9) {
                        alpha += 0.4;
                    }
                    if (alpha > 1) {
                        alpha = 1;
                    }
                }
                //cout << alpha << endl;
                bF[0] = alpha * res[0] + (1.0 - alpha) * intF[0];
                bF[1] = alpha * res[1] + (1.0 - alpha) * intF[1];
                bF[2] = alpha * res[2] + (1.0 - alpha) * intF[2];
                
                return bF;
                
               // return res;
            }
            */
            return res;
    }

    //==============================================================================//
    // Вычисление ускорения
    //==============================================================================//
    vector<double> InfluenceForce::extrapolator(double x, double y, double z) {
       // auto fs = std::chrono::high_resolution_clock::now();

        //auto s = std::chrono::high_resolution_clock::now();
        vector<double> V(3);
        double r = 1e6 * sqrt(x * x +
                              y * y +
                              z * z);

        double theta = acos(z / sqrt(x * x +
                                     y * y +
                                     z * z));
        double fi;
        if (x == 0 && y == 0)
            fi = 0;
        else {
            fi = atan2(y, x);
            if (fi < 0) {
                fi = 2 * pi + fi;
            }
        }
        x *= 1e6;
        y *= 1e6;
        z *= 1e6;
        
        //auto e = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> time = e - s;
        //cout << "Coord translation time: " << time.count() << endl;
        
       // theta -= pi / 2;
      //  fi += pi;
        //cout << "Input spherical coordinates: " << r << " " << theta << " " << fi << endl;
        
       // cout << "Transcripted cartesian coordinates: " << x << " " << y << " " << z << endl;
       // cout << "Step " << cc << endl;
       // cout << "Input cartesian coordinates: " << x << " " << y << " " << z << endl;
       // cout << "Mod spherical coordinates: " << r << " " << theta << " " << fi << endl;
       // theta += pi / 2;
        //For mapping uncomment:
        //fi -= pi;
        //cout << r << " " << theta << " " << fi << endl;
        int rIdx = max(int(floor((r - startRad) / stepRad)), 0);
       // cout << "R ind: "<< rIdx << endl;
       
       // s = std::chrono::high_resolution_clock::now();

        int zeroIdx = zeroSearcher(rIdx, x, y, z);
        
        // e = std::chrono::high_resolution_clock::now();
        //time = e - s;
        //cout << "Zero triangle search time: " << time.count() << endl;
        
        //cout << "Zero ind: " << zeroIdx << endl;
        //s = std::chrono::high_resolution_clock::now();

        int tr = searcher(zeroIdx, rIdx, x, y, z);
       // int tr = parSearcher(zeroIdx, rIdx, x, y, z);
        
         //e = std::chrono::high_resolution_clock::now();
        //time = e - s;
        //cout << "Final triangle search time: " << time.count() << endl;
        
        //this->curTr = tr;
        //cout << tr << endl;
        //V=layerCounterF(x, y, z, tr, rIdx);
        //cout << V[0] << " " << V[1] << " " << V[2] << endl;
        //V = layerCounter(theta, fi, tr, rIdx);
        //cout << V[0] << " " << V[1] << " " << V[2] << endl;
        //s = std::chrono::high_resolution_clock::now();

        V = counter(x, y, z, r, theta, fi, tr, rIdx);
        
        // e = std::chrono::high_resolution_clock::now();
        //time = e - s;
        //cout << "Counting time: " << time.count() << endl;
        
      //  cout << "\n";
        //cout << V[0] << " " << V[1] << " " << V[2] << endl;
        //cout << sin(theta) * cos(fi) * V[0] + cos(theta) * cos(fi) * V[1] - sin(fi)*V[2]<<" " << sin(theta) * sin(fi) * V[0] + cos(theta) * sin(fi) * V[1] + cos(fi) * V[2]
         //    << " " << cos(theta) * V[0] - sin(theta) * V[1] << endl;
        /*
        PotentialCounter* pc = new PotentialCounter;
        //V[0] += pc->U0(r) + pc->U2(r, theta, fi);
        V[0] += pc->accR0(r, theta, fi);
        V[1] += pc->accTh0(r, theta, fi);
        V[2] += pc->accFi0(r, theta, fi);
        delete pc;
        */
        
        //auto fe = std::chrono::high_resolution_clock::now();

        //time = fe - fs;
        //cout << "Full time: " << time.count() << endl;
        
        return V;
    }
}
/*
    Пример запуска:

    Extrapolation ex = {};
    ex.loader(100, 8, r0+1001000, 1000, r0+1000000);
    vector<double> eex = ex.extrapolator(r, theta, fi);
*/