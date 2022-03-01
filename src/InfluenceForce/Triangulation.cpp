//
//  Triangulation.cpp
//  Project
//
//  Created by Георгий on 09.12.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#include "InfluenceForce.h"
namespace Force
{
    //==============================================================================//
    // Радиус слоя
    //==============================================================================//
    void InfluenceForce::setRad(double r)
    {
        this->r = r;
    }

    //==============================================================================//
    // Преобразование декартовых координат в сферические
    //==============================================================================//
    void VerticeWithCompared::toSph()
    {
        this->r = sqrt(this->x * this->x +
            this->y * this->y +
            this->z * this->z);

        this->theta = acos(this->z / sqrt(this->x * this->x +
            this->y * this->y +
            this->z * this->z));
        if (this->x == 0 && this->y == 0)
            this->fi = 0;
        else {
            this->fi = atan2(this->y, this->x);

            if (this->fi < 0) {
                this->fi = 2 * pi + this->fi;
            }
        }
    }

    //==============================================================================//
    // Сравнение вершин по расстоянию
    //==============================================================================//
    bool VerticeWithCompared::operator<(const VerticeWithCompared& rht) const {
        return InfluenceForce::distCounter(*this, *compared) < InfluenceForce::distCounter(rht, *compared);
    }

    //==============================================================================//
    // Создание шаблона треугольника без вершин
    //==============================================================================//
    Triangle::Triangle() {
        index = 0;
        fatherInd = 0;
        V.clear();
    }

    //==============================================================================//
    // Создание шаблона треугольника с вершинами (по адресам)
    //==============================================================================//
    Triangle::Triangle(int index, VerticeWithCompared* V1, VerticeWithCompared* V2, VerticeWithCompared* V3, int fatherIndex) {
        fatherInd = fatherIndex;
        this->index = index;
        V.push_back(V1->index);
        V.push_back(V2->index);
        V.push_back(V3->index);
    }

    //==============================================================================//
    // Создание шаблона треугольника с вершинами
    //==============================================================================//
    Triangle::Triangle(int index, int V1, int V2, int V3, int fatherIndex) {
        fatherInd = fatherIndex;
        this->index = index;
        V.push_back(V1);
        V.push_back(V2);
        V.push_back(V3);
    }

    //==============================================================================//
    // Операция сравнения для треугольников
    //==============================================================================//
    bool Triangle::operator==(const Triangle& tr) const {
        bool eq = false;
        for (int i = 0; i < this->V.size(); i++) {
            for (int j = 0; j < tr.V.size(); j++) {
                if (this->V.at(i) == tr.V.at(j)) {
                    eq = true;
                    break;
                }
                else {
                    eq = false;
                }
            }
            if (eq == false) {
                break;
            }
        }
        return eq;
    }

    //==============================================================================//
    // Декартово расстояние между вершинами
    //==============================================================================//
    double InfluenceForce::distCounter(VerticeWithCompared V1, VerticeWithCompared V2) {
        double dist = sqrt((V1.x - V2.x) * (V1.x - V2.x) +
            (V1.y - V2.y) * (V1.y - V2.y) +
            (V1.z - V2.z) * (V1.z - V2.z));
        return dist;
    }

    double InfluenceForce::distCounter(Vertice V1, Vertice V2) {
        double dist = sqrt((V1.x - V2.x) * (V1.x - V2.x) +
            (V1.y - V2.y) * (V1.y - V2.y) +
            (V1.z - V2.z) * (V1.z - V2.z));
        return dist;
    }

    double InfluenceForce::distCounter(VerticeA V1, VerticeA V2) {
        double dist = sqrt((V1.x - V2.x) * (V1.x - V2.x) +
            (V1.y - V2.y) * (V1.y - V2.y) +
            (V1.z - V2.z) * (V1.z - V2.z));
        return dist;
    }

    //==============================================================================//
    // Нулевые треугольники
    //==============================================================================//
    void InfluenceForce::zero_triangles() {
        this->tr_arr_tr.clear();
        //1-4
        for (int i = 0; i < 4; i++) {
            Triangle tr = Triangle(i, &vert_arr_tr.at(0), &vert_arr_tr.at(i + 1), &vert_arr_tr.at(i + 2));
            this->tr_arr_tr.push_back(tr);
        }
        //5
        Triangle tr5 = Triangle(4, &vert_arr_tr.at(0), &vert_arr_tr.at(5), &vert_arr_tr.at(1));
        this->tr_arr_tr.push_back(tr5);
        //6,8,10,12
        for (int i = 1; i < 5; i++) {
            Triangle tr = Triangle(i + 4, &vert_arr_tr.at(i), &vert_arr_tr.at(i + 1), &vert_arr_tr.at(i + 5));
            this->tr_arr_tr.push_back(tr);
        }
        //7,9,11,13
        for (int i = 2; i < 6; i++) {
            Triangle tr = Triangle(i + 7, &vert_arr_tr.at(i), &vert_arr_tr.at(i + 4), &vert_arr_tr.at(i + 5));
            this->tr_arr_tr.push_back(tr);
        }
        //14-15
        for (int i = 5; i <= 6; i++) {
            Triangle tr = Triangle(i + 8, &vert_arr_tr.at(i), &vert_arr_tr.at(10), &vert_arr_tr.at(1));
            this->tr_arr_tr.push_back(tr);
        }
        //16-19
        for (int i = 6; i < 10; i++) {
            Triangle tr = Triangle(i + 9, &vert_arr_tr.at(11), &vert_arr_tr.at(i), &vert_arr_tr.at(i + 1));
            this->tr_arr_tr.push_back(tr);
        }
        //20
        Triangle tr20 = Triangle(19, &vert_arr_tr.at(11), &vert_arr_tr.at(10), &vert_arr_tr.at(6));
        this->tr_arr_tr.push_back(tr20);

        this->globalIndex = 19;
    }

    //==============================================================================//
    // Нулевая сетка
    //==============================================================================//
    void InfluenceForce::zeroMeshCreator(double r)
    {
        this->setRad(r);
        const float H_ANGLE = pi / 180 * 72;
        const float V_ANGLE = atan(1.0 / 2);
        this->vert_arr_tr.clear();
        VerticeWithCompared tre = {};
        this->vert_arr_tr.resize(12, tre);
        float z, xy;
        float hAngle1 = -pi / 2 - H_ANGLE / 2;
        float hAngle2 = -pi / 2;

        this->vert_arr_tr[0].x = 0;
        this->vert_arr_tr[0].y = 0;
        this->vert_arr_tr[0].z = this->r;
        this->vert_arr_tr[0].toSph();
        this->vert_arr_tr[0].index = 0;

        for (int i = 1; i <= 5; ++i)
        {
            z = this->r * sin(V_ANGLE);
            xy = this->r * cos(V_ANGLE);

            this->vert_arr_tr[i].x = xy * cos(hAngle1);

            this->vert_arr_tr[i + 5].x = xy * cos(hAngle2);

            this->vert_arr_tr[i].y = xy * sin(hAngle1);

            this->vert_arr_tr[i + 5].y = xy * sin(hAngle2);

            this->vert_arr_tr[i].z = z;

            this->vert_arr_tr[i + 5].z = -z;

            this->vert_arr_tr.at(i).toSph();
            this->vert_arr_tr[i + 5].toSph();
            this->vert_arr_tr[i].index = i;
            this->vert_arr_tr[i + 5].index = i + 5;

            hAngle1 += H_ANGLE;
            hAngle2 += H_ANGLE;
        }

        this->vert_arr_tr[11].x = 0;
        this->vert_arr_tr[11].y = 0;
        this->vert_arr_tr[11].z = -1.0 * this->r;
        this->vert_arr_tr[11].toSph();
        this->vert_arr_tr[11].index = 11;
        this->vIndex = 12;
        this->zero_triangles();
    }

    //==============================================================================//
    // Создание вершины
    //==============================================================================//
    VerticeWithCompared InfluenceForce::verticeCreator(VerticeWithCompared V1, VerticeWithCompared V2) {
        VerticeWithCompared newVertice;

        newVertice.x = (V1.x + V2.x) / 2.0;
        newVertice.y = (V1.y + V2.y) / 2.0;
        newVertice.z = (V1.z + V2.z) / 2.0;

        double norm = this->r / sqrt(newVertice.x * newVertice.x +
            newVertice.y * newVertice.y +
            newVertice.z * newVertice.z);
        newVertice.x *= norm;
        newVertice.y *= norm;
        newVertice.z *= norm;

        newVertice.toSph();

        return newVertice;
    }

    //==============================================================================//
    // Сравнение вершин
    //==============================================================================//
    bool InfluenceForce::vertDetector(VerticeWithCompared V) {
        for (int i = 0; i < this->vert_arr_tr.size(); i++) {
            if (this->vert_arr_tr.at(i).x == V.x &&
                this->vert_arr_tr.at(i).y == V.y &&
                this->vert_arr_tr.at(i).z == V.z) {
                this->localVIndex = i;
                return false;
            }
        }
        return true;
    }

    //==============================================================================//
    // Построение слоя сетки
    //==============================================================================//
    void InfluenceForce::mesher(double r, int degree) {
        this->zeroMeshCreator(r);
        vector<vector<Triangle>> all_triangles;
        all_triangles.push_back(tr_arr_tr);
        tr_arr_tr.clear();
        for (int iter = 1; iter <= degree; iter++)
        {
            all_triangles.push_back(vector<Triangle>());
            for (int tr_num = 0; tr_num < all_triangles.at(iter - 1).size(); tr_num++)
            {
                Triangle& t = all_triangles.at(iter - 1).at(tr_num);
                VerticeWithCompared Vs[3] = {};
                int Vptrs[3] = {};
                for (int idx = 0; idx < 3; idx++) {
                    Vs[idx] = verticeCreator(this->vert_arr_tr.at(t.V.at(idx % 3)), this->vert_arr_tr.at(t.V.at((idx + 1) % 3)));
                    if (vertDetector(Vs[idx])) {
                        Vs[idx].index = this->vIndex;
                        this->vIndex += 1;
                        this->vert_arr_tr.push_back(Vs[idx]);
                        Vptrs[idx] = Vs[idx].index;
                    }
                    else {
                        Vptrs[idx] = localVIndex;
                    }
                }

                Triangle tr1 = Triangle(this->globalIndex + 1, t.V[0], Vptrs[0], Vptrs[2], t.index);

                Triangle tr2 = Triangle(this->globalIndex + 2, t.V[1], Vptrs[0], Vptrs[1], t.index);

                Triangle tr3 = Triangle(this->globalIndex + 3, t.V[2], Vptrs[1], Vptrs[2], t.index);

                Triangle tr4 = Triangle(this->globalIndex + 4, Vptrs[0], Vptrs[1], Vptrs[2], t.index);

                t.childInd.clear();
                t.childInd.push_back(tr1.index);
                t.childInd.push_back(tr2.index);
                t.childInd.push_back(tr3.index);
                t.childInd.push_back(tr4.index);
                all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr1);
                all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr2);
                all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr3);
                all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr4);
                this->globalIndex = this->globalIndex + 4;
            }
            tr_arr_tr.insert(tr_arr_tr.end(), all_triangles.at(iter - 1).begin(), all_triangles.at(iter - 1).end());
            if (iter == degree) {
                tr_arr_tr.insert(tr_arr_tr.end(), all_triangles.at(iter).begin(), all_triangles.at(iter).end());
            }
        }
    }

    void Triangle::writeToFile(ofstream& file) {
        file.write(reinterpret_cast<const char*>(&index), sizeof(index));
        file.write(reinterpret_cast<const char*>(&fatherInd), sizeof(fatherInd));
        file.write(reinterpret_cast<const char*>(&childInd[0]), sizeof(childInd[0]));
        file.write(reinterpret_cast<const char*>(&childInd[1]), sizeof(childInd[1]));
        file.write(reinterpret_cast<const char*>(&childInd[2]), sizeof(childInd[2]));
        file.write(reinterpret_cast<const char*>(&childInd[3]), sizeof(childInd[3]));
        file.write(reinterpret_cast<const char*>(&V[0]), sizeof(V[0]));
        file.write(reinterpret_cast<const char*>(&V[1]), sizeof(V[1]));
        file.write(reinterpret_cast<const char*>(&V[2]), sizeof(V[2]));
    }

    //==============================================================================//
    // Создание файла треугольников
    //==============================================================================//
    void InfluenceForce::trianglesCreator(int degree) {
        string deg = to_string(degree);
        string fileName = "Triangles_" + deg;
        ofstream infile;
        infile.open("../maps_gr/" + fileName + "F.txt", std::ios_base::out | ios::binary | ios::trunc);
        cout << "Counting triangles with degree " << degree << endl;
        this->mesher(r0, degree);

        for (int j = 0; j < this->tr_arr_tr.size(); j++)
        {
            
            for (int i = 0; i < 3; i++) {
                this->tr_arr_tr[j].V[i] = this->vert_arr_tr.at(this->tr_arr_tr[j].V[i]).index;
            }
            this->tr_arr_tr[j].writeToFile(infile);
            

            /*
            infile << this->tr_arr_tr[j].index << " " << this->tr_arr_tr[j].fatherInd << " ";
            for (int k = 0; k < 4; k++) {
                infile << this->tr_arr_tr[j].childInd[k] << " ";
            }
            for (int i = 0; i < 3; i++) {
                infile << this->vert_arr_tr.at(this->tr_arr_tr[j].V[i]).index << " ";
            }
            infile << "\n";
            */
        }

       // infile << endl;
        infile.close();

        cout << "Triangles are done." << endl;
    }

    void VerticeWithCompared::writeToFile(ofstream& file) {
        file.write(reinterpret_cast<char*>(&index), sizeof(index));
        file.write(reinterpret_cast<char*>(&r), sizeof(r));
        file.write(reinterpret_cast<char*>(&theta), sizeof(theta));
        file.write(reinterpret_cast<char*>(&fi), sizeof(fi));
        file.write(reinterpret_cast<char*>(&x), sizeof(x));
        file.write(reinterpret_cast<char*>(&y), sizeof(y));
        file.write(reinterpret_cast<char*>(&z), sizeof(z));
        file.write(reinterpret_cast<char*>(&accR), sizeof(accR));
        file.write(reinterpret_cast<char*>(&accTh), sizeof(accTh));
        file.write(reinterpret_cast<char*>(&accFi), sizeof(accFi));
    }

    //==============================================================================//
    // Создание файла вершин
    //==============================================================================//
    void InfluenceForce::map(int degree, int polynomDegree, int startRad, int maxRad, int step) {
        //PotentialCounter pc = {};
        //pc.length = polynomDegree;
        //pc.LoadFromFile("../maps_gr/data.txt");
        //int c_w = 0;
        string deg = to_string(degree);
        string polD = to_string(polynomDegree);
        string st = to_string(step);
        string stR = to_string(startRad);
        string mxR = to_string(maxRad);
        ofstream infileV;
        
        infileV.open("../maps_gr/v_D" + deg + "_" + polD + "_" + st + "_" + stR + "-" + mxR + "F.txt", std::ios_base::out | ios::binary | ios::trunc);
        //cout << "Map in progress" << endl;
        time_t start, end;
        time(&start);

        // Если необходим файл с треугольниками, убрать комментирование
        //this->trianglesCreator(degree);
        
        for (int r = startRad; r <= maxRad; r += step) {
            time_t oneMapSt, oneMapF;
            time(&oneMapSt);
            this->mesher(r, degree);
            //cout << "Counting potentials on rad " << r << endl;

            for (int i = 0; i < this->vert_arr_tr.size(); i++) {
               // cout << vert_arr_tr.size() << endl;
                double x_m[3];
                double f_m[3];
                x_m[0] = this->vert_arr_tr[i].x/1e6;
                x_m[1] = this->vert_arr_tr[i].y/1e6;
                x_m[2] = this->vert_arr_tr[i].z/1e6;

                GetF_Harm_egm96(x_m, polynomDegree, f_m);

                this->vert_arr_tr[i].accR = f_m[0];
                this->vert_arr_tr[i].accTh = f_m[1];
                this->vert_arr_tr[i].accFi = f_m[2];

                this->vert_arr_tr[i].writeToFile(infileV);
                //c_w++;
                /*
                infileV << this->vert_arr_tr[i].index << " " <<
                    this->vert_arr_tr[i].r << " " <<
                    this->vert_arr_tr[i].theta << " " <<
                    this->vert_arr_tr[i].fi << " " <<
                    //this->vert_arr_tr[i].theta - pi / 2 << " " <<
                   // this->vert_arr_tr[i].fi + pi << " " <<
                    this->vert_arr_tr[i].x << " " <<
                    this->vert_arr_tr[i].y << " " <<
                    this->vert_arr_tr[i].z << " " <<
                    this->vert_arr_tr[i].accR << " " <<
                    this->vert_arr_tr[i].accTh << " " <<
                    this->vert_arr_tr[i].accFi << "\n";
                    */
            }
            time(&oneMapF);
                //cout << "Map with rad " << r << " is done in " << difftime(oneMapF, oneMapSt) / 60.0 << " min" << endl;
        }
        //cout << "Writes: " << c_w << endl;
        infileV.close();
        time(&end);
        cout << "Global map "<< degree << " is done in " << difftime(end, start) / 60.0 << " min" << endl;
           
    }


    void VerticeWithCompared::writeToFileA(ofstream& file) {
        file.write(reinterpret_cast<char*>(&index), sizeof(index));
        file.write(reinterpret_cast<char*>(&r), sizeof(r));
        file.write(reinterpret_cast<char*>(&theta), sizeof(theta));
        file.write(reinterpret_cast<char*>(&fi), sizeof(fi));
        file.write(reinterpret_cast<char*>(&x), sizeof(x));
        file.write(reinterpret_cast<char*>(&y), sizeof(y));
        file.write(reinterpret_cast<char*>(&z), sizeof(z));
        file.write(reinterpret_cast<char*>(&U), sizeof(U));
    }


    void InfluenceForce::mapAtm(int degree, double t, int data, int startRad, int maxRad, int step) {
        
        string deg = to_string(degree);
        string dat = to_string(data);
        string tt = to_string(t);
        string st = to_string(step);
        string stR = to_string(startRad);
        string mxR = to_string(maxRad);
        ofstream infileV;

        infileV.open("../maps_gr/v_A" + deg + "_" + tt + "_" + dat + "_" + st + "_" + stR + "-" + mxR + "F.txt", std::ios_base::out | ios::binary | ios::trunc);
        cout << "Map in progress" << endl;
        time_t start, end;
        time(&start);

        // Если необходим файл с треугольниками, убрать комментирование
        //this->trianglesCreator(degree);

        double ajd;
        double delt;
        double tim;
        
        jddeltt(data, t, &ajd, &delt, &tim);

        double F107 = 100;
        double F81 = 100;
        double aKp = 3;

        for (int i = 0; i < ListAtmIndex.size(); i++)
        {
            if (ListAtmIndex[i].data == data)
            {
                F107 = ListAtmIndex[i].F107;
                F81 = ListAtmIndex[i].F81;
                aKp = ListAtmIndex[i].aKp;
                break;
            }
        }

        for (int r = startRad; r <= maxRad; r += step) {
            time_t oneMapSt, oneMapF;
            time(&oneMapSt);
            this->mesher(r, degree);

            for (int i = 0; i < this->vert_arr_tr.size(); i++) {
                double x_m[3];
                
                x_m[0] = this->vert_arr_tr[i].x / 1e6;
                x_m[1] = this->vert_arr_tr[i].y / 1e6;
                x_m[2] = this->vert_arr_tr[i].z / 1e6;

                double f = Roa2004_2(tim, x_m, ajd, delt, F107, F81, aKp);

                this->vert_arr_tr[i].U = f;

                this->vert_arr_tr[i].writeToFileA(infileV);
                /*
                infileV << this->vert_arr_tr[i].index << " " <<
                    this->vert_arr_tr[i].r << " " <<
                    this->vert_arr_tr[i].theta << " " <<
                    this->vert_arr_tr[i].fi << " " <<
                    //this->vert_arr_tr[i].theta - pi / 2 << " " <<
                   // this->vert_arr_tr[i].fi + pi << " " <<
                    this->vert_arr_tr[i].x << " " <<
                    this->vert_arr_tr[i].y << " " <<
                    this->vert_arr_tr[i].z << " " <<
                    this->vert_arr_tr[i].U << "\n";
                    */
            }
            time(&oneMapF);
            cout << "Map with rad " << r << " is done in " << difftime(oneMapF, oneMapSt) / 60.0 << " min" << endl;
        }
        infileV.close();
        time(&end);
        cout << "Global map is done in " << difftime(end, start) / 60.0 << " min" << endl;

    }
    
}
/*
    Пример запуска:

    Triangulation tr = {};
    tr.map(6, 100, r0+20000000, r0+20001000, 1000);
*/
