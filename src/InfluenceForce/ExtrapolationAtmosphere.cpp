//
//  ExtrapolationAtmosphere.cpp
//  Project
//
//  Created by Георгий on 15.09.2021.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "InfluenceForce.h"
namespace Force
{
	void VerticeA::ReadFromStr(const char* str)
	{
		sscanf(str, "%d %lf %lf %lf %lf %lf %lf %lf", &index, &r, &theta, &fi, &x, &y, &z, &p);
	}

	void VerticeA::readFromFile(ifstream& file) {
		file.read(reinterpret_cast<char*>(&index), sizeof(index));
		file.read(reinterpret_cast<char*>(&r), sizeof(r));
		file.read(reinterpret_cast<char*>(&theta), sizeof(theta));
		file.read(reinterpret_cast<char*>(&fi), sizeof(fi));
		file.read(reinterpret_cast<char*>(&x), sizeof(x));
		file.read(reinterpret_cast<char*>(&y), sizeof(y));
		file.read(reinterpret_cast<char*>(&z), sizeof(z));
		file.read(reinterpret_cast<char*>(&p), sizeof(p));
	}

	void InfluenceForce::layer_num_set(double r_s, double r_e, double step) {
		this->layer_num = (r_e - r_s) / step;
	}

	void InfluenceForce::map_param_set(int deg)
	{
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
	}

	void InfluenceForce::single_am_loader(int data, int deg)
	{
		this->map_buff.clear();
		string degree = to_string(deg);
	
		string path = "../maps_gr/v_A" + degree + "_0.000000_" + to_string(data) + "_20000_6498136-7878136F.txt";
		ifstream inf;
		inf.open(path, ios::binary | ios::in);
		if (!inf)
			cout << "No such vertices map" << endl;
		else {
			//const char ch = '\n';
			//char mass[500] = {};


			for (int rad = 0; rad <= this->layer_num; rad++) {
				VerticeA vertex = {};
				this->map_buff.emplace_back(vertAm + 1, vertex);

				for (int i = 0; i <= this->vertAm; i++) {
					VerticeA& newVert = this->map_buff[rad][i];
					newVert.readFromFile(inf);
				}
				/*
				vector<VerticeA> v = {};
				this->map_buff.push_back(v);
				for (int i = 0; i <= this->vertAm; i++) {
					VerticeA newVert;
					memset(mass, 0, 500);
					in.getline(mass, 499, ch);
					newVert.ReadFromStr(mass);
					this->map_buff[rad].push_back(newVert);
				}
				*/
			}
		}

		if (this->map_buff.size() == 0) {
			cout << "No vertices file" << endl;
			exit;
		}
	}

	bool InfluenceForce::isInTrA(int rIdx, int idx, double x, double y, double z)
	{
		Tr& t = this->tr_arr[idx];
		VerticeA& V0 = this->cur_map->at(rIdx).at(t.V[0]);
		VerticeA& V1 = this->cur_map->at(rIdx).at(t.V[1]);
		VerticeA& V2 = this->cur_map->at(rIdx).at(t.V[2]);

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

	int Force::InfluenceForce::zeroSearcherA(int rIdx, double x, double y, double z)
	{
		for (int i = 0; i < 20; i++) {
			if (isInTrA(rIdx, i, x, y, z))
				return i;
		}
		return -1;
	}

	int InfluenceForce::searcherA(int fthrIdx, int rIdx, double x, double y, double z) {
		Tr& t = this->tr_arr.at(fthrIdx);
		
		if (t.childInd[0] != -1) {
			
			for (int j = 0; j < 4; j++) {
				if (isInTrA(rIdx, t.childInd[j], x, y, z)) {
					
					return searcherA(t.childInd[j], rIdx, x, y, z);
				}
			}
		}
		else
			return fthrIdx;
		return fthrIdx;
	}

	double InfluenceForce::layerCounterA(double theta, double fi, int trIdx, int rIdx)
	{
		double res;;
		Tr& t = this->tr_arr.at(trIdx);
		double delta = 1e-10;

		double d1 = sphDist(theta, fi, t.V[0], rIdx) + delta;
		double d2 = sphDist(theta, fi, t.V[1], rIdx) + delta;
		double d3 = sphDist(theta, fi, t.V[2], rIdx) + delta;

		res = (this->cur_map->at(rIdx).at(t.V[0]).p / d1 +
			   this->cur_map->at(rIdx).at(t.V[1]).p / d2 +
			   this->cur_map->at(rIdx).at(t.V[2]).p / d3) / (1 / d1 + 1 / d2 + 1 / d3);
		return res;
	}

	double InfluenceForce::funcFA(double l, double x, double y, VerticeA A, VerticeA B, VerticeA C, VerticeA D, VerticeA E, VerticeA F) {
		return x * x * (2 * B.p + 2 * A.p - 4 * D.p) / (l * l) + 2 * x * y * (2 * E.p - 2 * F.p - B.p + A.p) / (l * l) +
			   y * y * (A.p / 2.0 + B.p / 2.0 + 2 * C.p + D.p - 2 * E.p - 2 * F.p) / (l * l) +
			   x * (4 * D.p - B.p - 3 * A.p) / l + y * (-1.5 * A.p + B.p / 2.0 - C.p - 2 * D.p + 4 * F.p) / l + A.p;
	}

	double InfluenceForce::layerCounterFA(double x, double y, double z, int trIdx, int rIdx) {
		double p;
		vector<double> xt(3);
		xt = projectOnLayer(x, y, z, rIdx);

		x = xt[0];
		y = xt[1];
		z = xt[2];
		Tr& t = this->tr_arr.at(trIdx);

		VerticeA& A = this->cur_map->at(rIdx).at(t.V[0]);
		VerticeA& B = this->cur_map->at(rIdx).at(t.V[1]);
		VerticeA& C = this->cur_map->at(rIdx).at(t.V[2]);

		double AB = distCounter(A, B);
		double BC = distCounter(B, C);
		double CA = distCounter(C, A);
		double l = (AB + BC + CA) / 3.0;

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

		eu_x = eu_x / AB;
		eu_y = eu_y / AB;
		eu_z = eu_z / AB;

		double ev = sqrt(ev_x * ev_x + ev_y * ev_y + ev_z * ev_z);

		ev_z = ev_z / ev;
		ev_y = ev_y / ev;
		ev_x = ev_x / ev;

		double xn = (x - A.x) * eu_x + (y - A.y) * eu_y + (z - A.z) * eu_z;
		double yn = (x - A.x) * ev_x + (y - A.y) * ev_y + (z - A.z) * ev_z;

		double Bn_x = (B.x - A.x) * eu_x + (B.y - A.y) * eu_y + (B.z - A.z) * eu_z;

		double Cn_x = (C.x - A.x) * eu_x + (C.y - A.y) * eu_y + (C.z - A.z) * eu_z;
		double Cn_y = (C.x - A.x) * ev_x + (C.y - A.y) * ev_y + (C.z - A.z) * ev_z;

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

		double inC = sqrt((in_x - Cn_x) * (in_x - Cn_x) + Cn_y * Cn_y);
		double inX = sqrt((in_x - xn) * (in_x - xn) + yn * yn);

		double alpha = in_x / Bn_x;
		double beta = inX / inC;

		p = (A.p * (1 - alpha) + B.p * alpha) * (1 - beta) + C.p * beta;

		return p;
	}

	double InfluenceForce::counterA(double x, double y, double z, double r, double theta, double fi, int tr, int rIdx)
	{
		double p1, p;
		double delta = 1e-10;
		//p1 = layerCounterA(theta, fi, tr, rIdx);
		p1 = layerCounterFA(x, y, z, tr, rIdx);
		double w1 = 1 / (abs(r - (6498136 + (20000 * rIdx))) + delta);
		if (rIdx != this->layer_num) {
			double p2 = layerCounterFA(x, y, z, tr, rIdx + 1);
			double w2 = 1 / (abs(6498136 + (20000 * (rIdx + 1)) - r) + delta);

			p = (p1 * w1 + p2 * w2) / (w1 + w2);
		}
		else p = p1;

		return p;
	}

	vector<double> InfluenceForce::timeMove(double x, double y, double z, double time) {
		vector<double> res(6);
		double r = sqrt(x * x + y * y + z * z);

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

		double pace = 2 * pi / 86400;
		fi -= time * pace;
		if (fi < 0)
			fi = 2 * pi + fi;

		res[0] = r * sin(theta) * cos(fi);
		res[1] = r * sin(theta) * sin(fi);
		res[2] = r * cos(theta);
		res[3] = r;
		res[4] = theta;
		res[5] = fi;

		return res;
	}

	double InfluenceForce::extrapolatorA(double* vec, double time) {
		double x = vec[0]*1e6;
		double y = vec[1]*1e6;
		double z = vec[2]*1e6;
		
		vector<double> n_vec = timeMove(x, y, z, time);
	
		int rIdx = int(floor((n_vec[3] - 6498136) / 20000));

		if (n_vec[3] >= 6378136 + 1500000){//rIdx > this->layer_num) {
			cout << "Rad: " << n_vec[3] << " over 1500km" << endl;
			return 0;
		}

		if (rIdx < 0) {
			cout << "Rad: " << n_vec[3] << " under 120km" << endl;
			double h0, ast0, akst1, akst2;
			if (n_vec[3] < 0.0) {
				cout << "Collision with Earth" << endl;
				exit;
			}
			else if ((n_vec[3] > 0.0) && (n_vec[3] < 6378136 + 20000))
			{
				h0 = 0.0;
				ast0 = 1.228;
				akst1 = -9.0764E-2;
				akst2 = -2.0452E-3;
			}
			else if ((n_vec[3] >= 6378136 + 20000) && ( n_vec[3] < 6378136 + 60000))
			{
				h0 = 20.0;
				ast0 = 9.013E-2;
				akst1 = -0.16739;
				akst2 = 6.2669E-4;
			}
			else if ((n_vec[3] >= 6378136 + 60000) && (n_vec[3] < 6378136 + 100000))
			{
				
				h0 = 60.0;
				ast0 = 3.104E-4;
				akst1 = -0.137;
				akst2 = -7.8653E-4;
			}
			else if ((n_vec[3] >= 6378136 + 100000) && (n_vec[3] <= 6378136 + 120000))
			{
				h0 = 100.0;
				ast0 = 3.66E-7;
				akst1 = -0.18553;
				akst2 = 1.5397E-3;
			}

			double power = (n_vec[3]*1e-3 - 6378.136 - h0) * (akst1 + akst2 * (n_vec[3]*1e-3 - 6378.136 - h0));
			return ast0 * exp(power);
		}
		//cout << "Rad: " << n_vec[3] << " Rad Idx: " << rIdx << endl;
		int zeroIdx = zeroSearcherA(rIdx, n_vec[0], n_vec[1], n_vec[2]);
		//cout << zeroIdx << endl;
		
		int tr = searcherA(zeroIdx, rIdx, n_vec[0], n_vec[1], n_vec[2]);
		
		double p = counterA(n_vec[0], n_vec[1], n_vec[2], n_vec[3], n_vec[4], n_vec[5], tr, rIdx);

		return p;
	}
}