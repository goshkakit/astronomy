
#pragma once
#ifndef __ÑOORDINATE_H__
#define __ÑOORDINATE_H__

#include "math.h"
#include <vector>
#include <iostream>

namespace VecMath
{

	struct SCoordinate
	{
	public:
		double x;
		double y;
		SCoordinate(double _x = 0, double _y = 0) : x(_x), y(_y) {}
		void move(double _x, double _y);
		void print();
		bool operator==(SCoordinate& aRight) const;
		SCoordinate operator+(SCoordinate& r);
		SCoordinate operator-(SCoordinate& r);
		SCoordinate operator*(double r);
		SCoordinate operator/(double r);
		SCoordinate direction();
		double length();

		bool isZero();
	};

	struct S3DCoordinate
	{
		double x;
		double y;
		double z;
		S3DCoordinate(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z)
		{}
		void print();
		double norm();
		bool operator==(S3DCoordinate& aRight) const;
		//	S3DCoordinate operator=(S3DCoordinate &r);
		S3DCoordinate direction();
		S3DCoordinate operator=(const S3DCoordinate& r);
		S3DCoordinate operator+(S3DCoordinate& r);
		S3DCoordinate operator-(S3DCoordinate& r);
		S3DCoordinate operator*(double r);
		double operator*(S3DCoordinate& r);
		S3DCoordinate operator/(double r);
		S3DCoordinate operator^(S3DCoordinate& r);
		friend S3DCoordinate operator-(const S3DCoordinate& r);
	};

	struct SRect
	{
		int xLeft, xRight, yTop, yBottom;
		SRect(int _xLeft = 0, int _xRight = 0, int _yTop = 0, int _yBottom = 0)
			: xLeft(_xLeft), xRight(_xRight), yTop(_yTop), yBottom(_yBottom)
		{}
		SRect(const char* FileName);
		void scale(const double alfa) { xLeft *= alfa;  xRight *= alfa; yBottom *= alfa; yTop *= alfa; }

		bool isVoid();
		int width();
		int height();
	};

	struct SLine
	{
		int xLeft, xRight, yTop, yBottom;
		SLine(int _xLeft = 0, int _yBottom = 0, int _xRight = 0, int _yTop = 0)
			: xLeft(_xLeft), yBottom(_yBottom), xRight(_xRight), yTop(_yTop)
		{}
	};

	struct SSize
	{
		int width;
		int height;
		SSize(int aWidth = 0, int aHeight = 0)
			: width(aWidth), height(aHeight)
		{}
		int Area();
	};

	struct SRGB
	{
		int R;
		int G;
		int B;
		SRGB(int aR = 255, int aG = 255, int aB = 255) : R(aR), G(aG), B(aB)
		{};
	};

	struct STriangle
	{
		STriangle(int _n1 = 0, int _n2 = 0, int _n3 = 0);
		int n1, n2, n3;

	};

	struct SCTriangle
	{
		SCTriangle(S3DCoordinate _n1 = S3DCoordinate(), S3DCoordinate _n2 = S3DCoordinate(), S3DCoordinate _n3 = S3DCoordinate());
		S3DCoordinate n1, n2, n3;
		S3DCoordinate centre();
		double angle(S3DCoordinate aVertex);
	};

	double distance3d(S3DCoordinate aPoint1, S3DCoordinate aPoint2);

	double distance2d(SCoordinate aPoint1, SCoordinate aPoint2);

	struct SMatchingPointsPair
	{
		SCoordinate firstImagePoint;
		SCoordinate secondImagePoint;
		SMatchingPointsPair(SCoordinate first = SCoordinate(), SCoordinate second = SCoordinate()) : firstImagePoint(first), secondImagePoint(second)
		{}
		bool operator==(SMatchingPointsPair& aRight)
		{
			return firstImagePoint == aRight.firstImagePoint && secondImagePoint == aRight.secondImagePoint;
		}
	};

	struct SMatchingInfo
	{
		std::vector<SMatchingPointsPair> matchingPointsVector;
		SMatchingPointsPair bestPair;
		int dispertion;
		double discrepancy;
		double err3d;
		int meanBright;
		SMatchingInfo(int _dispertion = 0, int _discrepancy = 0, double _err3d = 0,
			int _meanBright = 0) : dispertion(_dispertion), discrepancy(_discrepancy),
			err3d(_err3d), meanBright(_meanBright)
		{}
		SMatchingPointsPair getBestPair();
	};

	class S3DMatrix
	{
	public:
		S3DMatrix();
		double a[3][3];

		S3DMatrix operator+(S3DMatrix r);
		S3DMatrix operator-(S3DMatrix r);
		S3DCoordinate operator*(S3DCoordinate r);
		S3DMatrix operator*(S3DMatrix r);
		S3DMatrix operator*(double alfa);
		S3DMatrix inverse();

		double GetDET()
		{
			double det = a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
			//printf( "det = %f\n", det );
			return det;
		}

		static S3DMatrix E() {
			S3DMatrix M;
			M.a[0][0] = 1;
			M.a[1][1] = 1;
			M.a[2][2] = 1;
			return M;
		}
		static S3DMatrix R12(float alfa) {
			S3DMatrix M;
			M.a[0][0] = cos(alfa);  M.a[0][1] = sin(alfa);
			M.a[1][0] = -sin(alfa); M.a[1][1] = cos(alfa);
			M.a[2][2] = 1;
			return M;
		}
		static S3DMatrix R13(float alfa) {
			S3DMatrix M;
			M.a[0][0] = cos(alfa);  M.a[0][2] = sin(alfa);
			M.a[2][0] = -sin(alfa); M.a[2][2] = cos(alfa);
			M.a[1][1] = 1;
			return M;
		}
		static S3DMatrix R23(float alfa) {
			S3DMatrix M;
			M.a[1][1] = cos(alfa);  M.a[1][2] = sin(alfa);
			M.a[2][1] = -sin(alfa); M.a[2][2] = cos(alfa);
			M.a[0][0] = 1;
			return M;
		}
		static S3DMatrix R123(float alfa12, float alfa23, float alfa13) {
			return R13(alfa13) * R23(alfa23) * R12(alfa12);
		}

		void print()
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					printf("%f ", a[i][j]);
				}
				printf("\n");
			}
		}

	};

	class CameraProjector
	{
	public:
		S3DMatrix R;
		S3DCoordinate t;
		CameraProjector() {
			centerCalculated = false;
		};
		CameraProjector(const double a[3][4]);
		CameraProjector(const S3DMatrix _R, const S3DCoordinate _t);
		CameraProjector(const char* FileName);
		CameraProjector inverse();
		S3DCoordinate getCenter();
		void MoveArea(double dx, double dy);

		SCoordinate operator*(S3DCoordinate r);
		CameraProjector operator*(double alfa);
	private:
		bool centerCalculated;
		S3DCoordinate center;
	};
}
#endif