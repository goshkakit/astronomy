#ifndef _ORBIT_PREDICTORBITTYPE_H_
#define _ORBIT_PREDICTORBITTYPE_H_
#include <math.h>
//==============================================================================//
// Andrianov N.G.
// opbit predict module
//==============================================================================//
// name
namespace Orbit
{
	struct PointDist
	{
		double d;
		double t;

		void init()
		{
			d = 0;
			t = 0;
		};
	};

	struct PointWithMinDist
	{
		double t;
		double d;
		double d_verify;

		int Nlist;
		int norad;
	};
	//================================================//
	// ������� ��������� �������� �� ������
	//================================================//
	struct SatStateVector
	{
		double posx;
		double posy;
		double posz;
		double posr;

		double vx;
		double vy;
		double vz;

		int flag;
		// ������� ������� �� �������
		void InitFromArray( double *x )
		{
			posx = x[0];
			posy = x[1];
			posz = x[2];
			vx = x[3];
			vy = x[4];
			vz = x[5];
			
			posr = posx*posx + posy*posy + posz*posz;
			posr = sqrt( posz );
		};
		// ����������� �������� ������� � ������
		void CopyToArray( double *x)
		{
			x[0] = posx;
			x[1] = posy;
			x[2] = posz;
			x[3] = vx;
			x[4] = vy;
			x[5] = vz;
		};
	};

	//================================================//
	// ������ ��� ��������������� ������
	//================================================//
	struct OrbitPoints
	{
		// ����� �����������
		int NY;
		int sizeMem;
		// ��������� ������ �������
		double TP0;
		double TPN;
		// ��������� ���������
		double *X0;
		// ����������������� ������
		// ��������� ����� ��� ���������� ���������������
		double *FP0;
		double *FP1;
		double *FP2;
		double *FP3;
		double *FP4;
		double *FP5;
		double *FP6;
		double *FP7;

		// ��������� ������
		void AllocMemory( int ny )
		{
			NY = ny;
			sizeMem = NY*sizeof(double);
			X0 = new double[ny];
			FP0 = new double[ny];
			FP1 = new double[ny];
			FP2 = new double[ny];
			FP3 = new double[ny];
			FP4 = new double[ny];
			FP5 = new double[ny];
			FP6 = new double[ny];
			FP7 = new double[ny];
		};
		// �������� ������
		void DeleteMemory()
		{
			delete X0;
			delete FP0;
			delete FP1;
			delete FP2;
			delete FP3;
			delete FP4;
			delete FP5;
			delete FP6;
			delete FP7;
		};
		// ������������� ��������� ����������
		void SetInitPosition( double *x0, double t0 )
		{
			TP0 = t0;
			for( int it = 0; it < NY; it++ )
			{
				X0[it] = x0[it];
				FP0[it] = x0[it];
			}
		};
	};

	//================================================//
	// ������ ���������
	//================================================//
	struct SatelliteArray
	{
		int NY;
		int Block;

		// �����
		double *TP0;
		// ������ ���������
		double *X0;
		double *Atm;
		double *Sun;

		// ����� ��������
		int *NoradID;

		void AllocMemory( int ny, int block )
		{
			NY = ny;
			Block = block;
			// ������ ��� ���� ��������
			TP0 = new double [block];
			X0 = new double [block*NY];
			NoradID = new int [block];

			Atm = new double [block];
			Sun= new double [block];
		};

		void DeleteMemory()
		{
			delete TP0;
			delete X0;
			delete NoradID;

			delete Atm;
			delete Sun;
		};

		// ���������� ������ ���������
		void SetSatellite( int it, double *inx, double inT )
		{
			// �����
			TP0[it] = inT; 
			// ������ ���������
			X0[it + 0*Block] = inx[0];
			X0[it + 1*Block] = inx[1];
			X0[it + 2*Block] = inx[2];
			X0[it + 3*Block] = inx[3];
			X0[it + 4*Block] = inx[4];
			X0[it + 5*Block] = inx[5];
		};
		// ��������� ������������� 
		void SetCoeffAtmSun( int it, double atm, double sun )
		{
			Atm[it] = atm;
			Sun[it] = sun;
		};
		// ����������� �����
		void SetSatelliteID( int it, double id )
		{
			NoradID[it] = id;
		};
		int GetSatelliteID( int it )
		{
			return NoradID[it];
		};
		// ��������� ��������
		void GetSatellite( int it, double *inx, double *inT )
		{
			// ������
			inT[0] = TP0[it]; 
			// ������ ���������
			inx[0] = X0[it + 0*Block];
			inx[1] = X0[it + 1*Block];
			inx[2] = X0[it + 2*Block];
			inx[3] = X0[it + 3*Block];
			inx[4] = X0[it + 4*Block];
			inx[5] = X0[it + 5*Block];
		};
		double GetCoeffAtm( int it )
		{
			return Atm[it];
		};
		double GetCoeffSun( int it )
		{
			return Sun[it];
		};
	};

	//================================================//
	// ������ �����
	//================================================//
	struct OrbitArrayPointsCPU
	{
		// ������ �������
		unsigned int mem_size;
		unsigned int mem_length;
		// ����� �����
		int N_hlist;
		// ����� �����
		int Length; 

		// ���������
		double **array_list; 

		void  InitArrayList( int N, int L )
		{
			// ������ ������ ��������
			mem_size = sizeof( double)*(L*4+1);
			mem_length = L*4+1;

			// ��������� ������ �� ����������
			N_hlist = N;
			Length = L;

			array_list = new double*[N_hlist];

			// ��������� ������� ������
			for ( int i = 0; i < N_hlist; i++ )
				array_list[i] = new double[ mem_length ];

		};

		void FreeArrayList()
		{
			// ������� ������� � �������
			for (int i = 0; i < N_hlist; i++)
			{
				delete array_list[i];
			}
			delete array_list;
		};
	};
	//================================================//
}
#endif