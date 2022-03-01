//==============================================================================//
// Andrianov N.G.
// opbit predict 
// Load data orbit Optics
//==============================================================================//

#include <vector>

#ifndef _SUMPLELOADER_H_
#define _SUMPLELOADER_H_

struct OpticPoint
{
	double jd;
	double Ra;
	double Dec;

	// data UTC YYYYMMDD HHMMSS.SSSS
	double dataUTC;
	double timeUTC;

	// data MDB  YYYYMMDD HHMMSS.SSSS UTC + 3:00
	double dataMDB;
	double timeMDB;
};

//==============================================================================//
//
//==============================================================================//

class SimpleDatLoader
{
private:
	bool loaddata;
	std::vector< OpticPoint > OpticList;

	double tx, ty, tz;
	double TelPos[3];

public:

	SimpleDatLoader();
	~SimpleDatLoader();
	void LoadData( const char *fname, bool AddFlag );
	void clear();

	double GetJDTimeFromFile( const char *fname );

	std::vector< OpticPoint > GetOpticList(){ return OpticList; };

	void GetTelPos( double *lteitrf )
	{
		lteitrf[0] = TelPos[0];
		lteitrf[1] = TelPos[1];
		lteitrf[2] = TelPos[2];
	}

	// дата первого измерения
	int GetDataMeg()
	{ 
		if( OpticList.size() > 0 )
		{
			return (int)OpticList[0].dataUTC;
		}
		return -1;
	};

	double GetJDMeg()
	{ 
		if( OpticList.size() > 0 )
		{
			return OpticList[0].jd;
		}
		return -1;
	}
};
//==============================================================================//
#endif