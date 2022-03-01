//==============================================================================//
// Andrianov N.G.
// opbit predict 
// Load data orbit TLE
//==============================================================================//
#include "SimpleDatLoader.h"
#include "DataConverter.h"

#include<iostream>
#include<fstream>
#include <string>

SimpleDatLoader::SimpleDatLoader()
{
	bool loaddata = false;
};
SimpleDatLoader::~SimpleDatLoader()
{
};
//==============================================================================//
// 
//==============================================================================//
double SimpleDatLoader::GetJDTimeFromFile( const char *fname )
{
	std::string strg;
	std::ifstream fin( fname );
	if(!fin)
	{
		std::cerr<<"File not found"; 
	}
	std::getline( fin, strg, '\n' );

	// comment
	if( strg[0] == '#' )
		std::getline( fin, strg, '\n' );

	// строка с параметрами
	std::getline(fin, strg, '\n');

	const char *buffl = strg.c_str();
	double ra, dec, jd;
	sscanf( buffl, "%lf %lf %lf", &ra, &dec, &jd );

	fin.close();

	return jd;
};
//==============================================================================//
// 
//==============================================================================//
void SimpleDatLoader::LoadData( const char *fname, bool AddFlag )
{
	if( loaddata == true && AddFlag == false )
		clear();

	DataConverter Dconv;

	FILE *fre;
	fre = fopen( fname, "r" );

	std::string strg;
	std::ifstream fin( fname );
	if(!fin)
	{
		std::cerr<<"File not found"; 
	}

	std::getline( fin, strg, '\n' );

	// comment
	if( strg[0] == '#' )
		std::getline( fin, strg, '\n' );

	const char *buffl = strg.c_str();
	float tmp, tmp2;
	//sscanf ( buffl,"%lf %lf %lf",&tx, &ty, &tz, &tmp ); // FIX read!!!!
	sscanf(buffl, "%lf %lf %lf %f %f", &tx, &ty, &tz, &tmp, &tmp2); // FIX read!!!! TODO!!
	// положение телескопа
	//fscanf( fre, "%lf %lf %lf", &tx, &ty, &tz );
	TelPos[0] = tx;
	TelPos[1] = ty;
	TelPos[2] = tz;

	double ra, dec, jd, px, py;
	double t1, t2, t3;
	char pathfits[512];
	char pathfits1[512];

	while( std::getline(fin, strg, '\n') )
	//while( fscanf( fre, "%lf %lf %lf %lf %lf %s", &ra, &dec, &jd, &px, &py, pathfits ) == 6 )
	//while( fscanf( fre, "%lf %lf %lf %s %s %lf %lf %lf %lf", &ra, &dec, &jd, pathfits, pathfits1, &px, &py, &t1, &t2 ) == 9 )
	//while( fscanf( fre, "%lf %lf %lf %s %s %lf %lf %lf %lf %lf", &ra, &dec, &jd, pathfits, pathfits1, &px, &py, &t1, &t2, &t3 ) == 10 )
	//while( fscanf( fre, "%lf %lf %lf %lf %lf %lf %lf %lf %s %s", &ra, &dec, &jd, &px, &py, &t1, &t2, &t3, pathfits, pathfits1 ) == 10 )
	{
		buffl = strg.c_str();
		sscanf( buffl, "%lf %lf %lf", &ra, &dec, &jd );
		//printf("%f %f %f\n", ra, dec, jd );

		//TMP !!!!!!
		//printf( "CORRECT TIME 1 sec\n" );
		//jd = jd + 1.0/86400.0;

		OpticPoint pt;
		pt.jd = jd;
		pt.Ra = ra;
		pt.Dec = dec;

		// UTC
		pt.jd = jd;
		pt.dataUTC = Dconv.JDtoYYYYMMDD( jd );
		pt.timeUTC = Dconv.SECtoHHMMSS( pt.dataUTC, jd );

		// MDB
		jd = jd + 0.125;
		pt.dataMDB = Dconv.JDtoYYYYMMDD( jd );
		pt.timeMDB = Dconv.SECtoHHMMSS( pt.dataMDB, jd );

		OpticList.push_back( pt );
	}
	fclose( fre );
	fin.close();

	printf( "OpticSize = %d\n", OpticList.size() );
	printf( "TelPosition = %f %f %f\n", tx, ty, tz );
	int it = OpticList.size()-1;
	printf( "end Ra = %f dec = %f jd = %f\n", OpticList[it].Ra, OpticList[it].Dec, OpticList[it].jd );
	loaddata = true;
};
//==============================================================================//
//
//==============================================================================//
void SimpleDatLoader::clear()
{
	OpticList.clear();
	loaddata = false;
};