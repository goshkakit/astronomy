//==============================================================================//
// Andrianov N.G.
// opbit predict 
// Load data orbit TLE
//==============================================================================//
#include<iostream>
#include<fstream>

#include "TLELoader.h"

TLELoader::TLELoader()
{
	bool loaddata = false;
};
TLELoader::~TLELoader()
{
};
//==============================================================================//
//
//==============================================================================//
void TLELoader::LoadData( const char *fname, int numbLine )
{
	if( loaddata == true )
		clear();

	std::ifstream input( fname );
	while(1)
	{
		string str1;
		string str2;
		string str3;

		// трехстрочный норад
		if( numbLine == 3 )
		{
			std::getline(input, str1);
			std::getline(input, str2);
			std::getline(input, str3);
		}
		// двухстрочный норад
		else if( numbLine == 2 )
		{
			str1 = "Sat NONE";
			std::getline(input, str2);
			std::getline(input, str3);
		}

		if( input.eof() || input.fail() )
			break;

		cTle tleSGP4(str1, str2, str3);

		// vector
		cOrbit *dataOrbit = new cOrbit( tleSGP4 );
		NORADList.push_back( dataOrbit );
	}

	input.close();

	printf( "NORAD Size = %zu\n", NORADList.size() );

	// Print TLE data
	printf("%s\t",  NORADList[0]->SatName().c_str() );
	printf("ID = %s\n", NORADList[0]->SatId().c_str() );

	int it = NORADList.size() - 1;
	printf("%s\t",  NORADList[it]->SatName().c_str() );
	printf("ID = %s\n", NORADList[it]->SatId().c_str() );

	loaddata = true;
};
//==============================================================================//
//
//==============================================================================//
void TLELoader::clear()
{
	for( int it = 0; it < NORADList.size(); it++ )
		delete NORADList[it];

	NORADList.clear();
	loaddata = false;
	printf( "clear TLE data\n" );
};
//==============================================================================//
//
//==============================================================================//
cOrbit* TLELoader::GetFromID( int id )
{
	for( int it = 0; it < NORADList.size(); it++ )
	{
		int lid = atoi( NORADList[it]->SatId().c_str() );
		if( id == lid )
			return NORADList[it];
	}

	return NULL;
}
//==============================================================================//
