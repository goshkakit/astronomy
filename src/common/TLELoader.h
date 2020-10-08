//==============================================================================//
// Andrianov N.G.
// opbit predict 
// Load data orbit TLE
//==============================================================================//

#include <vector>
#include "..\\norad\\std_add.h"
#include "..\\norad\\coreLib.h"
#include "..\\norad\\cOrbit.h"

#ifndef _TLELOADER_H_
#define _TLELOADER_H_

//==============================================================================//
//
//==============================================================================//
class TLELoader
{
private:
	bool loaddata;
	

public:

	std::vector< cOrbit* > NORADList;

	TLELoader();
	~TLELoader();
	void LoadData( const char *fname, int numbLine, int* SatIDList = NULL, int SatIDList_size = 0);
	void clear();

	cOrbit* GetFromID( int id );
};
//==============================================================================//
#endif