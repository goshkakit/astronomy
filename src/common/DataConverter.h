//==============================================================================//
// Andrianov N.G.
// opbit predict 
// Time Convert
//==============================================================================//

//==============================================================================//
//
//==============================================================================//
class DataConverter
{
private:
public:
	DataConverter();
	~DataConverter();

	double SECtoHHMMSS( double data, double jd );
	double YYYYMMDDtoJD( double dt );
	double JDtoYYYYMMDD( double ajd );
	int NINTC( double A );
};
//==============================================================================//