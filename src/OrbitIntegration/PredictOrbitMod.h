#include "PredictOrbitSat.h"
#include "IPredictOrbitMod.h"
#include "DataConverter.h"
//==============================================================================//
// windows
#ifdef WIN32

class PredictOrbitMod  : public IPredictOrbitMod
{
private:

	friend CPredictOrbitMod;

	Force::InfluenceForce *IForce;
	Orbit::PredictOrbitSat POSat;
	DataConverter DC;

	int IntegrationListNorad();

public:

	PredictOrbitMod()
	{ 
		printf("PredictOrbitMod()\n");
		Init();
	};
	~PredictOrbitMod()
	{ 
		printf("~PredictOrbitMod()\n");
	};

	virtual int _stdcall Init();
	virtual int _stdcall DeInit();
	virtual int _stdcall RunTest();
	virtual int _stdcall TestAllForce();
	virtual int _stdcall IntegrationList();
	virtual int _stdcall GetNewPosition( SatParamToPredict &sptr );
	virtual int _stdcall GetNewPositionForJD(SatParamToPredictJD &sptr);
};

//==============================================================================//
// Linux
#else

class PredictOrbitMod
{
private:

	Force::InfluenceForce *IForce;
	Orbit::PredictOrbitSat POSat;

	int IntegrationListNorad();

public:

	PredictOrbitMod(){ printf("PredictOrbitMod()\n"); };
	~PredictOrbitMod(){ printf("~PredictOrbitMod()\n"); };

	int Init();
	int DeInit();
	int RunTest();
	int TestAllForce();
	int IntegrationList();
	int GetNewPosition( SatParamToPredict &sptr );
};
//==============================================================================//
#endif