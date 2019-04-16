#include "InfluenceForce/InfluenceForce.h"

#include "common/TLELoader.h"
#include "common/SimpleDatLoader.h"
#include "common/StructTypes.h"
#include "common/mytypes.h"

#include "CPFLoader.h"

#ifndef _CorrespondenceData_H_
#define _CorrespondenceData_H_

class CorrespondenceData
{
private:
	// модуль воздействий
	Force::InfluenceForce *IForce;
	// загрузчии данных
	TLELoader tleload;
	SimpleDatLoader optLoad;
	CPFLoader cpfe;

	// орбита по NORAD
	cOrbit *orbitN;
	// измерения
	std::vector< OpticPoint > OpticArray;

	void RefractionCorrect( double inRa, double inDec, double *outRa, double *outDec, double *tel_icrf );
	double GetZenitAngle( double inRa, double inDec, double *tel_icrf);
	void CorrespondenceCPF( int it, bool Tcorr, double *Tpos, bool abbcorr, bool refcorr );
	void CalcError( double Ra, double Dec, double Ra1, double Dec1, double optRa, double optDec, double *Telicrf );
	void refractcorrect( double optRa, double optDec, double *Telicrf, double *outRa, double *outDec, double SatDist, double rDist );
	void ConvertXYZtoRADEC( double *resultPosition, double *inTelescopePosition, double *Ra, double *Dec );
	void ConvertTEMEtoICRF( double *inPTEME, double *outPICRF, TimePoint calcTime );

	// novas доработать
	void Ter2Cel( double jday, double jdmeg, double *vec1, double *vec2 ); 
	void EquToHor( double jd, double *pos, double *zd, double *az, double *rar, double *decr );
	//double TLEidentify( char *fname, double *Tel_ITRF, bool printrep );

	// rafraction
	double refractivity( double pMillibar, double tKelvin );
	void calculateF( double *F1p, double *F2p, double *F3p, double *F4p, double F0, double F1, double FF1, double FF2 );
	void calculateG( double *G1p, double *G2p, double *G3p, double *G4p, double G0, double G1, double GG1, double GG2 );
	double dE_Theta( double pMillibar, double tKelvin, double Rkm, double Dkm, double THETAradians );
	double dE_Elev( double pMillibar, double tKelvin, double Rkm, double Dkm, double Eradians );
	void TestCorrectRefract();

public:

	CorrespondenceData(){};
	~CorrespondenceData(){};

	void InitModyle();

	void RunCorrespondenceData( char *optfname, char *htsName, char *tleName, double *Tpos );

	void RunCorrespondenceDataRLS(  );
};

#endif