#include "CPFLoader.h"
#include "ephinit.h"
#include "ephinterp.h"
#include "testcpf.h"
#include <stdio.h>

struct eh ephihdr;
struct ei ephi;

CPFLoader::CPFLoader()
{
};
CPFLoader::~CPFLoader()
{
};
void CPFLoader::LoadCpffile( char *Fname )
{
	if ( cephinit( Fname, ephihdr, ephi ) != 0)
	{
		printf("File %s does not exist.\n", Fname);
		return;
	}
};
void CPFLoader::GetPos( double jd, double secinday, double *outvec )
{
	double backvec[6], aberout[3], aberback[3];
	double relcorrout, relcorrback;
	int fvel= 1;
	int ircode;
	//double outvec[6];

	// момент времени дл€ начала
	double jdint= jd;

	// секунды
	double seci= secinday;
	// интерпол€ци€
	cinterp_pvac (jdint, seci, outvec, backvec, aberout, aberback, &relcorrout, &relcorrback, fvel, &ircode, ephihdr, ephi );
	// out res
	//printf ("%9.1f %8.2f%20.3f%20.3f%20.3f\n", jdint-2400000.5e0, seci, outvec[0],outvec[1],outvec[2]);
	//printf ("                  %20.6f%20.6f%20.6f\n", outvec[3],outvec[4],outvec[5] );


	for( int it = 0; it < 6; it++ )
		outvec[it] = outvec[it]/1000.0;
};

int CPFLoader::GetNoradId()
{
	return ephihdr.norad;
}