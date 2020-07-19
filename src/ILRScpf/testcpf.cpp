#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ephinit.h"
#include "ephinterp.h"
#include "conio.h"
#include "testcpf.h"

struct eh ephihdr;
struct ei ephi;
//---------------------------------------------------------------------------------//
// ����
//---------------------------------------------------------------------------------//
int testcpf()
{
	double outvec[6], backvec[6], aberout[3], aberback[3];
	double relcorrout, relcorrback;
	int fvel= 1;
	int ircode;

	// �������� ������ ���������
	char str[256] = "tabl_y2005d080t1353_1500.0018";
	if ( cephinit(str, ephihdr, ephi ) != 0)
	{
		printf("File %s does not exist.\n",str);
		getch();
		return 0;
	}

	// ������ ������� ��� ������
	double jdint= 2400000.5e0+ 53450.0;
	for ( int i = 0; i <= 10; i++ )
	{
		// �������
		double seci= 50400 + i/20.0;
		// ������������
		cinterp_pvac (jdint, seci, outvec, backvec, aberout, aberback, &relcorrout, &relcorrback, fvel, &ircode, ephihdr, ephi );
		// out res
		printf ("%9.1f %8.2f%20.3f%20.3f%20.3f\n", jdint-2400000.5e0, seci, outvec[0],outvec[1],outvec[2]);
		printf ("                  %20.6f%20.6f%20.6f\n", outvec[3],outvec[4],outvec[5] );
	}

	getch();

	return 0;
}
//---------------------------------------------------------------------------------//
