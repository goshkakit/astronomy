//==============================================================================//
// Andrianov N.G.
// opbit predict module
// Integration motion
// RK method
//==============================================================================//
#include "PredictOrbitSat.h" 

#ifdef GPUCOMPILE

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "cutilNP.h"



// name
namespace Orbit
{

		// атмосфера
static const double hh_tick[42] = {  
	500.0, 600.0, 640.0, 1500.0, 600.0, 640.0,  
	500.0, 660.0, 700.0, 1500.0, 700.0, 660.0,  
	500.0, 760.0, 760.0, 1500.0, 780.0, 740.0,  
	500.0, 800.0, 820.0, 1500.0, 800.0, 800.0,  
	500.0, 860.0, 860.0, 1500.0, 800.0, 860.0,  
	500.0, 900.0, 920.0, 1500.0, 900.0, 900.0,  
	500.0,1000.0, 980.0, 1500.0, 760.0, 900.0 };
static const double h_aa0[14] = { 
		.26862900E+02, .27459800E+02, .28639500E+02, .29641800E+02,  
		.30167100E+02, .29757800E+02, .30785400E+02,  
		.17878100E+02,-.25490900E+01,-.13959900E+02,-.23307900E+02,  
		-.14726400E+02,-.49120000E+01,-.54095200E+01 };
static const double h_aa1[14] = { 
		-.45167400E+00,-.46366800E+00,-.49098700E+00,-.51495700E+00,  
		-.52783700E+00,-.51791500E+00,-.54569500E+00,  
		-.13202500E+00, .14006400E-01, .84495100E-01, .13514100E+00,  
		.71325600E-01, .10832600E-01, .55074900E-02 };
static const double h_aa2[14] = { 
		.29039700E-02, .29740000E-02, .32064900E-02, .34192600E-02,  
		.35321100E-02, .34269900E-02, .37032800E-02,  
		.22771700E-03,-.16946000E-03,-.32887500E-03,-.42080200E-03,  
		-.22801500E-03,-.81054600E-04,-.37885100E-04 };
static const double h_aa3[14] = { 
		-.10695300E-04,-.10753000E-04,-.11681000E-04,-.12578500E-04,  
		-.13022700E-04,-.12413700E-04,-.13707200E-04,  
		-.22543000E-06, .32719600E-06, .50591800E-06, .57371700E-06,  
		.28487000E-06, .11571200E-06, .24808000E-07 };
static const double h_aa4[14] = { 
		.22159800E-07, .21705900E-07, .23684700E-07, .25727000E-07,  
		.26645500E-07, .24820900E-07, .28061400E-07,  
		.13357400E-09,-.28763000E-09,-.39229900E-09,-.40323800E-09,  
		-.17438300E-09,-.81329600E-10, .49218300E-11 };
static const double h_aa5[14] = { 
		-.24294100E-10,-.23024900E-10,-.25180900E-10,-.27587400E-10,  
		-.28543200E-10,-.25841300E-10,-.30018400E-10,  
		-.45045800E-13, .12262500E-12, .15227900E-12, .14284600E-12,  
		.50807100E-13, .30491300E-13,-.86501100E-14 };
static const double h_aa6[14] = { 
		.10992600E-13, .10012300E-13, .10953600E-13, .12109100E-13,  
		.12500900E-13, .10938300E-13, .13114200E-13,  
		.67208600E-17,-.20573600E-16,-.23557600E-16,-.20172600E-16,  
		-.53495500E-17,-.49498900E-17, .19849000E-17 };
static const double h_bb0[14] = { 
		.68789400E-01, .15073000E+00, .47945100E-01, .22344800E-01,  
		.32639100E-02,-.51474900E-01,-.10725500E+00,  
		.23158400E+02, .33273200E+02, .39196100E+02, .43246900E+02,  
		.49573800E+02, .11278000E+02,-.52618400E+02 };
static const double h_bb1[14] = { 
		-.28407700E-02,-.40088900E-02,-.23945300E-02,-.19798000E-02,  
		-.15986900E-02,-.92105900E-03,-.17434300E-03,  
		-.80214700E-01,-.11109900E+00,-.12352000E+00,-.12697300E+00,  
		-.13861300E+00, .14347800E-02, .21468900E+00 };
static const double h_bb2[14] = { 
		.18392200E-04, .24393700E-04, .17033500E-04, .15410100E-04,  
		.14044300E-04, .11514700E-04, .90275900E-05,  
		.10582400E-03, .14142100E-03, .14901500E-03, .14263700E-03,  
		.14785100E-03,-.36984600E-04,-.29488200E-03 };
static const double h_bb3[14] = { 
		.91960500E-08,-.99277200E-08,-.13162600E-08,-.23543000E-08,  
		-.30228700E-08,-.12290100E-08,-.31651200E-09,  
		-.61503600E-07,-.79495200E-07,-.79705000E-07,-.70998500E-07,  
		-.69636100E-07, .35831800E-07, .17117100E-06 };
static const double h_bb4[14] = { 
		-.41687300E-10,-.18223900E-10,-.17403200E-10,-.12499400E-10,  
		-.92016000E-11,-.81310400E-11,-.61400000E-11,  
		.13245300E-10, .16583600E-10, .15877200E-10, .13164600E-10,  
		.12159500E-10,-.99122500E-11,-.36058200E-10 };
static const double h_cc0[14] = { 
		-.10482500E+01,-.93106000E+00,-.82086700E+00,-.74404700E+00,  
		-.72247100E+00,-.68748200E+00,-.73998400E+00,  
		.50503400E+02, .61624000E+02, .53262300E+02, .18223600E+02,  
		-.31844200E+02,-.48720800E+02,-.14785900E+03 };
static const double h_cc1[14] = { 
		.16630500E-01, .14153700E-01, .11991600E-01, .10474300E-01,  
		.98031700E-02, .91659400E-02, .95285400E-02,  
		-.17054100E+00,-.19296700E+00,-.14434200E+00,-.84002400E-02,  
		.16832700E+00, .22299600E+00, .53165200E+00 };
static const double h_cc2[14] = { 
		-.92426300E-04,-.72986200E-04,-.57983500E-04,-.47854400E-04,  
		-.42524500E-04,-.38093200E-04,-.36272700E-04,  
		.21723200E-03, .22806100E-03, .14659000E-03,-.38800000E-04,  
		-.26260300E-03,-.32188400E-03,-.67193700E-03 };
static const double h_cc3[14] = { 
		.27238200E-06, .20029400E-06, .15070700E-06, .11851300E-06,  
		.99554400E-07, .85127500E-07, .73887000E-07,  
		-.12190200E-06,-.11871500E-06,-.64644300E-07, .43138400E-07,  
		.16545400E-06, .19149500E-06, .36478700E-06 };
static const double h_cc4[14] = { 
		-.24135500E-09,-.16200600E-09,-.11302600E-09,-.83149800E-10,  
		-.65517500E-10,-.52997200E-10,-.42390700E-10,  
		.25403700E-10, .22963800E-10, .10422700E-10,-.12383200E-10,  
		-.36935500E-10,-.40806700E-10,-.72626800E-10 };
static const double h_dd0[14] = { 
		-.35189900E+00,-.47813000E-01, .20981000E+00, .26517400E+00,  
		.23047000E+00, .17007400E+00, .88141000E-01,  
		-.35189900E+00,-.47813000E-01, .20981000E+00, .26517400E+00,  
		.23047000E+00, .17007400E+00, .88141000E-01 };
static const double h_dd1[14] = { 
		.57705600E-02, .38081300E-02, .26288100E-02, .27583600E-02,  
		.33833100E-02, .40613100E-02, .46825300E-02,  
		.57705600E-02, .38081300E-02, .26288100E-02, .27583600E-02,  
		.33833100E-02, .40613100E-02, .46825300E-02 };
static const double h_dd2[14] = { 
		.99581900E-06, .42277100E-05, .42437900E-05, .20866800E-05,  
		-.55230500E-06,-.28211400E-05,-.42460900E-05,  
		.99581900E-06, .42277100E-05, .42437900E-05, .20866800E-05,  
		-.55230500E-06,-.28211400E-05,-.42460900E-05 };
static const double h_dd3[14] = { 
		-.72532400E-08,-.86682600E-08,-.66732800E-08,-.36954300E-08,  
		-.82360700E-09, .13836900E-08, .25350900E-08,  
		-.72532400E-08,-.86682600E-08,-.66732800E-08,-.36954300E-08,  
		-.82360700E-09, .13836900E-08, .25350900E-08 };
static const double h_dd4[14] = { 
		.29759000E-11, .30671200E-11, .21349600E-11, .11186200E-11,  
		.22134900E-12,-.42790800E-12,-.72903100E-12,  
		.29759000E-11, .30671200E-11, .21349600E-11, .11186200E-11,  
		.22134900E-12,-.42790800E-12,-.72903100E-12 };
static const double h_ee0[14] = { 
		-.73159600E+00,-.75217500E+00,-.57047600E+00,-.94957300E+00,  
		-.96759800E+00,-.10227800E+01,-.75790300E+00,  
		.38619900E+02, .51249000E+02, .68474600E+02, .58422000E+02,  
		.72018800E+01, .21594800E+02,-.88407600E+02 };
static const double h_ee1[14] = { 
		.59734500E-02, .56592500E-02, .29580200E-02, .81312100E-02,  
		.84199100E-02, .92363300E-02, .60606800E-02,  
		-.13214700E+00,-.16737300E+00,-.21565900E+00,-.16666400E+00,  
		.21610900E-01,-.20223900E-01, .33851800E+00 };
static const double h_ee2[14] = { 
		-.58203700E-05, .18082000E-05, .16889600E-04,-.38781300E-05,  
		-.35850000E-05,-.61012800E-05, .78529600E-05,  
		.17541100E-03, .21183200E-03, .26227300E-03, .18548600E-03,  
		-.65288200E-04,-.17202900E-04,-.44558100E-03 };
static const double h_ee3[14] = { 
		.68463400E-07, .33382200E-07,-.47475600E-08, .23769400E-07,  
		.17480100E-07, .17821100E-07,-.97489100E-08,  
		-.10241700E-06,-.11822100E-06,-.14097200E-06,-.91234500E-07,  
		.53707700E-07, .28301700E-07, .25172900E-06 };
static const double h_ee4[14] = { 
		-.95048300E-10,-.51396500E-10,-.17271100E-10,-.27746900E-10,  
		-.19622100E-10,-.17007300E-10, .15837700E-11,  
		.22144600E-10, .24505500E-10, .28228500E-10, .16711800E-10,  
		-.14095000E-10,-.89448600E-11,-.52030000E-10 };
static const double h_ff1[7] = { 
		.54110000E+00, .55150000E+00, .55850000E+00, .55850000E+00,  
		.55850000E+00, .55850000E+00, .55850000E+00 };
static const double h_ee5[7] = { 
		-.20670000E+00,-.16971000E+00,-.14671000E+00,-.13150000E+00,  
		-.12091600E+00,-.11363000E+00,-.10444000E+00 };
static const double h_ee6[7] = { 
		.97533000E-01, .79830000E-01, .68808000E-01, .61603000E-01,  
		.56538000E-01, .53178000E-01, .48551000E-01 };
static const double h_ee7[7] = { 
		-.11817000E-01,-.94393000E-02,-.79836000E-02,-.70866000E-02,  
		-.64324000E-02,-.60436000E-02,-.53567000E-02 };
static const double h_ee8[7] = { 
		.16145000E-02, .12622000E-02, .10535000E-02, .92813000E-03,  
		.83723000E-03, .77982000E-03, .68809000E-03 };
static const double h_eet5[14] = { 
		-.20610000E+00,-.16927900E+00,-.14637700E+00,-.13121000E+00,  
		-.12067000E+00,-.11339900E+00,-.10424300E+00,  
		-.20610000E+00,-.16927900E+00,-.14637700E+00,-.13121000E+00,  
		-.12067000E+00,-.11339900E+00,-.10424300E+00 };
static const double h_eet6[14] = { 
		.94449000E-01, .77599000E-01, .67052000E-01, .60105000E-01,  
		.55232000E-01, .51994000E-01, .47573000E-01,  
		.94449000E-01, .77599000E-01, .67052000E-01, .60105000E-01,  
		.55232000E-01, .51994000E-01, .47573000E-01 };
static const double h_eet7[14] = { 
		-.87953000E-02,-.71375000E-02,-.60951000E-02,-.54388000E-02,  
		-.49580000E-02,-.46876000E-02,-.41711000E-02,  
		-.87953000E-02,-.71375000E-02,-.60951000E-02,-.54388000E-02,  
		-.49580000E-02,-.46876000E-02,-.41711000E-02 };
static const double h_eet8[14] = { 
		.88385000E-03, .69025000E-03, .57456000E-03, .50585000E-03,  
		.45512000E-03, .42548000E-03, .37068000E-03,  
		.88385000E-03, .69025000E-03, .57456000E-03, .50585000E-03,  
		.45512000E-03, .42548000E-03, .37068000E-03 };
static const double h_aal0[14] = { 
		-.40776800E+00,-.90273900E+00,-.73303700E+00,-.13144400E+01,  
		-.12002600E+01,-.15215800E+01,-.16766400E+01,  
		.48653600E+02, .54486700E+02, .60126700E+02, .47099600E+02,  
		.50617400E+02, .80194200E+01,-.15572800E+02 };
static const double h_aal1[14] = { 
		.14850600E-02, .82680300E-02, .52339600E-02, .13312400E-01,  
		.11408700E-01, .15704000E-01, .17719400E-01,  
		-.17029100E+00,-.17829800E+00,-.18314400E+00,-.12526000E+00,  
		-.12904700E+00, .18530200E-01, .93670400E-01 };
static const double h_aal2[14] = { 
		.12535700E-04,-.12544800E-04, .63566700E-05,-.25558500E-04,  
		-.14732400E-04,-.30285900E-04,-.36949800E-04,  
		.22624200E-03, .22272500E-03, .21248100E-03, .12635200E-03,  
		.12484200E-03,-.61473300E-04,-.14903600E-03 };
static const double h_aal3[14] = { 
		.37731100E-07, .61285300E-07, .10906500E-07, .54398100E-07,  
		.27804000E-07, .45766800E-07, .50913400E-07,  
		-.13203200E-06,-.12270000E-06,-.10849700E-06,-.55158400E-07,  
		-.52499300E-07, .49767400E-07, .94215100E-07 };
static const double h_aal4[14] = { 
		-.77895300E-10,-.70796600E-10,-.26142700E-10,-.43378400E-10,  
		-.22632000E-10,-.28292600E-10,-.28287800E-10,  
		.28519300E-10, .25131600E-10, .20571000E-10, .87527200E-11,  
		.80827200E-11,-.12616200E-10,-.20961000E-10 };

__device__ __constant__ int KKR[12];

__device__ __constant__ double Dca[120]; 
__device__ __constant__ double Dca1[40];
__device__ __constant__ double Dcw[76];
__device__ __constant__ double Dcg[6];
__device__ __constant__ double Dcc[16];

__device__ __constant__ double CUtick[42];

__device__ __constant__ double CUaal0[14];
__device__ __constant__ double CUaal1[14];
__device__ __constant__ double CUaal2[14];
__device__ __constant__ double CUaal3[14];
__device__ __constant__ double CUaal4[14];

__device__ __constant__ double CUaa0[14];
__device__ __constant__ double CUaa1[14];
__device__ __constant__ double CUaa2[14];
__device__ __constant__ double CUaa3[14];
__device__ __constant__ double CUaa4[14];
__device__ __constant__ double CUaa5[14];
__device__ __constant__ double CUaa6[14];

__device__ __constant__ double CUbb0[14];
__device__ __constant__ double CUbb1[14];
__device__ __constant__ double CUbb2[14];
__device__ __constant__ double CUbb3[14];
__device__ __constant__ double CUbb4[14];

__device__ __constant__ double CUcc0[14];
__device__ __constant__ double CUcc1[14];
__device__ __constant__ double CUcc2[14];
__device__ __constant__ double CUcc3[14];
__device__ __constant__ double CUcc4[14];

__device__ __constant__ double CUdd0[14];
__device__ __constant__ double CUdd1[14];
__device__ __constant__ double CUdd2[14];
__device__ __constant__ double CUdd3[14];
__device__ __constant__ double CUdd4[14];

__device__ __constant__ double CUee0[14];
__device__ __constant__ double CUee1[14];
__device__ __constant__ double CUee2[14];
__device__ __constant__ double CUee3[14];
__device__ __constant__ double CUee4[14];
		
__device__ __constant__ double CUee5[7];
__device__ __constant__ double CUee6[7];
__device__ __constant__ double CUee7[7];
__device__ __constant__ double CUee8[7];
__device__ __constant__ double CUff1[7];

__device__ __constant__ double CUeet5[14];
__device__ __constant__ double CUeet6[14];
__device__ __constant__ double CUeet7[14];
__device__ __constant__ double CUeet8[14];

// вычисление воздействий на планеты
#include "cuInfluencePlanet.cu"
#include "cuInfluenceEarthRotation.cu"
#include "cuInfluenceEGM96.cu"
#include "cuInfluenceSun.cu"
#include "cuInfluenceAtmosphere.cu"

//=====================================================================//
// массивы с орбитами
//=====================================================================//
struct OrbitArrayPointsGPU
{
	// размер массива
	unsigned int mem_size;
	// число орбит
	int N_hlist;
	// число точек
	int Length; 

	// указатели
	double **h_array_list; 
	double **d_array_list;

	void InitArrayList( int N, int L )
	{
		// размер памяти массивов
		mem_size = sizeof(double)*(L*4+1);

		// выделение памяти на видеокарте
		N_hlist = N;
		Length = L;

		h_array_list = new double*[N_hlist];
		cutilSafeCall( cudaMalloc((void**)&d_array_list, N_hlist*sizeof( double* ) ));

		// выделение массива данных
		for ( int i = 0; i < N_hlist; i++ )
		{
			// allocate memory for one arrays on the device
			cutilSafeCall( cudaMalloc( (void**)&h_array_list[i], mem_size ) );
		}
		// копирование на устройство массива с указателями на буферы
		cutilSafeCall( cudaMemcpy( d_array_list, h_array_list, N_hlist*sizeof(double*), cudaMemcpyHostToDevice ));
	};

	void FreeArrayList()
	{
		// очистка массива и буферов
		for (int i = 0; i < N_hlist; i++)
		{
			cudaFree( h_array_list[i] );
		}
		cudaFree( d_array_list );

		delete h_array_list;
	};

	void CopyFromGPU( Orbit::OrbitArrayPointsCPU &outArr )
	{
		for (int i = 0; i < N_hlist; i++)
		{
			cutilSafeCall( cudaMemcpy( outArr.array_list[i], h_array_list[i], mem_size, cudaMemcpyDeviceToHost ) );
		}
	};
};
//=====================================================================//
// структура данных для GPU
//=====================================================================//
struct gpuOrbitPoint
{
	int NY;
	int Block;
	int sizeMem;
	int sizeMemOne;
	int sizeMemT;

	double *TP0;
	double *X0;
	double *FP0;
	double *FP1;
	double *FP2;
	double *FP3;
	double *FP4;
	double *FP5;
	double *FP6;
	double *FP7;

	double ajd0;
	double delt0;
	double Satm;
	double Ssun;

	CurrentIntegrateParam *Iparam;

	// указатель на эфемериды
	double *d_EF;
	// указатель на таблицы поправок времени
	double *dT_finals_tab;
	int *dT_finals_n;
	// массивы поправок для нутации
	double *dNUT_AMPL;
	double *dNUT_ARG;
	// для гармоник
	double *d_Garmonic;

	void AllocMemory( int ny, int block )
	{
		NY = ny;
		Block = block;
		sizeMem = block*NY*sizeof(double);
		sizeMemOne = NY*sizeof(double);
		sizeMemT = block*sizeof(double);

		// память под блок векторов
		cudaMalloc((void**)&TP0, sizeMemT );
	
		cudaMalloc((void**)&X0, sizeMem );

		cudaMalloc((void**)&FP0, sizeMem );
		cudaMalloc((void**)&FP1, sizeMem );
		cudaMalloc((void**)&FP2, sizeMem );
		cudaMalloc((void**)&FP3, sizeMem );
		cudaMalloc((void**)&FP4, sizeMem );
		cudaMalloc((void**)&FP5, sizeMem );
		cudaMalloc((void**)&FP6, sizeMem );
		cudaMalloc((void**)&FP7, sizeMem );

		cudaMalloc((void**)&Iparam, block*sizeof( CurrentIntegrateParam ) );
	};

	void DeleteMemory()
	{
		cudaFree(TP0);
		cudaFree(X0);
		cudaFree(FP0);
		cudaFree(FP1);
		cudaFree(FP2);
		cudaFree(FP3);
		cudaFree(FP4);
		cudaFree(FP5);
		cudaFree(FP6);
		cudaFree(FP7);
		cudaFree(Iparam);
	};

	void CopyToGPU( Orbit::SatelliteArray &inListSat )
	{
		cutilSafeCall( cudaMemcpy( TP0, inListSat.TP0, sizeMemT, cudaMemcpyHostToDevice ) );
		// задание векторов состояния объектов - положение, скорость
		// хранение последовательное х1 х2..... y1 y2.....
		cutilSafeCall( cudaMemcpy( X0, inListSat.X0, sizeMem, cudaMemcpyHostToDevice ) );
		cutilSafeCall( cudaMemcpy( FP0, inListSat.X0, sizeMem, cudaMemcpyHostToDevice ) );

		//cudaMemcpy( FP1, inOP.FP1, sizeMem, cudaMemcpyHostToDevice );
		//cudaMemcpy( FP2, inOP.FP2, sizeMem, cudaMemcpyHostToDevice );
		//cudaMemcpy( FP3, inOP.FP3, sizeMem, cudaMemcpyHostToDevice );
		//cudaMemcpy( FP4, inOP.FP4, sizeMem, cudaMemcpyHostToDevice );
		//cudaMemcpy( FP5, inOP.FP5, sizeMem, cudaMemcpyHostToDevice );
		//cudaMemcpy( FP6, inOP.FP6, sizeMem, cudaMemcpyHostToDevice );
		//cudaMemcpy( FP7, inOP.FP7, sizeMem, cudaMemcpyHostToDevice );
	};

	void CopyFromGPU( Orbit::SatelliteArray &outListSat )
	{
		// хранение последовательное х1 х2..... y1 y2.....
		cutilSafeCall( cudaMemcpy( outListSat.X0, X0, sizeMem, cudaMemcpyDeviceToHost ) );

		//cudaMemcpy( outOP.FP0, FP0, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP1, FP1, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP2, FP2, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP3, FP3, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP4, FP4, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP5, FP5, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP6, FP6, sizeMem, cudaMemcpyDeviceToHost );
		//cudaMemcpy( outOP.FP7, FP7, sizeMem, cudaMemcpyDeviceToHost);
	};
};

//=====================================================================//
// функция правой части на gpu
//=====================================================================//
// фуекция записывает в глобальную память и читает из глобальной
//__device__ void kernalFFxyzW(double t, double *xt, double *fx )
//{
//	double mm = 9.822*6370000.0*6370000.0;
//	int idt = blockDim.x*blockIdx.x + threadIdx.x;
//
//	double x = xt[idt + 0*CU_BlockXYZ];
//	double y = xt[idt + 1*CU_BlockXYZ];
//	double z = xt[idt + 2*CU_BlockXYZ];
//
//	double r = sqrt( x*x + y*y + z*z );
//	double IR3 = 1.0/(r*r*r);
//
//	fx[idt + 0*CU_BlockXYZ] = xt[idt + 3*CU_BlockXYZ];
//	fx[idt + 1*CU_BlockXYZ] = xt[idt + 4*CU_BlockXYZ];
//	fx[idt + 2*CU_BlockXYZ] = xt[idt + 5*CU_BlockXYZ];
//
//	fx[idt + 3*CU_BlockXYZ] = -x*mm*IR3;
//	fx[idt + 4*CU_BlockXYZ] = -y*mm*IR3;
//	fx[idt + 5*CU_BlockXYZ] = -z*mm*IR3;
//};

// фуекция записывает в глобальную память
//__device__ void kernalFFxyzRK(double t, double *xt, double *fx)
//{
//	int idt = blockDim.x*blockIdx.x + threadIdx.x;
//
//	double mm = 9.822*6370000.0*6370000.0;
//
//	double x = xt[0];
//	double y = xt[1];
//	double z = xt[2];
//
//	double vx = xt[3];
//	double vy = xt[4];
//	double vz = xt[5];
//
//	double r = sqrt( x*x + y*y + z*z );
//	double R3 = r*r*r;
//
//	double ffx = vx;
//	double ffy = vy;
//	double ffz = vz;
//	double ffvx = -x*mm/R3;
//	double ffvy = -y*mm/R3;
//	double ffvz = -z*mm/R3;
//
//	fx[idt + 0*CU_BlockXYZ] = ffx;
//	fx[idt + 1*CU_BlockXYZ] = ffy;
//	fx[idt + 2*CU_BlockXYZ] = ffz;
//
//	fx[idt + 3*CU_BlockXYZ] = ffvx;
//	fx[idt + 4*CU_BlockXYZ] = ffvy;
//	fx[idt + 5*CU_BlockXYZ] = ffvz;
//};
//// обычная функция
//__device__ void kernalFFxyz(double t, double *xt, double *fx )
//{
//	double mm = 9.822*6370000.0*6370000.0;
//
//	double x = xt[0];
//	double y = xt[1];
//	double z = xt[2];
//
//	double r = sqrt( x*x + y*y + z*z );
//	double IR3 = 1.0/(r*r*r);
//
//	fx[0] = xt[3];
//	fx[1] = xt[4];
//	fx[2] = xt[5];
//
//	fx[3] = -x*mm*IR3;
//	fx[4] = -y*mm*IR3;
//	fx[5] = -z*mm*IR3;
//};
//=====================================================================//
// вычисление возмущения
// фуекция записывает в глобальную память
//=====================================================================//
__device__ void  kernalFFxyzRK( double t, double *x, double *f, gpuOrbitPoint *OP )
{
	int idt = blockDim.x*blockIdx.x + threadIdx.x;

	double PLCOORD[11*3];

	double f_gr[3];
	double A_rot[9];

	double x_g[6];
	double f_hrm[3];
	double f_ah[3];
	

	//double S_ajd0 = 2456192.50;
	//double S_delt0 = -10.73281600;
	double S_ajd0 = OP->ajd0;
	double S_delt0 = OP->delt0;
	//double sp_q = 0.5E-05;
	//double sigma_up = 0.3E-2;
	double sp_q = OP->Ssun;
	double sigma_up =  OP->Satm;

	//A_rot[0] = 1;
	//A_rot[1] = 0;
	//A_rot[2] = 0;

	//A_rot[3] = 0;
	//A_rot[4] = 1;
	//A_rot[5] = 0;

	//A_rot[6] = 0;
	//A_rot[7] = 0;
	//A_rot[8] = 1;

	// матрица поворота
	kernaliers_mat( t, A_rot, S_ajd0, S_delt0, OP->dNUT_ARG, OP->dNUT_AMPL, OP->dT_finals_tab, OP->dT_finals_n[0] );
	//if( idt == SetOrig-1 )
	//{
	//	for( int j = 0; j < 9; j++ )
	//		printf( "%.10f ", A_rot[j] );
	//	printf( "%.10e %.10e %.10e\n", t, S_ajd0, S_delt0);
	//}
	// положение планет
	kernalplanets_update_geo( t, S_ajd0, S_delt0, OP->d_EF, PLCOORD );

	f[idt + 0*CU_BlockXYZ] = x[3];
	f[idt + 1*CU_BlockXYZ] = x[4];
	f[idt + 2*CU_BlockXYZ] = x[5];

	// гравитация планет
	kernalplanets_grav( x, f_gr, PLCOORD );

	// поворот вектора в систему координат чемли
	kernalstate_to_itrf( t, x, x_g, A_rot );

	// влияние гармоник
	kernalGetF_Harm_egm96( x_g, 75, f_hrm, OP->d_Garmonic );
	
	// солнечное давление
	double f_sp[3];
	double pln_coords11[3];
	kernalpln_coords( 10, pln_coords11, PLCOORD );
	kernalsp_cannonballForce( x, sp_q, f_sp, pln_coords11 );

	//--- влияние атмосферы ---//
	double f_atm[3];
	kernalAtm_drag( x_g, t, f_atm, sigma_up, S_ajd0, S_delt0 );
	f_hrm[0] = f_atm[0] + f_hrm[0];
	f_hrm[1] = f_atm[1] + f_hrm[1];
	f_hrm[2] = f_atm[2] + f_hrm[2];
	//------------------------//

	// перенос вектора воздействия гармоник и атмосферы
	double invA[9];
	kernaltranspose( A_rot, invA ); 
	kernalmatVecMul( invA, f_hrm, f_ah );
	
	//if( idt == SetOrig-1 )
	//{
	//	printf( "%.10e %.10e %.10e\n", f_ah[0], f_ah[1], f_ah[2] );
	//}
	// суммирование воздействий
	f[idt + 3*CU_BlockXYZ] = f_gr[0] + f_ah[0] + f_sp[0];
	f[idt + 4*CU_BlockXYZ] = f_gr[1] + f_ah[1] + f_sp[1];
	f[idt + 5*CU_BlockXYZ] = f_gr[2] + f_ah[2] + f_sp[2];
}
//=====================================================================//
// вычисление возмущения
//=====================================================================//
__device__ void  kernalFFxyz( double t, double *x, double *f, gpuOrbitPoint *OP )
{
	int idt = blockDim.x*blockIdx.x + threadIdx.x;

	double PLCOORD[11*3];
	double f_gr[3];
	double A_rot[9];

	double x_g[6];
	double f_hrm[3];
	double f_ah[3];

	//double S_ajd0 = 2456192.50;
	//double S_delt0 = -10.73281600;
	double S_ajd0 = OP->ajd0;
	double S_delt0 = OP->delt0;
	//double sp_q = 0.5E-05;
	//double sigma_up = 0.3E-2;
	double sp_q = OP->Ssun;
	double sigma_up =  OP->Satm;

	//A_rot[0] = 1;
	//A_rot[1] = 0;
	//A_rot[2] = 0;

	//A_rot[3] = 0;
	//A_rot[4] = 1;
	//A_rot[5] = 0;

	//A_rot[6] = 0;
	//A_rot[7] = 0;
	//A_rot[8] = 1;
	// матрица поворота
	kernaliers_mat( t, A_rot, S_ajd0, S_delt0, OP->dNUT_ARG, OP->dNUT_AMPL, OP->dT_finals_tab, OP->dT_finals_n[0] );

	// положение планет
	kernalplanets_update_geo( t, S_ajd0, S_delt0, OP->d_EF, PLCOORD );

	f[0] = x[3];
	f[1] = x[4];
	f[2] = x[5];

	// гравитация планет
	kernalplanets_grav( x, f_gr, PLCOORD );

	// поворот вектора в систему координат чемли
	kernalstate_to_itrf( t, x, x_g, A_rot );

	// влияние гармоник
	kernalGetF_Harm_egm96( x_g, 75, f_hrm, OP->d_Garmonic );
	
	// солнечное давление
	double f_sp[3];
	double pln_coords11[3];
	kernalpln_coords( 10, pln_coords11, PLCOORD );
	kernalsp_cannonballForce( x, sp_q, f_sp, pln_coords11 );

	//--- влияние атмосферы ---//
	double f_atm[3];
	kernalAtm_drag( x_g, t, f_atm, sigma_up, S_ajd0, S_delt0 );
	f_hrm[0] = f_atm[0] + f_hrm[0];
	f_hrm[1] = f_atm[1] + f_hrm[1];
	f_hrm[2] = f_atm[2] + f_hrm[2];
	//------------------------//

	// перенос вектора воздействия гармоник и атмосферы
	double invA[9];
	kernaltranspose( A_rot, invA ); 
	kernalmatVecMul( invA, f_hrm, f_ah );

	// суммирование воздействий
	f[3] = f_gr[0] + f_ah[0] + f_sp[0];
	f[4] = f_gr[1] + f_ah[1] + f_sp[1];
	f[5] = f_gr[2] + f_ah[2] + f_sp[2];
}

//=====================================================================//
// вычисление начальных точек методом Рунге-Кутты
//=====================================================================//
__global__ void kernalGetStartPoint( gpuOrbitPoint *OP, double h )
{
	int idt = blockDim.x*blockIdx.x + threadIdx.x;

	//int kkr[12] = { 3,4,5,6,1,2,4,3,2,1,6,5 };

	double g[6] = {	0.069431844202973712388026755553595247452137,
					0.330009478207571867598667120448377656399712,
					0.669990521792428132401332879551622343600287,
					0.930568155797026287611973244446404752547862,
					1.069431844202973712388026755553595247452137,
					1.330009478207571867598667120448377656399712 };
	int NY;
	double X0[6];
	double TP0;
	// загрузка значений
	// x1 x2 x3 ..... y1 y2 y3 ..... z1 z2 z3 .... vx . vy . vz
	// число значений на block
	X0[0] = OP->X0[idt + 0*CU_BlockXYZ];
	X0[1] = OP->X0[idt + 1*CU_BlockXYZ];
	X0[2] = OP->X0[idt + 2*CU_BlockXYZ];
	X0[3] = OP->X0[idt + 3*CU_BlockXYZ];
	X0[4] = OP->X0[idt + 4*CU_BlockXYZ];
	X0[5] = OP->X0[idt + 5*CU_BlockXYZ];
	TP0 = OP->TP0[idt];

	NY = OP->NY;
	
	//===============================================================//
	//if( kp == 1 || kp == 2 )
	// MAIN cycle of starting procedure
	// TR,HR - current time & current step
	// P  = F1+2*F2+2*F3+F4
	// fun(double t, double *xt, double *fx )
	double TPN;
	double XN[6];	// вектор состояния
	double FN1[6];	// значения функции
	double FN2[6];	// значения функции

	for (int j = 1; j <= 6; ++j)
	{
		//int kk = kkr[j - 1];				// nodes following: 3,4,5,6,1,2 
		int kk = KKR[j - 1];
		double hr = h * (1.0 - g[kk - 1]);	// step
		//kk = kkr[j + 5];
		kk = KKR[j + 5];
		double tr = TP0 + hr;				// next time step

		kernalFFxyz( TP0, X0, FN1, OP );
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		// F2=F(T0+1/2H,X0+1/2HF1) 
		for( int it = 0; it < NY; it++ )
			XN[it] = hr * FN1[it]*0.5 + X0[it];
		TPN = TP0 + hr*0.5;
		kernalFFxyz( TPN, XN, FN2, OP );
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		// F3=F(T0+1/2H,X0+1/2HF2)
		for( int it = 0; it < NY; it++ )
		{
			XN[it] = hr * FN2[it]*0.5 + X0[it];
			FN1[it] += FN2[it]*2.0;
		}
		TPN = TP0 + hr*0.5;
		kernalFFxyz( TPN, XN, FN2, OP );
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		// F4=F(T0+H,X0+HF3)
		for( int it = 0; it < NY; it++ )
		{
			XN[it] = hr * FN2[it] + X0[it];
			FN1[it] += FN2[it]*2.0;
		}
		TPN = TP0 + hr;
		kernalFFxyz( TPN, XN, FN2, OP );
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		// X1=X0+1/6H(F1+2F2+2F3+F4)
		for( int it = 0; it < NY; it++ )
			XN[it] = X0[it] +  hr / 6.0 *( FN1[it] + FN2[it] );

		// вычисление точек
		// 4,3,2,1,6,5 []
		if( j == 1)	kernalFFxyzRK(tr, XN, OP->FP4, OP );
		if( j == 2)	kernalFFxyzRK(tr, XN, OP->FP3, OP );
		if( j == 3)	kernalFFxyzRK(tr, XN, OP->FP2, OP );
		if( j == 4)	kernalFFxyzRK(tr, XN, OP->FP1, OP );
		if( j == 5)	kernalFFxyzRK(tr, XN, OP->FP6, OP );
		if( j == 6)	kernalFFxyzRK(tr, XN, OP->FP5, OP );
	}
	//=================================================================//
}
//=====================================================================//
// начальные параметры
//=====================================================================//
__global__ void kernalOrbitInitParam( gpuOrbitPoint *OP, double initStep )
{
	int idt = blockDim.x*blockIdx.x + threadIdx.x;
	
	OP->Iparam[idt].hh = initStep;
	OP->Iparam[idt].NY = 6;	
	OP->Iparam[idt].NP = 0;
	OP->Iparam[idt].IP9 = 2;

	OP->Iparam[idt].l = 0;
	OP->Iparam[idt].l1 = 0;
	OP->Iparam[idt].l2 = 0;
	OP->Iparam[idt].ltek = 0;

	OP->Iparam[idt].kkbeg = 1;
	OP->Iparam[idt].InvertStart = false;
}
//=====================================================================//
// прогноз орбиты
//=====================================================================//
__global__ void kernalOrbitPredict( gpuOrbitPoint *OP, double h, double t, double e, double **d_array, bool SavePoint ) //double *indc_a, double *indc_a1, double *indc_w, double *indc_g, double *indc_c
{
	int tx = threadIdx.x;
	int idt = blockDim.x*blockIdx.x + threadIdx.x;

	__shared__ int la[5];		// = { 0,24,48,72,96 };
	__shared__ int la1[5];		// = { 0,8,16,24,32 };
	__shared__ int la2[5];		// = { 0,24,24,24,48 };
	__shared__ double alim[2];	// = { 0.0403536069, 1.0 };
	__shared__ int isinv[5];	// = { 5,2,2,2,1 };
	__shared__ int isge[5];		// = { 4,5,5,5,3 };
	__shared__ int isle[5];		// = { 4,1,1,1,3 };
	__shared__ int iseq[5];		// = { 4,2,2,2,3 };

	// ptr to orbit for this thread
	
	
	//__shared__ double dc_a[120];
	//__shared__ double dc_a1[40];
	//__shared__ double dc_w[76];
	//__shared__ double dc_g[6];
	//__shared__ double dc_c[16];

	//__shared__ double *X0;
	//__shared__ double *FP0;
	//__shared__ double *FP1;
	//__shared__ double *FP2;
	//__shared__ double *FP3;
	//__shared__ double *FP4;
	//__shared__ double *FP5;
	//__shared__ double *FP6;
	//__shared__ double *FP7;

	double* s_array;
	if( SavePoint )
		s_array = d_array[idt];

	// инициализация констант
	if( tx == 0 )
	{
		la[0] = 0;
		la[1] = 24;
		la[2] = 48;
		la[3] = 72;
		la[4] = 96;

		la1[0] = 0;
		la1[1] = 8;
		la1[2] = 16;
		la1[3] = 24;
		la1[4] = 32;
	
		la2[0] = 0;
		la2[1] = 24;
		la2[2] = 24;
		la2[3] = 24;
		la2[4] = 48;

		alim[0] = 0.0403536069;
		alim[1] = 1.0;

		isinv[0] = 5;
		isinv[1] = 2;
		isinv[2] = 2;
		isinv[3] = 2;
		isinv[4] = 1;

		isge[0] = 4;
		isge[1] = 5;
		isge[2] = 5;
		isge[3] = 5;
		isge[4] = 3;

		isle[0] = 4;
		isle[1] = 1;
		isle[2] = 1;
		isle[3] = 1;
		isle[4] = 3;

		iseq[0] = 4;
		iseq[1] = 2;
		iseq[2] = 2;
		iseq[3] = 2;
		iseq[4] = 3;

		//X0 = OP->X0;
		//FP0 = OP->FP0;
		//FP1 = OP->FP1;
		//FP2 = OP->FP2;
		//FP3 = OP->FP3;
		//FP4 = OP->FP4;
		//FP5 = OP->FP5;
		//FP6 = OP->FP6;
		//FP7 = OP->FP7;
	}

	// загрузка коэффициентов
	//if( tx < 120 )
	//	dc_a[tx] = indc_a[tx];
	//if( tx < 40 )
	//	dc_a1[tx] = indc_a1[tx];
	//if( tx < 76 )
	//	dc_w[tx] = indc_w[tx];
	//if( tx < 6 )
	//	dc_g[tx] = indc_g[tx];
	//if( tx < 16 )
	//	dc_c[tx] = indc_c[tx];

	__syncthreads ();

	//// локальные переменные
	//// их слишком много
	double X0[6];
	double FP0[6];
	double FP1[6];
	double FP2[6];
	double FP3[6];
	double FP4[6];
	double FP5[6];
	double FP6[6];
	double FP7[6];
	
	//// загружаем данные для каждого потока
	//// в целом места всеравно нет, так что надо работать с глобальной памятью
	//// позже переделаем
	for( int it = 0; it < 6; it++ )
	{
		X0[it] = OP->X0[idt + it*CU_BlockXYZ];
		FP0[it] = OP->FP0[idt + it*CU_BlockXYZ];
		FP1[it] = OP->FP1[idt + it*CU_BlockXYZ];
		FP2[it] = OP->FP2[idt + it*CU_BlockXYZ];
		FP3[it] = OP->FP3[idt + it*CU_BlockXYZ];
		FP4[it] = OP->FP4[idt + it*CU_BlockXYZ];
		FP5[it] = OP->FP5[idt + it*CU_BlockXYZ];
		FP6[it] = OP->FP6[idt + it*CU_BlockXYZ];
		FP7[it] = OP->FP7[idt + it*CU_BlockXYZ];
	}

	double TP0 = OP->TP0[idt];
	double TPN = 0;
	//===============================================================//
	// параметры для сохранения
	//double hh = h;
	//int NY = 6;	
	//int NP = 0;
	//int IP9 = 2;

	//int l = 0;
	//int l1 = 0;
	//int l2 = 0;
	//int ltek = 0;

	//int kkbeg = 1;
	//double delt;
	//double rr1;
	//double rr2;
	//double rotn;
	//double dd1;
	//bool InvertStart = false;

	double hh = OP->Iparam[idt].hh;
	int NY = OP->Iparam[idt].NY;	
	int NP = OP->Iparam[idt].NP;
	int IP9 = OP->Iparam[idt].IP9;

	int l = OP->Iparam[idt].l;
	int l1 = OP->Iparam[idt].l1;
	int l2 = OP->Iparam[idt].l2;
	int ltek = OP->Iparam[idt].ltek;

	int kkbeg = OP->Iparam[idt].kkbeg;
	double delt = OP->Iparam[idt].delt;
	double rr1 = OP->Iparam[idt].rr1;
	double rr2 = OP->Iparam[idt].rr2;
	double rotn = OP->Iparam[idt].rotn;
	double dd1 = OP->Iparam[idt].dd1;
	bool InvertStart = OP->Iparam[idt].InvertStart;
	//===============================================================//

	//!!!!! чтобы работало
	idt = 0;
	int BlockXYZi = 1;
	//===============================================================//

	//------------------------------------------------//
	if( SavePoint )
	{
		// запись начального положения
		int itwrite = 1;
		double w_t = TP0;
		double w_x = X0[0];
		double w_y = X0[1];
		double w_z = X0[2];
		s_array[ itwrite ] = w_t;
		s_array[ itwrite + 1] = w_x;
		s_array[ itwrite + 2] = w_y;
		s_array[ itwrite + 3] = w_z;
	}
	//------------------------------------------------//

	if (hh * (t - TP0) < 0.0)
	{ 
		//printf("direct change\n");
		// direct change.

		for( int it = 0; it < NY; it++ )
		{
			int itt = idt + it*BlockXYZi;

			double tmp = FP1[itt];
			FP1[itt] = FP6[itt];
			FP6[itt] = tmp;

			tmp = FP2[itt];
			FP2[itt] = FP5[itt];
			FP5[itt] = tmp;

			tmp = FP3[itt];
			FP3[itt] = FP4[itt];
			FP4[itt] = tmp;
		}

		// step change
		hh = -hh;
		if (IP9 == 1) {	hh /= 0.7;	}
		if (IP9 == 5) {	hh *= 0.7;	}

		// case nuber cange 
		IP9 = isinv[IP9 - 1];

		//  2 steps of extrapolation
		kkbeg = 3;

		//goto L1400;
		InvertStart = true;
	}
	//===============================================================//

	int isstep = 0;
	// main cyrcle
	while( 1 )
	{
		if( InvertStart == false )
		{
			//===============================================================//
			//  S4: time overflow 
			if (hh * (t - TP0 - hh) <= 0.0 )
				break; 
			
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if( isstep > 2990 )
				break;

			// if overflow absent => step of implicit RK-method 
			// A5: step of implicit RK-method + error estimation
			TP0 +=hh;

			// DELT = maximum error
			// RR1,RR2 - working values
			// ROTH - error for i-th equation
			delt = 0.0;
			for (int in = 0; in < NY; in++ )
			{
				int itt = idt + in*BlockXYZi;
				//rr2 = FP3[itt]*dc_w[72] + FP4[itt]*dc_w[73] + FP5[itt]*dc_w[74] + FP6[itt]*dc_w[75];
				rr2 = FP3[itt]*Dcw[72] + FP4[itt]*Dcw[73] + FP5[itt]*Dcw[74] + FP6[itt]*Dcw[75];

				//step counter starting setting NP = 0;
				if (NP == 0 || kkbeg == 3 ){	
					FP0[itt] +=  hh * rr2;
					continue;
				}

				//rr1 = FP7[itt] + FP2[itt]*dc_a1[l1 + 4] + FP3[itt]*dc_a1[l1 + 5] + FP4[itt]*dc_a1[l1 + 6] + FP5[itt]*dc_a1[l1 + 7];
				rr1 = FP7[itt] + FP2[itt]*Dca1[l1 + 4] + FP3[itt]*Dca1[l1 + 5] + FP4[itt]*Dca1[l1 + 6] + FP5[itt]*Dca1[l1 + 7];
				dd1 = hh * (rr1 - rr2) / e;
				rotn = fabs(dd1);
				if (rotn >= delt){	delt = rotn;}

				FP0[itt] +=  hh * rr2;
			}
			NP++;
			//===============================================================//

			//===============================================================//
			//STEP: increase, decrease or do not change? 
			if ( NP == 1 || kkbeg == 3 ){ 
				IP9 = iseq[IP9 - 1];
				kkbeg = 1;
			}
			else if(alim[0] >= delt){	
				// A6:step increase 
				if (IP9 == 1) {	hh = hh; }
				if (IP9 == 2) {	hh /= 0.7;	}
				if (IP9 == 3) {	hh /= 0.7;	}
				if (IP9 == 4) {	hh /= 0.7;	}
				if (IP9 == 5) {	hh = hh;	}
				IP9 = isge[IP9 - 1];
			}
			else if (alim[1] <= delt) {	
				// A7:step decrease 
				if (IP9 == 1) {	hh = hh; }
				if (IP9 == 2) {	hh *= 0.7;	}
				if (IP9 == 3) {	hh *= 0.7;	}
				if (IP9 == 4) {	hh *= 0.7;	}
				if (IP9 == 5) {	hh = hh;	}
				IP9 = isle[IP9 - 1];
			}
			else
			{
				IP9 = iseq[IP9 - 1];
				kkbeg = 1;
			}
			//===============================================================//
		}

		InvertStart = false;
		//===============================================================//
		// A4: extrapolation 2|4 points + interpolation 
		l = la[IP9 - 1];
		l1 = la1[IP9 - 1];
		l2 = la2[IP9 - 1];
		//===============================================================//

		//===============================================================//
		for(int ii = 0; ii < NY; ii++ ) 
			FP7[idt + ii*BlockXYZi ] = 0.0;


		//  Main cycle of extrapolation 
		// FP0 + hh( a1*FP1 + a2*FP2 + a3*FP3 + a4*FP4 + a5*FP5  + a6*FP6 )
		for (int kk = kkbeg; kk <= 4; ++kk) 
		{
			ltek = l + kk * 6 - 6;
			// время в первой точке FP1
			//TPN = TP0 + dc_g[kk - 1] * hh;
			TPN = TP0 + Dcg[kk - 1] * hh;

			for(int in = 0; in < NY; in++ ) 
			{
				int itt = idt + in*BlockXYZi;
				//X0[itt] = FP0[itt] + hh * (dc_a[ltek] * FP1[itt] + dc_a[ltek + 1] * FP2[itt] + dc_a[ltek + 2] * FP3[itt] + dc_a[ltek + 3] * FP4[itt] + dc_a[ltek + 4]* FP5[itt] + dc_a[ltek + 5] *FP6[itt]);
				X0[itt] = FP0[itt] + hh * (Dca[ltek] * FP1[itt] + Dca[ltek + 1] * FP2[itt] + Dca[ltek + 2] * FP3[itt] + Dca[ltek + 3] * FP4[itt] + Dca[ltek + 4]* FP5[itt] + Dca[ltek + 5] *FP6[itt]);

				// accumulate values for error estimation
				//FP7[itt] += FP2[itt] * dc_a1[l1 + kk - 1];
				FP7[itt] += FP2[itt] * Dca1[l1 + kk - 1];
			}
			kernalFFxyz( TPN, X0, FP1, OP );

			// offset point
			for( int in = 0; in < NY; in++ )
			{
				int itt = idt + in*BlockXYZi;
				double tmp = FP1[itt];
				FP1[itt] = FP2[itt];
				FP2[itt] = FP3[itt];
				FP3[itt] = FP4[itt];
				FP4[itt] = FP5[itt];
				FP5[itt] = FP6[itt];
				FP6[itt] = tmp;
			}
		}
		//===============================================================//

		//===============================================================//
		// Main cycle of interpolation
		// X0 = FP0 + hh*( w1*FP1 + w2*FP2 + w3*FP3 + w4*FP4 + w5*FP5 + w6*FP6 )
		for (int kk = 1; kk <= 4; ++kk)
		{
			for ( int in = 0; in < NY; in++ )
			{
				int itt = idt + in*BlockXYZi;
				//X0[itt] =  FP0[itt] + hh * (dc_w[l2]*FP1[itt] + dc_w[l2 + 1]*FP2[itt] + dc_w[l2 + 2]*FP3[itt] + dc_w[l2 + 3]*FP4[itt] + dc_w[l2 + 4]*FP5[itt] + dc_w[l2 + 5]*FP6[itt] );
				X0[itt] =  FP0[itt] + hh * (Dcw[l2]*FP1[itt] + Dcw[l2 + 1]*FP2[itt] + Dcw[l2 + 2]*FP3[itt] + Dcw[l2 + 3]*FP4[itt] + Dcw[l2 + 4]*FP5[itt] + Dcw[l2 + 5]*FP6[itt] );
			}
			// new point for FP3 .... FP6
			//TPN = TP0 + dc_g[kk - 1] * hh;
			TPN = TP0 + Dcg[kk - 1] * hh;

			if( kk == 1) kernalFFxyz( TPN, X0, FP3, OP );
			if( kk == 2) kernalFFxyz( TPN, X0, FP4, OP );
			if( kk == 3) kernalFFxyz( TPN, X0, FP5, OP );
			if( kk == 4) kernalFFxyz( TPN, X0, FP6, OP );
			l2 += 6;
		}
		//===============================================================//

		//------------------------------------------------//
		if( SavePoint )
		{
			isstep++;
			// позиция для записи начиная со второго места плюс 1
			int itwrite = isstep*4 + 1;

			// запись значений
			double w_t = TP0;
			double w_x = FP0[idt + 0*BlockXYZi ];
			double w_y = FP0[idt + 1*BlockXYZi ];
			double w_z = FP0[idt + 2*BlockXYZi ];

			s_array[ itwrite ] = w_t;
			s_array[ itwrite + 1] = w_x;
			s_array[ itwrite + 2] = w_y;
			s_array[ itwrite + 3] = w_z;
		}
		//------------------------------------------------//
	}

	//===============================================================//
	// A8: interpolatin at the destination time 
	delt = t - TP0;
	for ( int in = 0; in < NY; in++ ) 
	{
		int itt = idt + in*BlockXYZi;
		rr1 = 0.0;
		for (int j = 1; j <= 13; j += 4) 
		{
			//rr1 = delt / hh * (dc_c[j - 1] * FP3[itt] + dc_c[j] * FP4[itt] + dc_c[j + 1] * FP5[itt] + dc_c[j + 2] * FP6[itt] + rr1);
			rr1 = delt / hh * (Dcc[j - 1] * FP3[itt] + Dcc[j] * FP4[itt] + Dcc[j + 1] * FP5[itt] + Dcc[j + 2] * FP6[itt] + rr1);

			X0[itt] = FP0[itt] + rr1 * hh;
		}
	}

	//===============================================================//
	idt = blockDim.x*blockIdx.x + threadIdx.x;
	// возвращаем шаг интегрирования
	OP->Iparam[idt].hh = hh;
	OP->Iparam[idt].NY = NY;	
	OP->Iparam[idt].NP = NP;
	OP->Iparam[idt].IP9 = IP9;

	OP->Iparam[idt].l = l;
	OP->Iparam[idt].l1 = l1;
	OP->Iparam[idt].l2 = l2;
	OP->Iparam[idt].ltek = ltek;

	OP->Iparam[idt].kkbeg = kkbeg;
	OP->Iparam[idt].delt = delt;
	OP->Iparam[idt].rr1 = rr1;
	OP->Iparam[idt].rr2 = rr2;
	OP->Iparam[idt].rotn = rotn;
	OP->Iparam[idt].dd1 = dd1;
	OP->Iparam[idt].InvertStart = InvertStart;

	// начальное время, это текущее законченное
	OP->TP0[idt] = TP0;
	for( int it = 0; it < 6; it++ )
	{
		OP->X0[idt + it*CU_BlockXYZ] = X0[it];
		OP->FP0[idt + it*CU_BlockXYZ] = FP0[it];
		OP->FP1[idt + it*CU_BlockXYZ] = FP1[it];
		OP->FP2[idt + it*CU_BlockXYZ] = FP2[it];
		OP->FP3[idt + it*CU_BlockXYZ] = FP3[it];
		OP->FP4[idt + it*CU_BlockXYZ] = FP4[it];
		OP->FP5[idt + it*CU_BlockXYZ] = FP5[it];
		OP->FP6[idt + it*CU_BlockXYZ] = FP6[it];
		OP->FP7[idt + it*CU_BlockXYZ] = FP7[it];
	}
	//===============================================================//

	//------------------------------------------------//
	if( SavePoint )
	{
		// write numb points in array
		s_array[ 0 ] = isstep + 1;
	}
	//------------------------------------------------//
	for( int it = 0; it < 6; it++ )
		OP->X0[idt + it*CU_BlockXYZ ] = X0[it];
}
//##############################################################################//
__global__ void kernalOrbitPredict_step( gpuOrbitPoint *OP, double h, double tstart_step, double t, double e, double **d_array, bool SavePoint ) //double *indc_a, double *indc_a1, double *indc_w, double *indc_g, double *indc_c
{
	int tx = threadIdx.x;
	int idt = blockDim.x*blockIdx.x + threadIdx.x;

	__shared__ int la[5];		// = { 0,24,48,72,96 };
	__shared__ int la1[5];		// = { 0,8,16,24,32 };
	__shared__ int la2[5];		// = { 0,24,24,24,48 };
	__shared__ double alim[2];	// = { 0.0403536069, 1.0 };
	__shared__ int isinv[5];	// = { 5,2,2,2,1 };
	__shared__ int isge[5];		// = { 4,5,5,5,3 };
	__shared__ int isle[5];		// = { 4,1,1,1,3 };
	__shared__ int iseq[5];		// = { 4,2,2,2,3 };

	// ptr to orbit for this thread
	
	
	//__shared__ double dc_a[120];
	//__shared__ double dc_a1[40];
	//__shared__ double dc_w[76];
	//__shared__ double dc_g[6];
	//__shared__ double dc_c[16];

	//__shared__ double *X0;
	//__shared__ double *FP0;
	//__shared__ double *FP1;
	//__shared__ double *FP2;
	//__shared__ double *FP3;
	//__shared__ double *FP4;
	//__shared__ double *FP5;
	//__shared__ double *FP6;
	//__shared__ double *FP7;

	double* s_array;
	if( SavePoint )
		s_array = d_array[idt];

	// инициализация констант
	if( tx == 0 )
	{
		la[0] = 0;
		la[1] = 24;
		la[2] = 48;
		la[3] = 72;
		la[4] = 96;

		la1[0] = 0;
		la1[1] = 8;
		la1[2] = 16;
		la1[3] = 24;
		la1[4] = 32;
	
		la2[0] = 0;
		la2[1] = 24;
		la2[2] = 24;
		la2[3] = 24;
		la2[4] = 48;

		alim[0] = 0.0403536069;
		alim[1] = 1.0;

		isinv[0] = 5;
		isinv[1] = 2;
		isinv[2] = 2;
		isinv[3] = 2;
		isinv[4] = 1;

		isge[0] = 4;
		isge[1] = 5;
		isge[2] = 5;
		isge[3] = 5;
		isge[4] = 3;

		isle[0] = 4;
		isle[1] = 1;
		isle[2] = 1;
		isle[3] = 1;
		isle[4] = 3;

		iseq[0] = 4;
		iseq[1] = 2;
		iseq[2] = 2;
		iseq[3] = 2;
		iseq[4] = 3;

		//X0 = OP->X0;
		//FP0 = OP->FP0;
		//FP1 = OP->FP1;
		//FP2 = OP->FP2;
		//FP3 = OP->FP3;
		//FP4 = OP->FP4;
		//FP5 = OP->FP5;
		//FP6 = OP->FP6;
		//FP7 = OP->FP7;
	}

	// загрузка коэффициентов
	//if( tx < 120 )
	//	dc_a[tx] = indc_a[tx];
	//if( tx < 40 )
	//	dc_a1[tx] = indc_a1[tx];
	//if( tx < 76 )
	//	dc_w[tx] = indc_w[tx];
	//if( tx < 6 )
	//	dc_g[tx] = indc_g[tx];
	//if( tx < 16 )
	//	dc_c[tx] = indc_c[tx];

	__syncthreads ();

	//// локальные переменные
	//// их слишком много
	double X0[6];
	double FP0[6];
	double FP1[6];
	double FP2[6];
	double FP3[6];
	double FP4[6];
	double FP5[6];
	double FP6[6];
	double FP7[6];

	double Xres[6];
	
	//// загружаем данные для каждого потока
	//// в целом места всеравно нет, так что надо работать с глобальной памятью
	//// позже переделаем
	for( int it = 0; it < 6; it++ )
	{
		X0[it] = OP->X0[idt + it*CU_BlockXYZ];
		FP0[it] = OP->FP0[idt + it*CU_BlockXYZ];
		FP1[it] = OP->FP1[idt + it*CU_BlockXYZ];
		FP2[it] = OP->FP2[idt + it*CU_BlockXYZ];
		FP3[it] = OP->FP3[idt + it*CU_BlockXYZ];
		FP4[it] = OP->FP4[idt + it*CU_BlockXYZ];
		FP5[it] = OP->FP5[idt + it*CU_BlockXYZ];
		FP6[it] = OP->FP6[idt + it*CU_BlockXYZ];
		FP7[it] = OP->FP7[idt + it*CU_BlockXYZ];
	}

	double TP0 = OP->TP0[idt];
	double TPN = 0;
	//===============================================================//
	// параметры для сохранения
	//double hh = h;
	//int NY = 6;	
	//int NP = 0;
	//int IP9 = 2;

	//int l = 0;
	//int l1 = 0;
	//int l2 = 0;
	//int ltek = 0;

	//int kkbeg = 1;
	//double delt;
	//double rr1;
	//double rr2;
	//double rotn;
	//double dd1;
	//bool InvertStart = false;

	double hh = OP->Iparam[idt].hh;
	int NY = OP->Iparam[idt].NY;	
	int NP = OP->Iparam[idt].NP;
	int IP9 = OP->Iparam[idt].IP9;

	int l = OP->Iparam[idt].l;
	int l1 = OP->Iparam[idt].l1;
	int l2 = OP->Iparam[idt].l2;
	int ltek = OP->Iparam[idt].ltek;

	int kkbeg = OP->Iparam[idt].kkbeg;
	double delt = OP->Iparam[idt].delt;
	double rr1 = OP->Iparam[idt].rr1;
	double rr2 = OP->Iparam[idt].rr2;
	double rotn = OP->Iparam[idt].rotn;
	double dd1 = OP->Iparam[idt].dd1;
	bool InvertStart = OP->Iparam[idt].InvertStart;

	// new t_steppt
	double t_steppt = tstart_step;//TP0; // время старта записи точек
	int itwrites = 0;
	double PointTStep = STEPDTIME/1000.0;
	//===============================================================//

	//!!!!! чтобы работало
	idt = 0;
	int BlockXYZi = 1;
	//===============================================================//

	//------------------------------------------------//
	//if( SavePoint )
	//{
	//	// запись начального положения
	//	int itwrite = 1;
	//	double w_t = TP0;
	//	double w_x = X0[0];
	//	double w_y = X0[1];
	//	double w_z = X0[2];
	//	s_array[ itwrite ] = w_t;
	//	s_array[ itwrite + 1] = w_x;
	//	s_array[ itwrite + 2] = w_y;
	//	s_array[ itwrite + 3] = w_z;
	//}
	//------------------------------------------------//

	if (hh * (t - TP0) < 0.0)
	{ 
		//printf("direct change\n");
		// direct change.

		for( int it = 0; it < NY; it++ )
		{
			int itt = idt + it*BlockXYZi;

			double tmp = FP1[itt];
			FP1[itt] = FP6[itt];
			FP6[itt] = tmp;

			tmp = FP2[itt];
			FP2[itt] = FP5[itt];
			FP5[itt] = tmp;

			tmp = FP3[itt];
			FP3[itt] = FP4[itt];
			FP4[itt] = tmp;
		}

		// step change
		hh = -hh;
		if (IP9 == 1) {	hh /= 0.7;	}
		if (IP9 == 5) {	hh *= 0.7;	}

		// case nuber cange 
		IP9 = isinv[IP9 - 1];

		//  2 steps of extrapolation
		kkbeg = 3;

		//goto L1400;
		InvertStart = true;
	}
	//===============================================================//

	//int isstep = 0;
	// main cyrcle
	while( 1 )
	{
		if( InvertStart == false )
		{
			//===============================================================//
			//  S4: time overflow 
			if (hh * (t - TP0 - hh) <= 0.0 )
				break; 

			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
			// Add new 
			if( SavePoint )
			{
				while(hh * ( t_steppt - TP0 - hh) <= 0.0 )
				{
					// A8: interpolatin at the destination time 
					double delts = t_steppt - TP0;
					
					//for ( int ii = 0; ii < NY; ii++ ) 
					//{
					//	rr1 = 0.0;
					//	for (int j = 1; j <= 13; j += 4) 
					//	{
					//		rr1 = delts / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
					//		Xres[ii] = FP0[ii] + rr1 * hh;
					//	}
					//}

					for ( int in = 0; in < NY; in++ ) 
					{
						int itt = idt + in*BlockXYZi;
						rr1 = 0.0;
						for (int j = 1; j <= 13; j += 4) 
						{
							//rr1 = delt / hh * (dc_c[j - 1] * FP3[itt] + dc_c[j] * FP4[itt] + dc_c[j + 1] * FP5[itt] + dc_c[j + 2] * FP6[itt] + rr1);
							rr1 = delts / hh * (Dcc[j - 1] * FP3[itt] + Dcc[j] * FP4[itt] + Dcc[j + 1] * FP5[itt] + Dcc[j + 2] * FP6[itt] + rr1);
							Xres[in] = FP0[itt] + rr1 * hh;
						}
					}
					//FILE *fre = fopen( "pt.log", "at" );
					//for ( int ii = 0; ii < NY; ii++ ) 
					//	fprintf( fre, "%f\t", Xres[ii] );
					//fprintf( fre, "%f\n", t_steppt );
					//fclose( fre );

					int iw = itwrites*4 + 1;
					s_array[ iw  ] = t_steppt;
					s_array[ iw  + 1] = Xres[0];
					s_array[ iw  + 2] = Xres[1];
					s_array[ iw  + 3] = Xres[2];
					itwrites++;
					
					// next pt
					t_steppt +=  PointTStep;
				}
			}
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
			
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//if( isstep > 2990 )
			//	break;

			// if overflow absent => step of implicit RK-method 
			// A5: step of implicit RK-method + error estimation
			TP0 +=hh;

			// DELT = maximum error
			// RR1,RR2 - working values
			// ROTH - error for i-th equation
			delt = 0.0;
			for (int in = 0; in < NY; in++ )
			{
				int itt = idt + in*BlockXYZi;
				//rr2 = FP3[itt]*dc_w[72] + FP4[itt]*dc_w[73] + FP5[itt]*dc_w[74] + FP6[itt]*dc_w[75];
				rr2 = FP3[itt]*Dcw[72] + FP4[itt]*Dcw[73] + FP5[itt]*Dcw[74] + FP6[itt]*Dcw[75];

				//step counter starting setting NP = 0;
				if (NP == 0 || kkbeg == 3 ){	
					FP0[itt] +=  hh * rr2;
					continue;
				}

				//rr1 = FP7[itt] + FP2[itt]*dc_a1[l1 + 4] + FP3[itt]*dc_a1[l1 + 5] + FP4[itt]*dc_a1[l1 + 6] + FP5[itt]*dc_a1[l1 + 7];
				rr1 = FP7[itt] + FP2[itt]*Dca1[l1 + 4] + FP3[itt]*Dca1[l1 + 5] + FP4[itt]*Dca1[l1 + 6] + FP5[itt]*Dca1[l1 + 7];
				dd1 = hh * (rr1 - rr2) / e;
				rotn = fabs(dd1);
				if (rotn >= delt){	delt = rotn;}

				FP0[itt] +=  hh * rr2;
			}
			NP++;
			//===============================================================//

			//===============================================================//
			//STEP: increase, decrease or do not change? 
			if ( NP == 1 || kkbeg == 3 ){ 
				IP9 = iseq[IP9 - 1];
				kkbeg = 1;
			}
			else if(alim[0] >= delt){	
				// A6:step increase 
				if (IP9 == 1) {	hh = hh; }
				if (IP9 == 2) {	hh /= 0.7;	}
				if (IP9 == 3) {	hh /= 0.7;	}
				if (IP9 == 4) {	hh /= 0.7;	}
				if (IP9 == 5) {	hh = hh;	}
				IP9 = isge[IP9 - 1];
			}
			else if (alim[1] <= delt) {	
				// A7:step decrease 
				if (IP9 == 1) {	hh = hh; }
				if (IP9 == 2) {	hh *= 0.7;	}
				if (IP9 == 3) {	hh *= 0.7;	}
				if (IP9 == 4) {	hh *= 0.7;	}
				if (IP9 == 5) {	hh = hh;	}
				IP9 = isle[IP9 - 1];
			}
			else
			{
				IP9 = iseq[IP9 - 1];
				kkbeg = 1;
			}
			//===============================================================//
		}

		InvertStart = false;
		//===============================================================//
		// A4: extrapolation 2|4 points + interpolation 
		l = la[IP9 - 1];
		l1 = la1[IP9 - 1];
		l2 = la2[IP9 - 1];
		//===============================================================//

		//===============================================================//
		for(int ii = 0; ii < NY; ii++ ) 
			FP7[idt + ii*BlockXYZi ] = 0.0;


		//  Main cycle of extrapolation 
		// FP0 + hh( a1*FP1 + a2*FP2 + a3*FP3 + a4*FP4 + a5*FP5  + a6*FP6 )
		for (int kk = kkbeg; kk <= 4; ++kk) 
		{
			ltek = l + kk * 6 - 6;
			// время в первой точке FP1
			//TPN = TP0 + dc_g[kk - 1] * hh;
			TPN = TP0 + Dcg[kk - 1] * hh;

			for(int in = 0; in < NY; in++ ) 
			{
				int itt = idt + in*BlockXYZi;
				//X0[itt] = FP0[itt] + hh * (dc_a[ltek] * FP1[itt] + dc_a[ltek + 1] * FP2[itt] + dc_a[ltek + 2] * FP3[itt] + dc_a[ltek + 3] * FP4[itt] + dc_a[ltek + 4]* FP5[itt] + dc_a[ltek + 5] *FP6[itt]);
				X0[itt] = FP0[itt] + hh * (Dca[ltek] * FP1[itt] + Dca[ltek + 1] * FP2[itt] + Dca[ltek + 2] * FP3[itt] + Dca[ltek + 3] * FP4[itt] + Dca[ltek + 4]* FP5[itt] + Dca[ltek + 5] *FP6[itt]);

				// accumulate values for error estimation
				//FP7[itt] += FP2[itt] * dc_a1[l1 + kk - 1];
				FP7[itt] += FP2[itt] * Dca1[l1 + kk - 1];
			}
			kernalFFxyz( TPN, X0, FP1, OP );

			// offset point
			for( int in = 0; in < NY; in++ )
			{
				int itt = idt + in*BlockXYZi;
				double tmp = FP1[itt];
				FP1[itt] = FP2[itt];
				FP2[itt] = FP3[itt];
				FP3[itt] = FP4[itt];
				FP4[itt] = FP5[itt];
				FP5[itt] = FP6[itt];
				FP6[itt] = tmp;
			}
		}
		//===============================================================//

		//===============================================================//
		// Main cycle of interpolation
		// X0 = FP0 + hh*( w1*FP1 + w2*FP2 + w3*FP3 + w4*FP4 + w5*FP5 + w6*FP6 )
		for (int kk = 1; kk <= 4; ++kk)
		{
			for ( int in = 0; in < NY; in++ )
			{
				int itt = idt + in*BlockXYZi;
				//X0[itt] =  FP0[itt] + hh * (dc_w[l2]*FP1[itt] + dc_w[l2 + 1]*FP2[itt] + dc_w[l2 + 2]*FP3[itt] + dc_w[l2 + 3]*FP4[itt] + dc_w[l2 + 4]*FP5[itt] + dc_w[l2 + 5]*FP6[itt] );
				X0[itt] =  FP0[itt] + hh * (Dcw[l2]*FP1[itt] + Dcw[l2 + 1]*FP2[itt] + Dcw[l2 + 2]*FP3[itt] + Dcw[l2 + 3]*FP4[itt] + Dcw[l2 + 4]*FP5[itt] + Dcw[l2 + 5]*FP6[itt] );
			}
			// new point for FP3 .... FP6
			//TPN = TP0 + dc_g[kk - 1] * hh;
			TPN = TP0 + Dcg[kk - 1] * hh;

			if( kk == 1) kernalFFxyz( TPN, X0, FP3, OP );
			if( kk == 2) kernalFFxyz( TPN, X0, FP4, OP );
			if( kk == 3) kernalFFxyz( TPN, X0, FP5, OP );
			if( kk == 4) kernalFFxyz( TPN, X0, FP6, OP );
			l2 += 6;
		}
		//===============================================================//

		//------------------------------------------------//
		//if( SavePoint )
		//{
		//	isstep++;
		//	// позиция для записи начиная со второго места плюс 1
		//	int itwrite = isstep*4 + 1;

		//	// запись значений
		//	double w_t = TP0;
		//	double w_x = FP0[idt + 0*BlockXYZi ];
		//	double w_y = FP0[idt + 1*BlockXYZi ];
		//	double w_z = FP0[idt + 2*BlockXYZi ];

		//	s_array[ itwrite ] = w_t;
		//	s_array[ itwrite + 1] = w_x;
		//	s_array[ itwrite + 2] = w_y;
		//	s_array[ itwrite + 3] = w_z;
		//}
		//------------------------------------------------//
	}

	//===============================================================//
	// A8: interpolatin at the destination time 
	delt = t - TP0;
	for ( int in = 0; in < NY; in++ ) 
	{
		int itt = idt + in*BlockXYZi;
		rr1 = 0.0;
		for (int j = 1; j <= 13; j += 4) 
		{
			//rr1 = delt / hh * (dc_c[j - 1] * FP3[itt] + dc_c[j] * FP4[itt] + dc_c[j + 1] * FP5[itt] + dc_c[j + 2] * FP6[itt] + rr1);
			rr1 = delt / hh * (Dcc[j - 1] * FP3[itt] + Dcc[j] * FP4[itt] + Dcc[j + 1] * FP5[itt] + Dcc[j + 2] * FP6[itt] + rr1);

			X0[itt] = FP0[itt] + rr1 * hh;
		}
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	// Add new 
	if( SavePoint )
	{
		while( hh*( t_steppt - t) <= 0.0 )
		{
			// A8: interpolatin at the destination time 
			double delts = t_steppt - TP0;

			//for ( int ii = 0; ii < NY; ii++ ) 
			//{
			//	rr1 = 0.0;
			//	for (int j = 1; j <= 13; j += 4) 
			//	{
			//		rr1 = delts / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
			//		Xres[ii] = FP0[ii] + rr1 * hh;
			//	}
			//}

			for ( int in = 0; in < NY; in++ ) 
			{
				int itt = idt + in*BlockXYZi;
				rr1 = 0.0;
				for (int j = 1; j <= 13; j += 4) 
				{
					//rr1 = delt / hh * (dc_c[j - 1] * FP3[itt] + dc_c[j] * FP4[itt] + dc_c[j + 1] * FP5[itt] + dc_c[j + 2] * FP6[itt] + rr1);
					rr1 = delts / hh * (Dcc[j - 1] * FP3[itt] + Dcc[j] * FP4[itt] + Dcc[j + 1] * FP5[itt] + Dcc[j + 2] * FP6[itt] + rr1);

					Xres[in] = FP0[itt] + rr1 * hh;
				}
			}

			//FILE *fre = fopen( "pt.log", "at" );
			//for ( int ii = 0; ii < NY; ii++ ) 
			//	fprintf( fre, "%f\t", Xres[ii] );
			//fprintf( fre, "%f\n", t_steppt );
			//fclose( fre );

			int iw = itwrites*4 + 1;
			s_array[ iw ] = t_steppt;
			s_array[ iw + 1] = Xres[0];
			s_array[ iw + 2] = Xres[1];
			s_array[ iw + 3] = Xres[2];
			itwrites++;

			// next pt
			t_steppt += PointTStep;
		}
	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	__syncthreads ();

	//===============================================================//
	idt = blockDim.x*blockIdx.x + threadIdx.x;
	// возвращаем шаг интегрирования
	OP->Iparam[idt].hh = hh;
	OP->Iparam[idt].NY = NY;	
	OP->Iparam[idt].NP = NP;
	OP->Iparam[idt].IP9 = IP9;

	OP->Iparam[idt].l = l;
	OP->Iparam[idt].l1 = l1;
	OP->Iparam[idt].l2 = l2;
	OP->Iparam[idt].ltek = ltek;

	OP->Iparam[idt].kkbeg = kkbeg;
	OP->Iparam[idt].delt = delt;
	OP->Iparam[idt].rr1 = rr1;
	OP->Iparam[idt].rr2 = rr2;
	OP->Iparam[idt].rotn = rotn;
	OP->Iparam[idt].dd1 = dd1;
	OP->Iparam[idt].InvertStart = InvertStart;

	// начальное время, это текущее законченное
	OP->TP0[idt] = TP0;
	for( int it = 0; it < 6; it++ )
	{
		OP->X0[idt + it*CU_BlockXYZ] = X0[it];
		OP->FP0[idt + it*CU_BlockXYZ] = FP0[it];
		OP->FP1[idt + it*CU_BlockXYZ] = FP1[it];
		OP->FP2[idt + it*CU_BlockXYZ] = FP2[it];
		OP->FP3[idt + it*CU_BlockXYZ] = FP3[it];
		OP->FP4[idt + it*CU_BlockXYZ] = FP4[it];
		OP->FP5[idt + it*CU_BlockXYZ] = FP5[it];
		OP->FP6[idt + it*CU_BlockXYZ] = FP6[it];
		OP->FP7[idt + it*CU_BlockXYZ] = FP7[it];
	}
	//===============================================================//

	//------------------------------------------------//
	if( SavePoint )
	{
		// write numb points in array
		s_array[ 0 ] = itwrites;
	}
	//------------------------------------------------//
	for( int it = 0; it < 6; it++ )
		OP->X0[idt + it*CU_BlockXYZ ] = X0[it];
}
//##############################################################################//

//==============================================================================//
// функция интерполяции по логранжу
//==============================================================================//
__device__ double cuGetInterpolatePoint( double *x, double *y, int np, double xf )
{
	double res = 0.0;
	// ссумируем функции лагранжа
	for( int it = 0; it < np; it++ )
	{
		double tmp = 1.0;
		// множители функции
		for( int j = 0; j < np; j++ )
		{
			// пропускаем множитель
			if( j == it )
				j++; 

			tmp = tmp*(xf-x[j])/(x[it]-x[j]);
		}

		// умножаем на значение в этой точке
		tmp = tmp*y[it];

		// ссумируем
		res = res + tmp;
	}
	return res;
};
//==============================================================================//
// поиск расстояния между точками
//==============================================================================//
__device__ double cuGetDistTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2 )
{
	double dist = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
	return dist;
}
//==============================================================================//
// проверка опасных сближений
//==============================================================================//
__global__ void cuFindSmallDistant( double **d_array, double *ResultFindMinDist, int *NumbResult, int size )
{
	// номер потока
	int ns = blockDim.x*blockIdx.x + threadIdx.x;
	if( ns < size-1 )
	{
		// указатели на массив
		double *OrbitPointsArray_S0 = d_array[0];
		double *OrbitPointsArray_S1 = d_array[ns+1];

		// Настройки поиска
		// число точек для интерполяции
		const int Nlograng = 5;
		const int Nlograng_half = 2;
		// позиция для старта, с третьей точки
		int Ns = 3;
		// отступаем от конца
		int Ne = 3;
		// шаг по времени в тыс секундах
		double dT = 1.0/1000.0;	
		double minDist = 10000000000.0;
		double minDistTime = 0;

		double Ts = OrbitPointsArray_S0[ Ns*4+1 ];		// Время - точка старта, отступаем от начала
		int Np0 = (int)OrbitPointsArray_S0[0];				// чисто точек
		int Np1 = (int)OrbitPointsArray_S1[0];				// чисто точек
		double Te0 = OrbitPointsArray_S0[ (Np0-1-Ne)*4+1 ];	// точка остановки, отступаем от конца
		double Te1 = OrbitPointsArray_S1[ (Np1-1-Ne)*4+1 ];	// точка остановки, отступаем от конца
		// выбираем минимальное время, чтобы не вылезти за границу массива точек
		double Te;
		if( Te0 < Te1 )
			Te = Te0;
		else 
			Te = Te1;

		// параметры
		//printf("%d\t %d\t %d\t %f\t %f\t %f\n", ns, Np0, Np1, Ts, Te, dT );

		// текущая позиция - начинае с позиции старта
		int crN0 = Ns;  // эту точку мы знаем
		int crN1 = 1;	// эту точку нужно уточнить
		
		// для поиска минимума
		double LastDist;
		int Direction = 0;
		int iterFind = 0;
		PointDist pt1; // текущая
		PointDist pt2; // средняя
		PointDist pt3; // два шага назад

		pt1.d = 0;
		pt1.t = 0;
		pt2.d = 0;
		pt2.t = 0;
		pt3.d = 0;
		pt3.t = 0;

		__syncthreads ();
		// цыкл по времени с шагом
		for( double iT = Ts; iT < Te; iT += dT )
		{
			// iT - текущее время
			// со старта совпадает с точкой старта
			// проверка не вышли ли за границу
			if( crN0 >= Np0-1-Ne || crN1 >= Np1-1-Ne )
				break;

			// проверяем какая точка сейчас текущая
			while( 1 ){
				// правильное положение точки и времени
				if( iT > OrbitPointsArray_S0[ (crN0-1)*4+1 ] && iT < OrbitPointsArray_S0[ (crN0+1)*4+1 ] )	break;
				// итератор времени меньше позиции, уменьшаем позицию
				if( iT <= OrbitPointsArray_S0[ (crN0-1)*4+1 ] )	crN0--;
				// итератор времени больше текущей позиуции, увеличиваем позицию
				if( iT >= OrbitPointsArray_S0[ (crN0+1)*4+1 ] )	crN0++;
			}
			while( 1 )	{
				// правильное положение точки и времени
				if( iT > OrbitPointsArray_S1[ (crN1-1)*4+1 ] && iT < OrbitPointsArray_S1[ (crN1+1)*4+1 ] )	break;
				// итератор времени меньше позиции, уменьшаем позицию
				if( iT <= OrbitPointsArray_S1[ (crN1-1)*4+1 ] )	crN1--;
				// итератор времени больше текущей позиуции, увеличиваем позицию
				if( iT >= OrbitPointsArray_S1[ (crN1+1)*4+1 ] )	crN1++;
			}

			// проверка не вышли ли за границу
			if( crN0 >= Np0-1-Ne || crN1 >= Np1-1-Ne )
				break;

			if( crN0 <= Nlograng_half || crN1 <= Nlograng_half )
				continue;
			// данные для интерполяции
			// crN0, crN1 - текущие центральные точки
			// iT - время
			// копируем точки для интерполяции, время Т где то по середине
					// массивы
			double PT0[Nlograng];
			double SX0[Nlograng];
			double SY0[Nlograng];
			double SZ0[Nlograng];
			int k = 0;
			for( int it = crN0 - Nlograng_half; it <= crN0 + Nlograng_half; it++)
			{
				PT0[k] = OrbitPointsArray_S0[ it*4+1 ];
				SX0[k] = OrbitPointsArray_S0[ it*4+1+1 ];
				SY0[k] = OrbitPointsArray_S0[ it*4+1+2 ];
				SZ0[k] = OrbitPointsArray_S0[ it*4+1+3 ];
				k++;
			}
			// берем соседние точки и интерполируем
			double x0 = cuGetInterpolatePoint( PT0, SX0, Nlograng, iT );
			double y0 = cuGetInterpolatePoint( PT0, SY0, Nlograng, iT );
			double z0 = cuGetInterpolatePoint( PT0, SZ0, Nlograng, iT );

			k = 0;
			for( int it = crN1 - Nlograng_half; it <= crN1 + Nlograng_half; it++)
			{
				PT0[k] = OrbitPointsArray_S1[ it*4+1 ];
				SX0[k] = OrbitPointsArray_S1[ it*4+1+1 ];
				SY0[k] = OrbitPointsArray_S1[ it*4+1+2 ];
				SZ0[k] = OrbitPointsArray_S1[ it*4+1+3 ];
				k++;
			}
			double x1 = cuGetInterpolatePoint( PT0, SX0, Nlograng, iT );
			double y1 = cuGetInterpolatePoint( PT0, SY0, Nlograng, iT );
			double z1 = cuGetInterpolatePoint( PT0, SZ0, Nlograng, iT );

			// получаем расстояние между точками
			double d =  cuGetDistTwoPoints( x0, y0, z0, x1, y1, z1 );
			d = d*1000.0; // км

			if( d < minDist )
			{
				minDist = d;
				minDistTime = iT;
			}

			//==========================================//
			// флаг прохождения минимума
			bool flagReverse = false;
			if( iterFind > 0 )
			{
				// направление на уменьшение положительное
				int nowDirection;
				if( d < LastDist )
					nowDirection = 1;
				else
					nowDirection = -1;

				// направление поменялось
				if( nowDirection < 0 && Direction > 0 )
				{
					flagReverse = true;
				}

				// сохраняем направление
				Direction = nowDirection;
			}
			// сохраняем предыдущую разницу
			LastDist = d;

			// сдвиг точек
			pt3 = pt2;
			pt2 = pt1;
			// точка следующая после минимальным расстоянием
			pt1.t = iT;
			pt1.d = d;

			// если было изменение направления
			// интерполируем расстояние
			double InterpolateDist = 100000000000;
			double InterpolateTime = 0;
			if( flagReverse == true && iterFind > 2 )
			{
				double x1 = pt3.t;
				double x2 = pt2.t;
				double x3 = pt1.t;

				double y1 = pt3.d;
				double y2 = pt2.d;
				double y3 = pt1.d;

				double a = y3 - ( x3*(y2-y1) + x2*y1 - x1*y2 )/(x2 - x1);
				a = a/( x3*( x3 - x1 - x2 ) + x1*x2 );

				double b = (y2 - y1)/(x2 - x1) - a*( x1 + x2 );

				double c = (x2*y1 - x1*y2)/(x2 - x1) + a*x1*x2;

				double discr = b*b - 4.0*a*c;

				double yv = 0;
				double xv = 0;
				if( a!= 0 )
				{
					yv = -discr/4.0/a;
					xv = -b/2.0/a;
				}
				
				InterpolateDist = yv;
				InterpolateTime = xv;
			}
			//==========================================//
			// коррекция шага
			if( d >= 100 && d < 300 )
				dT = 1.0/1000.0;
			else if( d < 100 )
				dT = 0.1/1000.0;
			else
				dT = 5.0/1000.0;

			if( flagReverse == true && InterpolateDist < 100 )
			{
				int npos = atomicAdd( &NumbResult[0], 1 );
				//int npos = (int)ResultFindMinDist[0];

				int indexpos = 3*npos + 1;
				ResultFindMinDist[ indexpos ] = ns+1;
				ResultFindMinDist[ indexpos + 1] = InterpolateTime;
				ResultFindMinDist[ indexpos + 2 ] = InterpolateDist;
			}
			iterFind ++;
		}
	}
};
//=====================================================================//

//#####################################################################//
//						функция прогноза на GPU
//#####################################################################//
	//==============================================================================//
	// инициализация интегратора на GPU
	//==============================================================================//
	void PredictOrbitSat::InitIntegrationGPU()
	{
		printf("Run Init Integration GPU 1\n");
		// constant memory for constant
		const double *hc_a = IC_a;
		const double *hc_a1 = IC_a1;
		const double *hc_w = IC_w;
		const double *hc_g = IC_g;
		const double *hc_c = IC_c;
		// Global memory for constant
		//double *dc_a;
		//double *dc_a1;
		//double *dc_w;
		//double *dc_g;
		//double *dc_c;
		//dc_a[120]; dc_a1[40]; dc_w[76]; dc_g[6]; dc_c[16];
		//cutilSafeCall( cudaMalloc( (void**)&dc_a, sizeof(double)*120 ));
		//cutilSafeCall( cudaMalloc( (void**)&dc_a1,sizeof(double)*40 ));
		//cutilSafeCall( cudaMalloc( (void**)&dc_w, sizeof(double)*76 ));
		//cutilSafeCall( cudaMalloc( (void**)&dc_g, sizeof(double)*6 ));
		//cutilSafeCall( cudaMalloc( (void**)&dc_c, sizeof(double)*16 ));
		//cutilSafeCall( cudaMemcpy( dc_a, hc_a, sizeof(double)*120, cudaMemcpyHostToDevice ));
		//cutilSafeCall( cudaMemcpy( dc_a1, hc_a1, sizeof(double)*40, cudaMemcpyHostToDevice ));
		//cutilSafeCall( cudaMemcpy( dc_w, hc_w, sizeof(double)*76, cudaMemcpyHostToDevice ));
		//cutilSafeCall( cudaMemcpy( dc_g, hc_g, sizeof(double)*6, cudaMemcpyHostToDevice ));
		//cutilSafeCall( cudaMemcpy( dc_c, hc_c, sizeof(double)*16, cudaMemcpyHostToDevice ));

		int inkkr[12] = { 3,4,5,6,1,2,4,3,2,1,6,5 };
		cutilSafeCall( cudaMemcpyToSymbol( KKR, inkkr, sizeof(int)*12, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol( Dca, hc_a, sizeof(double)*120, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol( Dca1, hc_a1, sizeof(double)*40, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol( Dcw, hc_w, sizeof(double)*76, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol( Dcg, hc_g, sizeof(double)*6, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol( Dcc, hc_c, sizeof(double)*16, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUtick, hh_tick, sizeof(double)*42, 0, cudaMemcpyHostToDevice ));

		cutilSafeCall( cudaMemcpyToSymbol(  CUaal0, h_aal0, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaal1, h_aal1, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaal2, h_aal2, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaal3, h_aal3, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaal4, h_aal4, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa0, h_aa0, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa1, h_aa1, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa2, h_aa2, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa3, h_aa3, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa4, h_aa4, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa5, h_aa5, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUaa6, h_aa6, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));

		cutilSafeCall( cudaMemcpyToSymbol(  CUbb0, h_bb0, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUbb1, h_bb1, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUbb2, h_bb2, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUbb3, h_bb3, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUbb4, h_bb4, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));

		cutilSafeCall( cudaMemcpyToSymbol(  CUcc0, h_cc0, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUcc1, h_cc1, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUcc2, h_cc2, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUcc3, h_cc3, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUcc4, h_cc4, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));

		cutilSafeCall( cudaMemcpyToSymbol(  CUdd0, h_dd0, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUdd1, h_dd1, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUdd2, h_dd2, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUdd3, h_dd3, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUdd4, h_dd4, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));

		cutilSafeCall( cudaMemcpyToSymbol(  CUee0, h_ee0, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee1, h_ee1, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee2, h_ee2, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee3, h_ee3, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee4, h_ee4, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		
		cutilSafeCall( cudaMemcpyToSymbol(  CUee5, h_ee5, sizeof(double)*7, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee6, h_ee6, sizeof(double)*7, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee7, h_ee7, sizeof(double)*7, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUee8, h_ee8, sizeof(double)*7, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUff1, h_ff1, sizeof(double)*7, 0, cudaMemcpyHostToDevice ));

		cutilSafeCall( cudaMemcpyToSymbol(  CUeet5, h_eet5, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUeet6, h_eet6, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUeet7, h_eet7, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));
		cutilSafeCall( cudaMemcpyToSymbol(  CUeet8, h_eet8, sizeof(double)*14, 0, cudaMemcpyHostToDevice ));

		cudaEventCreate ( &start );
		cudaEventCreate ( &stop );

		printf("Init Integration GPU OK\n");
	}
	//==============================================================================//
	// delete memory
	//==============================================================================//
	void PredictOrbitSat::DeleteIntegrationGPU()
	{
		// no Create momory today
		cudaEventDestroy ( start );
		cudaEventDestroy ( stop );
	};
	//==============================================================================//
	// вычисление прогноза для списка
	//==============================================================================//
	void PredictOrbitSat::RunIntegrationGpu( double t_s, double t, double h, double e, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu )
	{
		//---------------------------------------------//
		cudaEventRecord ( start, 0 );
		// таблица ядра
		dim3 threadsRK( CU_TableSize );
		dim3 gridRK( CU_GridSize );
		// выделение русурсов
		// указатель на структуру на видеокарте
		// копирование структуры
		gpuOrbitPoint d_OP;
		d_OP.AllocMemory( 6, CU_BlockXYZ );
		d_OP.CopyToGPU( ListSat );
		// задание параметров
		d_OP.d_EF = IF->d_FileMemDE403;
		d_OP.dNUT_AMPL =IF->d_AMPL;
		d_OP.dNUT_ARG = IF->d_ARG;
		d_OP.dT_finals_n = IF->d_finals_n;
		d_OP.dT_finals_tab = IF->d_finals_tab;
		d_OP.d_Garmonic = IF->d_egm96;

		d_OP.ajd0 = IF->S_ajd0;
		d_OP.delt0 = IF->S_delt0;
		d_OP.Satm = IF->SIGMA_ATM;
		d_OP.Ssun = IF->SIGMA_SUN;

		// копирование
		gpuOrbitPoint *dpointer_OP;
		cutilSafeCall( cudaMalloc( (void**)&dpointer_OP, sizeof( gpuOrbitPoint ) ) );
		cutilSafeCall( cudaMemcpy( dpointer_OP, &d_OP, sizeof( gpuOrbitPoint ), cudaMemcpyHostToDevice) );

		cudaEventRecord ( stop, 0 );
		cudaEventSynchronize ( stop );
		cudaEventElapsedTime ( &gpuTime, start, stop );
		printf("GPU prepare Time =  %.2f ms\n", gpuTime );
		//---------------------------------------------//

		cudaEventRecord ( start, 0 );
		
		kernalOrbitInitParam<<< gridRK, threadsRK>>>( dpointer_OP, h );
		kernalGetStartPoint<<< gridRK, threadsRK>>>( dpointer_OP, h );
		kernalOrbitPredict<<< gridRK, threadsRK>>>( dpointer_OP, h, t, e, NULL, false );

		cudaDeviceSynchronize();
		cudaEventRecord ( stop, 0 );
		cudaEventSynchronize ( stop );
		cudaEventElapsedTime ( &gpuTime, start, stop );
		printf("GPU cuGetStartPoint and OrbitPredict Time =  %.2f ms\n", gpuTime );

		// From GPU
		d_OP.CopyFromGPU( ListSatResGpu );
		d_OP.DeleteMemory();
		return;
	};
	//==============================================================================//
	// вычисление орбиты и поиск опасных сближений
	// t_s - время старта
	// t - время конца
	// h - шаг
	// е - ошибка интегрирования
	// ListSat - каталог
	// ListSatVerify - список для проверки
	//==============================================================================//
	void PredictOrbitSat::RunIntegrationGpu_Approach( double t_s, double t, double h, double e, SatelliteArray &ListSat, SatelliteArray &ListSatVerify, std::vector< PointWithMinDist > &result )
	{
		printf("GPU Version (Big Step Time)\n"); 

		//---------------------------------------------//
		cudaEventRecord ( start, 0 );
		// таблица ядра
		dim3 threadsRK( CU_TableSize );
		dim3 gridRK( CU_GridSize );

		// выделение русурсов
		// указатель на структуру на видеокарте
		// копирование структуры
		gpuOrbitPoint d_OP;
		d_OP.AllocMemory( 6, CU_BlockXYZ );
		d_OP.CopyToGPU( ListSat );

		// задание параметров
		d_OP.d_EF = IF->d_FileMemDE403;
		d_OP.dNUT_AMPL =IF->d_AMPL;
		d_OP.dNUT_ARG = IF->d_ARG;
		d_OP.dT_finals_n = IF->d_finals_n;
		d_OP.dT_finals_tab = IF->d_finals_tab;
		d_OP.d_Garmonic = IF->d_egm96;

		d_OP.ajd0 = IF->S_ajd0;
		d_OP.delt0 = IF->S_delt0;
		d_OP.Satm = IF->SIGMA_ATM;
		d_OP.Ssun = IF->SIGMA_SUN;

		// копирование на gpu
		gpuOrbitPoint *dpointer_OP;
		cutilSafeCall( cudaMalloc( (void**)&dpointer_OP, sizeof( gpuOrbitPoint ) ) );
		cutilSafeCall( cudaMemcpy( dpointer_OP, &d_OP, sizeof( gpuOrbitPoint ), cudaMemcpyHostToDevice) );

		// buffer constant
		double NumbdT = 8.0;
		const double TimeStep = 86.4/NumbdT;			// шаг большой части интегрировая
		const int MaxSizePoints = 21600 / NumbdT + 10;		//21600	// максимальное число точек на орбите17280

		// диапазон времен
		double dtt = t - t_s;
		// число шагов по полсуток
		int ntt = (int)(dtt/TimeStep);
		printf( "%f %d %f %f\n", dtt, ntt, t_s, t );

		// массив с контрольными точками
		double *TPR = new double[ntt+1];
		double *TPR_s = new double[ntt+1];
		for( int it = 0; it < ntt; it++ )
		{
			TPR[it] = t_s + ((float)it+1)*TimeStep;
			TPR_s[it] = t_s + ((float)it)*TimeStep;
		}

		// init orbit arr
		OrbitArrayPointsGPU LAD;
		LAD.InitArrayList( CU_BlockXYZ, MaxSizePoints );
		
		OrbitArrayPointsCPU LAC;
		LAC.InitArrayList( CU_BlockXYZ, MaxSizePoints );
		
		cudaEventRecord ( stop, 0 );
		cudaEventSynchronize ( stop );
		cudaEventElapsedTime ( &gpuTime, start, stop );
		printf("GPU prepare Time =  %.2f ms\n", gpuTime );
		//---------------------------------------------//

		//---------------------------------------------//
		// подготовка выполнение и подготовка точек
		kernalOrbitInitParam<<< gridRK, threadsRK>>>( dpointer_OP, h );
		kernalGetStartPoint<<< gridRK, threadsRK>>>( dpointer_OP, h );
		for( int it = 0; it < ntt; it++ )
		{
			cudaEventRecord ( start, 0 );
			kernalOrbitPredict_step<<< gridRK, threadsRK>>>( dpointer_OP, h, TPR_s[it], TPR[it], e, LAD.d_array_list, true );
			LAD.CopyFromGPU( LAC );
			cudaDeviceSynchronize();
			cudaEventRecord ( stop, 0 );
			cudaEventSynchronize ( stop );
			cudaEventElapsedTime ( &gpuTime, start, stop );
			printf("GPU cuGetStartPoint and OrbitPredict Time =  %.2f ms\n", gpuTime );

			int dt1 = clock();

			for( int ig = 1; ig < CU_BlockXYZ; ig++ )
			{
				FindSmallDistant_cpu( LAC.array_list[0], LAC.array_list[ig], ig, MaxSizePoints, result );

				//FILE *fres = fopen( "gpures.log", "w" );
				//int n = LAC.array_list[ig][ 0 ];
				//fprintf(fres, "%d\t", n );
				//for( int di = 0; di < n; di++ )
				//	fprintf( fres, "%f ", LAC.array_list[ig][ di*4+1 ] );
				//fprintf(fres, "\n");
				//fclose( fres );
			}

			for (int it = 0; it < result.size(); it++)
			{
				int iarr = result[it].Nlist;
				int idn = result[it].norad;
				double d = result[it].d;
				double dv = 0;// result[it].d_verify;
				double t = result[it].t;
				double dd = 0;// d - dv;
				printf("%d\t %d\t %f\t %f\t %f\t %.12f\n", iarr, idn, t, d, dv, dd);
			}

			int dt2 = clock();
			printf("TIME ALL GPU Find Close Approach %f  ms\n", (double)(dt2-dt1)/CLOCKS_PER_SEC*1000.0 );
		}
		//---------------------------------------------//

		
		// From GPU
		d_OP.DeleteMemory();
		LAD.FreeArrayList();
		LAC.FreeArrayList();
		delete TPR;
		return;
	};

	/*void PredictOrbitSat::RunIntegrationGpu( double t_s, double t, double h, double e, Orbit::OrbitArrayPointsCPU &outarr, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu )
	{
		printf("New Version GPU mod start .....\n"); 
	
		cudaEventRecord ( start, 0 );
		// таблица ядра
		dim3 threadsRK( CU_TableSize );
		dim3 gridRK( CU_GridSize );

		// выделение русурсов
		// указатель на структуру на видеокарте
		// копирование структуры
		gpuOrbitPoint d_OP;
		d_OP.AllocMemory( 6, CU_BlockXYZ );
		d_OP.CopyToGPU( ListSat );
		// задание параметров
		d_OP.d_EF = IF->d_FileMemDE403;
		d_OP.dNUT_AMPL =IF->d_AMPL;
		d_OP.dNUT_ARG = IF->d_ARG;
		d_OP.dT_finals_n = IF->d_finals_n;
		d_OP.dT_finals_tab = IF->d_finals_tab;
		d_OP.d_Garmonic = IF->d_egm96;
		// копирование
		gpuOrbitPoint *dpointer_OP;
		cutilSafeCall( cudaMalloc( (void**)&dpointer_OP, sizeof( gpuOrbitPoint ) ) );
		cutilSafeCall( cudaMemcpy( dpointer_OP, &d_OP, sizeof( gpuOrbitPoint ), cudaMemcpyHostToDevice) );

		cudaEventRecord ( stop, 0 );
		cudaEventSynchronize ( stop );
		cudaEventElapsedTime ( &gpuTime, start, stop );
		printf("GPU prepare Time =  %.2f ms\n", gpuTime );

		kernalOrbitInitParam<<< gridRK, threadsRK>>>( dpointer_OP, h );
		kernalGetStartPoint<<< gridRK, threadsRK>>>( dpointer_OP, h );

		// диапазон времен
		double dtt = t - t_s;
		// число шагов по полсуток
		int ntt = (int)(dtt/40.0);
		printf( "%f %d %f %f\n", dtt, ntt, t_s, t );
		// массив с контрольными точками
		double *TPR = new double[ntt+1];
		for( int it = 0; it < ntt; it++ )
			TPR[it] = t_s + (it+1)*40.0;
		TPR[ntt] = t;

		//-----------------------------------------//
		// число сближений
		int *h_numres = new int[1];
		int *d_numres;
		cutilSafeCall( cudaMalloc( (void**)&d_numres, sizeof( int ) ) );

		h_numres[0] = 0;
		cutilSafeCall( cudaMemcpy( d_numres, h_numres, sizeof( int ), cudaMemcpyHostToDevice) );

		int sizeresarr = 3*10000+1;
		double *h_ResultFindMinDist = new double[sizeresarr];
		double *d_ResultFindMinDist;
		cutilSafeCall( cudaMalloc( (void**)&d_ResultFindMinDist, sizeresarr*sizeof( double ) ) );
		//-----------------------------------------//

		OrbitArrayPointsGPU LAD;
		LAD.InitArrayList( CU_BlockXYZ, 10000 );

		for( int it = 0; it <= ntt; it++ )
		{

			cudaEventRecord ( start, 0 );
			kernalOrbitPredict<<< gridRK, threadsRK>>>( dpointer_OP, h, TPR[it], e, LAD.d_array_list );
			cudaDeviceSynchronize();
			cudaEventRecord ( stop, 0 );
			cudaEventSynchronize ( stop );
			cudaEventElapsedTime ( &gpuTime, start, stop );
			printf("GPU cuGetStartPoint and OrbitPredict Time =  %.2f ms\n", gpuTime );

			cudaEventRecord ( start, 0 );
			cuFindSmallDistant<<< gridRK, threadsRK>>>( LAD.d_array_list, d_ResultFindMinDist, d_numres, CU_BlockXYZ );
			cudaDeviceSynchronize();
			cudaEventRecord ( stop, 0 );
			cudaEventSynchronize ( stop );
			cudaEventElapsedTime ( &gpuTime, start, stop );
			printf("GPU cuFindSmallDistant Time =  %.2f ms\n", gpuTime );


			LAD.CopyFromGPU( outarr );
			// Nall, ni Ti Di, ....
			double *ResultFindMinDist = new double[3*10000+1]; 
			ResultFindMinDist[0] = 0;
				
			int dt1 = clock();
			for( int it = 1; it < CU_BlockXYZ; it++ )
				FindSmallDistant( outarr.array_list[0], outarr.array_list[it], it, ResultFindMinDist );

			int dt2 = clock();
			printf("TIME ALL GPU Find Close Approach %f  ms\n", (double)(dt2-dt1)/CLOCKS_PER_SEC*1000.0 );

			FILE *Fres = fopen( "DistInterpolate_gpu_ww.txt", "at" );
			int nf = ResultFindMinDist[0];
			for( int it = 0; it < nf; it++ )
			{
				int index = 3*it + 1;
				fprintf( Fres, "%.0f\t %f\t %f\n", ResultFindMinDist[index], ResultFindMinDist[index+1], ResultFindMinDist[index+2] );  
			}
			fclose( Fres );

			delete ResultFindMinDist;
		}
	
		//-----------------------------------------//
		// обработка результатов
		cutilSafeCall( cudaMemcpy( h_ResultFindMinDist, d_ResultFindMinDist, sizeresarr*sizeof( double ), cudaMemcpyDeviceToHost) );
		cutilSafeCall( cudaMemcpy( h_numres, d_numres, sizeof( int ), cudaMemcpyDeviceToHost) );


		printf( "h_numres = %d\n", h_numres[0] );
		h_ResultFindMinDist[0] = h_numres[0];

		FILE *Fres = fopen( "DistInterpolate_gpu_arr.txt", "w" );
		int nf = h_ResultFindMinDist[0];
		for( int it = 0; it < nf; it++ )
		{
			int index = 3*it + 1;
			fprintf( Fres, "%.0f\t %f\t %f\n", h_ResultFindMinDist[index], h_ResultFindMinDist[index+1], h_ResultFindMinDist[index+2] );  
		}
		fclose( Fres );
		//-----------------------------------------//

		// From GPU
		LAD.FreeArrayList();
		d_OP.CopyFromGPU( ListSatResGpu );
		d_OP.DeleteMemory();
		return;
	};*/


/*
printf("New Version GPU mod start .....\n"); 
		//---------------------------------------------//
		cudaEventRecord ( start, 0 );
		// таблица ядра
		dim3 threadsRK( CU_TableSize );
		dim3 gridRK( CU_GridSize );
		// выделение русурсов
		// указатель на структуру на видеокарте
		// копирование структуры
		gpuOrbitPoint d_OP;
		d_OP.AllocMemory( 6, CU_BlockXYZ );
		d_OP.CopyToGPU( ListSat );
		// задание параметров
		d_OP.d_EF = IF->d_FileMemDE403;
		d_OP.dNUT_AMPL =IF->d_AMPL;
		d_OP.dNUT_ARG = IF->d_ARG;
		d_OP.dT_finals_n = IF->d_finals_n;
		d_OP.dT_finals_tab = IF->d_finals_tab;
		d_OP.d_Garmonic = IF->d_egm96;

		d_OP.ajd0 = IF->S_ajd0;
		d_OP.delt0 = IF->S_delt0;
		d_OP.Satm = IF->SIGMA_ATM;
		d_OP.Ssun = IF->SIGMA_SUN;

		// копирование
		gpuOrbitPoint *dpointer_OP;
		cutilSafeCall( cudaMalloc( (void**)&dpointer_OP, sizeof( gpuOrbitPoint ) ) );
		cutilSafeCall( cudaMemcpy( dpointer_OP, &d_OP, sizeof( gpuOrbitPoint ), cudaMemcpyHostToDevice) );

		// buffer constant
		const double TimeStep = 45.0;		// шаг интегрировая
		const int MaxSizePoints = 3000;		// максимальное число точек на орбите
		const int MaxMinDist = 3000;		// максимальное число опасных сближений

		// диапазон времен
		double dtt = t - t_s;
		// число шагов по полсуток
		int ntt = (int)(dtt/TimeStep);
		printf( "%f %d %f %f\n", dtt, ntt, t_s, t );
		// массив с контрольными точками
		double *TPR = new double[ntt+1];
		for( int it = 0; it < ntt; it++ )
			TPR[it] = t_s + (it+1)*TimeStep;
		TPR[ntt] = t;

		// init orbit arr
		OrbitArrayPointsGPU LAD;
		LAD.InitArrayList( CU_BlockXYZ, MaxSizePoints );
		
		OrbitArrayPointsCPU LAC;
		LAC.InitArrayList( CU_BlockXYZ, MaxSizePoints );
		
		// сближения
		double *ResultFindMinDist = new double[3*MaxMinDist+1]; 
		ResultFindMinDist[0] = 0;

		cudaEventRecord ( stop, 0 );
		cudaEventSynchronize ( stop );
		cudaEventElapsedTime ( &gpuTime, start, stop );
		printf("GPU prepare Time =  %.2f ms\n", gpuTime );
		//---------------------------------------------//

		int maxNp1 = 0;

		// подготовка выполнение и подготовка точек
		kernalOrbitInitParam<<< gridRK, threadsRK>>>( dpointer_OP, h );
		kernalGetStartPoint<<< gridRK, threadsRK>>>( dpointer_OP, h );
		for( int it = 0; it <= ntt; it++ )
		{
			cudaEventRecord ( start, 0 );
			kernalOrbitPredict<<< gridRK, threadsRK>>>( dpointer_OP, h, TPR[it], e, LAD.d_array_list, true );
			LAD.CopyFromGPU( LAC );
			cudaDeviceSynchronize();
			cudaEventRecord ( stop, 0 );
			cudaEventSynchronize ( stop );
			cudaEventElapsedTime ( &gpuTime, start, stop );
			printf("GPU cuGetStartPoint and OrbitPredict Time =  %.2f ms\n", gpuTime );

			int dt1 = clock();
			for( int it = 1; it < CU_BlockXYZ; it++ )
			{
				int res = FindSmallDistant( LAC.array_list[0], LAC.array_list[it], it, ResultFindMinDist );
				if( res > maxNp1 )
					maxNp1 = res;
				printf( "%d\r", maxNp1 );
			}
			int dt2 = clock();
			printf("TIME ALL GPU Find Close Approach %f  ms\n", (double)(dt2-dt1)/CLOCKS_PER_SEC*1000.0 );
		}

		FILE *Fres = fopen( "DistInterpolate_gpul.txt", "w" );
		int nf = ResultFindMinDist[0];
		for( int it = 0; it < nf; it++ )
		{
			int index = 3*it + 1;
			int iarr = ResultFindMinDist[index];
			int idn = ListSat.GetSatelliteID( iarr );
			fprintf( Fres, "%d\t %d\t %f\t %f\n", iarr, idn, ResultFindMinDist[index+1], ResultFindMinDist[index+2] );  
		}
		fclose( Fres );

		// From GPU
		d_OP.CopyFromGPU( ListSatResGpu );
		d_OP.DeleteMemory();
		LAD.FreeArrayList();
		LAC.FreeArrayList();
		delete ResultFindMinDist;
		delete TPR;
*/
};
#endif
//#####################################################################//
