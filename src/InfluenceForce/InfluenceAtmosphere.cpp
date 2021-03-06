//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Earth Atmosphere
//==============================================================================//
//#include <math.h>
#include <cmath>
#include <stdio.h>

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iostream>

#include "InfluenceForce.h"

namespace Force
{
	//==============================================================================//
	// ???????
	//==============================================================================//
	int nint( double X )
	{
		double res;
		int iX = (int)X;
		double delta = X - (double)iX;

		if( delta > 0.5 )
			res = iX + 1;
		else
			res = iX;

		return iX;
	};
	double mmax( double X, double Y )
	{
		if( X > Y )
			return X;
		else
			return Y;
	};
	double mmin( double X, double Y )
	{
		if( X < Y )
			return X;
		else
			return Y;
	};
	double mod( double X, double Y )
	{
		int s = (int)(X/Y);
		double res = X - ((double)s)*Y;

		return res;
	};

	std::vector<std::string> split(std::string strToSplit, char delimeter)
	{
		std::stringstream ss(strToSplit);
		std::string item;
		std::vector<std::string> splittedStrings;
		while (std::getline(ss, item, delimeter))
		{
			splittedStrings.push_back(item);
		}
		return splittedStrings;
	}

	int InfluenceForce::InitAtmMaps(double st_date, double en_date){
		int m_num = en_date - st_date + 1;

		int st_datei = st_date;
		int en_datei = en_date;

		r_s = 6498136;
		r_e = 7878136;
		step = 20000;
		this->layer_num_set(r_s, r_e, step);

		int degree = this->degreeG;
		this->map_param_set(degree);
		this->triangleLoader("Triangles_" + to_string(degree) + "F.txt");

		for (int i = 0; i < m_num; i++) {
			single_am_loader(st_datei + i, degree);
			this->maps_atm.push_back(map_buff);
		}
		this->cur_map = &maps_atm[0];

		//infileV.open("../maps_gr/fl_6m.txt", std::ios_base::app);
		return 0;
	}

	void InfluenceForce::InitAtm()
	{
		char *fname = "data/atm.config";

		std::ifstream file(fname);
		std::string str;

		int count = 0;
		while (std::getline(file, str))
		{
			count++;
			// Process str
			if (count >= 3)
			{
				std::vector<std::string> arr;
				arr.push_back(str.substr(2, 4));
				arr.push_back(str.substr(7, 2));
				arr.push_back(str.substr(10, 2));
				arr.push_back(str.substr(16, 7));
				arr.push_back(str.substr(27, 7));
				arr.push_back(str.substr(39, 6));
				if (arr.size() == 6 )
				{
					//YYYY MM DD F_10_7   F81   KP
					atmIndex pt;
					pt.YYYY = std::stoi(arr[0]);
					pt.MM = std::stoi(arr[1]);
					pt.DD = std::stoi(arr[2]);

					pt.F107 = std::stod(arr[3]);
					pt.F81 = std::stod(arr[4]);
					pt.aKp = std::stod(arr[5]);

					pt.data = pt.DD + 100 * pt.MM + 10000* pt.YYYY;
					ListAtmIndex.push_back(pt);
				}
			}
		}
		if (ListAtmIndex.size() > 0)
		{
			printf("Load Atm index: %d from %d to %d \n", ListAtmIndex.size(), ListAtmIndex[0].data, ListAtmIndex[ListAtmIndex.size() - 1].data);
		}
		else
		{
			printf("Zero Load Atm index\n");
		}
	}
	void InfluenceForce::InitAtm(std::string path_atmconf)
	{
		const char *fname = path_atmconf.data();

		std::ifstream file(fname);
		std::string str;

		int count = 0;
		while (std::getline(file, str))
		{
			count++;
			// Process str
			if (count >= 3)
			{
				std::vector<std::string> arr;
				arr.push_back(str.substr(2, 4));
				arr.push_back(str.substr(7, 2));
				arr.push_back(str.substr(10, 2));
				arr.push_back(str.substr(16, 7));
				arr.push_back(str.substr(27, 7));
				arr.push_back(str.substr(39, 6));
				if (arr.size() == 6)
				{
					//YYYY MM DD F_10_7   F81   KP
					atmIndex pt;
					pt.YYYY = std::stoi(arr[0]);
					pt.MM = std::stoi(arr[1]);
					pt.DD = std::stoi(arr[2]);

					pt.F107 = std::stod(arr[3]);
					pt.F81 = std::stod(arr[4]);
					pt.aKp = std::stod(arr[5]);

					pt.data = pt.DD + 100 * pt.MM + 10000 * pt.YYYY;
					ListAtmIndex.push_back(pt);
				}
			}
		}
		if (ListAtmIndex.size() > 0)
		{
			printf("Load Atm index: %d from %d to %d \n", ListAtmIndex.size(), ListAtmIndex[0].data, ListAtmIndex[ListAtmIndex.size() - 1].data);
		}
		else
		{
			printf("Zero Load Atm index\n");
		}
	}
	//==============================================================================//
	// ??????? ?????? ????????? ??????? ????? ????????? ???????? ?????.
	//==============================================================================//
	double InfluenceForce::Roa2004_2(double time, double *x, double ajd0, double delt0, double F107, double F81, double aKp)
	{
		double f0t[7] = {  75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 200.0 };
		double Re = 6378.136;
		double om = 6.300388008;
		double alpha = 0.0033528037;
		double b1900 = 2415020.31352;
		double  btau = 365.2422;
		double  ro0 = 1.58868E-8;
		double   an0 =  .20580000E+01;
		double   an1 =  .58870000E-02;
		double  an2 = -.40120000E-05;

		//double saem[4];
		//saem[0] = 100.0;
		//saem[1] = 100.0;
		//saem[2] = 100.0;
		//saem[3] = 3.0;

		double h_tick[42] = {  
			 500.0, 600.0, 640.0, 1500.0, 600.0, 640.0,  
			 500.0, 660.0, 700.0, 1500.0, 700.0, 660.0,  
			 500.0, 760.0, 760.0, 1500.0, 780.0, 740.0,  
			 500.0, 800.0, 820.0, 1500.0, 800.0, 800.0,  
			 500.0, 860.0, 860.0, 1500.0, 800.0, 860.0,  
			 500.0, 900.0, 920.0, 1500.0, 900.0, 900.0,  
			 500.0,1000.0, 980.0, 1500.0, 760.0, 900.0 };
		double aa0[14] = { 
			  .26862900E+02, .27459800E+02, .28639500E+02, .29641800E+02,  
			  .30167100E+02, .29757800E+02, .30785400E+02,  
			  .17878100E+02,-.25490900E+01,-.13959900E+02,-.23307900E+02,  
			 -.14726400E+02,-.49120000E+01,-.54095200E+01 };
		double aa1[14] = { 
			 -.45167400E+00,-.46366800E+00,-.49098700E+00,-.51495700E+00,  
			 -.52783700E+00,-.51791500E+00,-.54569500E+00,  
			 -.13202500E+00, .14006400E-01, .84495100E-01, .13514100E+00,  
			  .71325600E-01, .10832600E-01, .55074900E-02 };
		double aa2[14] = { 
			  .29039700E-02, .29740000E-02, .32064900E-02, .34192600E-02,  
			  .35321100E-02, .34269900E-02, .37032800E-02,  
			  .22771700E-03,-.16946000E-03,-.32887500E-03,-.42080200E-03,  
			 -.22801500E-03,-.81054600E-04,-.37885100E-04 };
		double aa3[14] = { 
			 -.10695300E-04,-.10753000E-04,-.11681000E-04,-.12578500E-04,  
			 -.13022700E-04,-.12413700E-04,-.13707200E-04,  
			 -.22543000E-06, .32719600E-06, .50591800E-06, .57371700E-06,  
			  .28487000E-06, .11571200E-06, .24808000E-07 };
		double aa4[14] = { 
			  .22159800E-07, .21705900E-07, .23684700E-07, .25727000E-07,  
			  .26645500E-07, .24820900E-07, .28061400E-07,  
			  .13357400E-09,-.28763000E-09,-.39229900E-09,-.40323800E-09,  
			 -.17438300E-09,-.81329600E-10, .49218300E-11 };
		double aa5[14] = { 
			 -.24294100E-10,-.23024900E-10,-.25180900E-10,-.27587400E-10,  
			 -.28543200E-10,-.25841300E-10,-.30018400E-10,  
			 -.45045800E-13, .12262500E-12, .15227900E-12, .14284600E-12,  
			  .50807100E-13, .30491300E-13,-.86501100E-14 };
		double aa6[14] = { 
			  .10992600E-13, .10012300E-13, .10953600E-13, .12109100E-13,  
			  .12500900E-13, .10938300E-13, .13114200E-13,  
			  .67208600E-17,-.20573600E-16,-.23557600E-16,-.20172600E-16,  
			 -.53495500E-17,-.49498900E-17, .19849000E-17 };
		double bb0[14] = { 
			  .68789400E-01, .15073000E+00, .47945100E-01, .22344800E-01,  
			  .32639100E-02,-.51474900E-01,-.10725500E+00,  
			  .23158400E+02, .33273200E+02, .39196100E+02, .43246900E+02,  
			  .49573800E+02, .11278000E+02,-.52618400E+02 };
		double bb1[14] = { 
			 -.28407700E-02,-.40088900E-02,-.23945300E-02,-.19798000E-02,  
			 -.15986900E-02,-.92105900E-03,-.17434300E-03,  
			 -.80214700E-01,-.11109900E+00,-.12352000E+00,-.12697300E+00,  
			 -.13861300E+00, .14347800E-02, .21468900E+00 };
		double bb2[14] = { 
			  .18392200E-04, .24393700E-04, .17033500E-04, .15410100E-04,  
			  .14044300E-04, .11514700E-04, .90275900E-05,  
			  .10582400E-03, .14142100E-03, .14901500E-03, .14263700E-03,  
			  .14785100E-03,-.36984600E-04,-.29488200E-03 };
		double bb3[14] = { 
			  .91960500E-08,-.99277200E-08,-.13162600E-08,-.23543000E-08,  
			 -.30228700E-08,-.12290100E-08,-.31651200E-09,  
			 -.61503600E-07,-.79495200E-07,-.79705000E-07,-.70998500E-07,  
			 -.69636100E-07, .35831800E-07, .17117100E-06 };
		double bb4[14] = { 
			 -.41687300E-10,-.18223900E-10,-.17403200E-10,-.12499400E-10,  
			 -.92016000E-11,-.81310400E-11,-.61400000E-11,  
			  .13245300E-10, .16583600E-10, .15877200E-10, .13164600E-10,  
			  .12159500E-10,-.99122500E-11,-.36058200E-10 };
		double cc0[14] = { 
			 -.10482500E+01,-.93106000E+00,-.82086700E+00,-.74404700E+00,  
			 -.72247100E+00,-.68748200E+00,-.73998400E+00,  
			  .50503400E+02, .61624000E+02, .53262300E+02, .18223600E+02,  
			 -.31844200E+02,-.48720800E+02,-.14785900E+03 };
		double cc1[14] = { 
			  .16630500E-01, .14153700E-01, .11991600E-01, .10474300E-01,  
			  .98031700E-02, .91659400E-02, .95285400E-02,  
			 -.17054100E+00,-.19296700E+00,-.14434200E+00,-.84002400E-02,  
			  .16832700E+00, .22299600E+00, .53165200E+00 };
		double cc2[14] = { 
			 -.92426300E-04,-.72986200E-04,-.57983500E-04,-.47854400E-04,  
			 -.42524500E-04,-.38093200E-04,-.36272700E-04,  
			  .21723200E-03, .22806100E-03, .14659000E-03,-.38800000E-04,  
			 -.26260300E-03,-.32188400E-03,-.67193700E-03 };
		double cc3[14] = { 
			  .27238200E-06, .20029400E-06, .15070700E-06, .11851300E-06,  
			  .99554400E-07, .85127500E-07, .73887000E-07,  
			 -.12190200E-06,-.11871500E-06,-.64644300E-07, .43138400E-07,  
			  .16545400E-06, .19149500E-06, .36478700E-06 };
		double cc4[14] = { 
			 -.24135500E-09,-.16200600E-09,-.11302600E-09,-.83149800E-10,  
			 -.65517500E-10,-.52997200E-10,-.42390700E-10,  
			  .25403700E-10, .22963800E-10, .10422700E-10,-.12383200E-10,  
			 -.36935500E-10,-.40806700E-10,-.72626800E-10 };
		double dd0[14] = { 
			 -.35189900E+00,-.47813000E-01, .20981000E+00, .26517400E+00,  
			  .23047000E+00, .17007400E+00, .88141000E-01,  
			 -.35189900E+00,-.47813000E-01, .20981000E+00, .26517400E+00,  
			  .23047000E+00, .17007400E+00, .88141000E-01 };
		double dd1[14] = { 
			  .57705600E-02, .38081300E-02, .26288100E-02, .27583600E-02,  
			  .33833100E-02, .40613100E-02, .46825300E-02,  
			  .57705600E-02, .38081300E-02, .26288100E-02, .27583600E-02,  
			  .33833100E-02, .40613100E-02, .46825300E-02 };
		double dd2[14] = { 
			  .99581900E-06, .42277100E-05, .42437900E-05, .20866800E-05,  
			 -.55230500E-06,-.28211400E-05,-.42460900E-05,  
			  .99581900E-06, .42277100E-05, .42437900E-05, .20866800E-05,  
			 -.55230500E-06,-.28211400E-05,-.42460900E-05 };
		double dd3[14] = { 
			 -.72532400E-08,-.86682600E-08,-.66732800E-08,-.36954300E-08,  
			 -.82360700E-09, .13836900E-08, .25350900E-08,  
			 -.72532400E-08,-.86682600E-08,-.66732800E-08,-.36954300E-08,  
			 -.82360700E-09, .13836900E-08, .25350900E-08 };
		double dd4[14] = { 
			  .29759000E-11, .30671200E-11, .21349600E-11, .11186200E-11,  
			  .22134900E-12,-.42790800E-12,-.72903100E-12,  
			  .29759000E-11, .30671200E-11, .21349600E-11, .11186200E-11,  
			  .22134900E-12,-.42790800E-12,-.72903100E-12 };
		double ee0[14] = { 
			 -.73159600E+00,-.75217500E+00,-.57047600E+00,-.94957300E+00,  
			 -.96759800E+00,-.10227800E+01,-.75790300E+00,  
			  .38619900E+02, .51249000E+02, .68474600E+02, .58422000E+02,  
			  .72018800E+01, .21594800E+02,-.88407600E+02 };
		double ee1[14] = { 
			  .59734500E-02, .56592500E-02, .29580200E-02, .81312100E-02,  
			  .84199100E-02, .92363300E-02, .60606800E-02,  
			 -.13214700E+00,-.16737300E+00,-.21565900E+00,-.16666400E+00,  
			  .21610900E-01,-.20223900E-01, .33851800E+00 };
		double ee2[14] = { 
			 -.58203700E-05, .18082000E-05, .16889600E-04,-.38781300E-05,  
			 -.35850000E-05,-.61012800E-05, .78529600E-05,  
			  .17541100E-03, .21183200E-03, .26227300E-03, .18548600E-03,  
			 -.65288200E-04,-.17202900E-04,-.44558100E-03 };
		double ee3[14] = { 
			  .68463400E-07, .33382200E-07,-.47475600E-08, .23769400E-07,  
			  .17480100E-07, .17821100E-07,-.97489100E-08,  
			 -.10241700E-06,-.11822100E-06,-.14097200E-06,-.91234500E-07,  
			  .53707700E-07, .28301700E-07, .25172900E-06 };
		double ee4[14] = { 
			 -.95048300E-10,-.51396500E-10,-.17271100E-10,-.27746900E-10,  
			 -.19622100E-10,-.17007300E-10, .15837700E-11,  
			  .22144600E-10, .24505500E-10, .28228500E-10, .16711800E-10,  
			 -.14095000E-10,-.89448600E-11,-.52030000E-10 };
		double ff1[7] = { 
			  .54110000E+00, .55150000E+00, .55850000E+00, .55850000E+00,  
			  .55850000E+00, .55850000E+00, .55850000E+00 };
		double ee5[7] = { 
			 -.20670000E+00,-.16971000E+00,-.14671000E+00,-.13150000E+00,  
			 -.12091600E+00,-.11363000E+00,-.10444000E+00 };
		double ee6[7] = { 
			  .97533000E-01, .79830000E-01, .68808000E-01, .61603000E-01,  
			  .56538000E-01, .53178000E-01, .48551000E-01 };
		double ee7[7] = { 
			 -.11817000E-01,-.94393000E-02,-.79836000E-02,-.70866000E-02,  
			 -.64324000E-02,-.60436000E-02,-.53567000E-02 };
		double ee8[7] = { 
			  .16145000E-02, .12622000E-02, .10535000E-02, .92813000E-03,  
			  .83723000E-03, .77982000E-03, .68809000E-03 };
		double eet5[14] = { 
			 -.20610000E+00,-.16927900E+00,-.14637700E+00,-.13121000E+00,  
			 -.12067000E+00,-.11339900E+00,-.10424300E+00,  
			 -.20610000E+00,-.16927900E+00,-.14637700E+00,-.13121000E+00,  
			 -.12067000E+00,-.11339900E+00,-.10424300E+00 };
		double eet6[14] = { 
			  .94449000E-01, .77599000E-01, .67052000E-01, .60105000E-01,  
			  .55232000E-01, .51994000E-01, .47573000E-01,  
			  .94449000E-01, .77599000E-01, .67052000E-01, .60105000E-01,  
			  .55232000E-01, .51994000E-01, .47573000E-01 };
		double eet7[14] = { 
			 -.87953000E-02,-.71375000E-02,-.60951000E-02,-.54388000E-02,  
			 -.49580000E-02,-.46876000E-02,-.41711000E-02,  
			 -.87953000E-02,-.71375000E-02,-.60951000E-02,-.54388000E-02,  
			 -.49580000E-02,-.46876000E-02,-.41711000E-02 };
		double eet8[14] = { 
			  .88385000E-03, .69025000E-03, .57456000E-03, .50585000E-03,  
			  .45512000E-03, .42548000E-03, .37068000E-03,  
			  .88385000E-03, .69025000E-03, .57456000E-03, .50585000E-03,  
			  .45512000E-03, .42548000E-03, .37068000E-03 };
		double aal0[14] = { 
			 -.40776800E+00,-.90273900E+00,-.73303700E+00,-.13144400E+01,  
			 -.12002600E+01,-.15215800E+01,-.16766400E+01,  
			  .48653600E+02, .54486700E+02, .60126700E+02, .47099600E+02,  
			  .50617400E+02, .80194200E+01,-.15572800E+02 };
		double aal1[14] = { 
			  .14850600E-02, .82680300E-02, .52339600E-02, .13312400E-01,  
			  .11408700E-01, .15704000E-01, .17719400E-01,  
			 -.17029100E+00,-.17829800E+00,-.18314400E+00,-.12526000E+00,  
			 -.12904700E+00, .18530200E-01, .93670400E-01 };
		double aal2[14] = { 
			  .12535700E-04,-.12544800E-04, .63566700E-05,-.25558500E-04,  
			 -.14732400E-04,-.30285900E-04,-.36949800E-04,  
			  .22624200E-03, .22272500E-03, .21248100E-03, .12635200E-03,  
			  .12484200E-03,-.61473300E-04,-.14903600E-03 };
		double aal3[14] = { 
			  .37731100E-07, .61285300E-07, .10906500E-07, .54398100E-07,  
			  .27804000E-07, .45766800E-07, .50913400E-07,  
			 -.13203200E-06,-.12270000E-06,-.10849700E-06,-.55158400E-07,  
			 -.52499300E-07, .49767400E-07, .94215100E-07 };
		double aal4[14] = { 
			 -.77895300E-10,-.70796600E-10,-.26142700E-10,-.43378400E-10,  
			 -.22632000E-10,-.28292600E-10,-.28287800E-10,  
			  .28519300E-10, .25131600E-10, .20571000E-10, .87527200E-11,  
			  .80827200E-11,-.12616200E-10,-.20961000E-10 };

		//double F107, F81, aKp,
		double Rmod, sinF, h, F0;
		double e0, e1, e2, e3, e4, e5, e6, e7, e8;
		double a0, a1, a2, a3, a4, a5, a6;
		double al0, al1, al2, al3, al4;
		double c0, c1, c2, c3, c4;
		double d0, d1, d2, d3, d4;
		double b0, b1, b2, b3, b4;
		double fi1, power, aKs0, aKs1, aKs2, aKs3, aKs4, aKss1, aKss4, ron;
		double ct, sz, xss, yss, zss, rr;
		double cosFi1, sinFi1, app, s_app, c_app, xs, ys, zs, rs, cosFi, cos05, ak0, ak1, d, ad, ak2, ak3, ak4;
		double h0, akst1, akst2, ast0, t;
		int iq, iq6 ,ish_a, ish_b, ish_c, ish_d, ish_e, ish_l;
		
		t = ajd0+(delt0+time)/86.4;
		
		Rmod = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
		sinF = x[2]/Rmod;
		h = mmax( Rmod*1.0E+3-Re*(1.0-alpha*sinF*sinF), 0.0);
		//cout << time << " " << h << endl;
		double roa2004 = 0;

		if (h >= 1500.0) 
		{
			roa2004 = 0.0;
			return roa2004;
		}
		else if (h > 120.0) 
		{
			//F107 = saem[1];
			//F81  = saem[2];
			//aKp  = saem[3];
			iq = mmin(mmax(nint((F81-75.0)/25.0), 0), 7) +1;
			
			F0 = f0t[iq-1]; //!!!

			if( F81 > 200.0)
			{
				iq = 6;
			}
			else if (F81 > 225.0) 
			{
				F0 = 250.0;
				iq = 7;
			}

			iq6 = (iq-1)*6;

			iq6 = iq6 - 1; //!! ????????? ???????

			ish_a=0;
			ish_b=0;
			ish_c=0;
			ish_d=0;
			ish_e=0;
			ish_l=0;
			//??????????? ???????? ??????? ?? ????????
			if (h > h_tick[iq6+1]) ish_a=7;
			if (h > h_tick[iq6+2]) ish_b=7;
			if (h > h_tick[iq6+3]) ish_c=7;
			if (h > h_tick[iq6+4]) ish_d=7;
			if (h > h_tick[iq6+5]) ish_e=7;
			if (h > h_tick[iq6+6]) ish_l=7;

			// ------------------
			iq = iq - 1; //!! ????????? ???????
		   // a - coeffs
		   a0=aa0[ish_a+iq];
		   a1=aa1[ish_a+iq];
		   a2=aa2[ish_a+iq];
		   a3=aa3[ish_a+iq];
		   a4=aa4[ish_a+iq];
		   a5=aa5[ish_a+iq];
		   a6=aa6[ish_a+iq];
		   // al - coeffs
		   al0=aal0[ish_l+iq];
		   al1=aal1[ish_l+iq];
		   al2=aal2[ish_l+iq];
		   al3=aal3[ish_l+iq];
		   al4=aal4[ish_l+iq];
		   // c - coeffs
		   c0=cc0[ish_c+iq];
		   c1=cc1[ish_c+iq];
		   c2=cc2[ish_c+iq];
		   c3=cc3[ish_c+iq];
		   c4=cc4[ish_c+iq];
		   // d - coeffs
		   d0=dd0[ish_d+iq];
		   d1=dd1[ish_d+iq];
		   d2=dd2[ish_d+iq];
		   d3=dd3[ish_d+iq];
		   d4=dd4[ish_d+iq];
		   // b - coeffs
		   b0=bb0[ish_b+iq];
		   b1=bb1[ish_b+iq];
		   b2=bb2[ish_b+iq];
		   b3=bb3[ish_b+iq];
		   b4=bb4[ish_b+iq];
		   // e - coeffs
		   e0=ee0[ish_e+iq];
		   e1=ee1[ish_e+iq];
		   e2=ee2[ish_e+iq];
		   e3=ee3[ish_e+iq];
		   e4=ee4[ish_e+iq];
		   e5=ee5[iq];
		   e6=ee6[iq];
		   e7=ee7[iq];
		   e8=ee8[iq];
		   // ------------------

		   fi1 = ff1[iq];

		   power = a0+h*(a1+h*(a2+h*(a3+h*(a4+h*(a5+h*a6)))));
		   aKs0 = al0+h*(al1+h*(al2+h*(al3+h*al4)));
		   aKs1 = c0+h*(c1+h*(c2+h*(c3+h*c4)));
		   aKss1 = an0+h*(an1+h*an2);
		   aKs2 = d0+h*(d1+h*(d2+h*(d3+h*d4)));
		   aKs3 = b0+h*(b1+h*(b2+h*(b3+h*b4)));
		   aKs4 = e0+h*(e1+h*(e2+h*(e3+h*e4)));
		   aKss4 = e5+aKp*(e6+aKp*(e7+aKp*e8));
		   ron = ro0*exp(power);

		   ct  = (t-2415020.0)/36525.0;
		   sz  = 628.33195099*ct+1.739935890;
		   xss = -cos(sz);
		   sinF = -sin(sz);
		   yss = 0.91747*sinF;
		   zss = 0.397805*sinF;
		   rr = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );

		   cosFi1 = cos(fi1);
		   sinFi1 = sin(fi1);
		   app = sz+( mod(t, 1.0)-0.5)*om;
		   s_app = sin(app);
		   c_app = cos(app);
		   xs = xss*c_app + yss*s_app;
		   ys = yss*c_app - xss*s_app;
		   zs = zss;
		   rs = sqrt(xs*xs+ys*ys+zs*zs);
		   cosFi = (x[0]*(xs*cosFi1-ys*sinFi1)+ x[1]*(ys*cosFi1+xs*sinFi1)+ x[2]*zs)/(rs*Rmod);
		   cos05 = sqrt((1.0+cosFi)/2.0);
		   ak0 = 1.0 + aKs0*(F81-F0)/F0;
		   //aK1 = aKs1*(cos05**aKss1); //!!!
		   ak1 = aKs1*pow(cos05, aKss1);
		   d = mod(t-b1900, btau);
		   ad = -2.53418E-2 + d*(-2.44075E-3+d*(3.08389E-6 +d*( 
				 2.90115E-6 + d*(-4.99606E-8+d*(3.36327E-10+d*( 
				-1.0966E-12 + d*(1.73227E-15+d*(-1.06271E-18))))))));
		   ak2 = ad*aKs2;
		   ak3 = aKs3*(F107-F81)/(F81+fabs(F107-F81)); //!!!!!!linux
		   ak4 = aKs4*aKss4;
		   roa2004 = ron*ak0*(1.0 + ak1 + ak2 + ak3 + ak4);
		   return roa2004;
	   }
		else if ((h > 0.0) && (h < 20.0)) 
		{
			h0=0.0;
			ast0=1.228;
			akst1=-9.0764E-2;
			akst2=-2.0452E-3;
	   }
		else if ((h >= 20.0)&&(h < 60.0))
		{
			h0=20.0;
			ast0=9.013E-2;
			akst1=-0.16739;
			akst2=6.2669E-4;
		}
		else if ((h >= 60.0)&&(h < 100.0)) 
		{
			h0=60.0;
			ast0=3.104E-4;
			akst1=-0.137;
			akst2=-7.8653E-4;
		}
		else if ((h >= 100.0)&&(h <= 120.0)) 
		{
			h0=100.0;
			ast0=3.66E-7;
			akst1=-0.18553;
			akst2=1.5397E-3;
		}

		power = (h-h0)*(akst1+akst2*(h-h0));
		roa2004 = ast0*exp(power);
		return roa2004;
	}
	//==============================================================================//
	// ????????? ?????? ??????????, ?????????? ?????????????? ?????????.
	//==============================================================================//
	void InfluenceForce::Atm_drag( double *x, double t, double *f, double sigma_up, double ajd0, double delt0 )
	{
		double v;
		v = sqrt( x[3]*x[3] + x[4]*x[4] + x[5]*x[5] );

		double jd_current = ajd0 + (delt0 + t) / 86.40;
		double data_current = Ajd_dt(jd_current);
		int datai_current = data_current;
		/*
		if (datai_current > this->curDate) {
			if (this->dayC == 1)
				this->stTime = t;
			this->curDate = datai_current;
			this->dayC += 1;
			this->dC = true;
		}

		if (dayC == 1) {
			this->dayTime = t * 1000;
		}
		else {
			this->dayTime = 1000 * (t - this->stTime - (dayC - 2) * 86.4);
			if (this->dayTime < 0)
				this->dayTime = 0;
		}
		if (dayTime == 0 && dC == true) {
			this->cur_map = &this->maps_atm[dayC - 1];
			this->dC = false;
		}

		//cout << "Secs: " << t << " Data: " << datai_current << " Height: " << 100 * (sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) - 6.378136) << "km" << endl;
		//cout << "Day: " << dayC << " Day Time: " << dayTime << "\n" << endl;
		*/
		double F107 = 100;
		double F81 = 100;
		double aKp = 3;
		
		// find index
		for (int i = 0; i < ListAtmIndex.size(); i++)
		{
			if (ListAtmIndex[i].data == datai_current)
			{
				F107 = ListAtmIndex[i].F107;
				F81 = ListAtmIndex[i].F81;
				aKp = ListAtmIndex[i].aKp;
			}
		}

		double rc = Roa2004_2(t, x, ajd0, delt0, F107, F81, aKp );

		//cout << "P_c: " << rc  << " P_ex: " << extrapolatorA(x, dayTime) << endl;

		double coeff =  rc*sigma_up*v*1.0E+6;

		f[0] = -x[3]*coeff;
		f[1] = -x[4]*coeff;
		f[2] = -x[5]*coeff;

		//infileV << x[0] << " " << x[1] << " " << x[2] << " " << rc << "\n";

		this->cc++;
	}
	//==============================================================================//

	void InfluenceForce::Atm_extr(double* x, double t, double* f, double sigma_up, double ajd0, double delt0) {
		double v;
		v = sqrt(x[3] * x[3] + x[4] * x[4] + x[5] * x[5]);

		double jd_current = ajd0 + (delt0 + t) / 86.40;
		double data_current = Ajd_dt(jd_current);
		int datai_current = data_current;
		
		if (datai_current > this->curDate) {
			if (this->dayC == 1)
				this->stTime = t;
			this->curDate = datai_current;
			this->dayC += 1;
			this->dC = true;
		}

		if (dayC == 1) {
			this->dayTime = t * 1000;
		}
		else {
			this->dayTime = 1000 * (t - this->stTime - (dayC - 2) * 86.4);
			if (this->dayTime < 0)
				this->dayTime = 0;
		}
		if (dayTime == 0 && dC == true) {
			this->cur_map = &this->maps_atm[dayC - 1];
			this->dC = false;
		}

		double p = extrapolatorA(x, dayTime);
		double coeff = p * sigma_up * v * 1.0E+6;
		//cout << "Secs: " << t << " Data: " << datai_current << " Height: " << 1000 * (sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) - 6.378136) << "km" << endl;
		//cout << "Day: " << dayC << " Day Time: " << dayTime << " P: "<< p << "\n" << endl;
		f[0] = -x[3] * coeff;
		f[1] = -x[4] * coeff;
		f[2] = -x[5] * coeff;
		//infileV << x[0] << " " << x[1] << " " << x[2] << " " << p << "\n";
		this->cc++;
	}
};