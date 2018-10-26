//==============================================================================//
// Andrianov N.G.
// opbit predict module
// get constant
//==============================================================================//
#include "PredictOrbitSat.h" 
#include <string.h>
// name
namespace Orbit
{
	//==============================================================================//
	// константы для интегрирования
	//==============================================================================//
	double *PredictOrbitSat::InitParamA(  )
	{
		double a[120] = {
			-.019610491275593472272648676054287534656347,
			.196654108088638376428489472644939971840735,
			-.321928796454549357401512838287345619553191,
			.235672466002613193733835568578061677895444,
			-.164066973920507500617311204599663329905077,
			.142711531762372472517174433271890081830574,
			-.059059662809989025767068851444422824806559,
			.166192854270799103153093817522920008949859,
			-.260196351165171078142337121212043098781757,
			.44332221026179213842713490576279169819185,
			-.844248687165086253284249107803523489857774,
			.883999114815226983212093477622655362704093,
			-.032675285416564170270845712585139854060305,
			.180296891114417569043759375460402590930473,
			-.800237467888389222743369511892480107652254,
			2.907974028345013451886177840182695538682369,
			-3.023785871199941118021636472000107659609341,
			1.438418226837891622507247360386251835309345,
			-.017681213114125620098767835352200212821843,
			.240222971150109887390239971743070402581719,
			-1.662867944252536651104110457834383378395254,
			2.604729744651109519132085149384001678001547,
			-1.000602115694225667562024288969642989572428,
			.766766713056694819854550705475559252754122,
			-.023532181738490017855643312188788238280074,
			.235444756170251191549847055842371449076444,
			-.384772895059710101885959450842700235144721,
			.280277061548448647017143363431794558688285,
			-.191853363223592382328568171253865820953489,
			.153868466506066375891207270564783534065692,
			-.146670283911363133609247046084677679985881,
			.408576897315018840605999856033494803547051,
			-.622097822321896441665322664459324031190941,
			.978417029123618733697495555596512027730709,
			-1.468393563775878338561234828552904200973287,
			1.180177221778072207130976247915276737272062,
			-.118621923190321147172858215228057443861943,
			.617878268720229825361190580870947622440374,
			-2.371043542562823248107790279326182156140948,
			6.19652008105860789533941540398303330911738,
			-5.303214333080455591906418331155518924216553,
			1.648471970847190398887793720407399936261977,
			-.068057174159688090243430870215237526711285,
			.75374056490445931173113239275094491283215,
			-3.486664263728150950438696181446088443453779,
			4.248792491507461402502278320683982542109365,
			-1.310613085106969068034324784717157700118059,
			.793369622379913682095014367389960967889471,
			-.062146172814071375705457480209210212359075,
			.407024119471837170130712624475966640701218,
			-.547887773099124699838423039613992609652666,
			.316687942349131283304243594632524222261229,
			-.199454449272498849220920238818740441773095,
			.155208177567700183717871295087047648274526,
			-.194872334236218107555685235739360994948363,
			.473182071290933425061584904357486045715166,
			-.646082900331347695307710329488322284833587,
			.991054273272769239147018346486949532020165,
			-1.476181396692545906866774028900451861728154,
			1.182909764903980913120233463732077220174485,
			-.118621923190321147172858215228057443861943,
			.617878268720229825361190580870947622440374,
			-2.371043542562823248107790279326182156140948,
			6.19652008105860789533941540398303330911738,
			-5.303214333080455591906418331155518924216553,
			1.648471970847190398887793720407399936261977,
			-.068057174159688090243430870215237526711285,
			.75374056490445931173113239275094491283215,
			-3.486664263728150950438696181446088443453779,
			4.248792491507461402502278320683982542109365,
			-1.310613085106969068034324784717157700118059,
			.793369622379913682095014367389960967889471,
			-.008103146406158448677872072208745688456267,
			.13009003179182798473397244567571284716706,
			-.2650845056331891676153610217798865664421,
			.243338988773724019892167834343839306342222,
			-.183043266536408743485913219928597265245541,
			.152233742213178067541032789451272614086765,
			-.102269416077909346007225985559845256692367,
			.343971723339104256150414807709503561378935,
			-.591871891042586996587956344324486149817975,
			.961551115977580082998759727556905756896287,
			-1.457814619278325323796892806006000397382759,
			1.176442565289709194841567721072300142017591,
			-.118621923190321147172858215228057443861943,
			.617878268720229825361190580870947622440374,
			-2.371043542562823248107790279326182156140948,
			6.19652008105860789533941540398303330911738,
			-5.303214333080455591906418331155518924216553,
			1.648471970847190398887793720407399936261977,
			-.068057174159688090243430870215237526711285,
			.75374056490445931173113239275094491283215,
			-3.486664263728150950438696181446088443453779,
			4.248792491507461402502278320683982542109365,
			-1.310613085106969068034324784717157700118059,
			.793369622379913682095014367389960967889471,
			-.030051710115129077581681735881115797900829,
			.299676540603301018486578557069466576120501,
			-.488529849551401933578800022789160073393629,
			.353310966819570175441792096154015460263621,
			-.236211359621476871935801920973725530852068,
			.171237256068110401555939781974114613214541,
			-.400264386734502240991337548051234253112989,
			1.100823585493004667164009295371979941952903,
			-1.619305068434174827560528963865929927105323,
			2.320010646477567734538041871804541455997053,
			-2.734480172388982218522418556895501384824969,
			1.663224873794658752970901022084521823493036,
			-.448103670160934358919717664801989887053879,
			2.186876411257981880603851997798769137212916,
			-7.193268552543644111264887439099565795229543,
			13.752350993262369737650885015228183622265358,
			-9.525572375919678486453013116200217939991001,
			1.897707715896333470784214086626443206396437,
			-.252464785650279997657634910018216718925469,
			2.275939001252218329666290219271889746940691,
			-7.275674128748307197371479253318499407022339,
			7.030306062218702959287448619294249714263117,
			-1.666069721977970908140078497268661074773157,
			.81853172870266310182742706648564249206502
		};

		double *AA = new double[120];
		memcpy( AA, a, 120*sizeof(double) );

		return AA;
	}
	double *PredictOrbitSat::InitParamA1(  )
	{
		double a1[40] = { 
			-.129592686941370908617580445973016746178066932,
			.434136323735473134894250613858187912416716491,
			-.987914370363320998878515984232082126334485011,
			3.1222452323109401203806598302813841909101614,
			-11.03032354412456396898374694970409280198656786,
			12.76368350616289255974939552218638585723743236,
			-4.56091436337124178335169068409572155402704975,
			1.388679902591191844807228097678955267961859311,
			-.7740805706922659413173329423028026512407,
			2.505009402606041554544778841453833835061934265,
			-5.222985356593715655316251130614533764771955814,
			13.463657546628178761222984014590961416953608795,
			-31.84807348570252746214763575966166045109197442,
			28.13870829136599135075908637560894613014930483,
			-6.78771587418565072066559081301806752957967985,
			1.52548004657394811291996141394332301451948885,
			-1.056669062373016843993506431485552299573017393,
			2.916376810847920606166822406190804564140333076,
			-5.452942157674383292161753818181573760198833782,
			13.70963751693707076752315808873074767289984713,
			-32.18681999563676740631375241153515125987660265,
			28.36303357195290876979977892973359985450772774,
			-6.819990734463271135742857335722437872421255793,
			1.527374050409538534722110572269563100521801669,
			-.5196195265988806688329845368426604489663069435,
			2.093641994364162502922735276716863105983535453,
			-4.93319440974187701693999064566133582326125689,
			13.1353678684411969135962104639967391049311701,
			-31.3879223029948817906640628959332231373290313,
			27.83211555849954051734355785472428233861065579,
			-6.7432408733108706240946595465126679834382256,
			1.522851691341610166669194029512002843469459436,
			-4.78857569763568338324577989535990533025099254,
			14.91650929423222885715245220487518638283883385,
			-28.3234521159288792013528788769006072528307845,
			59.48238716289572615351242534673187751077764969,
			-97.4142209916184941716918363352549169067123911,
			65.49673575703264377442135561348774797085427757,
			-10.0430608226989869296095281272137303379515464,
			1.673677413721444900813790069634347963274953637
		};

		double *AA = new double[40];
		memcpy( AA, a1, 40*sizeof(double) );

		return AA;
	}
	double *PredictOrbitSat::InitParamW(  )
	{
		double w[76] = { 
			-2.03619749319694527721136401200837535635761238e-4,
			.0100588730542575556312123736026889173993056912416,
			.0635267293502763979216722340379978335150493864926,
			-.004666083656234159297488993888391763703407783342,
			8.63632719419451458267354680823348627426076789e-4,
			-1.4768751542583879791507647832225085060029975e-4,
			9.3456715666352400398108898909515928942148745e-4,
			-.021570688555602401495459603015036206487043960196,
			.2336385874372896694805220892121749517615044386168,
			.1240916708842203308325637720806805944341829659965,
			-.008341212505381763261099670439591711866157354564,
			.0012565537903825080381594436210548692678044878091,
			-.00179498739163126748342402051092903967453146308,
			.03198322351851511893852050358658088820002723034,
			.10359328404101594826646869360339682512391444826,
			.40503828844000138420163573787259110527571875982,
			.138740439675533464863588978471615156586316583827,
			-.0075697264910065163854570134716325919111576243,
			.00181670318116950208515521365054774692717531406,
			-.02957990120959961979838519123433117243910299208,
			.23481182227693165179695003710953013979232039041,
			.26874615788534257980758269140772038808131788241,
			.37335814367559341460921357803540634855572088786,
			.08141522998758875911145691547753130163043120713,
			-4.5758876186858558853051394587863228494118844e-4,
			.0143225171615434229596069649750566353520110506,
			.05839817163491159240220997214130914575555779108,
			-.00329383165004849786947143072920712668978123114,
			5.539031431352273508991805566398244755027041e-4,
			-9.132732469944686668741744432459915621181602e-5,
			.002584710659658300908234274292720338183975473851,
			-.03528481247685728716961278802145243734518128206,
			.24949402149343216611150872865380554046668362757,
			.11940199032303122724112521299750859561766632657,
			-.0072396663824022019151094922247575805582069888,
			.00105323459070966242252118475055320003477490804,
			-.0051504257102334129573932868178140055876185482,
			.056074748301911536232602646980515043511348952978,
			.076091256214193581461859431935854289567385794703,
			.413427767371522744170846882439004129647665714171,
			.136747656558863835933957038358194925848259534881,
			-.00720048094383015244053983334413203938675351363,
			.00526769465051677804665744664871693972981241263,
			-.05329870538411934208474757280637193276080573006,
			.261774975095319788296983768307621196774706715291,
			.260437663320203594374334870198153881930670086967,
			.375338730781660135697725786946595528372010911491,
			.081047797333445333281018945151689138501468293487,
			-9.4972951312378740276782794532906807180730842e-4,
			.01917310610670488920609264961888605181651483581,
			.053199525126511786394570373909455573328276581466,
			-.00228784804694189907437209106289120652802088288,
			3.53291483300181211659515204581442321149304183e-4,
			-5.650095347745794715586417110754541397521996e-5,
			.006575746031982001372117094286342339531502085268,
			-.05559507851799395929940501638939937965476064727,
			.269370385669247311182987906228059322753765522879,
			.115114539224667756657349391543924912928104828062,
			-.00635371142020822554531027077996046116932964874,
			8.97597219876983230928015559410922010429924923e-4,
			-.01348113093893672123006165466901291289056671362,
			.094333906798463505019466372046272743426194953342,
			.039445655746544959419225867626311783836176787199,
			.421535856115575090568554984019105552323486996195,
			.135059404216850630795600631975215319713160952304,
			-.00690317014606933217145332144627014280816504053,
			.01389568490940985339398614283356428038598821427,
			-.09179920150628932515238141674307889485714658257,
			.298411953860218858594007278831222541063834300469,
			.252268892238495417242873216192016836257142631281,
			.377043459334854189898798094531319338155858061075,
			.080747366960337293634689928801360651542186065294,
			.173927422568726928686531974610999703617674,
			.326072577431273071313468025389000296382326,
			.326072577431273071313468025389000296382326,
			.173927422568726928686531974610999703617674 
		};

		double *AA = new double[76];
		memcpy( AA, w, 76*sizeof(double) );

		return AA;
	}
	double *PredictOrbitSat::InitParamG(  )
	{
		double g[6] = {
			0.069431844202973712388026755553595247452137,
			0.330009478207571867598667120448377656399712,
			0.669990521792428132401332879551622343600287,
			0.930568155797026287611973244446404752547862,
			1.069431844202973712388026755553595247452137,
			1.330009478207571867598667120448377656399712 
		};

		double *AA = new double[6];
		memcpy( AA, g, 6*sizeof(double) );

		return AA;
	}

	double *PredictOrbitSat::InitParamC(  )
	{
		double c[16] = { 
			-1.855135017009736526300160598113958365545275481,
			4.698862351888765202815428907536820553842760509,
			-4.698862351888765202815428907536820553842760509,
			1.855135017009736526300160598113958365545275481,
			4.775286118057296050988873551673920880443110091,
			-10.4627407878153534040194101705196452570380801,
			8.33270861973970740724230545962763695833296188,
			-2.64525394998165005421176884078191258173799183,
			-4.273011803936099382986509257099496330211525014,
			6.90358346284478853307934737014000804337764127,
			-3.708535210731319537913690303801995595319963883,
			1.077963551822630387820852190761483882153847627,
			1.5267881254572667869843282781505335189313647521,
			-0.8136324494869272605618980817681830437999959804,
			0.4007615203116504048002817771001794872120881554,
			-0.1139171962819899312227119734825299623434569271 
		};
		double *AA = new double[16];
		memcpy( AA, c, 16*sizeof(double) );

		return AA;
	}
	//==============================================================================//
	// Get constant
	//==============================================================================//
	void PredictOrbitSat::InitIntegrationConstant()
	{
		printf("Get constant   ....");
		// получение констант
		IC_a = InitParamA();
		IC_a1 = InitParamA1();
		IC_w = InitParamW();
		IC_g = InitParamG();
		IC_c = InitParamC();
		printf("OK\n" );
	}
	//==============================================================================//
	// delete constant
	//==============================================================================//
	void PredictOrbitSat::DeleteIntegrationConstant()
	{
		delete IC_a;
		delete IC_a1;
		delete IC_w;
		delete IC_g;
		delete IC_c;
	};
	//==============================================================================//
}