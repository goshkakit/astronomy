#include "NucleoConnection.h"

#include "AudioFile.h"

uint16_t NucleoConnection::Crc_NextCrc16(uint16_t crc, uint8_t data)
{
	int i;
	crc = crc ^ ((uint16_t)data << 8);
	for (i = 0; i<8; i++)
	{
		if (crc & 0x8000)
			crc = (crc << 1) ^ 0x1021;
		else
			crc <<= 1;
	}

	return crc;
}
uint16_t NucleoConnection::GetCRC_16(const uint8_t *pcBlock, unsigned int len)
{
	unsigned int crc = 0;
	while (len--)
	{
		crc = Crc_NextCrc16(crc, *pcBlock++);
	}
	return crc;
}

StatusSuccess NucleoConnection::loadTask()
{
	StatusSuccess res;

	// verify input
	std::ifstream source("orbit_steps_az.txt", std::ios_base::in);
	int count = 0;
	double pt, az, N_az, Ni_az, stepT0, dstepTaz, maxErrMi_Az, stepTiaz_error;
	while (source >> pt >> az >> N_az >> Ni_az >> stepT0 >> dstepTaz >> maxErrMi_Az >> stepTiaz_error)
	{
		count++;
	}
	res.message = "TASK:\n";
	res.message += "Read points: " + std::to_string(count) + "\n";
	source.close();

	//read data
	taskData.Task_N.resize(count);
	taskData.Task_stepT.resize(count);
	taskData.Task_Angle.resize(count);
	taskData.Task_T.resize(count);
	taskData.count = count;

	std::ifstream ftask("orbit_steps_az.txt", std::ios_base::in);
	int index = 0;
	while (ftask >> pt >> az >> N_az >> Ni_az >> stepT0 >> dstepTaz >> maxErrMi_Az >> stepTiaz_error)
	{
		//Ni_az = 500;
		taskData.Task_T[index] = pt;
		taskData.Task_N[index] = abs(Ni_az);
		taskData.Task_stepT[index] = stepT0;
		//taskData.Task_stepT[index] = 1000000.0 / Ni_az;
		taskData.Task_Angle[index] = az;
		index++;
		if (index == count) break;
	}
	res.message += "Read task points: " + std::to_string(index) + "\n";
	ftask.close();

	taskData.crc1 = GetCRC_16((uint8_t*)(&taskData.Task_N[0]), taskData.Task_N.memorySize());
	taskData.crc2 = GetCRC_16((uint8_t*)(&taskData.Task_stepT[0]), taskData.Task_stepT.memorySize());

	res.message += "CRC-1: " + std::to_string(taskData.crc1) + "\n";
	res.message += "CRC-2: " + std::to_string(taskData.crc2) + "\n";
	res.message += "Count: " + std::to_string(taskData.count) + "\n";

	// make header;
	taskData.header[0] = 1;
	taskData.header[1] = taskData.count;
	taskData.header[2] = taskData.crc1;
	taskData.header[3] = taskData.crc2;

	res.success = 0;
	return res;
};

WavData NucleoConnection::LoadWavFile(std::string path)
{
	AudioFile<float> audioFile;
	WavData wd;
	if (audioFile.load(path) == false)
	{
		std::cout << "Error load WAV file" << std::endl;
		return wd;
	}

	std::cout << "BitDepth: " << audioFile.getBitDepth() << std::endl;
	std::cout << "LengthInSeconds: " << audioFile.getLengthInSeconds() << std::endl;
	std::cout << "NumChannels: " << audioFile.getNumChannels() << std::endl;
	std::cout << "SampleRate: " << audioFile.getSampleRate() << std::endl;
	std::cout << "Nsample: " << audioFile.samples[0].size() << std::endl;

	wd.BitDepth = audioFile.getBitDepth();
	wd.LengthInSeconds = audioFile.getLengthInSeconds();
	wd.NumChannels = audioFile.getNumChannels();
	wd.SampleRate = audioFile.getSampleRate();
	wd.Nsample = audioFile.samples[0].size();

	wd.signalA = new float[wd.Nsample];
	wd.signalB = new float[wd.Nsample];

	for (int i = 0; i < wd.Nsample; i++)
	{
		wd.signalA[i] = audioFile.samples[0][i];
		wd.signalB[i] = audioFile.samples[1][i];
	}
	return wd;
}

StatusSuccess NucleoConnection::runTest()
{
	StatusSuccess res;

	WavData wdata = LoadWavFile("dump-12.wav");
	

	float vmin = 0;
	float vmax = 0;
	for (int i = 0; i < wdata.Nsample; i++)
	{
		if (i == 0)
		{
			vmin = wdata.signalA[0];
			vmax = wdata.signalA[0];
		}

		if( vmin > wdata.signalA[i] ) vmin = wdata.signalA[i];
		if( vmax < wdata.signalA[i] ) vmax = wdata.signalA[i];
	}

	std::cout << "vmin: " << vmin << std::endl;
	std::cout << "vmax: " << vmax << std::endl;

	std::vector<double> Timpulse;
	bool PulseEnable = false;
	bool findStart = false;
	double iStart = 0;
	double impulseCount = 0;

	int countError = 0;
	for (int i = 1; i < wdata.Nsample; i++)
	{
		double val = wdata.signalA[i];
		if (val > 0.05 && !findStart)
		{
			findStart = true;
			iStart = i;

			// time start  ms
			double t = i;
			t = t / (double)wdata.SampleRate * 1000.0;
			std::cout << "T start: " << t << " ms." << std::endl;
		}

		// start impulse
		if (val > 0.05 && !PulseEnable)
		{
			PulseEnable = true;

			// time impulse start  ms
			double t = i - iStart;
			double T = t / (double)wdata.SampleRate;

			// verify interval
			if (Timpulse.size() > 2)
			{
				int id = Timpulse.size() - 1;
				double Tprev = Timpulse[id] - Timpulse[id-1];
				double Tcurr = T - Timpulse[id];
				double dT = (Tcurr - Tprev);
				if (abs(dT) > 0.2 * Tprev)
				{
					countError++;
					printf("Error impulse step N= %d T= %f prev= %f curr= %f dt= %f\n", (int)impulseCount, T* 1000.0, Tprev* 1000.0, Tcurr* 1000.0, dT* 1000.0);
				}
			}

			Timpulse.push_back(T);

			impulseCount +=1;

			//if (impulseCount > 40000)
			//{
			//	break;
			//}
		}

		if (val < 0.02 )
		{
			PulseEnable = false;
		}
	}
	printf("countError = %d\n", countError);

	// Freq Test
	
	//int nSec = 8;
	//int step = 500* nSec;
	//double Tsample = 1000000.0 / wdata.SampleRate;
	//printf("Tsample=%f\n", Tsample);
	//for (int i = step; i < Timpulse.size(); i += step )
	//{
	//	double T1 = Timpulse[i - step];
	//	double T2 = Timpulse[i];
	//	double dT = ((T2 - T1) - 1.0*nSec) * 1000000.0;
	//	printf("%d impulse = %f\n", step, dT);
	//}
	//return res;
	

	// make angle chart
	int _Nglobal = 700;
	int _StepperCount = 200;
	int _StepDivider = 16;
	double _StepSizeDeg = 360.0 / (double)_Nglobal / (double)_StepperCount / (double)_StepDivider;
	double _Angle0 = 192.157916891331;
	double _T0 = 23888.000000;
	double _Dir = -1.0;

	// идем по заданию
	int commonSteps = 0;

	
	std::vector<float> angleList;
	std::vector<float> timeList;
	int CountVerifyPoint = 0;
	for (int i = 0; i < taskData.Task_N.count(); i++)
	{
		// берем время для этого колличества шагов
		double _T = _T0 + Timpulse[commonSteps];
		// На какой угол мы повернулись
		double _commonSteps = (double)commonSteps;
		double _Anglei =  _Angle0 + (_Dir * _StepSizeDeg)*_commonSteps;
		angleList.push_back(_Anglei);
		timeList.push_back(_T);

		double dangle = taskData.Task_Angle[i] - _Anglei;
		dangle = dangle / _StepSizeDeg;

		double eT = ( taskData.Task_T[i] - _T)*1000000.0;

		printf("%lf %lf %lf\t  %.12lf %.12lf\t %lf\n", _T, taskData.Task_T[i], eT, _Anglei, taskData.Task_Angle[i], dangle);
		CountVerifyPoint++;

		// прибавляем число шагов на следующий интервал
		commonSteps += taskData.Task_N[i];
		if (commonSteps >= Timpulse.size())
		{
			break;
		}
	}

	printf("CountVerifyPoint= %d\n", CountVerifyPoint);

	//std::ofstream source("angle.txt", std::ios_base::out);
	//for (int i = 0; i < angleList.size(); i++)
	//{
	//	std::cout << timeList[i] << " " << angleList[i] << " " << taskData.Task_Angle[i] << std::endl;
	//}
	//source.close();
	

	return res;
}