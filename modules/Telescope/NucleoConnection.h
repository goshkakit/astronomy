#pragma once

#include "cpuptr.h"
#include <iostream>
#include <fstream>

struct WavData
{
	float *signalA = NULL;
	float *signalB = NULL;
	int BitDepth;
	int LengthInSeconds;
	int NumChannels;
	int SampleRate;
	int Nsample;
};

struct xyPoint
{
	double x;
	double y;
};

struct StatusSuccess
{
	int success = 0;
	std::string message;
};

class TaskData
{
public:
	unsigned count;
	unsigned int crc1, crc2;

	CPU_Ptr<unsigned int> header;
	CPU_Ptr<int> Task_N;
	CPU_Ptr<double> Task_stepT;
	CPU_Ptr<double> Task_Angle;
	CPU_Ptr<double> Task_T;

	TaskData()
	{
		header.resize(4);
	}

	~TaskData()
	{

	}
};


class NucleoConnection
{

public:
	NucleoConnection() {}
	~NucleoConnection() {}
	StatusSuccess loadTask();
	StatusSuccess runTest();

private:

	TaskData taskData;

	uint16_t Crc_NextCrc16(uint16_t crc, uint8_t data);
	uint16_t GetCRC_16(const uint8_t *pcBlock, unsigned int len);

	WavData LoadWavFile(std::string path);
};
