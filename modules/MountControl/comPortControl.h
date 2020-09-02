#pragma once
#include <windows.h>
#include <string>
#include <vector>

#define STATUS_MOVE 1
#define STATUS_ERROR 2
#define STATUS_COMPLETE 3


class ComPortControl
{

private:
	HANDLE hSerial;
	bool PortOpened = false;
	COMMTIMEOUTS CommTimeOuts;

	std::vector<std::string> split(const std::string& s, char delimiter);
public:
	ComPortControl();
	~ComPortControl();

	// Open
	int OpenCOM(LPCTSTR sPortName);

	// Send command
	bool WriteCommand(int pos1, int pos2);

	// Read responce
	std::string ReadResponce();

	// parse result string
	int parseResponce(std::string result);

	int getPosition();

	int runCommand(int pos1, int pos2);
};

