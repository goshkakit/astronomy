#include "comPortControl.h"

#include <sstream>
#include <iostream>
#include "stdio.h"

ComPortControl::ComPortControl() {}

ComPortControl::~ComPortControl() {}

std::vector<std::string> ComPortControl::split(const std::string& s, char delimiter)
{
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, delimiter))
	{
		tokens.push_back(token);
	}
	return tokens;
}

int ComPortControl::OpenCOM(LPCTSTR sPortName)
{
	hSerial = ::CreateFile(sPortName, GENERIC_READ | GENERIC_WRITE, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);

	if (hSerial == INVALID_HANDLE_VALUE)
	{
		if (GetLastError() == ERROR_FILE_NOT_FOUND)
		{
			std::cout << "serial port does not exist.\n";
			return -1;
		}
		std::cout << "some other error occurred.\n";
		return -1;
	}

	DCB dcbSerialParams = { 0 };
	dcbSerialParams.DCBlength = sizeof(dcbSerialParams);
	if (!GetCommState(hSerial, &dcbSerialParams))
	{
		std::cout << "getting state error\n";
	}
	dcbSerialParams.BaudRate = CBR_9600;
	dcbSerialParams.ByteSize = 8;
	dcbSerialParams.StopBits = ONESTOPBIT;
	dcbSerialParams.Parity = NOPARITY;
	if (!SetCommState(hSerial, &dcbSerialParams))
	{
		std::cout << "error setting serial port state\n";
		return -1;
	}
	PortOpened = true;

	CommTimeOuts.ReadIntervalTimeout = 200;
	CommTimeOuts.ReadTotalTimeoutMultiplier = 1;
	// значений этих тайм – аутов вполне хватает для уверенного приема 
	// даже на скорости 110 бод 
	CommTimeOuts.ReadTotalTimeoutConstant = 200;
	// используется в данном случае как время ожидания посылки 
	CommTimeOuts.WriteTotalTimeoutMultiplier = 0;
	CommTimeOuts.WriteTotalTimeoutConstant = 0;
	SetCommTimeouts(hSerial, &CommTimeOuts);

	return 0;
}

// Send command
bool ComPortControl::WriteCommand(int pos1, int pos2)
{
	if (PortOpened)
	{
		std::string command = std::to_string(pos1) + " " + std::to_string(pos2) + "\n";

		DWORD dwSize = command.size();	// размер этой строки
		DWORD dwBytesWritten;			// тут будет количество собственно переданных байт
		BOOL iRet = WriteFile(hSerial, command.c_str(), dwSize, &dwBytesWritten, NULL);
		if (dwSize != dwBytesWritten)
		{
			std::cout << "Send: " << dwSize << " Bytes in string. " << dwBytesWritten << " Bytes sended. " << std::endl;
			return false;
		}
	}
	else
	{
		std::cout << "PortOpened = false!" << std::endl;
		return false;
	}
}

// Read responce
std::string ComPortControl::ReadResponce()
{
	std::string result = "";
	if (PortOpened)
	{
		DWORD iSize;
		char sReceivedChar;
		//std::cout << "Recive: ";
		while (1)
		{
			// получаем 1 байт
			ReadFile(hSerial, &sReceivedChar, 1, &iSize, NULL);
			// если что-то принято, выводим
			if (iSize > 0)
			{
				//std::cout << sReceivedChar;
				result += std::string(1, sReceivedChar);
			}
			if (iSize > 0 && sReceivedChar == '\n')
			{
				result += std::string(1, sReceivedChar);
				break;
			}
			if (iSize == 0)
			{
				//std::cout << "timeout" << std::endl;
				result = "timeout:";
				break;
			}
		}
	}
	return result;
}

// parse result string
int ComPortControl::parseResponce(std::string result)
{
	if (result.size() > 0)
	{
		std::vector<std::string> arr = split(result, ' ');
		if (arr.size() == 3)
		{
			if (arr[0] == "move:")
			{
				int pos1 = std::atoi(arr[1].c_str());
				int pos2 = std::atoi(arr[2].c_str());
				printf("MOVE: %d %d\n", pos1, pos2);
				return STATUS_MOVE;
			}
			else if (arr[0] == "pos:")
			{
				int pos1 = std::atoi(arr[1].c_str());
				int pos2 = std::atoi(arr[2].c_str());
				printf("POS: %d %d\n", pos1, pos2);
				return STATUS_COMPLETE;
			}
			else
			{
				//printf("ERROR command\n");
				return STATUS_ERROR;
			}
		}
		else if (arr[0] == "timeout:")
		{
			//printf("ERROR timeout\n");
			return STATUS_ERROR;
		}
		else
		{
			//printf("ERROR count\n");
			return STATUS_ERROR;
		}
	}
	else
	{
		//printf("ERROR free\n");
		return STATUS_ERROR;
	}
}

int ComPortControl::getPosition()
{
	// service command, return position
	WriteCommand(1, 1);

	// wait recive command by mount
	std::string result = ReadResponce();
	//std::cout << "Recive: " << result;
	parseResponce(result);

	return 0;
}

int ComPortControl::runCommand(int pos1, int pos2)
{
	WriteCommand(pos1, pos2);

	// wait recive command by mount
	std::string result = ReadResponce();
	//std::cout << "Recive: " << result;
	parseResponce(result);

	// wait end run mount
	while (1) // оценка времени исполнения
	{
		result = ReadResponce();
		if (STATUS_COMPLETE == parseResponce(result))
		{
			printf("COMPLETE\n");
			break;
		}
		else
		{
			Sleep(100);
		}
	}

	return STATUS_COMPLETE;
}

int test_comPortControl(int argc, char** argv)
{
	printf("Connect\n");
	// Init
	ComPortControl CPC;
	int res = CPC.OpenCOM(L"COM4");
	if (res != 0)
	{
		return 0;
	}
	// Wait
	Sleep(2000);
	printf("Start stepper control!\n");
	// Read data after connect
	std::string result = CPC.ReadResponce();
	std::cout << "Recive: " << result << std::endl;

	CPC.getPosition();

	printf("\nSet 2 absolute position:\n");
	int k = 0;
	while (k < 100)
	{
		int pos1 = 0;
		int pos2 = 0;
		scanf("%d %d", &pos1, &pos2);
		CPC.runCommand(pos1, pos2);
		printf("\nSet 2 absolute position:\n");
		k++;
	}
}

