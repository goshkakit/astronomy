#pragma once
#include "NewConvertor.h"

struct CurrentDirection {
	Angs RaDec;
	Angs AzElev;
	Angs OwnAxes;
	double Jd;
};

struct CurrentPosition {
	double pos_ITRF[3];
	on_surface ang_pos;
};

struct MountSpec {
	double a_max = 2000.0;
	double V_max = 10000.0;
	double gear_ratio = 720.0;
	double divider = 16.0;
	double steps_in_circle_motor = 200.0;
};

class MountController {
public:
	MountController(std::unique_ptr<NewConvertor> convertor_, MountSpec mount_spec_, double* pos_ITRF_);
	MountController(std::unique_ptr<NewConvertor> convertor_, MountSpec mount_spec_, double Lat, double Lon, double Elev);
	~MountController();

	void SetCurrentDirectionRaDec(double Jd, Angs RaDec);
	void SetCurrentDirectionAzElev(double Jd, Angs AzElev);
	void SetCurentDirectionOwnAxes(double Jd, Angs OwnAxes);

	void SetCurrentPositionLatLonElev(double Lat, double Lon, double Elev);
	void SetCurrentPositionLatLonElev(double* pos);
	void SetCurrentPositionITRF(double x, double y, double z);
	void SetCurrentPositionITRF(double* pos);

	void SetMountSpec(double a_max, double V_max, double gear_ratio = 720.0, double divider = 16.0, double steps_in_circle_motor = 200.0);
	void SetMountSpec(MountSpec mount_spec_);

	void SetCalibrationMatrix(Matrix R);

	CurrentPosition GetCurrentPosition();
	CurrentDirection GetCurrentDirection();
	MountSpec GetMountSpec();
	Matrix GetCalibrationMatrix();

private:
	std::unique_ptr<NewConvertor> convertor;
	CurrentPosition cur_pos;
	CurrentDirection cur_dir;
	MountSpec mount_spec;
	Matrix R;
	Matrix InvR;

	Angs OwnAxes2AzElev(double Jd, Angs OwnAxes, on_surface pos, Matrix InvR);
	Angs AzElev2OwnAxes(double Jd, Angs AzElev, on_surface pos, Matrix R);
	Angs OwnAxes2RaDec(double Jd, Angs OwnAxes, on_surface pos, Matrix InvR);
	Angs RaDec2OwnAxes(double Jd, Angs RaDec, on_surface pos, Matrix R);
};
