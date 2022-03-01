#pragma once
#include <memory>
#include "Eigen/Dense"

#include "NewConvertor.h"
#include "MotionCalculator.h"

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

class MountController : public MotionCalculator {
public:
	MountController(std::shared_ptr<NewConvertor> convertor_, const MountSpecification& mount_spec_, double* pos_ITRF_);
	MountController(std::shared_ptr<NewConvertor> convertor_, const MountSpecification& mount_spec_, double Lat, double Lon, double Elev);
	~MountController();

	void SetCurrentDirectionRaDec(double Jd, const Angs& RaDec);
	void SetCurrentDirectionAzElev(double Jd, const Angs& AzElev);
	void SetCurrentDirectionOwnAxes(double Jd, const Angs& OwnAxes);

	void SetCurrentPositionLatLonElev(double Lat, double Lon, double Elev);
	void SetCurrentPositionLatLonElev(double* pos);
	void SetCurrentPositionITRF(double x, double y, double z);
	void SetCurrentPositionITRF(double* pos);

	void SetCalibrationMatrix(const Eigen::Matrix3d& R);

	CurrentPosition GetCurrentPosition() const;
	CurrentDirection GetCurrentDirection() const;
	Eigen::Matrix3d GetCalibrationMatrix() const;

private:
	std::shared_ptr<NewConvertor> convertor;
	CurrentPosition cur_pos;
	CurrentDirection cur_dir;
	Eigen::Matrix3d R;
	Eigen::Matrix3d InvR;

	Angs OwnAxes2AzElev(double Jd, const Angs& OwnAxes, const on_surface& pos, const Eigen::Matrix3d& InvR) const;
	Angs AzElev2OwnAxes(double Jd, const Angs& AzElev, const on_surface& pos, const Eigen::Matrix3d& R) const;
	Angs OwnAxes2RaDec(double Jd, const Angs& OwnAxes, const on_surface& pos, const Eigen::Matrix3d& InvR) const;
	Angs RaDec2OwnAxes(double Jd, const Angs& RaDec, const on_surface& pos, const Eigen::Matrix3d& R) const;
};
