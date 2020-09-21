#pragma once

#include "NewConvertor.h"

struct MountSpecification {
	double acceleration_max = 2000.0;
	double velocity_max = 10000.0;
	double gear_ratio = 720.0;
	double divider = 16.0;
	double steps_in_circle_motor = 200.0;

	MountSpecification operator=(const MountSpecification& mount_spec_);
};

struct Trajectory {
	double moving_time_mot1;	// sec
	double moving_time_mot2;	// sec
	double velocity_mot1;		// rad
	double velocity_mot2;		// rad
	double acceleration_mot1;	// rad
	double acceleration_mot2;	// rad
	Angs startOwnAxes;			// rad
	Angs endOwnAxes;			// rad
};

class MotionCalculator {
public:
	MotionCalculator(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor);
	MotionCalculator(const MountSpecification& mount_spec_);
	~MotionCalculator();

	void SetMountSpecification(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor);
	void SetMountSpecification(const MountSpecification& mount_spec_);

	MountSpecification GetMountSpecification() const;

	Trajectory CalculateTrajectP2P(const Angs& startOwnAxes, const Angs& endOwnAxes) const;

private:
	MountSpecification mount_spec;
};