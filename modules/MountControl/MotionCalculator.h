#pragma once

#include "NewConvertor.h"

struct MountSpecification {
	double acceleration_max = 2000.0;
	double velocity_max = 10000.0;
	double gear_ratio = 700.0;
	double divider = 16.0;
	double steps_in_circle_motor = 200.0;

	MountSpecification operator=(const MountSpecification& mount_spec_);
};

struct TaskPoint {
	double ang;		// deg
	double time;	// mks (e-6)
};

struct MountTask {
	double moving_time_mot1;	// s
	double moving_time_mot2;	// s
	double velocity_mot1;		// rad
	double velocity_mot2;		// rad
	double acceleration_mot1;	// rad
	double acceleration_mot2;	// rad
	Angs startOwnAxes;			// rad
	Angs endOwnAxes;			// rad

	std::vector<TaskPoint> task_mot1;
	std::vector<TaskPoint> task_mot2;

	double expected_steps_mot1;
	double expected_steps_mot2;
};

class MotionCalculator {
public:
	MotionCalculator(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor);
	MotionCalculator(const MountSpecification& mount_spec_);
	~MotionCalculator();

	void SetMountSpecification(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor);
	void SetMountSpecification(const MountSpecification& mount_spec_);

	MountSpecification GetMountSpecification() const;

	double Steps2Rad(double steps) const;
	double Rad2Steps(double rad) const;

	MountTask CalculatePoint2PointTask(const Angs& startOwnAxes, const Angs& endOwnAxes) const;

private:
	MountSpecification mount_spec;
	double steps_in_rad;
	double rads_in_step;

	std::vector<TaskPoint> FillPoint2PointDriveTask(const TaskPoint& start_point, const TaskPoint& end_point) const;
};