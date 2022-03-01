#include <algorithm>
#include <cmath>

#include "MotionCalculator.h"

MountSpecification MountSpecification::operator=(const MountSpecification& mount_spec_) {
	acceleration_max = mount_spec_.acceleration_max;
	velocity_max = mount_spec_.velocity_max;
	gear_ratio = mount_spec_.gear_ratio;
	divider = mount_spec_.divider;
	steps_in_circle_motor = mount_spec_.steps_in_circle_motor;

	return *this;
}

MotionCalculator::MotionCalculator(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor) {
	SetMountSpecification(a_max, V_max, gear_ratio, divider, steps_in_circle_motor);

	rads_in_step = 2 * M_PI / (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio);
	steps_in_rad = (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio) / (2 * M_PI);
}

MotionCalculator::MotionCalculator(const MountSpecification& mount_spec_) {
	SetMountSpecification(mount_spec_);

	rads_in_step = 2 * M_PI / (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio);
	steps_in_rad = (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio) / (2 * M_PI);
}

MotionCalculator::~MotionCalculator() {
}

void MotionCalculator::SetMountSpecification(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor) {
	mount_spec.acceleration_max = a_max;
	mount_spec.velocity_max = V_max;
	mount_spec.gear_ratio = gear_ratio;
	mount_spec.divider = divider;
	mount_spec.steps_in_circle_motor = steps_in_circle_motor;

	rads_in_step = 2 * M_PI / (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio);
	steps_in_rad = (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio) / (2 * M_PI);
}

void MotionCalculator::SetMountSpecification(const MountSpecification& mount_spec_) {
	mount_spec = mount_spec_;

	rads_in_step = 2 * M_PI / (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio);
	steps_in_rad = (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio) / (2 * M_PI);
}

MountSpecification MotionCalculator::GetMountSpecification() const {
	return mount_spec;
}

double MotionCalculator::Steps2Rad(double steps) const {
	return steps * rads_in_step;
}

double MotionCalculator::Rad2Steps(double rad) const {
	return rad * steps_in_rad;
}



MountTask MotionCalculator::CalculatePoint2PointTask(const Angs& startOwnAxes, const Angs& endOwnAxes) const {
	MountTask task;

	double steps_mot1 = Rad2Steps(endOwnAxes.ang1 - startOwnAxes.ang1);
	task.expected_steps_mot1 = steps_mot1;
	double steps_mot2 = Rad2Steps(endOwnAxes.ang2 - startOwnAxes.ang2);
	task.expected_steps_mot2 = steps_mot2;


	double sign_mot1 = sign(steps_mot1);
	double sign_mot2 = sign(steps_mot2);
	steps_mot1 = abs(steps_mot1);
	steps_mot2 = abs(steps_mot2);

	task.startOwnAxes = startOwnAxes;
	task.endOwnAxes = endOwnAxes;

	task.acceleration_mot1 = Steps2Rad(sign_mot1 * mount_spec.acceleration_max);
	task.acceleration_mot2 = Steps2Rad(sign_mot2 * mount_spec.acceleration_max);

	double time_for_acceleration_and_slowing = 2 * (mount_spec.velocity_max / mount_spec.acceleration_max);
	double steps_for_acceleration_and_slowing = 2 * (sqr(mount_spec.velocity_max) / (2 * mount_spec.acceleration_max));

	double moving_time_mot1;
	if (steps_mot1 > steps_for_acceleration_and_slowing) {
		task.velocity_mot1 = Steps2Rad(sign_mot1 * mount_spec.velocity_max);
		moving_time_mot1 = time_for_acceleration_and_slowing + (steps_mot1 - steps_for_acceleration_and_slowing) / mount_spec.velocity_max;

		task.task_mot1 = FillPoint2PointDriveTask({ Steps2Rad(sign_mot1 * steps_for_acceleration_and_slowing / 2), 
													time_for_acceleration_and_slowing / 2 },
												  { Steps2Rad(sign_mot1 * steps_mot1 - sign_mot1 * steps_for_acceleration_and_slowing / 2),
													moving_time_mot1 - time_for_acceleration_and_slowing / 2 });
	}
	else {
		double velocity = sqrt(2 * (steps_mot1 / 2) * mount_spec.acceleration_max);
		task.velocity_mot1 = Steps2Rad(sign_mot1 * velocity);
		moving_time_mot1 = 2 * velocity / mount_spec.acceleration_max;
	}

	double moving_time_mot2;
	if (steps_mot2 > steps_for_acceleration_and_slowing) {
		task.velocity_mot2 = Steps2Rad(sign_mot2 * mount_spec.velocity_max);
		moving_time_mot2 = time_for_acceleration_and_slowing + (steps_mot2 - steps_for_acceleration_and_slowing) / mount_spec.velocity_max;

		task.task_mot2 = FillPoint2PointDriveTask({ Steps2Rad(sign_mot2 * steps_for_acceleration_and_slowing / 2),
													time_for_acceleration_and_slowing / 2 },
												  { Steps2Rad(sign_mot2 * steps_mot2 - sign_mot2 * steps_for_acceleration_and_slowing / 2),
													moving_time_mot2 - time_for_acceleration_and_slowing / 2 });
	}
	else {
		double velocity = sqrt(2 * (steps_mot2 / 2) * mount_spec.acceleration_max);
		task.velocity_mot2 = Steps2Rad(sign_mot2 * velocity);
		moving_time_mot2 = 2 * velocity / mount_spec.acceleration_max;
	}

	task.moving_time_mot1 = moving_time_mot1;
	task.moving_time_mot2 = moving_time_mot2;

	return task;
}

std::vector<TaskPoint> MotionCalculator::FillPoint2PointDriveTask(const TaskPoint& start_point, const TaskPoint& end_point) const {
	vector<TaskPoint> result;
	result.push_back({ start_point.ang * 180 / M_PI, start_point.time * 1000000 });
	double velosity = (end_point.ang - start_point.ang) / (end_point.time - start_point.time);
	
	size_t n = std::floor(end_point.time - start_point.time);
	for (size_t i = 0; i < n; ++i) {
		result.push_back({ (start_point.ang + velosity * (i + 1)) * 180 / M_PI,
						   (start_point.time + (i + 1)) * 1000000 });
	}

	if (end_point.time - start_point.time > std::floor(end_point.time - start_point.time) + 0.0000001) {
		result.push_back({ end_point.ang * 180 / M_PI, end_point.time * 1000000 });
	}
	
	return result;
}