#include <algorithm>

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
}

MotionCalculator::MotionCalculator(const MountSpecification& mount_spec_) {
	SetMountSpecification(mount_spec_);
}

MotionCalculator::~MotionCalculator() {
}

void MotionCalculator::SetMountSpecification(double a_max, double V_max, double gear_ratio, double divider, double steps_in_circle_motor) {
	mount_spec.acceleration_max = a_max;
	mount_spec.velocity_max = V_max;
	mount_spec.gear_ratio = gear_ratio;
	mount_spec.divider = divider;
	mount_spec.steps_in_circle_motor = steps_in_circle_motor;
}

void MotionCalculator::SetMountSpecification(const MountSpecification& mount_spec_) {
	mount_spec = mount_spec_;
}

MountSpecification MotionCalculator::GetMountSpecification() const {
	return mount_spec;
}

Trajectory MotionCalculator::CalculateTrajectP2P(const Angs& startOwnAxes, const Angs& endOwnAxes) const {
	Trajectory trajectory;

	double step_in_rad = 2 * M_PI / (mount_spec.steps_in_circle_motor * mount_spec.divider * mount_spec.gear_ratio);



	double steps_mot1 = (endOwnAxes.ang1 - startOwnAxes.ang1) / step_in_rad;
	double steps_mot2 = (endOwnAxes.ang2 - startOwnAxes.ang2) / step_in_rad;

	double sign_mot1 = sign(steps_mot1);
	double sign_mot2 = sign(steps_mot2);
	steps_mot1 = abs(steps_mot1);
	steps_mot2 = abs(steps_mot2);

	trajectory.startOwnAxes = startOwnAxes;
	trajectory.endOwnAxes = endOwnAxes;

	trajectory.acceleration_mot1 = (sign_mot1 * mount_spec.acceleration_max) * step_in_rad;
	trajectory.acceleration_mot2 = (sign_mot2 * mount_spec.acceleration_max) * step_in_rad;

	double time_for_acceleration_and_slowing = 2 * (mount_spec.velocity_max / mount_spec.acceleration_max);
	double steps_for_acceleration_and_slowing = 2 * (sqr(mount_spec.velocity_max) / (2 * mount_spec.acceleration_max));

	double moving_time_mot1;
	if (steps_mot1 > steps_for_acceleration_and_slowing) {
		trajectory.velocity_mot1 = (sign_mot1 * mount_spec.velocity_max) * step_in_rad;
		moving_time_mot1 = time_for_acceleration_and_slowing + (steps_mot1 - steps_for_acceleration_and_slowing) / mount_spec.velocity_max;
	}
	else {
		double velocity = sqrt(2 * (steps_mot1 / 2) * mount_spec.acceleration_max); 
		trajectory.velocity_mot1 = (sign_mot1 * velocity) * step_in_rad;
		moving_time_mot1 = 2 * velocity / mount_spec.acceleration_max;
	}

	double moving_time_mot2;
	if (steps_mot2 > steps_for_acceleration_and_slowing) {
		trajectory.velocity_mot2 = (sign_mot2 * mount_spec.velocity_max) * step_in_rad;
		moving_time_mot2 = time_for_acceleration_and_slowing + (steps_mot2 - steps_for_acceleration_and_slowing) / mount_spec.velocity_max;
	}
	else {
		double velocity = sqrt(2 * (steps_mot2 / 2) * mount_spec.acceleration_max);
		trajectory.velocity_mot2 = (sign_mot2 * velocity) * step_in_rad;
		moving_time_mot2 = 2 * velocity / mount_spec.acceleration_max;
	}

	trajectory.moving_time_mot1 = moving_time_mot1 * step_in_rad;
	trajectory.moving_time_mot2 = moving_time_mot2 * step_in_rad;

	return trajectory;
}