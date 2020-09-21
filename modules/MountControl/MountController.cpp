#include "MountController.h"

MountController::MountController(std::unique_ptr<NewConvertor> convertor_, const MountSpecification& mount_spec_, double* pos_ITRF_) :
	MotionCalculator(mount_spec_),
	convertor(std::move(convertor_))
{
	SetCurrentPositionITRF(pos_ITRF_);
	cur_pos.ang_pos.pressure = 1000;
	cur_pos.ang_pos.temperature = 0;

	R.fill(0);
	R(0, 2) = 1;
	R(1, 1) = -1;
	R(2, 0) = 1;

	InvR = R.inverse();
}

MountController::MountController(std::unique_ptr<NewConvertor> convertor_, const MountSpecification& mount_spec_, double Lat, double Lon, double Elev) :
	MotionCalculator(mount_spec_),
	convertor(std::move(convertor_))
{
	double pos_LatLonElev[3] = { Lat, Lon, Elev };
	SetCurrentPositionLatLonElev(pos_LatLonElev);
	cur_pos.ang_pos.pressure = 1000;
	cur_pos.ang_pos.temperature = 0;

	R.fill(0);
	R(0, 2) = 1;
	R(1, 1) = -1;
	R(2, 0) = 1;

	InvR = R.inverse();
}

MountController::~MountController() {
}

void MountController::SetCurrentDirectionRaDec(double Jd, const Angs& RaDec) {
	Angs AzElev = convertor->RaDec2AzElev(Jd, RaDec, cur_pos.ang_pos);
	
	cur_dir.AzElev = AzElev;
	cur_dir.Jd = Jd;
	cur_dir.OwnAxes = AzElev2OwnAxes(Jd, AzElev, cur_pos.ang_pos, R);
	cur_dir.RaDec = RaDec;
}

void MountController::SetCurrentDirectionAzElev(double Jd, const Angs& AzElev) {
	cur_dir.AzElev = AzElev;
	cur_dir.Jd = Jd;
	cur_dir.OwnAxes = AzElev2OwnAxes(Jd, AzElev, cur_pos.ang_pos, R);
	cur_dir.RaDec = convertor->AzElev2RaDec(Jd, AzElev, cur_pos.ang_pos);
}

void MountController::SetCurrentDirectionOwnAxes(double Jd, const Angs& OwnAxes) {
	Angs AzElev = OwnAxes2AzElev(Jd, OwnAxes, cur_pos.ang_pos, InvR);

	cur_dir.AzElev = AzElev;
	cur_dir.Jd = Jd;
	cur_dir.OwnAxes = OwnAxes;
	cur_dir.RaDec = convertor->AzElev2RaDec(Jd, AzElev, cur_pos.ang_pos);
}

void MountController::SetCurrentPositionLatLonElev(double Lat, double Lon, double Elev) {
	cur_pos.ang_pos.height = Elev * 1000;
	cur_pos.ang_pos.latitude = Lat;
	cur_pos.ang_pos.longitude = Lon;

	convertor->WGS84_XYZ(Elev * 1000, Lat, Lon, cur_pos.pos_ITRF[0], cur_pos.pos_ITRF[1], cur_pos.pos_ITRF[2]);
}

void MountController::SetCurrentPositionLatLonElev(double* pos) {
	SetCurrentPositionLatLonElev(pos[0], pos[1], pos[2]);
}

void MountController::SetCurrentPositionITRF(double x, double y, double z) {
	cur_pos.pos_ITRF[0] = x;
	cur_pos.pos_ITRF[1] = y;
	cur_pos.pos_ITRF[2] = z;

	convertor->XYZ_WGS84(x, y, z, cur_pos.ang_pos.height, cur_pos.ang_pos.latitude, cur_pos.ang_pos.longitude);
}

void MountController::SetCurrentPositionITRF(double* pos) {
	SetCurrentPositionITRF(pos[0], pos[1], pos[2]);
}

void MountController::SetCalibrationMatrix(const Eigen::Matrix3d& R_) {
	R = R_;
	InvR = R_.inverse();
}

CurrentPosition MountController::GetCurrentPosition() const {
	return cur_pos;
}

CurrentDirection MountController::GetCurrentDirection() const {
	return cur_dir;
}

Eigen::Matrix3d MountController::GetCalibrationMatrix() const {
	return R;
}

Angs MountController::OwnAxes2AzElev(double Jd, const Angs& OwnAxes, const on_surface& pos, const Eigen::Matrix3d& InvR) const {
	Angs tmpOwnAxes = OwnAxes;
	if (tmpOwnAxes.ang1 > M_PI / 2) {
		tmpOwnAxes.ang1 -= M_PI;
		tmpOwnAxes.ang2 += M_PI / 2;
	}
	else {
		tmpOwnAxes.ang2 = M_PI / 2 - tmpOwnAxes.ang2;
	}

	std::vector<double> xyz = convertor->SphereCoord2XYZ({ 1.0, tmpOwnAxes.ang1, tmpOwnAxes.ang2 });
	Eigen::Vector3d v(xyz[0], xyz[1], xyz[2]);
	Eigen::Vector3d xyzAzElev = InvR * v;

	SphereCoord sphere = convertor->XYZ2SphereCoord({ xyzAzElev[0], xyzAzElev[1], xyzAzElev[2] });
	if (OwnAxes.ang2 < -M_PI / 2 || OwnAxes.ang2 > M_PI / 2) {
		sphere.ang1 += M_PI;
	}
	else if (OwnAxes.ang2 < M_PI / 2 && OwnAxes.ang1 < M_PI / 2) {
		sphere.ang1 += 2 * M_PI;
	}

	return { sphere.ang1, sphere.ang2 };
}

Angs MountController::AzElev2OwnAxes(double Jd, const Angs& AzElev, const on_surface& pos, const Eigen::Matrix3d& R) const {
	std::vector<double> xyz = convertor->SphereCoord2XYZ({ 1.0, AzElev.ang1, AzElev.ang2 });
	Eigen::Vector3d v(xyz[0], xyz[1], xyz[2]);
	Eigen::Vector3d xyz1 = R * v;

	SphereCoord sphere;
	if (xyz1[2] != 1) {
		sphere = convertor->XYZ2SphereCoord({ xyz1[0], xyz1[1], xyz1[2] });
		if (sphere.ang1 < 0) {
			sphere.ang1 += M_PI;
			sphere.ang2 -= M_PI / 2;
		}
		else {
			sphere.ang2 = M_PI / 2 - sphere.ang2;
		}
	}
	else {
		sphere.ang1 = M_PI / 2;
		sphere.ang2 = 0;
	}

	return { sphere.ang1, sphere.ang2 };
}

Angs MountController::OwnAxes2RaDec(double Jd, const Angs& OwnAxes, const on_surface& pos, const Eigen::Matrix3d& InvR) const {
	Angs AzElev = OwnAxes2AzElev(Jd, OwnAxes, pos, InvR);

	return convertor->AzElev2RaDec(Jd, AzElev, pos);
}

Angs MountController::RaDec2OwnAxes(double Jd, const Angs& RaDec, const on_surface& pos, const Eigen::Matrix3d& R) const {
	Angs AzElev = convertor->RaDec2AzElev(Jd, RaDec, pos);
	
	return AzElev2OwnAxes(Jd, AzElev, pos, R);
}