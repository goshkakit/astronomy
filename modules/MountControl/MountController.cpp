#include "MountController.h"

MountController::MountController(std::unique_ptr<NewConvertor> convertor_, MountSpec mount_spec_, double* pos_ITRF_) :
	convertor(std::move(convertor_))
{
	SetMountSpec(mount_spec_);
	SetCurrentPositionITRF(pos_ITRF_);
}

MountController::MountController(std::unique_ptr<NewConvertor> convertor_, MountSpec mount_spec_, double Lat, double Lon, double Elev) :
	convertor(std::move(convertor_))
{
	double pos_LatLonElev[3] = { Lat, Lon, Elev };
	SetMountSpec(mount_spec_);
	SetCurrentPositionLatLonElev(pos_LatLonElev);
}

MountController::~MountController() {
}

void MountController::SetCurrentDirectionRaDec(double Jd, Angs RaDec) {
	Angs AzElev = convertor->RaDec2AzElev(Jd, RaDec, cur_pos.ang_pos);
	
	cur_dir.AzElev = AzElev;
	cur_dir.Jd = Jd;
	cur_dir.OwnAxes = AzElev2OwnAxes(Jd, AzElev, cur_pos.ang_pos, R);
	cur_dir.RaDec = RaDec;
}

void MountController::SetCurrentDirectionAzElev(double Jd, Angs AzElev) {
	cur_dir.AzElev = AzElev;
	cur_dir.Jd = Jd;
	cur_dir.OwnAxes = AzElev2OwnAxes(Jd, AzElev, cur_pos.ang_pos, R);
	cur_dir.RaDec = convertor->AzElev2RaDec(Jd, AzElev, cur_pos.ang_pos);
}

void MountController::SetCurentDirectionOwnAxes(double Jd, Angs OwnAxes) {
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

void MountController::SetMountSpec(double a_max = 2000.0, double V_max = 10000.0, double gear_ratio = 720.0, double divider = 16.0, double steps_in_circle_motor = 200.0) {
	mount_spec = {
		a_max,
		V_max,
		gear_ratio,
		divider,
		steps_in_circle_motor
	};
}

void MountController::SetMountSpec(MountSpec mount_spec_) {
	mount_spec = mount_spec_;
}

void MountController::SetCalibrationMatrix(Matrix R_) {
	R = R_;
	InvR = R.Inversion();
}

CurrentPosition MountController::GetCurrentPosition() {
	return cur_pos;
}

CurrentDirection MountController::GetCurrentDirection() {
	return cur_dir;
}

MountSpec MountController::GetMountSpec() {
	return mount_spec;
}

Matrix MountController::GetCalibrationMatrix() {
	return R;
}

Angs MountController::OwnAxes2AzElev(double Jd, Angs OwnAxes, on_surface pos, Matrix InvR) {
	Angs tmpOwnAxes = OwnAxes;
	if (tmpOwnAxes.ang1 > M_PI / 2) {
		tmpOwnAxes.ang1 -= M_PI;
		tmpOwnAxes.ang2 += M_PI / 2;
	}
	else {
		tmpOwnAxes.ang2 = M_PI / 2 - tmpOwnAxes.ang2;
	}

	std::vector<double> xyz = convertor->SphereCoord2XYZ({ 1.0, tmpOwnAxes.ang1, tmpOwnAxes.ang2 });
	std::vector<double> xyzAzElev = Mat3x3XStolb3x1(InvR, xyz);

	SphereCoord sphere = convertor->XYZ2SphereCoord(xyzAzElev);
	if (OwnAxes.ang2 < -M_PI / 2 || OwnAxes.ang2 > M_PI / 2) {
		sphere.ang1 += M_PI;
	}
	else if (OwnAxes.ang2 < M_PI / 2 && OwnAxes.ang1 < M_PI / 2) {
		sphere.ang1 += 2 * M_PI;
	}

	return { sphere.ang1, sphere.ang2 };
}

Angs MountController::AzElev2OwnAxes(double Jd, Angs AzElev, on_surface pos, Matrix R) {
	std::vector<double> xyz = convertor->AzElevR2ITRF(AzElev.ang1, AzElev.ang2, 1.0, pos);
	std::vector<double> xyz1 = Mat3x3XStolb3x1(R, xyz);

	SphereCoord sphere;
	if (xyz1[2] != 1) {
		sphere = convertor->XYZ2SphereCoord(xyz1);
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

Angs MountController::OwnAxes2RaDec(double Jd, Angs OwnAxes, on_surface pos, Matrix InvR) {
	Angs AzElev = OwnAxes2AzElev(Jd, OwnAxes, pos, InvR);

	return convertor->AzElev2RaDec(Jd, AzElev, pos);
}

Angs MountController::RaDec2OwnAxes(double Jd, Angs RaDec, on_surface pos, Matrix R) {
	Angs AzElev = convertor->RaDec2AzElev(Jd, RaDec, pos);
	
	return AzElev2OwnAxes(Jd, AzElev, pos, R);
}