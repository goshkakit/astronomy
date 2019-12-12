#pragma once

// System 
#include "InfluenceForce\InfluenceForce.h"
#include "AstroDataManager.h"
#include "AstroTypes.h"

class TrackCalculation
{
private:

	DataConverter Dconv_;
	AstroConstant const_;
	Force::InfluenceForce *IForce;

public:
	TrackCalculation() {}

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	void Init()
	{
		// init rotation
		IForce = new Force::InfluenceForce();
		IForce->Init_CPU();
	}

	TelescopeViewParams ViewTEMESatFromTelescope(double JD_UTC, double* tel_pos, double *sat_pos)
	{
		// MDB
		double jdm = JD_UTC + 0.125;
		double date1 = Dconv_.JDtoYYYYMMDD(jdm);
		double time1 = Dconv_.SECtoHHMMSS(date1, jdm);

		// установка времени
		double int1, ajd1, delt1;
		IForce->set_time(date1, time1, &ajd1, &delt1, &int1);

		// матрица перевода в земную систему
		double Arot[9];
		IForce->iers_update_matrix(int1, Arot, ajd1, delt1);

		// матрица перехода из земной в нормальную систему
		//double invArot[9];
		//IForce->transpose(Arot, invArot);

		// матрица перехода из ICRF в TEME
		double A_Teme[9];
		IForce->GetTemeMatrix(int1, A_Teme, ajd1, delt1);

		// обратная матрица перехода из TEME в ICRF
		double invA_Teme[9];
		IForce->transpose(A_Teme, invA_Teme);


		// перевод вектора состояния в ICRF 
		double PICRF[3];
		IForce->matVecMul(invA_Teme, sat_pos, PICRF);
		// переход в земную
		double PITRF[3];
		IForce->matVecMul(Arot, PICRF, PITRF);

		// перевод координат телескопа в ICRF
		//double TICRF[3];
		//IForce->matVecMul(invArot, tel_pos, TICRF);

		// вектор телескопа
		S3DCoordinate Rtel;
		Rtel.x = tel_pos[0];
		Rtel.y = tel_pos[1];
		Rtel.z = tel_pos[2];

		// вектор наблюдения
		S3DCoordinate Rn;
		Rn.x = PITRF[0];
		Rn.y = PITRF[1];
		Rn.z = PITRF[2];
		S3DCoordinate Rmeg = Rn - Rtel;

		double cosa = (Rtel*Rmeg) / (Rtel.norm()*Rmeg.norm());
		double zenit = acos(cosa);

		TelescopeViewParams vp;

		vp.az = 0;
		vp.el = 90.0 - zenit*const_.RG;
		vp.range = Rmeg.norm();

		return vp;
	}

	void ProcessTrack_v1(TelescopeTrack &telescopeTrack, cOrbit* orbit, cSite siteView)
	{
		double stepSec = telescopeTrack.stepSec_;
		int Nglobal = 700;
		int StepperCount = 200;
		int StepDivider = 16;
		double StepSizeDeg = 360.0 / Nglobal / StepperCount / StepDivider;
		double StepSizeRad = StepSizeDeg * 3600.0;
		printf("stepSec = %f, StepSizeDeg = %f, StepSizeRad = %f\n", stepSec, StepSizeDeg, StepSizeRad);

		double Sum_N_az = 0;
		double Sum_N_el = 0;
		double Sum_Ni_az = 0;
		double Sum_Ni_el = 0;

		double Nd_az = 0;
		double Nd_el = 0;

		double maxErr_Az = 0;
		double maxErr_El = 0;

		double maxErrM_Az = 0;
		double maxErrM_El = 0;

		double a0_az = 0;
		double a0_el = 0;

		int countPoint = telescopeTrack.track.size();
		int idPointStart = 2;
		int idPointEnd = countPoint-2;
		double Position_Az = telescopeTrack.track[idPointStart].az;
		double Position_El = telescopeTrack.track[idPointStart].el;
		double Position_TimeAz = telescopeTrack.track[idPointStart].TimeAfterEpohe;
		double Position_TimeEl = telescopeTrack.track[idPointStart].TimeAfterEpohe;
		// find speed
		if (1)
		{
			TelescopeViewParams pt0 = telescopeTrack.track[idPointStart];
			double h = 0.01;
			cEciTime eci = orbit->GetPosition(pt0.TimeAfterEpohe + h/60.0);
			cTopo topoLook = siteView.GetLookAngle(eci);
			double el = topoLook.ElevationDeg();
			double az = topoLook.AzimuthDeg();

			// Deg/sec
			a0_az = (az - pt0.az) / h;
			a0_el = (el - pt0.el) / h;

			// Step/sec
			a0_az = a0_az / StepSizeDeg;
			a0_el = a0_el / StepSizeDeg;

			// Time step sec
			a0_az = 1.0 / a0_az;
			a0_el = 1.0 / a0_el;
		}
		printf("a0_az = %f a0_el = %f\n", a0_az, a0_el);


		// processing track
		for (int i = idPointStart; i < idPointEnd; i++)
		{
			TelescopeViewParams pt0 = telescopeTrack.track[i];
			TelescopeViewParams pt1 = telescopeTrack.track[i+1];

			// Speed
			double v_az = (pt1.az - pt0.az) / stepSec;
			double v_el = (pt1.el - pt0.el) / stepSec;

			// Delta
			double d_az = (pt1.az - pt0.az);
			double d_el = (pt1.el - pt0.el);

			// Count Step
			double N_az = d_az / StepSizeDeg + Nd_az;
			double N_el = d_el / StepSizeDeg + Nd_el;
			int Ni_az = N_az;
			int Ni_el = N_el;
			Nd_az = N_az - Ni_az;
			Nd_el = N_el - Ni_el;

			// Summary
			Sum_N_az += N_az;
			Sum_N_el += N_el;
			Sum_Ni_az += Ni_az;
			Sum_Ni_el += Ni_el;

			double q_az = (stepSec - abs(a0_az*Ni_az)) * 2.0 / (abs(Ni_az)*(abs(Ni_az) - 1.0));
			double q_el = (stepSec - abs(a0_el*Ni_el)) * 2.0 / (abs(Ni_el)*(abs(Ni_el) - 1.0));
			double a_az = a0_az;
			double a_el = a0_el;

			printf("q_az = %.2f q_el = %.2f ", q_az*1000000, q_el * 1000000);
			printf("a0_az = %.0f a0_el = %.0f ", a0_az * 1000000, a0_el * 1000000);

			// Emulation
			double maxErrMi_Az = 0;
			double maxErrMi_El = 0;
			for (int j = 0; j < abs(Ni_az); j++)
			{
				Position_Az += StepSizeDeg*sgn(Ni_az);
				a_az = a_az + q_az;
				Position_TimeAz += a_az/60.0;
				cEciTime eci = orbit->GetPosition(Position_TimeAz);
				cTopo topoLook = siteView.GetLookAngle(eci);
				double az = topoLook.AzimuthDeg();
				double dt_az = (az - Position_Az)*const_.DtoSec;
				if (maxErrM_Az < dt_az) maxErrM_Az = dt_az;
				if (maxErrMi_Az < dt_az) maxErrMi_Az = dt_az;
			}
			for (int j = 0; j < abs(Ni_el); j++)
			{
				Position_El += StepSizeDeg*sgn(Ni_el);
				a_el = a_el + q_el;
				Position_TimeEl += a_el/60.0;
				cEciTime eci = orbit->GetPosition(Position_TimeEl);
				cTopo topoLook = siteView.GetLookAngle(eci);
				double el = topoLook.ElevationDeg();
				double dt_el = (el - Position_El)*const_.DtoSec;
				if (maxErrM_El < dt_el) maxErrM_El = dt_el;
				if (maxErrMi_El < dt_el) maxErrMi_El = dt_el;
			}
			a0_az = a_az;
			a0_el = a_el;

			double Err_az = abs(Position_Az - pt1.az)*const_.DtoSec;
			double Err_el = abs(Position_El - pt1.el)*const_.DtoSec;
			if (maxErr_Az < Err_az) maxErr_Az = Err_az;
			if (maxErr_El < Err_el) maxErr_El = Err_el;


			printf("AZ: %5.2f %5.2f %5.2f %d [err:%5.4f, %5.4f]\t EL:  %5.2f %5.2f %5.2f %d [err:%5.4f %5.4f]\n", pt0.az, d_az, N_az, Ni_az, Err_az, maxErrMi_Az, pt0.el,  d_el, N_el, Ni_el, Err_el, maxErrMi_El);
		}

		double SumDist_az = telescopeTrack.track[idPointEnd].az - telescopeTrack.track[idPointStart].az;
		double SumDist_el = telescopeTrack.track[idPointEnd].el - telescopeTrack.track[idPointStart].el;

		
		printf("Max Err AZ: %f EL %f\n", maxErr_Az, maxErr_El);
		printf("Max ErrM AZ: %f EL %f\n", maxErrM_Az, maxErrM_El);
		printf("Common: %f %f %f\t %f %f, %f\n", SumDist_az, Sum_N_az*StepSizeDeg, Sum_Ni_az*StepSizeDeg, SumDist_el, Sum_N_el*StepSizeDeg, Sum_Ni_el*StepSizeDeg);
	}

	SCoordinate getSpeedAz(cOrbit* orbit, cSite siteView, double time_ae, double StepSizeDeg, SCoordinate &stepV, SCoordinate &stepT )
	{
		double h = 0.1;
		double h2 = 2.0*h;

		cEciTime eci = orbit->GetPosition(time_ae - h / 60.0);
		cTopo topoLook = siteView.GetLookAngle(eci);
		double el1 = topoLook.ElevationDeg();
		double az1 = topoLook.AzimuthDeg();

		eci = orbit->GetPosition(time_ae + h / 60.0);
		topoLook = siteView.GetLookAngle(eci);
		double el2 = topoLook.ElevationDeg();
		double az2 = topoLook.AzimuthDeg();

		// Deg/sec
		double a0_az = (az2 - az1) / h2;
		double a0_el = (el2 - el1) / h2;

		SCoordinate pt;
		pt.x = a0_az;
		pt.y = a0_el;

		// Step/sec
		a0_az = a0_az / StepSizeDeg;
		a0_el = a0_el / StepSizeDeg;
		stepV.x = a0_az;
		stepV.y = a0_el;

		// Time step sec
		a0_az = 1.0 / a0_az;
		a0_el = 1.0 / a0_el;

		stepT.x = abs(a0_az);
		stepT.y = abs(a0_el);

		return pt;
	}
	void ProcessTrack(TelescopeTrack &telescopeTrack, cOrbit* orbit, cSite siteView)
	{
		double stepSec = telescopeTrack.stepSec_;
		int Nglobal = 700;
		int StepperCount = 200;
		int StepDivider = 16;
		double StepSizeDeg = 360.0 / Nglobal / StepperCount / StepDivider;
		double StepSizeRad = StepSizeDeg * 3600.0;
		printf("stepSec = %f, StepSizeDeg = %f, StepSizeRad = %f\n", stepSec, StepSizeDeg, StepSizeRad);

		double Sum_N_az = 0;
		double Sum_N_el = 0;
		double Sum_Ni_az = 0;
		double Sum_Ni_el = 0;

		double Nd_az = 0;
		double Nd_el = 0;

		double maxErr_Az = 0;
		double maxErr_El = 0;

		double maxErrM_Az = 0;
		double maxErrM_El = 0;

		int countPoint = telescopeTrack.track.size();
		int idPointStart = 2;
		int idPointEnd = countPoint - 2;
		double Position_Az = telescopeTrack.track[idPointStart].az;
		double Position_El = telescopeTrack.track[idPointStart].el;
		double Position_TimeAz = telescopeTrack.track[idPointStart].TimeAfterEpohe;
		double Position_TimeEl = telescopeTrack.track[idPointStart].TimeAfterEpohe;

		// prev speed
		double stepTaz_prev = 0;
		double stepTel_prev = 0;

		// processing track
		for (int i = idPointStart; i < idPointEnd; i++)
		{
			double MCS = 1000000.0;

			// теущая скорость в точке
			SCoordinate stepV0;
			SCoordinate stepT0;
			TelescopeViewParams pt0 = telescopeTrack.track[i];
			double ptt0 = telescopeTrack.track[i].TimeAfterEpohe;
			SCoordinate Vsat0 = getSpeedAz(orbit, siteView, ptt0, StepSizeDeg, stepV0, stepT0);

			//if (i == idPointStart)
			{
				stepTaz_prev = stepT0.x;
				stepTel_prev = stepT0.y;
			}

			// следуюущая тока
			SCoordinate stepV1;
			SCoordinate stepT1;
			TelescopeViewParams pt1 = telescopeTrack.track[i + 1];
			double ptt1 = telescopeTrack.track[i+1].TimeAfterEpohe;
			SCoordinate Vsat1 = getSpeedAz(orbit, siteView, ptt1, StepSizeDeg, stepV1, stepT1);
			
			//printf("%f %f %f %f AZ_T: %.2f %.2f EL_T: %.2f %.2f\n", Vsat0.x, Vsat1.x, Vsat0.y, Vsat1.y, stepT0.x*MCS, stepT1.x*MCS, stepT0.y*MCS, stepT1.y*MCS);

			// Delta
			double d_az = (pt1.az - pt0.az);
			double d_el = (pt1.el - pt0.el);
			
			// Count Step
			double N_az = d_az / StepSizeDeg + Nd_az;
			double N_el = d_el / StepSizeDeg + Nd_el;
			int Ni_az = N_az;
			int Ni_el = N_el;
			Nd_az = N_az - Ni_az;
			Nd_el = N_el - Ni_el;

			// Summary
			Sum_N_az += N_az;
			Sum_N_el += N_el;
			Sum_Ni_az += Ni_az;
			Sum_Ni_el += Ni_el;

			double dstepTaz = 0;
			double dstepTel = 0;

			if (abs(Ni_az) > 1)
			{
				dstepTaz = (stepSec - abs(/*stepT0.x*/stepTaz_prev*Ni_az)) * 2.0 / (abs(Ni_az)*(abs(Ni_az) - 1.0));
			}
			if (abs(Ni_az) == 1)
			{
				dstepTaz = stepT1.x - stepTaz_prev;// stepT0.x;
			}

			if (abs(Ni_el) > 1)
			{
				dstepTel = (stepSec - abs(/*stepT0.y*/stepTel_prev*Ni_el)) * 2.0 / (abs(Ni_el)*(abs(Ni_el) - 1.0));
			}
			if (abs(Ni_el) == 1)
			{
				dstepTel = stepT1.y - stepTel_prev;// stepT0.y;
			}

			bool reverAz = false;
			if (sgn(stepV1.x) != sgn(stepV0.x))
			{
				//printf("Revert motion Az\n");
				reverAz = true;
			}

			bool reverEl = false;
			if (sgn(stepV1.y) != sgn(stepV0.y))
			{
				//printf("Revert motion El\n");
				reverEl = true;
			}

			//double a_az = a0_az;
			//double a_el = a0_el;
			//printf("q_az = %.2f q_el = %.2f ", q_az * 1000000, q_el * 1000000);
			//printf("a0_az = %.0f a0_el = %.0f ", a0_az * 1000000, a0_el * 1000000);

			// Emulation
			double maxErrMi_Az = 0;
			double maxErrMi_El = 0;
			double stepTaz_curr = stepTaz_prev;// stepT0.x;
			double stepTel_curr = stepTel_prev;// stepT0.y;
			for (int j = 0; j < abs(Ni_az); j++)
			{
				if (!reverAz)
				{
					cEciTime eci = orbit->GetPosition(Position_TimeAz);
					cTopo topoLook = siteView.GetLookAngle(eci);
					double az = topoLook.AzimuthDeg();
					double dt_az = abs(az - Position_Az)*const_.DtoSec;
					if (maxErrM_Az < dt_az) maxErrM_Az = dt_az;
					if (maxErrMi_Az < dt_az) maxErrMi_Az = dt_az;
					Position_TimeAz += stepTaz_curr / 60.0;
					stepTaz_curr += dstepTaz;
				}
				Position_Az += StepSizeDeg*sgn(Ni_az);
			}
			if (abs(Ni_az) == 0 || reverAz )
			{
				Position_TimeAz += stepSec / 60.0;
			}

			//printf("AZ: %f %f %f %e %e\t|\t T: %f %f %e AZ: %f %f\n", stepT0.x, stepT1.x, stepTaz_curr, (stepTaz_curr- stepT1.x)*MCS, dstepTaz, Position_TimeAz, ptt1, (Position_TimeAz-ptt1)*60.0, Position_Az, pt1.az);
			stepTaz_prev = stepTaz_curr;

			for (int j = 0; j < abs(Ni_el); j++)
			{
				if (!reverEl)
				{
					cEciTime eci = orbit->GetPosition(Position_TimeEl);
					cTopo topoLook = siteView.GetLookAngle(eci);
					double el = topoLook.ElevationDeg();
					double dt_el = abs(el - Position_El)*const_.DtoSec;
					if (maxErrM_El < dt_el) maxErrM_El = dt_el;
					if (maxErrMi_El < dt_el) maxErrMi_El = dt_el;
					Position_TimeEl += stepTel_curr / 60.0;
					stepTel_curr += dstepTel;
				}

				Position_El += StepSizeDeg*sgn(Ni_el);
			}
			if (abs(Ni_el) == 0 || reverEl )
			{
				Position_TimeEl += stepSec / 60.0;
			}

			//printf("EL: %f %f %f %e %e\t|\t T: %f %f %e EL: %f %f\n", stepT0.y, stepT1.y, stepTel_curr, (stepTel_curr - stepT1.y)*MCS, dstepTel, Position_TimeEl, ptt1, (Position_TimeEl-ptt1)*60.0, Position_El, pt1.el);
			stepTel_prev = stepTel_curr;

			double Err_az = abs(Position_Az - pt1.az)*const_.DtoSec;
			double Err_el = abs(Position_El - pt1.el)*const_.DtoSec;
			if (maxErr_Az < Err_az) maxErr_Az = Err_az;
			if (maxErr_El < Err_el) maxErr_El = Err_el;


			printf("AZ: %5.2f %5.2f %5.2f %d [err:%5.4f, %5.4f]\t EL:  %5.2f %5.2f %5.2f %d [err:%5.4f %5.4f]\n", pt0.az, d_az, N_az, Ni_az, Err_az, maxErrMi_Az, pt0.el, d_el, N_el, Ni_el, Err_el, maxErrMi_El);
			printf("\n");
		}

		double SumDist_az = telescopeTrack.track[idPointEnd].az - telescopeTrack.track[idPointStart].az;
		double SumDist_el = telescopeTrack.track[idPointEnd].el - telescopeTrack.track[idPointStart].el;


		printf("Max Err AZ: %f EL %f\n", maxErr_Az, maxErr_El);
		printf("Max ErrM AZ: %f EL %f\n", maxErrM_Az, maxErrM_El);
		printf("Common: %f %f %f\t %f %f, %f\n", SumDist_az, Sum_N_az*StepSizeDeg, Sum_Ni_az*StepSizeDeg, SumDist_el, Sum_N_el*StepSizeDeg, Sum_Ni_el*StepSizeDeg);
	}

	void PrepareTracks(TelescopObject tel, cOrbit* orbit)
	{
		double SAT_elv = 10.0;
		double jd_delta = 1.0;

		// Telescope point
		// Latitude  in degrees (negative south)
		// Longitude in degrees (negative west)
		// Altitude  in km
		cSite siteView(tel.lat, tel.lon, tel.height);

		double tel_pos[3];
		double tel_pos_wgs72[3];

		double lon = tel.lon*const_.GR;
		double lat = tel.lat*const_.GR;
		double alt = tel.height;

		double h = 6378.135 + tel.height;

		double theta = lon;
		double c = 1.0 / sqrt(1.0 + Zeptomoby::OrbitTools::F * (Zeptomoby::OrbitTools::F - 2.0) * sqr(sin(lat)));
		double s = sqr(1.0 - Zeptomoby::OrbitTools::F) * c;
		double achcp = (XKMPER_WGS72 * c + alt) * cos(lat);

		tel_pos_wgs72[0] = achcp * cos(theta);         // km
		tel_pos_wgs72[1] = achcp * sin(theta);         // km
		tel_pos_wgs72[2] = (XKMPER_WGS72 * s + alt) * sin(lat);   // km

		tel_pos[0] = h*cos(lat)*cos(lon);
		tel_pos[1] = h*cos(lat)*sin(lon);
		tel_pos[2] = h*sin(lat);

		// view position 
		for (int ti = 0; ti < 3; ti++)
		{
			printf("telesope D72 %f m for axes %d\n", (tel_pos_wgs72[ti] - tel_pos[ti])*1000.0, ti);
		}

		double Ap_km = orbit->Apogee();
		std::string idSat = orbit->SatId();
	
		double JDtle = orbit->Epoch().Date();
		double jds = (int)JDtle;
		double jde = (int)JDtle + jd_delta;

		double dayS = Dconv_.JDtoYYYYMMDD(jds);
		double dayE = Dconv_.JDtoYYYYMMDD(jde);
		//printf("JD start: %f ,  Day start: %f\n", jds, dayS);

		bool SatVisible = false;
		double stepSec = 1;
		double stepMin = stepSec / 60.0;

		double mpe_start = (jds - JDtle)*24.0*60.0;
		double mpe_stop = (jde - JDtle)*24.0*60.0;
		double JDstart = 0;

		TelescopeTrack telescopeTrack;
		telescopeTrack.stepSec_ = stepSec;

		// Calculate position, velocity
		// mpe = "minutes past epoch
		double fixAzValue = 0;
		for (double mpe = mpe_start; mpe <= mpe_stop; mpe += stepMin)
		{
			double currenJd = JDtle + mpe / 24.0 / 60.0;

			// Get the position of the satellite at time "mpe"
			cEciTime eci = orbit->GetPosition(mpe);

			// Now get the "look angle" from the site to the satellite. 
			// Note that the ECI object "eciSDP4" contains a time associated
			// with the coordinates it contains; this is the time at which
			// the look angle is valid.
			cTopo topoLook = siteView.GetLookAngle(eci);

			double sat_pos[3];
			sat_pos[0] = eci.Position().m_x;
			sat_pos[1] = eci.Position().m_y;
			sat_pos[2] = eci.Position().m_z;
			TelescopeViewParams vp = ViewTEMESatFromTelescope(currenJd, tel_pos, sat_pos);

			//double current_el = topoLook.ElevationDeg();
			double current_el = vp.el;
			//printf("el %f R %f\n", vp.el, vp.range);

			if (current_el >= SAT_elv && !SatVisible)
			{

				SatVisible = true;
				printf("==========================================================================================\n");
				printf("%s AP %f ", idSat.c_str(), Ap_km);
				printf("EL: %3.5f [%3.5f] R: %5.5f [%5.5f] ", vp.el, topoLook.ElevationDeg(), vp.range, vp.range - topoLook.RangeKm());
				printf("AZ: %3.3f S: %.7f\n", topoLook.AzimuthDeg(), currenJd);
				JDstart = currenJd;

				telescopeTrack.track.clear();
				fixAzValue = 0;
			}
			if (current_el < SAT_elv && SatVisible)
			{
				SatVisible = false;
				printf("==========================================================================================\n");
				double duration = (currenJd - JDstart)*86400.0;
				printf("EL: %3.5f [%3.5f] R: %5.5f [%5.5f] ", vp.el, topoLook.ElevationDeg(), vp.range, vp.range - topoLook.RangeKm());
				printf("AZ: %3.3f E: %.7f D: %.5f\n", topoLook.AzimuthDeg(), currenJd, duration);

				ProcessTrack(telescopeTrack, orbit, siteView);
				//break;
			}

			if (SatVisible)
			{
				double el = topoLook.ElevationDeg();
				double az = topoLook.AzimuthDeg();

				TelescopeViewParams point;
				point.az = az;
				point.el = el;
				point.TimeAfterEpohe = mpe;

				if (telescopeTrack.track.size() > 0)
				{
					int last = telescopeTrack.track.size() - 1;
					TelescopeViewParams pl = telescopeTrack.track[last];

					if (abs(pl.az - point.az) > 180.0)
					{
						fixAzValue = 360.0*sgn(pl.az - point.az);
					}
				}
				if (fixAzValue == 0)
				{
					telescopeTrack.track.push_back(point);
				}

				//printf("%.3f\t %.3f\t %.3f\t %.3f\n", vp.el, vp.az, el, az);
			}
		}
	}

};