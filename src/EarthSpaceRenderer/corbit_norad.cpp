#include "corbit_norad.h"
#include "Norad/cEci.h"

COrbit_Norad::COrbit_Norad(Zeptomoby::OrbitTools::cOrbit * _pOrbit)
	: pOrbit(_pOrbit)
{

}

struct Scene::vec3d_t COrbit_Norad::GetNewPosition(double time)
{
	Zeptomoby::OrbitTools::cEciTime eci1 = pOrbit->GetPosition(time);

	const Zeptomoby::OrbitTools::cVector & pos = eci1.Position();

	return Scene::vec3d_t(
		  pos.m_x * 1e-6
		, pos.m_y * 1e-6
		, pos.m_z * 1e-6
	);
}
