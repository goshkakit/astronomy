#include "space_force.h"

namespace Space
{
	class CStaticForceInitializer
	{
	public:
		~CStaticForceInitializer();
		Force::InfluenceForce * operator()() const;

	private:
		static Force::InfluenceForce *f;
	};
}

Force::InfluenceForce * Space::CStaticForceInitializer::f = (Force::InfluenceForce *)0;

Space::CStaticForceInitializer::~CStaticForceInitializer()
{
	if (f) {
		f->DeInit_CPU();
		delete f;
		f = (Force::InfluenceForce *)0;
	}
}

Force::InfluenceForce * Space::CStaticForceInitializer::operator()() const
{
	if (!f) {
		f = new Force::InfluenceForce();
		f->Init_CPU();
	}
	return f;
}

Space::CStaticForceInitializer StaticForceInitializer;

Force::InfluenceForce * Space::Force()
{
	return StaticForceInitializer();
}
