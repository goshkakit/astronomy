#include "space_force.h"

namespace Space
{
	class CStaticForceInitializer
	{
	public:
		CStaticForceInitializer();
		~CStaticForceInitializer();
		Force::InfluenceForce * operator()() const;

	private:
		Force::InfluenceForce *f;
	};
}

Space::CStaticForceInitializer::CStaticForceInitializer()
{
	f = new Force::InfluenceForce();
	f->Init_CPU();
}

Space::CStaticForceInitializer::~CStaticForceInitializer()
{
	f->DeInit_CPU();
	delete f;
}

Force::InfluenceForce * Space::CStaticForceInitializer::operator()() const
{
	return f;
}

Space::CStaticForceInitializer StaticForceInitializer;

Force::InfluenceForce *const Space::Force = StaticForceInitializer();
