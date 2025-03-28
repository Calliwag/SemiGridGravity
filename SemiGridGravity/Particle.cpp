#include "Particle.hpp"

void Particle::Update(double timestep)
{
	if (!active)
		return;
	vel += acc * (timestep * timestep);
	pos += vel;
	acc = { 0,0 };
}
