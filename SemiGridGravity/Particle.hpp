#pragma once

#include "Types.hpp"

class Particle
{
public:
	Vec2d pos = { 0,0 };
	Vec2d vel = { 0,0 };
	Vec2d acc = { 0,0 };
	double mass = 1;
	bool active = true;

	void Update(double timestep);
};