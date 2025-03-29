#include "Simulation.hpp"
//#define RAYMATH_DISABLE_CPP_OPERATORS
#include <raylib-cpp.hpp>

void main()
{
	printf(" > Entering main\n");
	raylib::Window window(800, 800, "Simulation");
	Simulation sim;
	sim.boundary = { {0,0},{100,100} };
	sim.maxBoundary = { { -50,-50 }, { 150,150 } };
	//sim.maxBoundary = sim.boundary;
	sim.cellSize = 1;
	sim.gridSize = { 100,100 };
	sim.cells = Grid<Cell>(101, 101);
	sim.rand.seed(0);

	//sim.FillArea(1000000);
	sim.FillCircle(1000000, { 50,50 }, 15);
	sim.gravityStrength = 0.5;

	for (int i = 0; i < sim.particles.size(); i++)
	{
		Particle& particle = sim.particles[i];
		Vec2d axis = { particle.pos.x - 50,particle.pos.y - 50 };
		double dist = Mag(axis);
		axis = Norm(axis);
		double mult = 4 * sqrt(dist) * dist;
		particle.vel = 0.01 * Vec2d{ -mult * axis.y,mult * axis.x };
	}

	Grid<double> kernel(3, 3);
	kernel[{0, 0}] = 1.0 / 16.0;
	kernel[{0, 1}] = 2.0 / 16.0;
	kernel[{0, 2}] = 1.0 / 16.0;
	kernel[{1, 0}] = 2.0 / 16.0;
	kernel[{1, 1}] = 4.0 / 16.0;
	kernel[{1, 2}] = 2.0 / 16.0;
	kernel[{2, 0}] = 1.0 / 16.0;
	kernel[{2, 1}] = 2.0 / 16.0;
	kernel[{2, 2}] = 1.0 / 16.0;

	int frame = 0;
	while (!window.ShouldClose())
	{
		printf("%i \n", frame++);
		sim.Step();
		BeginDrawing();
		window.ClearBackground(BLACK);

		for(int x = 0; x < 800; x++)
			for (int y = 0; y < 800; y++)
			{
				double fracX = x / 800.0;
				double fracY = y / 800.0;
				Vec2i gridPos = { sim.cells.X * fracX - sim.gridOffset.x,sim.cells.Y * fracY - sim.gridOffset.y };
				raylib::Vector2 screenPos(x, y);
				raylib::Color color = raylib::Color::Red();
				color.a = 255.0 * std::min(1.0, sim.cells[gridPos].mass / 250.0);
				screenPos.DrawPixel(color);
			}

		if (raylib::Mouse::IsButtonPressed(MOUSE_BUTTON_LEFT))
		{
			Vector2 mPosRL = raylib::Mouse::GetPosition();
			Vec2d mPos = { mPosRL.x / 8.0,mPosRL.y / 8.0 };
			for (int i = 0; i < sim.particles.size(); i++)
			{
				Particle& particle = sim.particles[i];
				Vec2d axis = (particle.pos - mPos);
				double dist = Mag(axis);
				axis = Norm(axis);
				double force = 10000.0 / (dist * dist);
				Vec2d push = force * axis;
				particle.acc += push;
			}
		}

		Grid<double> screen(800, 800);
		for (int i = 0; i < sim.particles.size(); i++)
		{
			if (!sim.particles[i].active)
				continue;
			if (!screen.CheckBound({ (int)(sim.particles[i].pos.x * 8), (int)(sim.particles[i].pos.y * 8) }))
				continue;
			screen[{(int)(sim.particles[i].pos.x * 8), (int)(sim.particles[i].pos.y * 8)}] += sim.particles[i].mass;
		}
		for(int i = 0; i < 0; i++)
			screen = Convolution(screen, kernel);
		for(int x = 0; x < 800; x++)
			for (int y = 0; y < 800; y++)
			{
				raylib::Vector2 screenPos(x, y);
				double frac = screen[{x, y}] / 75;
				//frac = frac * frac;
				frac = sqrt(frac);
				if (frac > 1) frac = 1;
				Color color1 = { 0,0,0,255 };
				Color color2 = { 192,64,0,255 };
				Color color3 = { 255,255,0,255 };
				Color color = BLACK;
				if (frac < 0.5)
				{
					double lFrac = 2 * frac;
					color.r = (1 - lFrac) * color1.r + lFrac * color2.r;
					color.g = (1 - lFrac) * color1.g + lFrac * color2.g;
					color.b = (1 - lFrac) * color1.b + lFrac * color2.b;
				}
				else
				{
					double lFrac = 2 * (frac-0.5);
					color.r = (1 - lFrac) * color2.r + lFrac * color3.r;
					color.g = (1 - lFrac) * color2.g + lFrac * color3.g;
					color.b = (1 - lFrac) * color2.b + lFrac * color3.b;
				}
				screenPos.DrawPixel({ color.r,color.g,color.b,unsigned char(255/* * frac*/) });
			}

		EndDrawing();
	}
}
