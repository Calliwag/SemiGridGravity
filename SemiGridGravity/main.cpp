#include "Simulation.hpp"
//#define RAYMATH_DISABLE_CPP_OPERATORS
#include <raylib-cpp.hpp>

struct ScreenData
{
	double mass = 0;
	Vec2i pos = { 0,0 };
};

int main()
{
	SetTraceLogLevel(LOG_FATAL);

	printf(" > Entering main\n");
	raylib::Window window(800, 800, "Simulation");
	Simulation sim;
	sim.boundary = { {0,0},{100,100} };
	sim.maxBoundary = { { -15,-15 }, { 115,115 } };
	//sim.maxBoundary = sim.boundary;
	sim.cellSize = 1;
	sim.gridSize = { 100,100 };
	sim.cells = Grid<Cell>(101, 101);
	sim.rand.seed(10);

	//sim.FillArea(1000000);
	sim.FillCircle(1000000, { 50,50 }, 15);
	sim.gravityStrength = 0.5;

	for (int i = 0; i < sim.particles.size(); i++)
	{
		Particle& particle = sim.particles[i];
		Vec2d axis = { particle.pos.x - 50,particle.pos.y - 50 };
		double dist = Mag(axis);
		axis = Norm(axis);
		double mult = 4.5 * sqrt(dist) * dist;
		particle.vel = 0.01 * Vec2d{ -mult * axis.y,mult * axis.x };
	}

	float kernel[9] = {
		1.0 / 16.0,
		2.0 / 16.0,
		1.0 / 16.0,
		2.0 / 16.0,
		4.0 / 16.0,
		2.0 / 16.0,
		1.0 / 16.0,
		2.0 / 16.0,
		1.0 / 16.0
	};

	int frame = 0;
	while (!window.ShouldClose())
	{
		int count = 0;
		for (int i = 0; i < sim.particles.size(); i++)
		{
			if (sim.particles[i].active)
				count++;
		}

		printf("%i: %i \n", frame++, count);
		sim.Step();
		BeginDrawing();
		window.ClearBackground(BLACK);

		Grid<ScreenData> screen(800, 800);
		for (int x = 0; x < screen.X; x++)
		{
			for (int y = 0; y < screen.Y; y++)
			{
				screen[{x, y}].mass = 0;
				screen[{x, y}].pos = { x,y };
			}
		}
		for (int i = 0; i < sim.particles.size(); i++)
		{
			if (!sim.particles[i].active)
				continue;
			if (!screen.CheckBound({ (int)(sim.particles[i].pos.x * 8), (int)(sim.particles[i].pos.y * 8) }))
				continue;
			screen[{(int)(sim.particles[i].pos.x * 8), (int)(sim.particles[i].pos.y * 8)}].mass += sim.particles[i].mass;
		}
		mutex drawMutex;
		raylib::Image screenImg(800, 800);
		std::for_each(std::execution::par_unseq, screen.arr, screen.arr + screen.X * screen.Y, [&](const ScreenData& data) {
			int x = data.pos.x;
			int y = data.pos.y;
			double frac = data.mass / 75;
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
				color.r = (unsigned char)((1 - lFrac) * color1.r + lFrac * color2.r);
				color.g = (unsigned char)((1 - lFrac) * color1.g + lFrac * color2.g);
				color.b = (unsigned char)((1 - lFrac) * color1.b + lFrac * color2.b);
			}
			else
			{
				double lFrac = 2 * (frac - 0.5);
				color.r = (unsigned char)((1 - lFrac) * color2.r + lFrac * color3.r);
				color.g = (unsigned char)((1 - lFrac) * color2.g + lFrac * color3.g);
				color.b = (unsigned char)((1 - lFrac) * color2.b + lFrac * color3.b);
			}
			color.a = 255;
			raylib::Vector2 screenPos((float)x, (float)y);
			{
				std::lock_guard<mutex> lock(drawMutex);
				screenImg.DrawPixel(x, y, color);
			}
			});
		screenImg.KernelConvolution(kernel, 9);
		raylib::Texture screenTex(screenImg);
		screenTex.Draw();

		EndDrawing();
	}

	return 1;
}
