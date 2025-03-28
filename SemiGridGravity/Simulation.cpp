#include "Simulation.hpp"
#include <cmath>
#include <random>

void Simulation::Step()
{
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> dist(0, 1);
	gridOffset = { dist(engine),dist(engine) };
	UpdateGridSize();
	UpdateCellMasses();
	MakeQuadTree();
	UpdateCellForces();
	UpdateParticles();
}

Vec2i Simulation::GetGridPos(Vec2d pos) const
{
	Vec2d offset = pos - (boundary.c1 - gridOffset);
	double fracX = offset.x * (1 / (boundary.c2.x - boundary.c1.x));
	double fracY = offset.y * (1 / (boundary.c2.y - boundary.c1.y));
	Vec2i ret = { (int)(fracX * gridSize.x),(int)(fracY * gridSize.y) };
	return ret;
}

void Simulation::MakeQuadTree()
{
	qTree = QuadTree({ 0,0 }, gridSize + Vec2i{ 1,1 }, cells);
}

void Simulation::UpdateGridSize()
{
	RectD currentBound = { 0.5 * (boundary.c1 + boundary.c2),0.5 * (boundary.c1 + boundary.c2) };
	for (int i = 0; i < particles.size(); i++)
	{
		if (!particles[i].active)
			continue;
		Vec2d pos = particles[i].pos;

		if (pos.x < currentBound.c1.x)
		{
			currentBound.c1.x = pos.x;
		}
		else if (pos.x > currentBound.c2.x)
		{
			currentBound.c2.x = pos.x;
		}
		if (pos.y < currentBound.c1.y)
		{
			currentBound.c1.y = pos.y;
		}
		else if (pos.y > currentBound.c2.y)
		{
			currentBound.c2.y = pos.y;
		}
	}
	boundary = currentBound;
	int sizeX = 1 + (int)(boundary.c2.x - boundary.c1.x) / (cellSize);
	int sizeY = 1 + (int)(boundary.c2.y - boundary.c1.y) / (cellSize);
	gridSize = { sizeX,sizeY };
	cells = Grid<Cell>(sizeX + 1, sizeY + 1);
}

void Simulation::UpdateCellMasses()
{
	for (int i = 0; i < cells.X * cells.Y; i++)
	{
		cells.arr[i] = Cell();
		cells.arr[i].particles.reserve(particles.size() / double(cells.X * cells.Y));
	}
	for (int i = 0; i < particles.size(); i++)
	{
		Particle& particle = particles[i];
		if (!particle.active)
			continue;
		Vec2i gridPos = GetGridPos(particle.pos);
		if (!cells.CheckBound((double)gridPos.x, (double)gridPos.y))
			continue;
		cells[gridPos].particles.push_back(&particle);
	}
	for (int i = 0; i < cells.X * cells.Y; i++)
	{
		Cell& cell = cells.arr[i];
		cell.center = { 0,0 };
		cell.mass = 0;
		for (int j = 0; j < cell.particles.size(); j++)
		{
			cell.center += cell.particles[j]->mass * cell.particles[j]->pos;
			cell.mass += cell.particles[j]->mass;
		}
		if(cell.particles.size() > 0)
			cell.center *= (1 / cell.mass);
	}
}

void Simulation::UpdateCellForces()
{
	for (int i = 0; i < cells.X * cells.Y; i++)
	{
		//if (cells.arr[i].mass == 0)
		//	continue;
		//for (int j = i + 1; j < cells.X * cells.Y; j++)
		//{
		//	if (cells.arr[j].mass == 0)
		//		continue;
		//	Vec2d axis = cells.arr[i].center - cells.arr[j].center;
		//	double dist = Mag(axis);
		//	axis *= (1.0 / dist);
		//	cells.arr[i].acc -= axis * (cells.arr[j].mass / (dist * dist));
		//	cells.arr[j].acc += axis * (cells.arr[i].mass / (dist * dist));
		//}
		qTree.CalculateForce(*this, cells.arr[i]);
		cells.arr[i].acc *= gravityStrength;
	}
}

void Simulation::UpdateParticles()
{
	for (int i = 0; i < particles.size(); i++)
	{
		Particle& particle = particles[i];
		if (!particle.active)
			continue;
		Vec2i gridPos = GetGridPos(particle.pos);
		if (!cells.CheckBound((double)gridPos.x,(double)gridPos.y)) continue;

		particle.acc += cells[gridPos].acc;
		//if (cells[gridPos].particles.size() > 10)
		//{
		//	Vec2d center = cells[gridPos].center * cells[gridPos].mass;
		//	center -= particle.pos;
		//	double mass = cells[gridPos].mass - particle.mass;
		//	center = center * (1 / mass);
		//	Vec2d axis = center - particle.pos;
		//	double dist = Mag(axis);
		//	axis *= (1.0 / dist);
		//	particle.acc += 0.01 * axis * (mass / (dist * dist + 0.001));
		//}
		particle.Update(0.01);

		if (!maxBoundary.CheckBound(particle.pos))
			particle.active = false;
	}
}

void Simulation::FillArea(int count)
{
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> distX(boundary.c1.x, boundary.c2.x);
	std::uniform_real_distribution<double> distY(boundary.c1.y, boundary.c2.y);
	for (int i = 0; i < count; i++)
	{
		Particle newParticle;
		newParticle.pos = { distX(engine),distY(engine) };
		newParticle.vel = { 0,0 };
		particles.push_back(newParticle);
	}
}

void Simulation::FillCircle(int count, Vec2d center, double radius)
{
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> distX(boundary.c1.x, boundary.c2.x);
	std::uniform_real_distribution<double> distY(boundary.c1.y, boundary.c2.y);
	for (int i = 0; i < count; i++)
	{
		Particle newParticle;
		while(Mag(newParticle.pos - center) > radius)
			newParticle.pos = { distX(engine),distY(engine) };
		newParticle.vel = { 0,0 };
		particles.push_back(newParticle);
	}
}

QuadTree::QuadTree(Vec2i GridPos, Vec2i Size, Grid<Cell>& cells)
{
	gridPos = GridPos;
	size = Size;
	children = {};
	if (size.x > 1 || size.y > 1)
	{
		Subdivide(cells);
		mass = 0;
		center = { 0,0 };
		for (int i = 0; i < children.size(); i++)
		{
			mass += children[i].mass;
			center += children[i].center * children[i].mass;
		}
		if (mass > 0)
			center *= (1 / mass);
		else
			center = { 0,0 };
	}
	else
	{
		mass = cells[gridPos].mass;
		center = cells[gridPos].center;
	}
}

void QuadTree::Subdivide(Grid<Cell>& cells)
{
	if (size.x > 1 && size.y > 1)
	{
		Vec2i halfSize = 0.5 * size;
		children.push_back(QuadTree(gridPos, halfSize, cells));
		children.push_back(QuadTree(gridPos + Vec2i{ halfSize.x,0 }, { size.x - halfSize.x, halfSize.y },cells));
		children.push_back(QuadTree(gridPos + Vec2i{ 0,halfSize.y }, { halfSize.x,size.y - halfSize.y }, cells));
		children.push_back(QuadTree(gridPos + halfSize, size - halfSize, cells));
		return;
	}
	if (size.x > 1 && size.y <= 1)
	{
		Vec2i halfSize = { 0.5 * size.x,size.y };
		children.push_back(QuadTree(gridPos, halfSize, cells));
		children.push_back(QuadTree(gridPos + Vec2i{ halfSize.x,0 }, { size.x - halfSize.x,size.y }, cells));
		return;
	}
	if (size.x <= 1 && size.y > 1)
	{
		Vec2i halfSize = { size.x,0.5 * size.y };
		children.push_back(QuadTree(gridPos, halfSize, cells));
		children.push_back(QuadTree(gridPos + Vec2i{ 0,halfSize.y }, { size.x,size.y - halfSize.y }, cells));
		return;
	}
}

void QuadTree::CalculateForce(Simulation& sim, Cell& cell)
{
	if (mass == 0)
		return;
	Vec2d axis = center - cell.center;
	double dist = Mag(axis);
	if (children.size() == 0 || Mag(size * sim.cellSize) / dist < ratio)
	{
		if (Mag(cell.center - center) < 0.0001)
			return;
		axis *= (1.0 / dist);
		cell.acc += axis * (mass / (dist * dist));
	}
	else
	{
		for (int i = 0; i < children.size(); i++)
		{
			children[i].CalculateForce(sim, cell);
		}
	}
}
