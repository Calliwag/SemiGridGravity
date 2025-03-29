#include "Simulation.hpp"
#include <cmath>
#include <execution>

void Simulation::Step()
{
	std::uniform_real_distribution<double> dist(0, 1);
	gridOffset = { dist(rand),dist(rand) };
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
	if (qTree)
		delete qTree;
	qTree = new QuadTree({ 0,0 }, gridSize + Vec2i{ 1,1 }, cells);
	qTree->Prune();
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
	boundary.c1 -= {1, 1};
	boundary.c2 += {1, 1};
	int sizeX = 1 + (int)(boundary.c2.x - boundary.c1.x) / (cellSize);
	int sizeY = 1 + (int)(boundary.c2.y - boundary.c1.y) / (cellSize);
	gridSize = { sizeX,sizeY };
	cells = Grid<Cell>(sizeX + 1, sizeY + 1);
}

void Simulation::UpdateCellMasses()
{
	for (int i = 0; i < particles.size(); i++)
	{
		Particle& particle = particles[i];
		if (!particle.active)
			continue;
		Vec2i gridPos = GetGridPos(particle.pos);
		//if (!cells.CheckBound((double)gridPos.x, (double)gridPos.y))
		//	continue;
		cells[gridPos].count++;
		cells[gridPos].mass += particle.mass;
		cells[gridPos].center += particle.pos;
	}
	for (int i = 0; i < cells.X * cells.Y; i++)
	{
		if(cells.arr[i].count > 0)
			cells.arr[i].center *= (1 / cells.arr[i].mass);
	}
}

void Simulation::UpdateCellForces()
{
	//for (int i = 0; i < cells.X * cells.Y; i++)
	//{
	//	if (cells.arr[i].count == 0)
	//		continue;
	//	qTree->CalculateForce(*this, cells.arr[i]);
	//	cells.arr[i].acc *= gravityStrength;
	//}
	std::for_each(std::execution::par, cells.arr, cells.arr + cells.X * cells.Y, [&](auto&& cell)
		{
			if (cell.count == 0)
				return;
			qTree->CalculateForce(*this, cell);
			cell.acc *= gravityStrength;
		});
}

void Simulation::UpdateParticles()
{
	for (int i = 0; i < particles.size(); i++)
	{
		Particle& particle = particles[i];
		if (!particle.active)
			continue;

		Vec2d offset = particle.pos - (boundary.c1 - gridOffset);
		double fracX = offset.x * (1 / (boundary.c2.x - boundary.c1.x));
		double fracY = offset.y * (1 / (boundary.c2.y - boundary.c1.y));
		Vec2i gridPos = { (int)(fracX * gridSize.x),(int)(fracY * gridSize.y) };

		if (!cells.CheckBound((double)gridPos.x, (double)gridPos.y)) continue;

		particle.acc += cells[gridPos].acc;
		particle.Update(0.01);

		if (!maxBoundary.CheckBound(particle.pos))
			particle.active = false;
	}
}

void Simulation::FillArea(int count)
{
	std::uniform_real_distribution<double> distX(boundary.c1.x, boundary.c2.x);
	std::uniform_real_distribution<double> distY(boundary.c1.y, boundary.c2.y);
	for (int i = 0; i < count; i++)
	{
		Particle newParticle;
		newParticle.pos = { distX(rand),distY(rand) };
		newParticle.vel = { 0,0 };
		particles.push_back(newParticle);
	}
}

void Simulation::FillCircle(int count, Vec2d center, double radius)
{
	std::uniform_real_distribution<double> distX(boundary.c1.x, boundary.c2.x);
	std::uniform_real_distribution<double> distY(boundary.c1.y, boundary.c2.y);
	for (int i = 0; i < count; i++)
	{
		Particle newParticle;
		while(Mag(newParticle.pos - center) > radius)
			newParticle.pos = { distX(rand),distY(rand) };
		newParticle.vel = { 0,0 };
		particles.push_back(newParticle);
	}
}

QuadTree::QuadTree(Vec2i GridPos, Vec2i Size, Grid<Cell>& cells) : gridPos(GridPos), size(Size)
{
	if (size.x > 1 || size.y > 1)
	{
		Subdivide(cells);
		for (int i = 0; i < children.size(); i++)
		{
			mass += children[i].mass;
			center += children[i].center * children[i].mass;
		}
		if (mass > 0)
			center *= (1 / mass);
		else
		{
			center = { 0,0 };
			active = false;
		}
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
		children.resize(4);
		Vec2i halfSize = 0.5 * size;
		children[0] = QuadTree(gridPos, halfSize, cells);
		children[1] = QuadTree(gridPos + Vec2i{ halfSize.x,0 }, { size.x - halfSize.x, halfSize.y }, cells);
		children[2] = QuadTree(gridPos + Vec2i{ 0,halfSize.y }, { halfSize.x,size.y - halfSize.y }, cells);
		children[3] = QuadTree(gridPos + halfSize, size - halfSize, cells);
		return;
	}
	if (size.x > 1 && size.y <= 1)
	{
		children.resize(2);
		Vec2i halfSize = { 0.5 * size.x,size.y };
		children[0] = QuadTree(gridPos, halfSize, cells);
		children[1] = QuadTree(gridPos + Vec2i{ halfSize.x,0 }, { size.x - halfSize.x,size.y }, cells);
		return;
	}
	if (size.x <= 1 && size.y > 1)
	{
		children.resize(2);
		Vec2i halfSize = { size.x,0.5 * size.y };
		children[0] = QuadTree(gridPos, halfSize, cells);
		children[1] = QuadTree(gridPos + Vec2i{0,halfSize.y}, {size.x,size.y - halfSize.y}, cells);
		return;
	}
}

void QuadTree::Prune()
{
	//if (childCount == 1)
	//{
	//	*this = children[0];
	//}

	//for (int i = 0; i < children.size(); i++)
	//{
	//	if (!children[i].active)
	//	{
	//		children.erase(children.begin() + i);
	//		i--;
	//	}
	//}

	//for (int i = 0; i < children.size(); i++)
	//{
	//	children[i].Prune();
	//}
}

void QuadTree::CalculateForce(Simulation& sim, Cell& cell) const
{
	Vec2d axis = center - cell.center;
	double dist = Mag(axis);
	if (children.size() == 0 || Mag(size * sim.cellSize) / dist < ratio)
	{
		if (dist < 0.0001)
			return;
		axis *= (1.0 / dist);
		cell.acc += axis * (mass / (dist * dist));
	}
	else
	{
		for (int i = 0; i < children.size(); i++)
		{
			if(children[i].active)
				children[i].CalculateForce(sim, cell);
		}
	}
}
