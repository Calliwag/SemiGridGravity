#pragma once

#include "Particle.hpp"
#include <vector>
#include <random>

using std::vector;

class Simulation;

struct Cell
{
	double mass = 0;
	Vec2d center = { 0,0 };
	Vec2d acc = { 0,0 };
	int count = 0;
};

class QuadTree
{
public:
	bool active = true;
	Vec2i gridPos;
	Vec2i size;
	vector<QuadTree> children = {};
	double mass = 0;
	Vec2d center = { 0,0 };
	double ratio = 1;

	QuadTree()
	{
		gridPos = { 0,0 };
		size = { 0,0 };
		children = {};
		mass = 0;
		center = { 0,0 };
	}
	QuadTree(Vec2i GridPos, Vec2i Size, Grid<Cell>& cells);
	void Subdivide(Grid<Cell>& cells);
	void Prune();
	void CalculateForce(Simulation& sim, Cell& cell) const;
};

class Simulation
{
public:
	vector<Particle> particles = {};
	RectD boundary;
	RectD maxBoundary;

	std::mt19937 rand;
	Vec2i gridSize;
	Grid<Cell> cells;
	Vec2d gridOffset = { 0,0 };
	double cellSize;
	double gravityStrength;

	QuadTree* qTree = nullptr;

	void Step();

	Vec2i GetGridPos(Vec2d pos) const;
	void MakeQuadTree();
	void UpdateGridSize();
	void UpdateCellMasses();
	void UpdateCellForces();
	void UpdateParticles();
	void FillArea(int count);
	void FillCircle(int count, Vec2d center, double radius);
};