#pragma once

#include "Particle.hpp"
#include <vector>
#include <random>
#include <mutex>

using std::vector;
using std::mutex;

class Simulation;

struct Cell
{
	double mass;
	Vec2d center;
	Vec2d acc;
	int count;
	mutex Mutex;

	Cell& operator =(const Cell& cell)
	{
		mass = cell.mass;
		center = cell.center;
		acc = cell.acc;
		count = cell.count;
		return *this;
	}
	Cell()
	{
		mass = 0;
		center = { 0,0 };
		acc = { 0,0 };
		count = 0;
	}
	Cell(const Cell& other)
	{
		mass = other.mass;
		center = other.center;
		acc = other.acc;
		count = other.count;
	}
};

struct QTreeSpec
{
	int index;
	Vec2i gridPos;
	Vec2i size;
};

struct NoInit {};
const NoInit noInit;

class QuadTree
{
public:
	bool active = true;
	Vec2i gridPos;
	Vec2i size;
	vector<QuadTree> children;
	double mass;
	Vec2d center;
	double ratio;

	QuadTree()
	{
		children = {};
		gridPos = { 0,0 };
		size = { 0,0 };
		mass = 0;
		center = { 0,0 };
		ratio = 1;
	}
	QuadTree(NoInit) {};
	QuadTree(Vec2i GridPos, Vec2i Size, const Grid<Cell>& cells);
	void Subdivide(const Grid<Cell>& cells);
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