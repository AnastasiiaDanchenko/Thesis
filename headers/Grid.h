#pragma once
#include "Particle.h"
#include <list>
#include <numeric>

struct GridCell2D {
	Eigen::Vector2i cellNumber;
	std::vector<Particle2D*> cellParticles;
};

struct GridCell {
	Eigen::Vector3i cellNumber;
	std::vector<Particle*> cellParticles;
};

class Grid2D {
protected:
	int gridWidth;
	int gridHeight;
	double cellSize;
	std::vector<GridCell2D> cells2D;

public:
	Grid2D(double size);

	void virtual initializeGrid();
	void virtual updateGrid();
	void neighborSearch(std::vector<Particle2D>& particles);
};

class Grid : public Grid2D {
protected:
	int gridDepth;
	std::vector<GridCell> cells;

public:
	Grid(double size);

	void initializeGrid() override;
	void updateGrid() override;
	void neighborSearch(std::vector<Particle>& particles);
};

void InitFluid();
void InitGhostFluid();

void InitBoundaries();
void MovingBoundary();

void InitFluid2D();
void InitFluidForBoundaryTest2D();
void InitFluidForRotatingTest2D();

void InitBoundaries2D();
void InitMovingThroughBoundaries2D();
void MovingBoundary2D();
void RotatingBoundary2D();
