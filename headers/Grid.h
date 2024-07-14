#pragma once
#include "Particle.h"
#include <list>
#include <numeric>

struct GridCell {
    Eigen::Vector3i cellNumber;
    Eigen::Vector3d minBounds;
    Eigen::Vector3d maxBounds;
    std::vector<Particle*> cellParticles;
};

struct GridCell2D {
	Eigen::Vector2i cellNumber;
	Eigen::Vector2d minBounds;
	Eigen::Vector2d maxBounds;
	std::vector<Particle2D*> cellParticles;
};

// Particle container
extern std::vector<Particle> particles;
extern std::vector<Particle2D> particles2D;
extern std::vector<Particle> ghostParticles;

extern std::vector<GridCell> grid;
extern std::vector<GridCell2D> grid2D;
extern std::vector<std::list<Particle*>> linearGrid;
extern std::vector<size_t> particleIndices;

extern int GRID_WIDTH;
extern int GRID_HEIGHT;
extern int GRID_DEPTH;

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

void UniformGrid();
void UniformGrid2D();
void GridUpdate();
void GridUpdate2D();
