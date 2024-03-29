#pragma once
#include "Particle.h"
#include <list>
#include <numeric>

struct GridCell {
    Eigen::Vector3i cellNumber;
    Eigen::Vector3f minBounds;
    Eigen::Vector3f maxBounds;
    std::vector<Particle*> cellParticles;
};

// Particle container
extern std::vector<Particle> particles;
extern std::vector<GridCell> grid;
extern std::vector<std::list<Particle*>> linearGrid;
extern std::vector<size_t> particleIndices;

extern int GRID_WIDTH;
extern int GRID_HEIGHT;
extern int GRID_DEPTH;

void InitFluid(const int lenght);
void InitBoundaries();

void UniformGrid();
void GridUpdate();
