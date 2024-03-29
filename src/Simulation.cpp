#include "..\headers\Simulation.h"

void Simulation() {
    NSUniformGrid();
    //QuadraticSearch();
    ComputeDensityPressure();
    ComputeAcceleration();
    UpdateParticles();
}

void Initialization(const int l) {
    InitBoundaries();
    InitFluid(l);
    UniformGrid();
}
