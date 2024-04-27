#include "..\headers\Simulation.h"

void Simulation() {
    NSUniformGrid();
    ComputeDensityPressure();
    ComputeAcceleration();
    UpdateParticles();
}

void SimulationIISPH() {
	NSUniformGrid();
	ComputeDensity();
    PredictVelocity();
    ComputeDensityError();
    ComputeLaplacian();
    CompressionConvergence();
    UpdateParticles();
}

void Initialization(const int l) {
    InitBoundaries();
    InitFluid(l);
    UniformGrid();
}
