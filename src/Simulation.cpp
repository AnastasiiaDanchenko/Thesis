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

void Initialization() {
    InitBoundaries();
    InitFluid();
    UniformGrid();
}

void Simulation2D() {
	NSUniformGrid2D();
	ComputeDensityPressure2D();
	ComputeAcceleration2D();
	UpdateParticles2D();
}

void Initialization2D() {
	InitBoundaries2D();
	InitFluid2D();
	UniformGrid2D();
}

void SimulationIISPH2D() {
	NSUniformGrid2D();

	ComputeDensity2D();
	ComputeSurface2D();

	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();

	UpdateParticles2D();

	//NSUniformGrid2D();
	//CalculateDensity2D();
	//PredictVelocityAdvection2D();
	//ComputeDii();
	//ComputeAii();
	//PredictDensity2D();
	//SolvePressure2D();
	//ComputePressureAcceleration2D();
	//AdvectParticles2D();
}
