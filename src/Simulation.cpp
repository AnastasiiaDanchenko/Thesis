#include "..\headers\Simulation.h"

void Simulation() {
    NSUniformGrid();
    ComputeDensityPressure();
    ComputeAcceleration();
    UpdateParticles();
}

void SimulationIISPH() {
	NSUniformGrid();
	BoundaryMassUpdate();
	ComputeDensity();
    PredictVelocity();
    ComputeDensityError();
    ComputeLaplacian();
    CompressionConvergence();
    UpdateParticles();

	NSghostsGrid();
	UpdateGhosts();
}

void Initialization() {
    InitBoundaries();
    InitFluid();
    UniformGrid();

	InitGhostFluid();
}

void Simulation2D() {
	NSUniformGrid2D();
	ComputeDensityPressure2D();
	ComputeAcceleration2D();
	Update2D();
}

void Initialization2D() {
	InitBoundaries2D();
	InitFluid2D();
	UniformGrid2D();
}

void MovingBoundaryInitialization() {
	InitMovingThroughBoundaries2D();
	InitFluidForBoundaryTest2D();
	UniformGrid2D();
}

void RotatingBoundaryInitialization() {
	InitBoundaries2D();
	RotatingBoundary2D();
	InitFluid2D();
	UniformGrid2D();
}

void RotatingBoundaryIISPH2D() {
	NSUniformGrid2D();
	BoundaryMassUpdate2D();
	ComputeDensity2D();
	ComputeSurface2D();
	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();
	RotateBoundary2D();
}

void SimulationIISPH2D() {
	NSUniformGrid2D();
	BoundaryMassUpdate2D();

	ComputeDensity2D();
	ComputeSurface2D();
	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();
	UpdateParticles2D();
}

void MovingBoundaryIISPH2D() {
	NSUniformGrid2D();

	BoundaryMassUpdate2D();
	ComputeDensity2D();
	ComputeSurface2D();
	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();
	UpdateParticles2D();
}
