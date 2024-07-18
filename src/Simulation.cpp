#include "..\headers\Simulation.h"

void Simulation(Grid& grid) {
	grid.updateGrid();
	grid.neighborSearch(particles);
    ComputeDensityPressure();
    ComputeAcceleration();
    UpdateParticles();
}

void SimulationIISPH(Grid& grid) {
	grid.updateGrid();
	grid.neighborSearch(particles);
	BoundaryMassUpdate();
	ComputeDensity();
    PredictVelocity();
    ComputeDensityError();
    ComputeLaplacian();
    CompressionConvergence();
    UpdateParticles();

	grid.neighborSearch(ghostParticles);
	UpdateGhosts();
}

void Initialization() {
    InitBoundaries();
    InitFluid();

	InitGhostFluid();
}

void Simulation2D(Grid2D& grid2D) {
	grid2D.updateGrid();
	grid2D.neighborSearch(particles2D);
	ComputeDensityPressure2D();
	ComputeAcceleration2D();
	Update2D();
}

void Initialization2D() {
	InitBoundaries2D();
	InitFluid2D();
}

void MovingBoundaryInitialization() {
	InitMovingThroughBoundaries2D();
	InitFluidForBoundaryTest2D();
}

void RotatingBoundaryInitialization() {
	InitBoundaries2D();
	RotatingBoundary2D();
	InitFluid2D();
}

void RotatingBoundaryIISPH2D(Grid2D& grid2D) {
	grid2D.updateGrid();
	grid2D.neighborSearch(particles2D);
	BoundaryMassUpdate2D();

	ComputeDensity2D();
	ComputeSurface2D();
	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();
	RotateBoundary2D();
}

void SimulationIISPH2D(Grid2D& grid2D) {
	grid2D.updateGrid();
	grid2D.neighborSearch(particles2D);
	BoundaryMassUpdate2D();

	ComputeDensity2D();
	ComputeSurface2D();
	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();
	UpdateParticles2D();
}

void MovingBoundaryIISPH2D(Grid2D& grid2D) {
	grid2D.updateGrid();
	grid2D.neighborSearch(particles2D);
	BoundaryMassUpdate2D();

	ComputeDensity2D();
	ComputeSurface2D();
	PredictVelocity2D();
	ComputeDensityError2D();
	ComputeLaplacian2D();
	CompressionConvergence2D();
	UpdateParticles2D();
}
