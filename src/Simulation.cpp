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

void Simulation2D(Solver2D& solver) {
	solver.neighborSearch();
	solver.computeDensityPressure();
	solver.computeAcceleration();
	solver.updateParticles();
}

void Initialization2D(Solver2D& solver) {
	solver.initBoundaries();
	solver.initFluid();
}

void MovingBoundaryInitialization() {
	InitMovingThroughBoundaries2D();
	InitFluidForBoundaryTest2D();
}

void RotatingBoundaryInitialization(Solver2D& solver) {
	solver.initBoundaries();
	RotatingBoundary2D();
	solver.initFluid();
}

void RotatingBoundaryIISPH2D(Solver2D& solver) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.computeSurface();
	solver.predictVelocity();
	solver.computeDensityError();
	solver.computeLaplacian();
	solver.compressionConvergence();

	RotateBoundary2D();
}

void SimulationIISPH2D(Solver2D& solver) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.computeSurface();
	solver.predictVelocity();
	solver.computeDensityError();
	solver.computeLaplacian();
	solver.compressionConvergence();
	solver.advectParticles();
}

void MovingBoundaryIISPH2D(Solver2D& solver) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.computeSurface();
	solver.predictVelocity();
	solver.computeDensityError();
	solver.computeLaplacian();
	solver.compressionConvergence();
	solver.advectParticles();
}
