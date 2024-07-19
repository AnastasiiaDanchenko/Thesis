#include "..\headers\Simulation.h"

void Simulation(Solver& solver) {
	solver.neighborSearch();
	solver.computeDensityPressure();
	solver.computeAcceleration();
	solver.updateParticles();
}

void SimulationIISPH(Solver& solver) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.predictVelocity();
	solver.computeDensityError();
	solver.computeLaplacian();
	solver.compressionConvergence();
	solver.updateParticles();

	solver.neighborSearchGhosts();
	solver.updateGhosts();
}

void Initialization(Solver& solver) {
	solver.initBoundaries();
	solver.initFluid();
	solver.initGhostFluid();
	solver.initGhostBoundary();
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

void MovingBoundaryInitialization(Solver2D& solver) {
	solver.initMovingBoundary();
	solver.initMovingFluid();
}

void RotatingBoundaryInitialization(Solver2D& solver) {
	solver.initBoundaries();
	solver.initRotatingBoundary();
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
	solver.rotateBoundary();
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
