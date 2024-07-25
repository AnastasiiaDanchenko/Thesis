#include "..\headers\Simulation.h"

void Simulation(Solver& solver) {
	solver.neighborSearch();
	solver.computeDensityPressure();
	solver.computeAcceleration();
	solver.updateParticles();
}

void SimulationIISPH(Solver& solver, int simulationCode) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.predictVelocity();
	solver.computeDensityError();
	solver.computeLaplacian();
	solver.compressionConvergence();
	solver.updateParticles();

	if (simulationCode == 0) {
		solver.neighborSearchGhosts();
		solver.updateGhosts();
	}
	else if (simulationCode == 1) {
		for (auto& body : solver.getRigidBodies()) {
			body.computeParticleQuantities();
			body.updateBodyQuantities();
			body.updateParticles();
		}
	}
}

void Initialization(Solver& solver, int simulationCode) {
	if (simulationCode == 0) {
		solver.initBoundaries();
		solver.initFluid();
		solver.initGhostFluid();
		solver.initGhostBoundary();
	}
	else if (simulationCode == 1) {
		solver.initBoundaries();
		solver.initFluid();
		solver.initRigidCube();
	}
	
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
