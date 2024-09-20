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
	solver.computeArtificialDensity();
	solver.predictVelocity();
	solver.computeSourceTerm();
	solver.computeDiagonalElement();
	solver.compressionConvergence();
	solver.updateParticles();

	if (parameters.simulationType == 0) {
		solver.neighborSearchGhosts();
		solver.updateGhosts();
	}
	else if (parameters.simulationType == 1) {
		for (auto& body : solver.getRigidBodies()) {
			body.updateBodyQuantities();
			body.updateParticles();
		}
	}
}

void Initialization(Solver& solver) {
	if (parameters.simulationType == 0) {
		solver.initBoundaries();
		solver.initFluid();
		solver.initGhostFluid();
		solver.initGhostBoundary();
	}
	else if (parameters.simulationType == 1) {
		//solver.initBoundaries();
		//solver.initFluid();
		if (parameters.rigidBodyType == "cubes") {
			solver.initRigidCube();

			for (auto& body : solver.getRigidBodies()) {
				for (auto& p : body.getOuterParticles()) {
					p.parentBody = &body;
				}
			}
		}
		else if (parameters.rigidBodyType == "cylinder") {
			solver.initRigidCylinder();
		}
		else if (parameters.rigidBodyType == "cuboid") {
			solver.initRigidCuboid();
		} 
		else if (parameters.rigidBodyType == "cubes_fall") {
			solver.initRigidCubesFalling();

			for (auto& body : solver.getRigidBodies()) {
				for (auto& p : body.getOuterParticles()) {
					p.parentBody = &body;
				}
			}
		}
		else {
			solver.addRigidBody(solver.sampleOBJ());
		}
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
	solver.computeSourceTerm();
	solver.computeDiagonalElement();
	solver.compressionConvergence();
	solver.rotateBoundary();
}

void SimulationIISPH2D(Solver2D& solver) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.computeSurface();
	solver.predictVelocity();
	solver.computeSourceTerm();
	solver.computeDiagonalElement();
	solver.compressionConvergence();
	solver.advectParticles();
}

void MovingBoundaryIISPH2D(Solver2D& solver) {
	solver.neighborSearch();
	solver.boundaryMassUpdate();
	solver.computeDensity();
	solver.computeSurface();
	solver.predictVelocity();
	solver.computeSourceTerm();
	solver.computeDiagonalElement();
	solver.compressionConvergence();
	solver.advectParticles();
}
