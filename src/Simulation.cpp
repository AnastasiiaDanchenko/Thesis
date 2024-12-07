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

	if (parameters.simulationType == 1) {
		for (auto& body : solver.getRigidBodies()) {
			if (!body.getBoundary()) {
				body.updateBodyQuantities();
				body.updateParticles();
			}
		}
	}
}

void Initialization(Solver& solver) {
	solver.initBoundaries();
	solver.initFluid();

	if (parameters.simulationType == 1) {
		if (parameters.rigidBodyType == "cubes") {
			solver.initRigidCube();
		}
		else if (parameters.rigidBodyType == "cubes_bounce") {
			solver.initRigidCubesBounce();
		}
		else if (parameters.rigidBodyType == "cylinder") {
			solver.initRigidCylinder();
		}
		else if (parameters.rigidBodyType == "cuboid") {
			solver.initRigidCuboid();
		} 
		else if (parameters.rigidBodyType == "cubes_fall") {
			solver.initRigidCubesFalling();
		}
		else if (parameters.rigidBodyType == "bunny") {
			solver.addRigidBody(solver.sampleOBJ());
		}
		else if (parameters.rigidBodyType == "duck") {
			for (int i = 1; i < 5; i++) {
				solver.addRigidBody(solver.sampleOBJ());
				// move the duck to the right
				for (auto& p : solver.getRigidBodies()[i].getOuterParticles()) {
					p.relativePosition += Eigen::Vector3d(0, 20, 10);
					p.position.x() += 150 * i - 400;
				}
				solver.getRigidBodies()[i].setCenteOfMass(solver.getRigidBodies()[i].getPositionCM() + Eigen::Vector3d(150 * i - 400, -20, -20));
			}
		}
		else if (parameters.rigidBodyType == "final") {
			for (int i = -1; i < 2; i++) {
				solver.addRigidBody(solver.sampleOBJ());
				// move the armadillo to the right
				for (auto& p : solver.getRigidBodies()[i + 2].getOuterParticles()) {
					p.position.x() += 150 * i;
				}
				solver.getRigidBodies()[i + 2].setCenteOfMass(solver.getRigidBodies()[i + 2].getPositionCM() + Eigen::Vector3d(150 * i, 0, 0));
			}
		}

		for (auto& body : solver.getRigidBodies()) {
			for (auto& p : body.getOuterParticles()) {
				p.parentBody = &body;
			}
			std::cout << "Rigid body with " << body.getOuterParticles().size() << " particles initialized." << std::endl;
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
