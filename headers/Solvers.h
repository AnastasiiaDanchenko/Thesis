#pragma once
#include "Kernel.h"
#include "..\utilities\OBJ_Loader.h"
#include <omp.h>
#include <limits>

class Solver2D {
private:
	Grid2D grid2D;

public:
	Solver2D();

	void virtual initBoundaries();
	void virtual initFluid();
	void virtual neighborSearch();

	// SPH functions
	void virtual computeDensityPressure();
	void virtual computeAcceleration();
	void virtual updateParticles();

	// IISPH functions
	void virtual boundaryMassUpdate();
	void virtual computeDensity();
	void virtual computeSurface();
	void virtual predictVelocity();
	void virtual computeDensityError();
	void virtual computeLaplacian();
	void virtual compressionConvergence();
	void advectParticles();

	void virtual initMovingBoundary();
	void initMovingFluid();
	void moveBoundary();

	void initRotatingBoundary();
	void rotateBoundary();
};

class Solver : public Solver2D {
private:
	Grid grid;
	std::vector<RigidBody> rigidBodies;

public:
	Solver();

	void initBoundaries() override;
	void initFluid() override;
	void neighborSearch() override;

	// SPH functions
	void computeDensityPressure() override;
	void computeAcceleration() override;
	void updateParticles() override;

	// IISPH functions
	void boundaryMassUpdate() override;
	void computeDensity() override;
	void predictVelocity() override;
	void computeDensityError() override;
	void computeLaplacian() override;
	void compressionConvergence() override;

	void initMovingBoundary() override;

	// Ghost particles
	void initGhostFluid();
	void initGhostBoundary();
	void neighborSearchGhosts();
	void updateGhosts();

	// Rigid bodies kinematics
	void initRigidCube();
	void addRigidBody(std::vector < std::vector<Particle>> body);
	void boundaryToRigidForces();
	std::vector <std::vector<Particle>> sampleOBJ(std::string fileName);
	std::vector<RigidBody>& getRigidBodies() { return rigidBodies; }
};
