#pragma once
#include "Kernel.h"
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
	void virtual advectParticles();
};

class Solver : public Solver2D {

};

void ComputeDensityPressure();
void ComputeAcceleration();
void UpdateParticles();

// IISPH functions
void ComputeDensity();
void PredictVelocity();
void ComputeDensityError();
void ComputeLaplacian();
void CompressionConvergence();

void BoundaryMassUpdate();

void RotateBoundary2D();

void UpdateGhosts();
