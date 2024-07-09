#pragma once
#include "Kernel.h"
#include <omp.h>
#include <limits>

void ComputeDensityPressure();
void ComputeAcceleration();
void UpdateParticles();

// IISPH functions
void ComputeDensity();
void PredictVelocity();
void ComputeDensityError();
void ComputeLaplacian();
void CompressionConvergence();

// 2D
void ComputeDensityPressure2D();
void ComputeAcceleration2D();
void Update2D();

// IISPH functions 2D
void ComputeDensity2D();
void ComputeSurface2D();
void PredictVelocity2D();
void ComputeDensityError2D();
void ComputeLaplacian2D();
void CompressionConvergence2D();
void UpdateParticles2D();

void BoundaryMassUpdate();
void BoundaryMassUpdate2D();

void RotateBoundary2D();
