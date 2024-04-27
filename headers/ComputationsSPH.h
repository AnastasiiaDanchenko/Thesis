#pragma once
#include "Kernel.h"

void ComputeDensityPressure();
void ComputeAcceleration();
void UpdateParticles();

// IISPH functions
void ComputeDensity();
void PredictVelocity();
void ComputeDensityError();
void ComputeLaplacian();
void CompressionConvergence();
