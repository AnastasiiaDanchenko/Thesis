#pragma once
#include "Kernel.h"
#include <omp.h>

void CalculateDensity2D();
void PredictVelocityAdvection2D();
void ComputeDii();
void ComputeAii();
void PredictDensity2D();
void SolvePressure2D();
void ComputePressureAcceleration2D();
void AdvectParticles2D();
