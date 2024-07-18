#pragma once
#include <chrono>

#include "ComputationsSPH.h"

void Initialization();
void Initialization2D();

void MovingBoundaryInitialization();
void RotatingBoundaryInitialization();

void Simulation(Grid& grid);
void SimulationIISPH(Grid& grid);

void Simulation2D(Grid2D& grid2d);
void SimulationIISPH2D(Grid2D& grid2d);
void MovingBoundaryIISPH2D(Grid2D& grid2d);
void RotatingBoundaryIISPH2D(Grid2D& grid2d);
