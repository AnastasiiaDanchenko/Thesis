#pragma once
#include <chrono>

#include "ComputationsSPH.h"

void Initialization();
void Initialization2D(Solver2D& solver);

void MovingBoundaryInitialization(Solver2D& solver);
void RotatingBoundaryInitialization(Solver2D& solver);

void Simulation(Grid& grid);
void SimulationIISPH(Grid& grid);

void Simulation2D(Solver2D& solver);
void SimulationIISPH2D(Solver2D& solver);
void MovingBoundaryIISPH2D(Solver2D& solver);
void RotatingBoundaryIISPH2D(Solver2D& solver);
