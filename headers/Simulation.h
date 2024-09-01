#pragma once
#include <chrono>

#include "Solvers.h"

// 3D
void Initialization(Solver& solver);

void MovingBoundaryInitialization(Solver2D& solver);
void RotatingBoundaryInitialization(Solver2D& solver);

void Simulation(Solver& solver);
void SimulationIISPH(Solver& solver);

// 2D
void Initialization2D(Solver2D& solver);
void Simulation2D(Solver2D& solver);
void SimulationIISPH2D(Solver2D& solver);
void MovingBoundaryIISPH2D(Solver2D& solver);
void RotatingBoundaryIISPH2D(Solver2D& solver);
