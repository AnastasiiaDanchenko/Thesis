#pragma once
#include <chrono>

#include "Solvers.h"

void Initialization(Solver& solver);
void Initialization2D(Solver2D& solver);

void MovingBoundaryInitialization(Solver2D& solver);
void RotatingBoundaryInitialization(Solver2D& solver);

void Simulation(Solver& solver);
void SimulationIISPH(Solver& solver);

void Simulation2D(Solver2D& solver);
void SimulationIISPH2D(Solver2D& solver);
void MovingBoundaryIISPH2D(Solver2D& solver);
void RotatingBoundaryIISPH2D(Solver2D& solver);
