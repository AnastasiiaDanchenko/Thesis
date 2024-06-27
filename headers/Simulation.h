#pragma once
#include <chrono>

#include "ComputationsSPH.h"
#include "ComputationsIISPH.h"
#include "NeighborSearch.h"
//#include "matplotlibcpp.h"

//namespace plt = matplotlibcpp;

void Initialization();
void Initialization2D();
void MovingBoundaryInitialization();

void Simulation();
void SimulationIISPH();
void Simulation2D();
void SimulationIISPH2D();
