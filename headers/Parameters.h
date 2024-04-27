#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>

// Read parameters for the SPH simulation from a .csv file
const std::string fileName = "parameters.csv";

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Window parameters
extern int WINDOW_WIDTH;
extern int WINDOW_HEIGHT;
extern int SCENE_DEPTH;
extern int PARTICLE_NEIGHBORS; // Visualized neighbors for a given particle

// Initial grid parameters
extern int PARTICLES_PER_DIMENSION;
extern float SPACING;
extern int CELL_SIZE;

// SPH parameters
extern float SUPPORT;
extern float REST_DENSITY;
extern float TIME_STEP;
extern float STIFFNESS;
extern float VISCOSITY;
extern Eigen::Vector3f GRAVITY;

// IISPH parameters
extern float GAMMA;
extern float OMEGA;

extern std::string NS_METHOD;
extern std::string SIMULATION;
extern bool VISUALIZATION;

void readParameters();
