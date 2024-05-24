#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>

// Read parameters for the SPH simulation from a .csv file
//const std::string fileName = "parameters.csv";
//const std::string fileName = "parametersIISPH2D.csv";
const std::string fileName = "parametersIISPH.csv";

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
extern int PARTICLES_X;
extern int PARTICLES_Y;
extern int PARTICLES_Z;
extern float SPACING;
extern int CELL_SIZE;

// SPH parameters
extern float SUPPORT;
extern float REST_DENSITY;
extern float TIME_STEP;
extern float STIFFNESS;
extern float VISCOSITY;
extern float COHESION;
extern Eigen::Vector3f GRAVITY;
extern Eigen::Vector2f GRAVITY2D;

// IISPH parameters
extern float GAMMA;
extern float OMEGA;
extern float AVG_DENSITY;
extern float DENSITY_ERR;
extern float ERR_THRESHOLD;
extern int NB_ITERATIONS;

extern std::string NS_METHOD;
extern std::string SIMULATION;
extern bool VISUALIZATION;
extern bool SURFACE_TENSION;
extern int DIMENSIONS;

void readParameters();
