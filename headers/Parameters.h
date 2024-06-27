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
extern int NB_FLUID_PARTICLES;
extern double SPACING;
extern double CELL_SIZE;

// SPH parameters
extern double SUPPORT;
extern double REST_DENSITY;
extern double TIME_STEP;
extern double STIFFNESS;
extern double VISCOSITY;
extern double COHESION;
extern Eigen::Vector3d GRAVITY;
extern Eigen::Vector2d GRAVITY2D;
extern double MAX_TIME_STEP;
extern int ITERATIONS_COUNT;

// IISPH parameters
extern double GAMMA;
extern double OMEGA;
extern double AVG_DENSITY;
extern double DENSITY_ERR;
extern double FIRST_ERR;
extern double ERR_THRESHOLD;
extern int NB_ITERATIONS;

extern std::string NS_METHOD;
extern std::string SIMULATION;
extern bool VISUALIZATION;
extern bool SURFACE_TENSION;
extern int DIMENSIONS;
extern bool FIRST_STEP;
extern double FIRST_STEP_CORRECTION;

void readParameters();
