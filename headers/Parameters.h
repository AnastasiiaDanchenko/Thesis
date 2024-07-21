#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

const std::string fileName = "parameters.json";

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct WindowSize {
	int width;
	int height;
	int depth;
};

struct ParticlesPerDimension {
	int x;
	int y;
	int z;
};

struct Parameters {
	WindowSize windowSize;
	int dimensions;
	int visualizeNeighbors;
	ParticlesPerDimension particlesPerDimension;

	double spacing;
	double support;

	double restDensity;
	double timeStep;
	double stiffness;
	double viscosity;
	
	bool surfaceTension;
	double cohesion;

	double maxTimeStep;
	int iterationsCount;

	Eigen::Vector3d gravity;
	Eigen::Vector2d gravity2D;

	double gamma;
	double omega;
	double avgDensity;
	double densityErr;
	double firstErr;
	double errThreshold;
	int nbIterations;

	int boundaryTestID;

	int simulationType;
	float slicingPlane;

	Parameters();
	void readParameters();
};

extern Parameters parameters;
