#pragma once
#include <vector>

#include "Parameters.h"

class Particle {
public:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d acceleration;

    double mass;
    double density;
    double pressure;

    bool isFluid;

    int ID;

    std::vector<Particle*> neighbors;

    Particle();
    Eigen::Vector3i getCellNumber();

    // IISPH
    Eigen::Vector3d predictedVelocity;
    double sourceTerm;
    double diagonal;
    Eigen::Vector3d pressureAcceleration;
};

class Particle2D {
public:
    Eigen::Vector2d position;
    Eigen::Vector2d velocity;
    Eigen::Vector2d acceleration;
    Eigen::Vector2d normal;

    double mass;
    double density;
    double pressure;

    bool isFluid;
    bool isSurface;

    int ID;

    std::vector<Particle2D*> neighbors;

    Particle2D();
    Eigen::Vector2i getCellNumber();

    // IISPH
    Eigen::Vector2d predictedVelocity;
    double sourceTerm;
    double diagonal;
    Eigen::Vector2d pressureAcceleration;
    Eigen::Vector2d dii;
    double aii;
    double predictedDensity;
    double predictedPressure;
    Eigen::Vector2d ci;
};
