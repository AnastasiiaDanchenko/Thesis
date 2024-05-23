#pragma once
#include <vector>

#include "Parameters.h"

class Particle {
public:
    Eigen::Vector3f position;
    Eigen::Vector3f velocity;
    Eigen::Vector3f acceleration;

    float mass;
    float density;
    float pressure;

    bool isFluid;

    int ID;

    std::vector<Particle*> neighbors;

    Particle();
    Eigen::Vector3i getCellNumber();

    // IISPH
    Eigen::Vector3f predictedVelocity;
    float sourceTerm;
    float diagonal;
    Eigen::Vector3f pressureAcceleration;
};

class Particle2D {
public:
    Eigen::Vector2f position;
    Eigen::Vector2f velocity;
    Eigen::Vector2f acceleration;

    float mass;
    float density;
    float pressure;

    bool isFluid;
    bool isSurface;

    int ID;

    std::vector<Particle2D*> neighbors;

    Particle2D();
    Eigen::Vector2i getCellNumber();

    // IISPH
    Eigen::Vector2f predictedVelocity;
    float sourceTerm;
    float diagonal;
    Eigen::Vector2f pressureAcceleration;
    Eigen::Vector2f dii;
    float aii;
    float predictedDensity;
    float predictedPressure;
    Eigen::Vector2f ci;
};
