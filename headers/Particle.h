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
    float coefficient;
    float diagonal;
    Eigen::Vector3f pressureAcceleration;
};
