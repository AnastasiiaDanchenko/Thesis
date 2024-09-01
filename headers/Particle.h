#pragma once
#include <vector>

#include "Parameters.h"

class RigidBody; // Forward declaration

class Particle {
public:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d acceleration;

    double mass;
    double density;
    double pressure;

    bool isFluid;
	bool isRigid;

    int ID;

    std::vector<Particle*> neighbors;

    Particle();

    // IISPH
    Eigen::Vector3d predictedVelocity;
    double sourceTerm;
    double diagonal;
    Eigen::Vector3d pressureAcceleration;

	// Rigid body
    RigidBody* parentBody;
	Eigen::Vector3d relativePosition;
    double artificialVolume;
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

    // IISPH
    Eigen::Vector2d predictedVelocity;
    double sourceTerm;
    double diagonal;
    Eigen::Vector2d pressureAcceleration;
};

extern std::vector<Particle> particles;
extern std::vector<Particle2D> particles2D;
extern std::vector<Particle> ghostParticles;
