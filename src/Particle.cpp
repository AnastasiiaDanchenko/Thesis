#include "..\headers\Particle.h"

std::vector<Particle> particles;
std::vector<Particle2D> particles2D;
std::vector<Particle> ghostParticles;

Particle::Particle() :
    position(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),
    velocity(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),
    acceleration(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),

    mass(parameters.spacing * parameters.spacing * parameters.spacing * parameters.restDensity),
    density(parameters.restDensity),
    pressure(0.0f),
    isFluid(true),
    ID(0),
    predictedVelocity(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),
    sourceTerm(0.0f),
    diagonal(0.0f),
    pressureAcceleration(Eigen::Vector3d(0.0f, 0.0f, 0.0f)) {}

Particle2D::Particle2D() :
	position(Eigen::Vector2d(0.0, 0.0)),
	velocity(Eigen::Vector2d(0.0, 0.0)),
	acceleration(Eigen::Vector2d(0.0, 0.0)),
	normal(Eigen::Vector2d(0.0, 0.0)),

	mass(parameters.spacing * parameters.spacing * parameters.restDensity),
	density(parameters.restDensity),
	pressure(0.0),
	isFluid(true),
	isSurface(false),
	ID(0),
	predictedVelocity(Eigen::Vector2d(0.0, 0.0)),
	sourceTerm(0.0),
	diagonal(0.0),
	pressureAcceleration(Eigen::Vector2d(0.0, 0.0)) {}
