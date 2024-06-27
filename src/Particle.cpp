#include "..\headers\Particle.h"

Particle::Particle() :
    position(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),
    velocity(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),
    acceleration(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),

    mass(SPACING * SPACING * SPACING * REST_DENSITY),
    density(REST_DENSITY),
    pressure(0.0f),
    isFluid(true),
    ID(0),
    predictedVelocity(Eigen::Vector3d(0.0f, 0.0f, 0.0f)),
    sourceTerm(0.0f),
    diagonal(0.0f),
    pressureAcceleration(Eigen::Vector3d(0.0f, 0.0f, 0.0f)) {}

Eigen::Vector3i Particle::getCellNumber() {
	int x = std::floor(position.x() / CELL_SIZE);
	int y = std::floor(position.y() / CELL_SIZE);
    int z = std::floor(position.z() / CELL_SIZE);
	return Eigen::Vector3i(x, y, z);
}

Particle2D::Particle2D() :
	position(Eigen::Vector2d(0.0, 0.0)),
	velocity(Eigen::Vector2d(0.0, 0.0)),
	acceleration(Eigen::Vector2d(0.0, 0.0)),
	normal(Eigen::Vector2d(0.0, 0.0)),

	mass(SPACING * SPACING * REST_DENSITY),
	density(REST_DENSITY),
	pressure(0.0),
	isFluid(true),
	isSurface(false),
	ID(0),
	predictedVelocity(Eigen::Vector2d(0.0, 0.0)),
	sourceTerm(0.0),
	diagonal(0.0),
	pressureAcceleration(Eigen::Vector2d(0.0, 0.0)),
	dii(Eigen::Vector2d(0.0f, 0.0f)),
	aii(0.0f),
	predictedDensity(0.0),
	predictedPressure(0.0),
	ci(Eigen::Vector2d(0.0f, 0.0f)) {}

Eigen::Vector2i Particle2D::getCellNumber() {
	int x = std::floor(position.x() / CELL_SIZE);
	int y = std::floor(position.y() / CELL_SIZE);
	return Eigen::Vector2i(x, y);
}